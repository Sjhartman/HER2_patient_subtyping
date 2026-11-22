# Step 3.0: Unsupervised subsets of primary tumors (PAM on top variable genes), compared with PAM50,
# HER2 status, and processing batch
# Run from repo root: Rscript analyses/03.0_pam_clusters.R
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(matrixStats); library(cluster)
  library(ConsensusClusterPlus); library(genefu); library(mclust); library(patchwork)
  library(tidyverse)   # tidyverse last: genefu masks select, mclust masks map
})
set.seed(1)
out <- "results/03.0_pam_clusters"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Data: uncorrected VST, primary tumors only (no plate/era correction; see Step 1.3) ----
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
vsd <- vsd[, vsd$sample_type == "Primary Tumor"]
vmat <- assay(vsd); colnames(vmat) <- vsd$sample_id
lab <- read_csv("results/02.0_her2_labels/her2_labels.csv", show_col_types = FALSE) |>
  mutate(Hist = case_when(Histology == "Breast Invasive Ductal Carcinoma" ~ "Ductal",
                          Histology == "Breast Invasive Lobular Carcinoma" ~ "Lobular", .default = "Other"))
stopifnot(identical(lab$sample_id, colnames(vmat)))

# Top variable genes, re-ranked on tumors only; genes mean-centered; distance = 1 - Pearson between tumors
top <- \(n) { m <- vmat[order(-rowVars(vmat))[1:n], ]; m - rowMeans(m) }
pam_cl <- \(m, k) pam(as.dist(1 - cor(m)), k, diss = TRUE)
x <- top(2000)

# ---- 1. Choose k: average silhouette (PAM on all tumors) + consensus stability (PAC) ----
ks <- 2:8
sil <- map_dbl(ks, \(k) pam_cl(x, k)$silinfo$avg.width)
# Consensus clustering: 200 resamples of 80% tumors and 80% genes, PAM on 1 - Pearson
cc_file <- file.path(out, "consensus.rds")   # ~8 min; cached (gitignored)
if (file.exists(cc_file)) cc <- readRDS(cc_file) else {
  cc <- ConsensusClusterPlus(x, maxK = max(ks), reps = 200, pItem = 0.8, pFeature = 0.8, clusterAlg = "pam",
                             distance = "pearson", seed = 1, plot = NULL, verbose = FALSE, title = tempdir())
  saveRDS(cc, cc_file)
}
# PAC = share of tumor pairs with ambiguous consensus (between 0.1 and 0.9); lower = more stable
pac <- map_dbl(ks, \(k) { cm <- cc[[k]]$consensusMatrix; v <- cm[lower.tri(cm)]; mean(v > 0.1 & v < 0.9) })
ksel <- tibble(k = ks, silhouette = sil, PAC = pac)
write_csv(ksel, file.path(out, "k_selection.csv")); print(ksel)

k <- 5   # chosen from k_selection.csv (see report Section 3.1)

# ---- 2. Final clusters; relabel C1..Ck by size ----
fit <- pam_cl(x, k)
relab <- \(cl) { r <- rank(-table(cl), ties.method = "first"); paste0("C", r[as.character(cl)]) }
lab$cluster <- relab(fit$clustering)
lab$silhouette <- fit$silinfo$widths[colnames(vmat), "sil_width"]
lab$consensus <- relab(cc[[k]]$consensusClass[colnames(vmat)])

# ---- 3. Robustness: top 1000 / 5000 genes, consensus clusters, k = 6 / 7 ----
# ARI vs. main (same k only); per-cluster recovery = % of a main cluster's tumors in its best-matching cluster
robust <- tibble(version = c("Top 1000 genes", "Top 5000 genes", "Consensus (top 2000)", "k = 6", "k = 7"),
                 cl = list(pam_cl(top(1000), k)$clustering, pam_cl(top(5000), k)$clustering, lab$consensus,
                           pam_cl(x, 6)$clustering, pam_cl(x, 7)$clustering)) |>
  mutate(ARI = map_dbl(cl, \(c) if (n_distinct(c) == k) mclust::adjustedRandIndex(c, lab$cluster) else NA),
         recovery = map(cl, \(c) { tb <- table(lab$cluster, c); round(100 * apply(tb, 1, max) / rowSums(tb)) |> as.list() })) |>
  unnest_wider(recovery) |> select(-cl)
write_csv(robust, file.path(out, "robustness.csv")); print(robust)

# ---- 4. PAM50 reference labels (genefu, robust centroids; aliases CDCA1 = NUF2, KNTC2 = NDC80) ----
data(pam50.robust, package = "genefu")
p50 <- t(vmat[intersect(c(pam50.robust$centroids.map$probe, "NUF2", "NDC80"), rownames(vmat)), ])
colnames(p50) <- recode(colnames(p50), NUF2 = "CDCA1", NDC80 = "KNTC2")
stopifnot(all(pam50.robust$centroids.map$probe %in% colnames(p50)))
lab$PAM50 <- intrinsic.cluster.predict(pam50.robust, p50, data.frame(probe = colnames(p50)), do.mapping = FALSE)$subtype |>
  as.character() |> factor(c("LumA", "LumB", "Her2", "Basal", "Normal"))
ari50 <- mclust::adjustedRandIndex(lab$cluster, lab$PAM50)
write_csv(tibble(comparison = "PAM cluster vs PAM50", ARI = ari50), file.path(out, "ari_pam50.csv"))
cat("ARI cluster vs PAM50:", round(ari50, 2), "\n")

# ---- 5. Cluster vs. biology, HER2 and batch ----
lab <- lab |> mutate(
  HER2 = factor(case_when(her2_status == "Positive" ~ "Positive", her2_low == "HER2-low (2+/FISH-)" ~ "HER2-low 2+",
                          her2_low == "HER2-low (1+)" ~ "HER2-low 1+", her2_low == "HER2-0" ~ "HER2-0",
                          her2_status == "Negative" ~ "Negative, no score"),
                c("HER2-0", "HER2-low 1+", "HER2-low 2+", "Negative, no score", "Positive")),
  GISTIC = factor(gistic), ERBB2_mRNA = ERBB2)
cat_vars <- c("PAM50", "her2_status", "HER2", "GISTIC", "ER", "PR", "Hist", "Era")
# Bias-corrected Cramer's V (Bergsma 2013): plain V is inflated for sparse tables with many levels (e.g. plate)
cramer <- \(tb) {
  tb <- tb[rowSums(tb) > 0, colSums(tb) > 0, drop = FALSE]
  n <- sum(tb); r <- nrow(tb); k <- ncol(tb)
  phi2 <- max(0, unname(suppressWarnings(chisq.test(tb))$statistic) / n - (k - 1) * (r - 1) / (n - 1))
  rt <- r - (r - 1)^2 / (n - 1); kt <- k - (k - 1)^2 / (n - 1)
  sqrt(phi2 / min(kt - 1, rt - 1))
}
assoc <- bind_rows(
  map_dfr(cat_vars, \(v) {
    # Drop non-calls; "Negative, no score" tracks label source/site rather than biology
    a <- lab |> filter(!is.na(.data[[v]]), !.data[[v]] %in% c("Missing", "Conflict", "Unresolved", "Negative, no score"))
    tb <- table(a$cluster, as.character(a[[v]]))   # as.character drops empty factor levels
    tibble(variable = v, n = sum(tb), test = "Chi-square (simulated p)", effect = "Cramer's V", value = cramer(tb),
           p = chisq.test(tb, simulate.p.value = TRUE, B = 1e4)$p.value)
  }),
  map_dfr(c("ERBB2_mRNA", "estimate_score", "Age"), \(v) {
    a <- lab |> filter(!is.na(.data[[v]]))
    kw <- kruskal.test(a[[v]], a$cluster)
    tibble(variable = v, n = nrow(a), test = "Kruskal-Wallis", effect = "epsilon^2",
           value = unname(kw$statistic) / (nrow(a) - 1), p = kw$p.value)
  })) |> mutate(padj = p.adjust(p, "BH"))
write_csv(assoc, file.path(out, "cluster_assoc.csv")); print(assoc)

# Per-cluster summary
summ <- lab |> group_by(cluster) |>
  summarise(n = n(), silhouette = median(silhouette),
            pct_HER2pos = 100 * mean(her2_status == "Positive") / mean(her2_status %in% c("Positive", "Negative")),
            pct_amplified = 100 * mean(gistic == 2, na.rm = TRUE), ERBB2_median = median(ERBB2),
            pct_ERpos = 100 * mean(ER == "Positive", na.rm = TRUE), pct_lobular = 100 * mean(Hist == "Lobular"),
            estimate_median = median(estimate_score), top_PAM50 = names(which.max(table(PAM50))),
            pct_top_PAM50 = 100 * max(table(PAM50)) / n())
write_csv(summ, file.path(out, "cluster_summary.csv")); print(summ)
write_csv(count(lab, cluster, PAM50) |> pivot_wider(names_from = PAM50, values_from = n, values_fill = 0),
          file.path(out, "cluster_vs_pam50.csv"))
write_csv(count(lab, cluster, HER2) |> pivot_wider(names_from = HER2, values_from = n, values_fill = 0),
          file.path(out, "cluster_vs_her2.csv"))

# Batch: era vs. cluster within histology (generalized CMH, ductal + lobular strata), then also within ER status.
# Late-era ductal tumors are more often ER- by IHC (independent of RNA): case mix. Era tests within ER strata
# show how much remains after accounting for it. Plate within ductal only.
b <- lab |> filter(Hist %in% c("Ductal", "Lobular"))
be <- b |> filter(!is.na(ER))
cmh <- mantelhaen.test(table(b$cluster, b$Era, b$Hist))
cmh_er <- mantelhaen.test(table(be$cluster, be$Era, interaction(be$Hist, be$ER, drop = TRUE)))
duct <- lab |> filter(Hist == "Ductal")
erp <- duct |> filter(ER == "Positive"); ern <- duct |> filter(ER == "Negative")
sim_p <- \(tb) chisq.test(tb, simulate.p.value = TRUE, B = 1e4)$p.value
batch <- tibble(test = c("Era vs cluster | histology (CMH)", "Era vs cluster | histology x ER (CMH)",
                         "Era vs cluster, ductal only (chi-square, simulated p)",
                         "Era vs cluster, ER+ ductal (chi-square, simulated p)",
                         "Era vs cluster, ER- ductal (chi-square, simulated p)",
                         "Plate vs cluster, ductal only (chi-square, simulated p)"),
                p = c(cmh$p.value, cmh_er$p.value, sim_p(table(duct$cluster, duct$Era)),
                      sim_p(table(erp$cluster, erp$Era)), sim_p(table(ern$cluster, ern$Era)),
                      sim_p(table(duct$cluster, duct$Plate))),
                cramers_v = c(NA, NA, cramer(table(duct$cluster, duct$Era)), cramer(table(erp$cluster, erp$Era)),
                              cramer(table(ern$cluster, ern$Era)), cramer(table(duct$cluster, duct$Plate))))
write_csv(erp |> count(Era, cluster) |> group_by(Era) |> mutate(pct = round(100 * n / sum(n), 1)),
          file.path(out, "cluster_by_era_er_pos_ductal.csv"))
write_csv(duct |> filter(!is.na(ER)) |> count(Era, ER) |> group_by(Era) |> mutate(pct = round(100 * n / sum(n), 1)),
          file.path(out, "er_by_era_ductal.csv"))
write_csv(batch, file.path(out, "batch_tests.csv")); print(batch)
write_csv(count(lab, Hist, Era, cluster) |> group_by(Hist, Era) |> mutate(pct = round(100 * n / sum(n), 1)),
          file.path(out, "cluster_by_era_histology.csv"))

# ---- 6. ADC candidates: where do HER2-low and amplified IHC-negative tumors fall? ----
# % HER2-low among HER2-negative tumors with a numeric IHC score (HER2-low needs the score; "no score" excluded)
scored_neg <- lab |> filter(label_source == "Numeric IHC score", HER2 %in% c("HER2-0", "HER2-low 1+", "HER2-low 2+"))
low_pct <- \(g) g |> summarise(n_scored_negative = n(), `HER2-0` = sum(HER2 == "HER2-0"),
                               `HER2-low 1+` = sum(HER2 == "HER2-low 1+"), `HER2-low 2+` = sum(HER2 == "HER2-low 2+"),
                               pct_HER2_low = round(100 * mean(HER2 != "HER2-0"), 1))
lowtab <- bind_rows(scored_neg |> group_by(group = cluster) |> low_pct(), scored_neg |> low_pct() |> mutate(group = "All"))
write_csv(lowtab, file.path(out, "her2_low_by_cluster.csv")); print(lowtab)
low_p <- chisq.test(table(scored_neg$cluster, scored_neg$HER2 != "HER2-0"), simulate.p.value = TRUE, B = 1e4)$p.value
cat("HER2-low vs HER2-0 across clusters (chi-square, simulated p):", signif(low_p, 2), "\n")
write_csv(tibble(test = "HER2-low vs HER2-0 across clusters (chi-square, simulated p)", n = nrow(scored_neg), p = low_p),
          file.path(out, "her2_low_test.csv"))
# Selection check: HER2- tumors without a numeric score (reported Negative or FISH only) vs. scored 0 / 1+ / 2+.
# If unscored negatives sit at HER2-0-like mRNA, the scored subset over-represents HER2-low.
lab |> filter(her2_status == "Negative") |> group_by(her2_low) |>
  summarise(n = n(), ERBB2_median = median(ERBB2), ERBB2_q1 = quantile(ERBB2, 0.25), ERBB2_q3 = quantile(ERBB2, 0.75)) |>
  write_csv(file.path(out, "her2_negative_mrna.csv"))
lab |> filter(sample_id %in% c("TCGA-A2-A0YJ-01", "TCGA-A8-A08J-01")) |>
  select(sample_id, cluster, PAM50, ERBB2, gistic, score) |> print()

# ---- 7. C5 (HER2-enriched) ERBB2 mRNA is bimodal: split at the density valley and compare metadata ----
c5 <- lab |> filter(cluster == "C5")
dn <- density(c5$ERBB2); w <- dn$x > 13 & dn$x < 16
valley <- dn$x[w][which.min(dn$y[w])]
c5 <- c5 |> mutate(mode = if_else(ERBB2 > valley, "ERBB2-high", "ERBB2-low"))
cat("C5 ERBB2 valley:", round(valley, 2), "\n")
# Split checks: valley under other bandwidths; 2- vs 1-component Gaussian mixture (BIC, higher = better in mclust)
bw_valley <- map_dbl(c(0.5, 0.75, 1, 1.25, 1.5), \(a) { d <- density(c5$ERBB2, adjust = a); w <- d$x > 13 & d$x < 16
  d$x[w][which.min(d$y[w])] })
gm <- mclust::Mclust(c5$ERBB2, G = 1:2, modelNames = "V", verbose = FALSE)
split_chk <- tibble(check = c(paste("Valley, bandwidth adjust", c(0.5, 0.75, 1, 1.25, 1.5)),
                              "Gaussian mixture BIC, 1 component", "Gaussian mixture BIC, 2 components",
                              "Gaussian mixture, tumors in lower component"),
                    value = c(bw_valley, gm$BIC[1:2], sum(gm$classification == which.min(gm$parameters$mean))))
write_csv(split_chk, file.path(out, "c5_split_checks.csv")); print(split_chk)
# How well each mode fits C5: PAM silhouette, nearest other cluster, consensus agreement; ER and purity contrasts
nb <- relab(fit$clustering)[match(fit$silinfo$widths[, "neighbor"], fit$clustering)]   # neighbor cluster, relabelled
names(nb) <- rownames(fit$silinfo$widths)
c5 <- c5 |> mutate(neighbor = nb[sample_id])
c5_fit <- c5 |> group_by(mode) |>
  summarise(n = n(), median_silhouette = median(silhouette), pct_negative_silhouette = 100 * mean(silhouette < 0),
            pct_consensus_C5 = 100 * mean(consensus == "C5"),
            neighbor = paste(names(sort(table(neighbor), decreasing = TRUE)), sort(table(neighbor), decreasing = TRUE),
                             sep = ": ", collapse = ", "),
            pct_ERpos = 100 * mean(ER == "Positive", na.rm = TRUE), estimate_median = median(estimate_score))
c5_tests <- tibble(test = c("ER+ by mode (Fisher)", "ESTIMATE by mode (Wilcoxon)"),
                   p = c(fisher.test(table(c5$mode, c5$ER))$p.value, wilcox.test(estimate_score ~ mode, c5)$p.value))
write_csv(c5_fit, file.path(out, "c5_mode_fit.csv")); write_csv(c5_tests, file.path(out, "c5_mode_tests.csv"))
print(c5_fit); print(c5_tests)
c5_vars <- c("gistic", "her2_status", "HER2", "PAM50", "ER", "PR", "Hist", "Era", "label_source")
c5_tab <- map_dfr(c5_vars, \(v) c5 |> count(mode, value = as.character(.data[[v]])) |> mutate(variable = v)) |>
  pivot_wider(names_from = mode, values_from = n, values_fill = 0) |> relocate(variable)
write_csv(c5_tab, file.path(out, "c5_erbb2_modes.csv"))
c5 |> group_by(mode) |> summarise(n = n(), across(c(ERBB2, GRB7, STARD3, ESR1, PGR, estimate_score, Age), \(x) median(x, na.rm = TRUE))) |>
  write_csv(file.path(out, "c5_erbb2_modes_continuous.csv"))
# Clinically HER2+ tumors in the ERBB2-low mode (label-noise check against Step 2 site flags)
c5 |> filter(mode == "ERBB2-low", her2_status == "Positive") |>
  select(sample_id, TSS, ERBB2, gistic, score, fish, label_source, rule, PAM50) |>
  write_csv(file.path(out, "c5_low_her2pos.csv"))
# Purity and infiltration of the two C5 groups (descriptive; nothing regressed out). Stroma / immune scores:
# mean z (across all tumors) of COL1A1, FAP, COL3A1 / PTPRC, CD8A, CD3E. Spearman rho: ERBB2 mRNA vs ESTIMATE.
mz <- t(scale(t(vmat[c("COL1A1", "FAP", "COL3A1", "PTPRC", "CD8A", "CD3E"), ])))
infil <- tibble(sample_id = colnames(vmat), stroma = colMeans(mz[1:3, ]), immune = colMeans(mz[4:6, ]))
c5 |> left_join(infil, by = "sample_id") |> group_by(mode) |>
  summarise(n = n(), estimate_median = median(estimate_score), stroma_median = median(stroma),
            immune_median = median(immune), rho_ERBB2_ESTIMATE = cor(ERBB2, estimate_score, method = "spearman")) |>
  write_csv(file.path(out, "c5_purity.csv"))
lab <- lab |> left_join(c5 |> select(sample_id, c5_mode = mode), by = "sample_id")
cat("Scored HER2-negative C5 tumors by mode:\n"); print(count(filter(lab, sample_id %in% scored_neg$sample_id, cluster == "C5"), c5_mode))

write_csv(lab |> select(sample_id, cluster, silhouette, consensus, PAM50, HER2, c5_mode), file.path(out, "clusters.csv"))

# ---- 8. Figures ----
# Cluster display labels: cluster + its most common PAM50 subtype, e.g. "C1 (LumB)"
cl_lab <- lab |> count(cluster, PAM50) |> slice_max(n, by = cluster, with_ties = FALSE) |>
  mutate(label = paste0(cluster, " (", PAM50, ")")) |> select(cluster, label) |> deframe()
lab <- lab |> mutate(cl_lab = factor(cl_lab[cluster], cl_lab[sort(names(cl_lab))]))
pal <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")[1:k],
                levels(lab$cl_lab))
# One palette per feature, shared across Figures 6-8 and distinct between features:
# clusters (Dark2 hues), PAM50 (browns/yellow/plum/aqua), HER2 category (grey ramp to black), GISTIC (blue-red)
pal50 <- c(LumA = "#8c510a", LumB = "#dfc27d", Her2 = "#ffd92f", Basal = "#40004b", Normal = "#c7eae5")
pal_her2 <- c(`HER2-0` = "#f0f0f0", `HER2-low 1+` = "#bdbdbd", `HER2-low 2+` = "#737373",
              `Negative, no score` = "#d9f0d3", Positive = "#000000")
pal_cn <- c(`-2` = "#2166ac", `-1` = "#92c5de", `0` = "#d9d9d9", `1` = "#f4a582", `2` = "#b2182b")

# k selection
p_k <- ksel |> pivot_longer(-k) |> ggplot(aes(k, value)) + geom_line() + geom_point() +
  geom_vline(xintercept = k, linetype = 2, color = "grey50") +
  facet_wrap(~name, scales = "free_y", labeller = as_labeller(c(PAC = "PAC (lower = more stable)",
                                                                 silhouette = "Average silhouette width"))) +
  labs(x = "Number of clusters (k)", y = NULL) + theme_bw()
ggsave(file.path(out, "k_selection.png"), p_k, width = 7, height = 3)

# PCA (tumor PC1/PC2 from Step 1) by cluster (labelled with majority PAM50)
pcs <- read_csv("results/01.0_qc_vst_pca/tumor_pcs_features.csv", show_col_types = FALSE) |>
  select(sample_id, PC1, PC2) |> inner_join(lab, by = "sample_id")
p_pca <- ggplot(pcs, aes(PC1, PC2, color = cl_lab)) + geom_point(size = 0.7) + scale_color_manual(values = pal) +
  labs(color = "PAM cluster\n(majority PAM50)") + theme_bw() + guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(out, "pca_clusters.png"), p_pca, width = 5.5, height = 4)

# Marker heatmap: z-scored VST, tumors ordered by cluster then ERBB2; annotation strips
markers <- list(HER2 = c("ERBB2", "GRB7", "STARD3"), Luminal = c("ESR1", "PGR", "FOXA1", "GATA3"),
                Basal = c("KRT5", "KRT14", "KRT17", "FOXC1"), Proliferation = c("MKI67", "CCNB1", "AURKA"),
                Immune = c("PTPRC", "CD8A"), Stroma = c("COL1A1", "FAP"))
ord <- lab |> arrange(cluster, -ERBB2) |> pull(sample_id)
# Heatmap cluster labels name the dominant PAM50 type, e.g. "C1 (LumB dominant)"
dom <- \(x) factor(str_replace(x, "\\)", " dominant)"), str_replace(levels(lab$cl_lab), "\\)", " dominant)"))
hm <- t(scale(t(vmat[unlist(markers), ord]))) |> pmin(3) |> pmax(-3) |> as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "z") |>
  left_join(enframe(markers, "set", "gene") |> unnest(gene), by = "gene") |>
  left_join(lab |> transmute(sample_id, cluster = dom(cl_lab)), by = "sample_id") |>
  mutate(sample_id = factor(sample_id, ord), gene = factor(gene, rev(unlist(markers))), set = fct_inorder(set))
p_hm <- ggplot(hm, aes(sample_id, gene, fill = z)) + geom_raster() +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", name = "z") +
  facet_grid(set ~ cluster, scales = "free", space = "free") +
  theme_minimal() + theme(axis.text.x = element_blank(), panel.spacing = unit(1, "pt"), strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = NULL)
# Annotation tracks: one fill scale per track so each keeps its own palette and legend
ann <- lab |> mutate(sample_id = factor(sample_id, ord), cluster = dom(cl_lab), cl_dom = cluster)
track <- \(var, name, values, title) {
  ggplot(ann, aes(sample_id, name, fill = .data[[var]])) + geom_raster() +
    facet_grid(~cluster, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = values, na.value = "white", name = title, drop = FALSE) +
    theme_minimal() + theme(axis.text.x = element_blank(), panel.spacing = unit(1, "pt"), strip.text = element_blank(),
                            legend.key.size = unit(8, "pt"), legend.text = element_text(size = 7),
                            legend.title = element_text(size = 8)) +
    labs(x = NULL, y = NULL)
}
ann <- ann |> mutate(gistic_f = factor(gistic, names(pal_cn)))
p_ann <- (track("cl_dom", "PAM cluster", setNames(pal, levels(ann$cl_dom)), "PAM cluster") + theme(strip.text = element_text(size = 8.8))) /
  track("PAM50", "PAM50", pal50, "PAM50") /
  track("HER2", "HER2 (clinical)", pal_her2, "HER2 (clinical)\nwhite = no call") /
  track("gistic_f", "ERBB2 GISTIC", pal_cn, "ERBB2 GISTIC\nwhite = NA")
ggsave(file.path(out, "marker_heatmap.png"),
       p_ann / p_hm + plot_layout(heights = c(0.3, 0.25, 0.25, 0.25, 6), guides = "collect"), width = 11, height = 8)

# HER2 categories and ERBB2 mRNA by cluster
# Right panel: points coloured by ERBB2 GISTIC (explains the bimodal C5 distribution); dashed line = C5 valley
# Left: HER2 category (same palette as Figure 7), with a cluster-colour strip under each bar (same as Figures 6-7)
# Right: points filled by ERBB2 GISTIC (same palette as Figure 7; NA = white), explaining the bimodal C5
# Left fill: Positive = cluster colour, HER2-low 2+ = lighter cluster colour, then greys; stacked top to bottom
# Positive, 2+, 1+, HER2-0, Negative no score. Legend keys use C5's colours as the example.
light <- \(col) colorRampPalette(c(col, "white"))(5)[3]
her2_lv <- c("Positive", "HER2-low 2+", "HER2-low 1+", "HER2-0", "Negative, no score")
bar_fill <- c(setNames(pal, paste(names(pal), "Positive")), setNames(map_chr(pal, light), paste(names(pal), "HER2-low 2+")),
              `HER2-low 1+` = "#636363", `HER2-0` = "#bdbdbd", `Negative, no score` = "white")
c5_lab <- levels(lab$cl_lab)[5]
bars <- lab |> filter(!is.na(HER2)) |>
  mutate(HER2 = factor(HER2, her2_lv),
         key = if_else(HER2 %in% c("Positive", "HER2-low 2+"), paste(cl_lab, HER2), as.character(HER2)),
         key = factor(key, c(as.vector(outer(levels(lab$cl_lab), her2_lv[1:2], paste)), her2_lv[3:5])))
bars <- bars |> mutate(key = fct_relevel(key, unlist(map(her2_lv, \(h) grep(h, levels(key), value = TRUE, fixed = TRUE)))))
ymin <- min(lab$ERBB2) - 0.6
p_her2 <- (ggplot(bars, aes(cl_lab, fill = key)) + geom_bar(position = "fill", color = "grey40", linewidth = 0.2) +
             scale_fill_manual(values = bar_fill, breaks = c(paste(c5_lab, her2_lv[1:2]), her2_lv[3:5]),
                               labels = her2_lv, name = "HER2 (clinical)\nPositive / 2+: cluster colour") +
             labs(x = NULL, y = "Fraction of tumors")) +
  (ggplot(lab, aes(cl_lab, ERBB2)) + geom_boxplot(outlier.shape = NA, color = "black") +
     geom_jitter(aes(fill = factor(gistic, names(pal_cn))), shape = 21, stroke = 0.15, color = "grey40",
                 width = 0.2, height = 0, size = 0.9) +
     geom_hline(yintercept = valley, linetype = 2, color = "grey50") +
     scale_fill_manual(values = pal_cn, na.value = "white", name = "ERBB2 GISTIC\nwhite = NA",
                       guide = guide_legend(override.aes = list(size = 2.5))) +
     ggnewscale::new_scale_fill() +
     geom_tile(data = distinct(lab, cl_lab), aes(y = ymin, fill = cl_lab), height = 0.3, show.legend = FALSE) +
     scale_fill_manual(values = pal) +
     labs(x = NULL, y = "ERBB2 mRNA (VST)")) &
  theme_bw() & theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "her2_by_cluster.png"), p_her2, width = 10, height = 4)
