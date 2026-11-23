# Step 5.0: Which pathways and immune cells separate the C5 ERBB2-low group from the luminal clusters (C1-C3 pooled)
# and from C5 ERBB2-high? limma on VST (all genes; TREAT for |log2 FC| > 1) -> GSEA (fgsea; Hallmark, Reactome) on the
# moderated t, checked with camera (sample-level test allowing for inter-gene correlation); with and without
# adjustment for ESTIMATE stromal and immune scores; immune/stromal cell abundance from MCP-counter markers.
# Run from repo root: Rscript analyses/05.0_c5_low_gsea.R
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(limma); library(fgsea); library(msigdbr); library(tidyestimate); library(tidyverse)
})
source("analyses/helpers.R")
out <- "results/05.0_c5_low_gsea"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Data: uncorrected VST, primary tumors, all 16,453 genes; current HGNC symbols (needed to match gene sets) ----
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
vsd <- vsd[, vsd$sample_type == "Primary Tumor"]
v <- assay(vsd); colnames(v) <- vsd$sample_id
rownames(v) <- suppressWarnings(current_symbols(rownames(v)))

# ESTIMATE stromal and immune scores (ssGSEA, rank-based) as purity covariates
est <- v |> filter_common_genes(id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE) |>
  estimate_score(is_affymetrix = FALSE) |> transmute(sample_id = sample, Stromal = stromal / 1000, Immune = immune / 1000)
lab <- read_csv("results/03.0_pam_clusters/clusters.csv", show_col_types = FALSE) |>
  left_join(est, by = "sample_id") |>
  mutate(grp = factor(if_else(cluster == "C5", paste0("C5_", c5_mode), cluster),
                      c("C1", "C2", "C3", "C4", "C5_ERBB2-high", "C5_ERBB2-low")))
stopifnot(identical(lab$sample_id, colnames(v)), !anyNA(lab$grp))
levels(lab$grp) <- make.names(levels(lab$grp))
n_lum <- table(lab$grp)[c("C1", "C2", "C3")]

# ---- 1. limma (trend, robust to genes with outlying variances) with each cluster as its own group, so C1-C3
# differences do not inflate the variance. C5-low vs C1-C3: pooled mean of all 718 C1-C3 tumors (clusters weighted by
# size). Two models: unadjusted, and adjusted for ESTIMATE stromal + immune scores (what differs beyond infiltration).
# Moderated t (eBayes) ranks genes for GSEA; significance uses TREAT, which tests |log2 FC| > 1 directly.
designs <- list(Unadjusted = model.matrix(~ 0 + grp, lab), `Purity-adjusted` = model.matrix(~ 0 + grp + Stromal + Immune, lab))
designs <- map(designs, \(X) { colnames(X) <- sub("^grp", "", colnames(X)); X })
w <- n_lum / sum(n_lum)
contr <- \(X) { cm <- makeContrasts(contrasts = c(sprintf("C5_ERBB2.low - (%f*C1 + %f*C2 + %f*C3)", w[1], w[2], w[3]),
                                                  "C5_ERBB2.low - C5_ERBB2.high"), levels = X)
  colnames(cm) <- c("C1-C3", "C5 ERBB2-high"); cm }
de <- imap_dfr(designs, \(X, adj) {
  cm <- contr(X); fit <- lmFit(v, X) |> contrasts.fit(cm)
  eb <- eBayes(fit, trend = TRUE, robust = TRUE); tr <- treat(fit, lfc = 1, trend = TRUE, robust = TRUE)
  map_dfr(colnames(cm), \(k) topTable(eb, coef = k, number = Inf, sort.by = "none") |> as_tibble(rownames = "gene") |>
            mutate(treat_p = tr$p.value[, k], treat_padj = p.adjust(treat_p, "BH"), comparator = k, adjustment = adj, .before = 1))
})
write_csv(de, file.path(out, "limma.csv"))
de |> group_by(comparator, adjustment) |> summarise(up = sum(treat_padj < 0.05 & logFC > 0), down = sum(treat_padj < 0.05 & logFC < 0),
                                                   .groups = "drop") |> print()

# Volcano plots (unadjusted model): significant = TREAT BH p < 0.05 (|log2 FC| > 1); labels = top 5 other genes each way
# by t, plus the 8 Step 4 compact-model genes (7 vs. C1-C3, STARD3 vs. C5 ERBB2-high) and ERBB2 highlighted
key <- read_csv("results/04.0_c5_low_enet/compact_model.csv", show_col_types = FALSE) |>
  filter(!term %in% c("(Intercept)", "Purity")) |>
  transmute(gene = symbol, key_lab = paste("Step 4 gene,", if_else(beta > 0, "higher", "lower"), "in C5 ERBB2-low vs.", comparator)) |>
  bind_rows(tibble(gene = "ERBB2")) |>     # ERBB2 shown with STARD3 (17q12 amplicon)
  mutate(key_lab = if_else(gene %in% c("ERBB2", "STARD3"), "ERBB2, and STARD3 (Step 4 gene, lower in C5 ERBB2-low vs. C5 ERBB2-high)", key_lab))
vd <- de |> filter(adjustment == "Unadjusted") |>
  mutate(dir = case_when(treat_padj >= 0.05 ~ "n.s.", logFC > 0 ~ "Higher in C5 ERBB2-low", TRUE ~ "Lower in C5 ERBB2-low")) |>
  group_by(comparator) |>
  mutate(comparator = sprintf("C5 ERBB2-low vs. %s (%d higher, %d lower)", comparator,
                              sum(dir == "Higher in C5 ERBB2-low"), sum(dir == "Lower in C5 ERBB2-low"))) |> ungroup()
vlab <- vd |> filter(dir != "n.s.", !gene %in% key$gene) |> group_by(comparator, dir) |> slice_max(abs(t), n = 5)
vkey <- vd |> inner_join(key, by = "gene")
p_vol <- ggplot(vd, aes(logFC, -log10(treat_padj))) + geom_point(aes(color = dir), size = 0.5, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3) +
  ggrepel::geom_text_repel(data = vlab, aes(label = gene, color = dir), size = 2.5, max.overlaps = Inf, show.legend = FALSE, segment.size = 0.2) +
  geom_point(data = vkey, aes(fill = key_lab), shape = 21, size = 2.6, color = "black") +
  ggrepel::geom_label_repel(data = vkey, aes(label = gene, fill = key_lab), size = 2.8, fontface = "bold", label.size = 0.2,
                            max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, show.legend = FALSE, box.padding = 0.9, force = 3,
                            ylim = c(5, NA)) +   # boxed labels kept off the crowded baseline (both panels reach well above 5)
  scale_color_manual(values = c(`Higher in C5 ERBB2-low` = "#b2182b", `Lower in C5 ERBB2-low` = "#2166ac", n.s. = "grey75"), name = NULL) +
  scale_fill_manual(values = c(`Step 4 gene, higher in C5 ERBB2-low vs. C1-C3` = "#fdb863",
                               `Step 4 gene, lower in C5 ERBB2-low vs. C1-C3` = "#b2abd2",
                               `ERBB2, and STARD3 (Step 4 gene, lower in C5 ERBB2-low vs. C5 ERBB2-high)` = "#66a61e"), name = NULL) +
  facet_wrap(~ comparator, scales = "free") + labs(x = "log2 fold change (VST)", y = "-log10 BH-adjusted p (TREAT, |log2 FC| > 1)") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)), fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out, "volcano.png"), p_vol, width = 10, height = 6, dpi = 200)

# ---- 2. GSEA: fgsea (multilevel, exact small p) on the moderated t; Hallmark (main) and Reactome (supplement) ----
sets <- list(Hallmark = msigdbr(species = "Homo sapiens", collection = "H"),
             Reactome = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")) |>
  map(\(m) split(m$gene_symbol, m$gs_name))
set.seed(1)
gsea <- de |> group_by(comparator, adjustment) |> group_modify(\(d, key) imap_dfr(sets, \(s, coll) {
  stat <- setNames(d$t, d$gene)
  fgsea(s, stat, minSize = 15, maxSize = 500, eps = 0) |> as_tibble() |> mutate(collection = coll, .before = 1)
})) |> ungroup() |>
  mutate(pathway = str_remove(pathway, "^(HALLMARK|REACTOME)_"), leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))
write_csv(gsea, file.path(out, "gsea.csv"))
gsea |> filter(collection == "Hallmark", padj < 0.05) |> select(comparator, adjustment, pathway, NES, padj) |>
  arrange(comparator, adjustment, -NES) |> print(n = Inf)

# Check: fgsea permutes genes, so with ~800 tumors nearly any coherent set is significant. camera tests each Hallmark set
# with the same design and contrasts, allowing for inter-gene correlation: (a) limma's default preset (0.01), and
# (b) each set's correlation estimated from the residuals (inter.gene.cor = NA; strict, since co-expressed programs
# such as proliferation have high correlation and so few effective genes).
idx <- ids2indices(sets$Hallmark, rownames(v))
cam <- imap_dfr(designs, \(X, adj) { cm <- contr(X)
  map_dfr(colnames(cm), \(k) map_dfr(c(Default = 0.01, Estimated = NA), .id = "correlation", \(rho)
    camera(v, idx, X, contrast = cm[, k], inter.gene.cor = rho, trend.var = TRUE) |>
      as_tibble(rownames = "pathway") |> select(pathway, any_of("Correlation"), Direction, PValue, FDR)) |>
    mutate(comparator = k, adjustment = adj, .before = 1)) }) |>
  mutate(pathway = str_remove(pathway, "^HALLMARK_"))
write_csv(cam, file.path(out, "camera_hallmark.csv"))
gsea |> filter(collection == "Hallmark") |> select(comparator, adjustment, pathway, NES, padj) |>
  left_join(cam |> pivot_wider(id_cols = c(comparator, adjustment, pathway), names_from = correlation, values_from = FDR),
            by = c("comparator", "adjustment", "pathway")) |>
  group_by(comparator, adjustment) |>
  summarise(fgsea = sum(padj < 0.05), camera_default = sum(Default < 0.05), camera_estimated = sum(Estimated < 0.05), .groups = "drop") |> print()

# Step 4 stable genes in each Hallmark leading edge (which programs the elastic net genes belong to)
stab <- read_csv("results/04.0_c5_low_enet/features.csv", show_col_types = FALSE) |> filter(sel_freq >= 0.7)
le_stable <- gsea |> filter(collection == "Hallmark", adjustment == "Unadjusted", padj < 0.05) |>
  mutate(stable_genes = map2_chr(str_split(leadingEdge, ";"), comparator,
                                 \(le, k) paste(intersect(le, stab$symbol[stab$comparator == k]), collapse = ";"))) |>
  filter(stable_genes != "") |> select(comparator, pathway, NES, padj, stable_genes)
write_csv(le_stable, file.path(out, "hallmark_stable_genes.csv")); print(le_stable)

# Figure: Hallmark NES, pathways with padj < 0.05 in any contrast; both adjustments side by side
hm <- gsea |> filter(collection == "Hallmark") |> group_by(pathway) |> filter(any(padj < 0.05)) |> ungroup() |>
  mutate(comparator = factor(paste("vs.", comparator), c("vs. C1-C3", "vs. C5 ERBB2-high")),
         adjustment = factor(adjustment, names(designs)), sig = padj < 0.05)
ord <- hm |> filter(comparator == "vs. C1-C3", adjustment == "Unadjusted") |> arrange(NES) |> pull(pathway)
p_hm <- ggplot(hm, aes(adjustment, factor(pathway, ord))) +
  geom_point(aes(size = -log10(padj), fill = NES, alpha = sig), shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "NES\n(+ = higher in\nC5 ERBB2-low)") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
  scale_size_continuous(range = c(1, 6), name = "-log10 adj. p") +
  facet_wrap(~ comparator) + labs(x = NULL, y = NULL) + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "hallmark_nes.png"), p_hm, width = 7.5, height = 0.17 * length(ord) + 1.8, dpi = 200)

# ---- 3. Immune and stromal cells: MCP-counter (Becht 2016) score = mean VST of each population's marker genes.
# Scores compare abundance between tumors (not fractions within a tumor). Group differences by linear model on the
# two groups compared, with HC3 robust SEs (groups of very different size and spread). The adjusted model adds the
# ESTIMATE immune + stromal scores: each population relative to the ESTIMATE signal, not a true composition.
mcp_genes <- mcp_markers() |> filter(gene %in% rownames(v))
mcp <- mcp_genes |> group_by(population) |> group_modify(\(g, k) tibble(sample_id = colnames(v), score = colMeans(v[g$gene, , drop = FALSE]))) |>
  ungroup() |> left_join(lab, by = "sample_id")
write_csv(mcp_genes |> count(population, name = "n_markers"), file.path(out, "mcp_counter_markers.csv"))
cmp_grp <- c(`C1-C3` = "C1|C2|C3", `C5 ERBB2-high` = "C5_ERBB2.high")
mcp_test <- expand_grid(population = unique(mcp$population), comparator = names(cmp_grp), adjustment = names(designs)) |>
  mutate(fit = pmap(list(population, comparator, adjustment), \(p, k, a) {
    d <- mcp |> filter(population == p, grp == "C5_ERBB2.low" | str_detect(grp, paste0("^(", cmp_grp[k], ")$"))) |>
      mutate(low = as.integer(grp == "C5_ERBB2.low"))
    f <- if (a == "Unadjusted") score ~ low else score ~ low + Stromal + Immune
    m <- lm(f, d); co <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))["low", ]
    tibble(diff = co[1], se = co[2], p = co[4])
  })) |> unnest(fit) |> group_by(adjustment) |> mutate(padj = p.adjust(p, "BH")) |> ungroup()
write_csv(mcp_test, file.path(out, "mcp_counter_tests.csv")); print(mcp_test, n = Inf)
write_csv(mcp |> group_by(population, grp) |> summarise(n = n(), median = median(score), .groups = "drop"),
          file.path(out, "mcp_counter_by_group.csv"))

# Figure: MCP-counter scores by group (z-scored across tumors), cluster colors from Step 3 (C5-low lighter green)
pal <- c(C1 = "#1b9e77", C2 = "#d95f02", C3 = "#7570b3", C4 = "#e7298a", `C5 ERBB2-high` = "#66a61e", `C5 ERBB2-low` = "#b8dc8f")
p_mcp <- mcp |> group_by(population) |> mutate(z = (score - mean(score)) / sd(score)) |> ungroup() |>
  mutate(grp = factor(str_replace(str_replace(grp, "_", " "), "\\.", "-"), names(pal))) |>
  ggplot(aes(grp, z, fill = grp)) + geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = pal, guide = "none") + facet_wrap(~ population, nrow = 2) +
  labs(x = NULL, y = "MCP-counter score (z)") + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out, "mcp_counter.png"), p_mcp, width = 10, height = 5.5, dpi = 200)

# ---- 4. Is low ER signalling explained by ER− tumors or weak ER staining? IHC ER+ tumors only ----
# ESR1 and an early-estrogen-response score (mean z of the Hallmark genes across primary tumors) in IHC ER+ tumors;
# ER % positive where recorded (file bins: <10% = ER-low range; >= 80% = strong).
er_genes <- intersect(sets$Hallmark$HALLMARK_ESTROGEN_RESPONSE_EARLY, rownames(v))
er_score <- colMeans((v[er_genes, ] - rowMeans(v[er_genes, ])) / matrixStats::rowSds(v[er_genes, ]))
er <- read_csv("../data_dir/brca_tcga_clinical_data.csv", show_col_types = FALSE) |>
  select(sample_id = `Sample ID`, er = `ER Status By IHC`, pct = `ER Status IHC Percent Positive`) |>
  inner_join(lab |> select(sample_id, grp, Stromal, Immune), by = "sample_id") |>
  filter(er == "Positive", grp %in% c("C1", "C2", "C3", "C5_ERBB2.low")) |>
  mutate(group = if_else(grp == "C5_ERBB2.low", "C5 ERBB2-low", "C1-C3"), ESR1 = v["ESR1", sample_id], er_score = er_score[sample_id],
         recorded = str_detect(coalesce(pct, ""), "%"), low = recorded & pct == "<10%", high = recorded & pct %in% c("80-89%", "90-99%"))
er_out <- er |> group_by(group) |>
  summarise(n_er_pos = n(), ESR1_median = median(ESR1), er_score_median = median(er_score),
            n_pct_recorded = sum(recorded), below_10 = sum(low), at_least_80 = sum(high)) |>
  mutate(wilcox_p_ESR1 = wilcox.test(ESR1 ~ group, er)$p.value, wilcox_p_er_score = wilcox.test(er_score ~ group, er)$p.value,
         fisher_p_below_10 = fisher.test(table(er$group[er$recorded], er$low[er$recorded]))$p.value,
         fisher_p_at_least_80 = fisher.test(table(er$group[er$recorded], er$high[er$recorded]))$p.value)
write_csv(er_out, file.path(out, "er_percent.csv")); print(er_out)
# Two follow-ups within IHC ER+ tumors: (a) adjusted for ESTIMATE stromal + immune (linear model, C5-low coefficient);
# (b) only tumors staining >= 80% (Wilcoxon)
er_adj <- map_dfr(c(ESR1 = "ESR1", er_score = "er_score"), .id = "measure", \(y) {
  co <- summary(lm(reformulate(c("group", "Stromal", "Immune"), y), er))$coefficients["groupC5 ERBB2-low", ]
  s80 <- filter(er, high)
  tibble(adj_diff = co[1], adj_p = co[4], n80_low = sum(s80$group == "C5 ERBB2-low"), n80_lum = sum(s80$group == "C1-C3"),
         median80_low = median(s80[[y]][s80$group == "C5 ERBB2-low"]), median80_lum = median(s80[[y]][s80$group == "C1-C3"]),
         wilcox_p80 = wilcox.test(reformulate("group", y), s80)$p.value) })
write_csv(er_adj, file.path(out, "er_pos_followup.csv")); print(er_adj)

writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
