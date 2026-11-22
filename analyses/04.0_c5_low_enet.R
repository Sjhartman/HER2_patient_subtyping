# Step 4.0: What separates the C5 ERBB2-low group from C5 ERBB2-high and from the luminal clusters (C1-C3 pooled)?
# Elastic net logistic regression per contrast, ESTIMATE purity unpenalized (genes selected beyond infiltration);
# selection frequency across CV training sets, then a compact logistic model on the top genes. All from the same
# repeated stratified CV.
# Run from repo root: Rscript analyses/04.0_c5_low_enet.R
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(matrixStats); library(glmnet); library(pROC); library(logistf)
  library(parallel); library(tidyverse)
})
out <- "results/04.0_c5_low_enet"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
auc_ <- \(y, x) as.numeric(auc(y, x, direction = "<", quiet = TRUE))

# ---- Data: uncorrected VST, primary tumors; top 5000 variable genes (superset of the 2000 used for PAM) ----
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
vsd <- vsd[, vsd$sample_type == "Primary Tumor"]
vall <- assay(vsd); colnames(vall) <- vsd$sample_id
vmat <- vall[order(-rowVars(vall))[1:5000], ]
lab <- read_csv("results/02.0_her2_labels/her2_labels.csv", show_col_types = FALSE) |>
  left_join(read_csv("results/03.0_pam_clusters/clusters.csv", show_col_types = FALSE), by = "sample_id") |>
  mutate(group = if_else(cluster == "C5", paste("C5", c5_mode), cluster),
         Purity = estimate_score / 1000,                    # ESTIMATE per 1000; higher = lower purity
         Lobular = as.integer(Histology == "Breast Invasive Lobular Carcinoma"),
         LateEra = as.integer(Era == "Plate A213 and after"))
stopifnot(identical(lab$sample_id, colnames(vmat)))
lab <- lab |> filter(!is.na(Histology))   # 1 tumor; dropped so main and sensitivity models use the same tumors
# Current HGNC symbols for display (TCGA uses older ones, e.g. C17orf37 = MIEN1)
sym <- \(g) coalesce(limma::alias2SymbolTable(g, "Hs"), g)

# ---- Contrasts: C5 ERBB2-low (y = 1) vs. C5 ERBB2-high, and vs. luminal clusters C1-C3 pooled (y = 0) ----
# ERBB2 is dropped from the C5-high contrast: the two groups were defined by ERBB2 mRNA.
comps <- c(`C5 ERBB2-high` = "C5 ERBB2-high", `C1-C3` = "C1|C2|C3")
alpha <- 0.5                       # elastic net mixing: keeps correlated genes together, still sparse
k <- 5; reps <- 20                 # outer CV: 5-fold x 20 repeats = 100 training sets (80% of tumors)
m_max <- 20                        # compact models: top 1..20 genes
cores <- 7

xy <- function(comp, covars = "Purity") {
  d <- lab |> filter(group == "C5 ERBB2-low" | str_detect(group, paste0("^(", comps[comp], ")$")))
  g <- setdiff(rownames(vmat), if (comp == "C5 ERBB2-high") "ERBB2")
  x <- cbind(as.matrix(d[covars]), t(vmat[g, d$sample_id]))
  list(x = x, y = as.integer(d$group == "C5 ERBB2-low"), ids = d$sample_id, erbb2 = vall["ERBB2", d$sample_id],
       pf = c(rep(0, length(covars)), rep(1, length(g))))
}
make_folds <- \(y, k) ave(seq_along(y), y, FUN = \(i) sample(rep_len(1:k, length(i))))   # stratified
# Inner 10-fold CV (fold ids passed in) picks lambda by the 1-SE rule
fit_enet <- \(x, y, pf, foldid) cv.glmnet(x, y, family = "binomial", alpha = alpha, penalty.factor = pf, foldid = foldid)
# Selected genes (covariates excluded) with standardized coefficient (beta x SD of the gene in the fitting data)
selected <- \(fit, x) { b <- coef(fit, s = "lambda.1se")[-1, 1]; b <- b[b != 0 & !names(b) %in% c("Purity", "LateEra", "Lobular")]
  b * colSds(x[, names(b), drop = FALSE]) }
glm_auc <- \(f, tr, te) auc_(te$y, predict(glm(f, binomial, tr), te))

# All random folds are drawn here, in the parent process, so results do not depend on the number of cores
set.seed(1)
plan <- map(set_names(names(comps)), \(comp) {
  y <- xy(comp)$y
  jobs <- expand_grid(r = 1:reps, j = 1:k)
  outer <- map(1:reps, \(r) make_folds(y, k))
  jobs$test <- map2(jobs$r, jobs$j, \(r, j) outer[[r]] == j)
  jobs$inner <- map(jobs$test, \(t) make_folds(y[!t], 10))
  list(jobs = jobs, full_inner = make_folds(y, 10))
})

# ---- 1. Repeated CV: per training set, record selected genes, gene ranking and test-fold AUC of each model ----
# Compact model (m genes): in each training set, the m genes with the largest |standardized coefficient| from that
# set's elastic net, refit as plain logistic regression with purity. CV evaluates this procedure, not a fixed gene set.
cv_contrast <- function(comp) {
  d <- xy(comp); jobs <- plan[[comp]]$jobs
  res <- mclapply(seq_len(nrow(jobs)), \(i) {
    test <- jobs$test[[i]]
    fit <- fit_enet(d$x[!test, ], d$y[!test], d$pf, jobs$inner[[i]])
    b <- selected(fit, d$x[!test, ])
    df <- data.frame(y = d$y, ERBB2_mRNA = d$erbb2, d$x[, c("Purity", names(b))], check.names = FALSE)
    top_g <- names(sort(-abs(b)))
    # Near-separation warnings (C5-high contrast) do not affect AUC, which uses only the ranking of predictions
    compact <- map_dbl(seq_len(min(m_max, length(b))), \(m)
      suppressWarnings(glm_auc(reformulate(c("Purity", sprintf("`%s`", top_g[1:m])), "y"), df[!test, ], df[test, ])))
    list(sel = names(b), top = top_g[seq_len(min(m_max, length(top_g)))],
         auc = c(`Purity only` = glm_auc(y ~ Purity, df[!test, ], df[test, ]),
                 `ERBB2 mRNA + purity` = if (comp == "C1-C3") glm_auc(y ~ ERBB2_mRNA + Purity, df[!test, ], df[test, ]) else NA,
                 set_names(compact, paste0("Compact (top ", seq_along(compact), ") + purity")),
                 `Elastic net (all genes) + purity` = auc_(d$y[test], predict(fit, d$x[test, ], s = "lambda.1se")[, 1])))
  }, mc.cores = cores)
  list(sel = map(res, "sel"), top = map(res, "top"),
       auc = map_dfr(res, \(r) enframe(r$auc, "model", "auc"), .id = "fold") |> mutate(comparator = comp, .before = 1))
}
cv <- map(set_names(names(comps)), cv_contrast)

# Corrected repeated k-fold t-test (Nadeau & Bengio 2003): mean, 95% CI and p (vs 0) over the 100 fold values
corrected <- \(x, k) {
  se <- sqrt((1 / length(x) + 1 / (k - 1)) * var(x)); df <- length(x) - 1
  tibble(est = mean(x), lo = mean(x) - qt(0.975, df) * se, hi = mean(x) + qt(0.975, df) * se, p = 2 * pt(-abs(mean(x) / se), df))
}
auc_long <- map_dfr(cv, "auc") |> filter(!is.na(auc))
auc <- auc_long |> group_by(comparator, model) |> summarise(corrected(auc, k), .groups = "drop") |>
  select(comparator, model, auc = est, lo, hi)
write_csv(auc, file.path(out, "cv_auc.csv")); print(auc, n = Inf)

# Compact size (parsimony choice): smallest m within 0.01 mean AUC of the best compact model (m = 1..20)
m_sel <- auc |> filter(str_detect(model, "Compact")) |> group_by(comparator) |>
  summarise(m = min(parse_number(model[auc >= max(auc) - 0.01])), best_m = parse_number(model[which.max(auc)]))
write_csv(m_sel, file.path(out, "compact_size.csv")); print(m_sel)
# % of training sets in which each gene was in the in-fold top m, and ranked first (how fixed the compact gene set is)
top_freq <- map_dfr(names(comps), \(comp) { m <- m_sel$m[m_sel$comparator == comp]; n_sets <- length(cv[[comp]]$top)
  tibble(comparator = comp, gene = unlist(map(cv[[comp]]$top, \(t) head(t, m)))) |> count(comparator, gene) |>
    mutate(pct_in_top_m = 100 * n / n_sets,
           rank1_pct = 100 * map_int(gene, \(g) sum(map_chr(cv[[comp]]$top, 1) == g)) / n_sets) |> select(-n) })

# Paired comparisons of the selected compact model (corrected t-test on per-fold AUC differences, BH)
cmp <- map_dfr(names(comps), \(comp) {
  w <- auc_long |> filter(comparator == comp) |> pivot_wider(names_from = model, values_from = auc)
  cm <- paste0("Compact (top ", m_sel$m[m_sel$comparator == comp], ") + purity")
  others <- intersect(c("Purity only", "ERBB2 mRNA + purity", "Elastic net (all genes) + purity"), names(w))
  map_dfr(set_names(others), \(o) corrected(w[[cm]] - w[[o]], k), .id = "vs") |> mutate(comparator = comp, model = cm, .before = 1)
}) |> rename(delta_auc = est) |> mutate(padj = p.adjust(p, "BH"))
write_csv(cmp, file.path(out, "cv_compare.csv")); print(cmp)

# ---- 2. Full-data elastic net (one fit per contrast, fixed inner folds) + selection frequency ----
# Selection frequency: share of the 100 CV training sets in which the gene had a non-zero coefficient (descriptive;
# overlapping training sets, no error control). diff = mean VST (C5-low) - mean VST (comparator)
full <- map(set_names(names(comps)), \(comp) { d <- xy(comp)
  selected(fit_enet(d$x, d$y, d$pf, plan[[comp]]$full_inner), d$x) })
features <- map_dfr(names(comps), \(comp) {
  d <- xy(comp); b <- full[[comp]]
  freq <- table(unlist(cv[[comp]]$sel)) / length(cv[[comp]]$sel); g <- names(freq)
  tibble(comparator = comp, gene = g, symbol = sym(g), sel_freq = as.numeric(freq[g]), coef_std = replace_na(unname(b[g]), 0),
         mean_c5low = colMeans(d$x[d$y == 1, g, drop = FALSE]), mean_comp = colMeans(d$x[d$y == 0, g, drop = FALSE]),
         diff = mean_c5low - mean_comp)
}) |> left_join(top_freq, by = c("comparator", "gene")) |> mutate(across(c(pct_in_top_m, rank1_pct), \(x) replace_na(x, 0))) |>
  arrange(comparator, -sel_freq, -abs(coef_std))
write_csv(features, file.path(out, "features.csv"))
stable <- features |> filter(sel_freq >= 0.7)    # "stable" = selected in >= 70% of training sets
print(count(stable, comparator))

# ---- 3. Final compact models (full data): the m genes most often in the in-fold top m (what CV evaluated) ----
# Firth's penalized logistic regression (handles the near-separation of the C5-high contrast); profile-likelihood
# 95% CI and p. Genes z-scored within the contrast, so each odds ratio is per 1 SD; purity per 1000 ESTIMATE units.
compact <- map_dfr(names(comps), \(comp) {
  d <- xy(comp); m <- m_sel$m[m_sel$comparator == comp]
  g <- features |> filter(comparator == comp) |> slice_max(pct_in_top_m, n = m, with_ties = FALSE) |> pull(gene)
  df <- data.frame(y = d$y, Purity = d$x[, "Purity"], scale(d$x[, g, drop = FALSE]), check.names = FALSE)
  fit <- logistf(reformulate(c("Purity", sprintf("`%s`", g)), "y"), df)
  tibble(comparator = comp, term = str_remove_all(names(coef(fit)), "`"), beta = coef(fit),
         lo = fit$ci.lower, hi = fit$ci.upper, p = fit$prob, apparent_auc = auc_(d$y, predict(fit))) |>
    mutate(symbol = if_else(term %in% c("(Intercept)", "Purity"), term, sym(term)), .after = term)
})
write_csv(compact, file.path(out, "compact_model.csv")); print(compact, n = Inf)

# ---- 4. Sensitivity (full-data fits, same inner folds): drop purity; add era + lobular histology (unpenalized) ----
# Overlap of selected genes with the main full-data model (Jaccard)
sens <- map_dfr(names(comps), \(comp) {
  main <- names(full[[comp]])
  map_dfr(list(`No purity covariate` = character(0), `Purity + era + lobular` = c("Purity", "LateEra", "Lobular")),
          \(cv_) { d <- xy(comp, covars = cv_); s <- names(selected(fit_enet(d$x, d$y, d$pf, plan[[comp]]$full_inner), d$x))
            tibble(n_main = length(main), n_sens = length(s), n_shared = length(intersect(main, s)),
                   jaccard = length(intersect(main, s)) / length(union(main, s))) }, .id = "model") |>
    mutate(comparator = comp, .before = 1)
})
write_csv(sens, file.path(out, "sensitivity.csv")); print(sens)

# ---- 5. Infiltration check: do stable genes track stroma / immune content? ----
# (a) Spearman rho within C1-C3 only (pooling both groups would make any C5-low marker track purity via group
#     membership) of each stable gene with ESTIMATE and with stroma / immune marker scores (mean z across all tumors
#     of COL1A1, FAP, COL3A1 / PTPRC, CD8A, CD3E; as in Step 3)
# (b) Purity-matched difference: each C5-low tumor matched to its 3 nearest C1-C3 tumors by ESTIMATE (with
#     replacement); diff_matched = mean VST (C5-low) - mean VST (matched C1-C3)
mz <- t(scale(t(vall[c("COL1A1", "FAP", "COL3A1", "PTPRC", "CD8A", "CD3E"), ])))
infil <- tibble(sample_id = colnames(vall), Stroma = colMeans(mz[1:3, ]), Immune = colMeans(mz[4:6, ]))
dl <- xy("C1-C3"); ii <- match(dl$ids, infil$sample_id); lu <- dl$y == 0
pur <- dl$x[, "Purity"]; ref <- which(lu)
matched <- unlist(map(which(!lu), \(i) ref[order(abs(pur[ref] - pur[i]))[1:3]]))
lum <- stable |> filter(comparator == "C1-C3")
infil_cor <- lum |> select(gene, symbol, diff) |>
  mutate(rho_ESTIMATE = map_dbl(gene, \(g) cor(dl$x[lu, g], pur[lu], method = "spearman")),
         rho_stroma = map_dbl(gene, \(g) cor(dl$x[lu, g], infil$Stroma[ii][lu], method = "spearman")),
         rho_immune = map_dbl(gene, \(g) cor(dl$x[lu, g], infil$Immune[ii][lu], method = "spearman")),
         diff_matched = map_dbl(gene, \(g) mean(dl$x[!lu, g]) - mean(dl$x[matched, g])))
write_csv(tibble(group = c("C5 ERBB2-low", "Matched C1-C3", "All C1-C3"),
                 median_ESTIMATE = c(median(pur[!lu]), median(pur[matched]), median(pur[lu])) * 1000),
          file.path(out, "purity_matching.csv"))
write_csv(infil_cor, file.path(out, "infiltration_corr.csv"))
print(infil_cor |> filter(pmax(abs(rho_ESTIMATE), abs(rho_stroma), abs(rho_immune)) >= 0.4), n = Inf)

# ---- 6. Figures ----
# (a) Top 20 genes per contrast by selection frequency (ties by |coefficient|); colour = direction
top <- features |> group_by(comparator) |> slice_head(n = 20) |> ungroup() |>
  mutate(comparator = factor(paste("C5 ERBB2-low vs", comparator), paste("C5 ERBB2-low vs", names(comps))),
         gene_lab = fct_reorder(paste(symbol, comparator, sep = "___"), sel_freq + abs(coef_std) / 100))  # unique per panel
p_top <- ggplot(top, aes(sel_freq, gene_lab, fill = diff > 0)) + geom_col() +
  geom_vline(xintercept = 0.7, linetype = 2, color = "grey50") +
  scale_y_discrete(labels = \(x) sub("___.*", "", x)) + facet_wrap(~comparator, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), labels = c("Lower in C5-low", "Higher in C5-low"),
                    name = NULL) +
  labs(x = "Selection frequency (100 CV training sets)", y = NULL) + theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "top_features.png"), p_top, width = 8, height = 5)

# (b) CV AUC vs compact model size, with purity-only, ERBB2 + purity and full elastic net as reference lines
ref <- auc |> filter(!str_detect(model, "Compact"))
p_auc <- auc |> filter(str_detect(model, "Compact")) |> mutate(m = parse_number(model)) |>
  ggplot(aes(m, auc)) + geom_hline(data = ref, aes(yintercept = auc, linetype = model), color = "grey40") +
  geom_pointrange(aes(ymin = lo, ymax = pmin(hi, 1)), size = 0.3) +
  geom_vline(data = m_sel, aes(xintercept = m), color = "#b2182b", linetype = 3) +
  facet_wrap(~paste("C5 ERBB2-low vs", comparator)) + scale_x_continuous(breaks = seq(2, m_max, 2)) +
  labs(x = "Compact model: number of top genes (+ purity)", y = "CV AUC (mean of test folds; corrected 95% CI)",
       linetype = NULL) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "compact_auc.png"), p_auc, width = 8, height = 4)

# (c) Group-mean z-score (across all tumors) of stable genes; all six groups shown
grp_lv <- c("C5 ERBB2-low", "C5 ERBB2-high", "C1", "C2", "C3", "C4")
hm_genes <- unique(stable$gene)
z <- t(scale(t(vmat[hm_genes, ])))
z_wide <- sapply(grp_lv, \(g) rowMeans(z[, lab$sample_id[lab$group == g], drop = FALSE]))
# Which groups share the C5-low vs C1-C3 profile? Median group z of its up / down genes, and Pearson r of each group's
# profile (mean z over these genes) with the C5-low profile. Genes were chosen to separate C5-low from C1-C3, so any
# non-luminal group correlates by construction; compare groups with each other, not with 0.
zl <- z_wide[lum$gene, ]
bind_rows(as_tibble(apply(zl[lum$diff > 0, ], 2, median)) |> mutate(row = "Median z, genes higher in C5-low"),
          as_tibble(apply(zl[lum$diff < 0, ], 2, median)) |> mutate(row = "Median z, genes lower in C5-low"),
          as_tibble(cor(zl)[, "C5 ERBB2-low"]) |> mutate(row = "r with C5 ERBB2-low profile")) |>
  mutate(group = rep(grp_lv, 3)) |> pivot_wider(names_from = group, values_from = value) |>
  write_csv(file.path(out, "luminal_profile_by_group.csv"))
blocks <- stable |> group_by(gene) |>
  summarise(contrast = if (n() > 1) "Both" else paste("vs", comparator))
hm <- as_tibble(z_wide, rownames = "gene") |> pivot_longer(-gene, names_to = "group", values_to = "z") |>
  left_join(blocks, by = "gene") |>
  mutate(group = factor(group, grp_lv), gene = factor(gene, rownames(z_wide)[hclust(dist(z_wide))$order]),
         contrast = factor(contrast, c("vs C5 ERBB2-high", "Both", "vs C1-C3")))
p_hm <- ggplot(hm, aes(group, gene, fill = pmax(pmin(z, 1.5), -1.5))) + geom_tile() +
  facet_grid(contrast ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = sym) +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", name = "Mean z") +
  labs(x = NULL, y = NULL) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6),
        strip.text.y = element_text(angle = 0))
ggsave(file.path(out, "stable_genes_heatmap.png"), p_hm, width = 5.5, height = max(4, 0.11 * length(hm_genes) + 1))
