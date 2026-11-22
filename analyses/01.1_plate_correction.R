# Step 1.1: Plate / processing-era assessment. Tests plate correction on primary tumors, but NO correction
# is carried forward: era is confounded with sample type and histology (see era_balance.csv).
# Run from repo root after 01.0: Rscript analyses/01.1_plate_correction.R
suppressPackageStartupMessages({library(tidyverse); library(SummarizedExperiment); library(limma); library(matrixStats)})
r0 <- "results/01.0_qc_vst_pca"
out <- "results/01.1_plate_correction"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
pid <- \(p) strtoi(substr(p, 2, 4), 36)          # plate IDs: sequential base-36 codes, all start with "A"
eras <- c("Plates before A213", "Plate A213 and after")

qc <- read_csv(file.path(r0, "sample_qc.csv"), show_col_types = FALSE)
tum <- read_csv(file.path(r0, "tumor_pcs_features.csv"), show_col_types = FALSE)

# ---- Plate composition: normal / tumor / metastasis per plate, counts and % of plate ----
comp <- qc |> dplyr::count(Plate = plate, Era = era, sample_type) |>
  pivot_wider(names_from = sample_type, values_from = n, values_fill = 0) |>
  transmute(Plate, Era, Normal = `Adjacent Normal`, Tumor = `Primary Tumor`, Metastasis,
            Total = Normal + Tumor + Metastasis,
            across(c(Normal, Tumor, Metastasis), \(x) round(100 * x / Total, 1), .names = "% {.col}")) |>
  arrange(pid(Plate))
write_csv(comp, file.path(out, "plate_composition.csv"))

# ---- Era balance: sample types and tumor biology by processing era (Fisher test, feature vs rest) ----
row <- \(label, x, era) {
  t <- table(factor(era, eras), x)
  tibble(Feature = label, !!eras[1] := sprintf("%.1f%%", 100 * t[1, "TRUE"] / sum(t[1, ])),
         !!eras[2] := sprintf("%.1f%%", 100 * t[2, "TRUE"] / sum(t[2, ])), p = signif(fisher.test(t)$p.value, 2))
}
counts <- qc |> dplyr::count(era, sample_type) |> pivot_wider(names_from = era, values_from = n, values_fill = 0) |>
  transmute(Feature = paste0(sample_type, " (n)"), !!eras[1] := as.character(.data[[eras[1]]]),
            !!eras[2] := as.character(.data[[eras[2]]]))
nt <- filter(qc, sample_type != "Metastasis")
balance <- bind_rows(
  counts,
  row("Adjacent normal (% of normal + tumor)", nt$sample_type == "Adjacent Normal", nt$era),
  with(filter(tum, !is.na(Histology)), row("Lobular histology", Histology == "Breast Invasive Lobular Carcinoma", Era)),
  with(filter(tum, !is.na(`ERBB2 CN (GISTIC)`)), row("ERBB2 amplified (GISTIC 2)", `ERBB2 CN (GISTIC)` == 2, Era)),
  with(filter(tum, !is.na(ER)), row("ER positive", ER == "Positive", Era)),
  with(filter(tum, !is.na(PR)), row("PR positive", PR == "Positive", Era)),
  with(filter(tum, !is.na(`HER2 IHC`)), row("HER2 IHC positive", `HER2 IHC` == "Positive", Era)))
write_csv(balance, file.path(out, "era_balance.csv"))
print(balance)

# ---- Plate drift: mean PC3/PC4 per plate in processing order ----
p <- tum |> group_by(Plate, Era) |> summarise(n = n(), PC3 = mean(PC3), PC4 = mean(PC4), .groups = "drop") |>
  mutate(order = rank(pid(Plate)), `Tumors on plate` = if_else(n < 10, "<10", "≥10"), Era = factor(Era, eras)) |>
  pivot_longer(c(PC3, PC4)) |>
  ggplot(aes(order, value, size = n, color = Era, shape = `Tumors on plate`)) + geom_point(alpha = 0.7) +
  facet_wrap(~name) + scale_shape_manual(values = c(`<10` = 1, `≥10` = 16)) +
  labs(x = "Plate (processing order)", y = "Mean PC score", color = NULL) + theme_bw()
ggsave(file.path(out, "plate_drift.png"), p, width = 9, height = 3.5)

# ---- Correction test (tumors only): pool plates with <10 tumors into nearest large plate by order ----
vsd <- readRDS(file.path(r0, "vsd.rds"))
vsd <- vsd[, vsd$sample_type == "Primary Tumor"]
ann <- tibble(sample_id = vsd$sample_id) |>
  left_join(select(tum, sample_id, Plate, ER, PR, `HER2 IHC`), by = "sample_id") |>
  mutate(across(c(ER, PR, `HER2 IHC`), \(x) replace_na(x, "Unknown")))
vmat <- assay(vsd)

plates <- ann |> dplyr::count(Plate) |> mutate(order = pid(Plate))
big <- plates |> filter(n >= 10)
plates <- plates |> mutate(pooled = map_chr(Plate, \(p) if (p %in% big$Plate) p else
  big$Plate[which.min(abs(big$order - pid(p)))]))
ann <- ann |> left_join(plates |> select(Plate, pooled), by = "Plate")

bio_design <- model.matrix(~ ER + PR + `HER2 IHC`, ann)
versions <- list(
  Uncorrected = vmat,
  `Plate (41)` = removeBatchEffect(vmat, batch = ann$Plate),
  `Plate (pooled)` = removeBatchEffect(vmat, batch = ann$pooled),
  `Plate (pooled, ER/PR/HER2 kept)` = removeBatchEffect(vmat, batch = ann$pooled, design = bio_design))
lvl <- names(versions)

# Adjusted R^2 of every row of m on one factor (fast, via shared design)
adj_r2 <- function(m, x) {
  ok <- !is.na(x) & x != "Unknown"; m <- m[, ok, drop = FALSE]; X <- model.matrix(~ factor(x[ok]))
  rss <- rowSums((m - t(qr.fitted(qr(X), t(m))))^2); tss <- rowSums((m - rowMeans(m))^2)
  n <- ncol(m); k <- ncol(X) - 1
  1 - (rss / (n - k - 1)) / (tss / (n - 1))
}
feats <- c("ER", "PR", "HER2 IHC")   # plate is ~0 after correction by construction, so not shown

# PCA of each version: PC1-5 vs biology
pc_assoc <- imap_dfr(versions, \(m, v) {
  pc <- prcomp(t(m[order(-rowVars(m))[1:2000], ]))
  pve <- round(100 * pc$sdev[1:5]^2 / sum(pc$sdev^2), 1)
  map_dfr(feats, \(f) tibble(version = v, PC = paste0("PC", 1:5, " (", pve, "%)"),
                             feature = f, adj_r2 = adj_r2(t(pc$x[, 1:5]), ann[[f]])))
})
write_csv(pc_assoc, file.path(out, "pc_assoc_by_version.csv"))
p <- pc_assoc |> mutate(version = factor(version, lvl), feature = fct_rev(factor(feature, feats)),
                        adj_r2 = pmax(adj_r2, 0)) |>
  ggplot(aes(PC, feature, fill = adj_r2)) + geom_tile(color = "white") +
  geom_text(aes(label = round(adj_r2, 2)), size = 2.8) +
  scale_fill_gradient(low = "white", high = "firebrick", name = expression(Adj.~R^2)) +
  facet_wrap(~version, ncol = 2, scales = "free_x") + labs(x = NULL, y = NULL) + theme_minimal()
ggsave(file.path(out, "pc_assoc_by_version.png"), p, width = 10, height = 3.5)

# Marker genes vs their clinical label
markers <- tribble(~gene, ~feature, "ESR1", "ER", "PGR", "PR", "ERBB2", "HER2 IHC", "GRB7", "HER2 IHC")
marker_r2 <- imap_dfr(versions, \(m, v) markers |>
  mutate(version = v, adj_r2 = map2_dbl(gene, feature, \(g, f) adj_r2(m[g, , drop = FALSE], ann[[f]]))))
write_csv(marker_r2, file.path(out, "marker_r2_by_version.csv"))
p <- ggplot(marker_r2, aes(paste0(gene, " ~ ", feature), adj_r2, fill = factor(version, rev(lvl)))) +
  geom_col(position = "dodge") + coord_flip() + guides(fill = guide_legend(reverse = TRUE, nrow = 2)) +
  labs(x = NULL, y = expression(Adj.~R^2), fill = NULL) + theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "marker_r2.png"), p, width = 6, height = 4)

# Small-plate tumors: ERBB2 and ESR1 under each version
small <- ann |> mutate(i = row_number()) |> filter(Plate != pooled)
small_tab <- imap_dfr(versions[c("Uncorrected", "Plate (41)", "Plate (pooled)")], \(m, v) small |>
  transmute(sample_id, Plate, `Pooled into` = pooled, ER, `HER2 IHC`, version = v,
            ERBB2 = m["ERBB2", i], ESR1 = m["ESR1", i])) |>
  pivot_wider(names_from = version, values_from = c(ERBB2, ESR1)) |> arrange(Plate)
write_csv(small_tab, file.path(out, "small_plate_samples.csv"))
print(marker_r2, n = 20)
