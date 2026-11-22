# Step 1.1: Test removing the plate effect from primary-tumor VST, and check ER/PR/HER2 biology is kept
# Run from repo root after 01.0: Rscript analyses/01.1_plate_correction.R
suppressPackageStartupMessages({library(tidyverse); library(SummarizedExperiment); library(limma); library(matrixStats)})
r0 <- "results/01.0_qc_vst_pca"
out <- "results/01.1_plate_correction"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

vsd <- readRDS(file.path(r0, "vsd.rds"))
vsd <- vsd[, vsd$sample_type == "Primary Tumor"]
ann <- tibble(sample_id = vsd$sample_id) |>
  left_join(read_csv(file.path(r0, "tumor_pcs_features.csv"), show_col_types = FALSE) |>
              select(sample_id, Plate, ER, PR, `HER2 IHC`, PC3, PC4), by = "sample_id") |>
  mutate(across(c(ER, PR, `HER2 IHC`), \(x) replace_na(x, "Unknown")))
vmat <- assay(vsd)

# ---- Plate sizes: pool plates with <10 tumors into the nearest large plate by processing order ----
# Plate IDs are sequential base-36 codes, so numeric distance approximates processing time.
pid <- \(p) strtoi(substr(p, 2, 4), 36)
plates <- ann |> dplyr::count(Plate) |> mutate(order = pid(Plate))
big <- plates |> filter(n >= 10)
plates <- plates |> mutate(pooled = map_chr(Plate, \(p) if (p %in% big$Plate) p else
  big$Plate[which.min(abs(big$order - pid(p)))]))
ann <- ann |> left_join(plates |> select(Plate, pooled), by = "Plate")
write_csv(plates, file.path(out, "plate_sizes.csv"))
print(filter(plates, Plate != pooled))

# Plate means of PC3/PC4 by processing order: shows the drift and the small plates
p <- ann |> group_by(Plate) |> summarise(n = n(), PC3 = mean(PC3), PC4 = mean(PC4)) |>
  mutate(order = rank(pid(Plate)), small = n < 10) |> pivot_longer(c(PC3, PC4)) |>
  ggplot(aes(order, value, size = n, color = small)) + geom_point(alpha = 0.7) + facet_wrap(~name) +
  scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "firebrick"), labels = c("≥10", "<10"),
                     name = "Tumors on plate") +
  labs(x = "Plate (processing order)", y = "Mean PC score", size = "n") + theme_bw()
ggsave(file.path(out, "plate_drift.png"), p, width = 8, height = 3.5)

# ---- Versions: uncorrected, all 41 plates, pooled plates, pooled with ER/PR/HER2 protected ----
bio_design <- model.matrix(~ ER + PR + `HER2 IHC`, ann)
versions <- list(
  Uncorrected = vmat,
  `Plate (41)` = removeBatchEffect(vmat, batch = ann$Plate),
  `Plate (pooled)` = removeBatchEffect(vmat, batch = ann$pooled),
  `Plate (pooled, ER/PR/HER2 kept)` = removeBatchEffect(vmat, batch = ann$pooled, design = bio_design))
lvl <- names(versions)
saveRDS(versions[["Plate (pooled)"]], file.path(out, "vst_tumor_plate_corrected.rds"))   # step 3 input; unprotected avoids circularity

# Adjusted R^2 of every row of m on one factor (fast, via shared design)
adj_r2 <- function(m, x) {
  ok <- !is.na(x) & x != "Unknown"; m <- m[, ok, drop = FALSE]; X <- model.matrix(~ factor(x[ok]))
  rss <- rowSums((m - t(qr.fitted(qr(X), t(m))))^2); tss <- rowSums((m - rowMeans(m))^2)
  n <- ncol(m); k <- ncol(X) - 1
  1 - (rss / (n - k - 1)) / (tss / (n - 1))
}
feats <- c("Plate", "ER", "PR", "HER2 IHC")

# ---- PCA of each version: PC1-5 vs plate and biology ----
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
ggsave(file.path(out, "pc_assoc_by_version.png"), p, width = 10, height = 4.5)

# ---- Gene level: marker genes vs their clinical label, and vs plate ----
markers <- tribble(~gene, ~feature, "ESR1", "ER", "PGR", "PR", "ERBB2", "HER2 IHC", "GRB7", "HER2 IHC")
marker_r2 <- imap_dfr(versions, \(m, v) markers |>
  mutate(version = v,
         `Clinical label` = map2_dbl(gene, feature, \(g, f) adj_r2(m[g, , drop = FALSE], ann[[f]])),
         Plate = map_dbl(gene, \(g) adj_r2(m[g, , drop = FALSE], ann$Plate))))
write_csv(marker_r2, file.path(out, "marker_r2_by_version.csv"))
p <- marker_r2 |> pivot_longer(c(`Clinical label`, Plate), names_to = "vs") |>
  ggplot(aes(paste0(gene, " ~ ", feature), value, fill = factor(version, rev(lvl)))) +
  geom_col(position = "dodge") + facet_wrap(~vs) + coord_flip() +
  guides(fill = guide_legend(reverse = TRUE, nrow = 2)) +
  labs(x = NULL, y = expression(Adj.~R^2), fill = NULL) + theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "marker_r2.png"), p, width = 7, height = 4)

# ---- Small-plate samples: ERBB2 and ESR1 under each version ----
small <- ann |> mutate(i = row_number()) |> filter(Plate != pooled)
small_tab <- imap_dfr(versions[1:3], \(m, v) small |>
  transmute(sample_id, Plate, `Pooled into` = pooled, ER, `HER2 IHC`, version = v,
            ERBB2 = m["ERBB2", i], ESR1 = m["ESR1", i])) |>
  pivot_wider(names_from = version, values_from = c(ERBB2, ESR1)) |> arrange(Plate)
write_csv(small_tab, file.path(out, "small_plate_samples.csv"))
print(small_tab, width = 160)

# Reference: median ERBB2 / ESR1 (uncorrected) by clinical label
ref <- ann |> mutate(ERBB2 = vmat["ERBB2", ], ESR1 = vmat["ESR1", ])
print(ref |> group_by(`HER2 IHC`) |> summarise(ERBB2 = median(ERBB2)))
print(ref |> group_by(ER) |> summarise(ESR1 = median(ESR1)))
print(marker_r2, width = 120)
