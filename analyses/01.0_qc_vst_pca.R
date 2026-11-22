# Step 1.0: QC, DESeq2 VST, PCA, and clinical/technical associations of PC1-5
# Run from repo root: Rscript analyses/01.0_qc_vst_pca.R
suppressPackageStartupMessages({library(tidyverse); library(DESeq2); library(matrixStats)})
set.seed(1)
out <- "results/01.0_qc_vst_pca"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
rsem <- read_csv("../data_dir/tcga.brca.rsem.csv", show_col_types = FALSE)
meta <- rsem |> select(1:3) |>
  mutate(sample_id = substr(bcr_patient_barcode, 1, 15),
         plate = substr(bcr_patient_barcode, 22, 25),        # sequencing plate (batch)
         tss = substr(bcr_patient_barcode, 6, 7))            # tissue source site
uq <- t(as.matrix(rsem[, -(1:3)]))                           # genes x samples
colnames(uq) <- meta$bcr_patient_barcode

# ---- QC ----
# Input is TCGA upper-quartile (UQ) scaled RSEM: UQ of non-zero values = 1000 in every sample,
# so read depth is already scaled out upstream. Rounded values are count-like, not true counts.
meta <- meta |> mutate(uq_total = colSums(uq), genes_detected = colSums(uq > 0),
                       upper_quartile = apply(uq, 2, \(x) quantile(x[x > 0], 0.75)))
counts <- round(uq)
keep <- rowSums(counts >= 10) >= 0.1 * ncol(counts)          # expressed in >=10% of samples
counts <- counts[keep, ]
cat("Genes kept:", sum(keep), "of", length(keep), "\n")

# ---- VST: approximate variance stabilizer on UQ-scaled data ----
# Size factors (~1) only correct residual composition; UQ scaling already handled depth.
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~1)
dds <- estimateSizeFactors(dds)
meta$size_factor <- sizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(out, "vsd.rds"))

qc <- meta |> select(sample_type, uq_total, upper_quartile, genes_detected, size_factor) |>
  mutate(uq_total = uq_total / 1e6) |>
  pivot_longer(-sample_type) |>
  mutate(name = recode_values(name, "uq_total" ~ "Total (M, UQ-scaled)", "upper_quartile" ~ "Upper quartile",
                              "genes_detected" ~ "Genes detected", "size_factor" ~ "DESeq2 size factor") |>
           fct_inorder())
p <- ggplot(qc, aes(sample_type, value, fill = sample_type)) + geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~name, scales = "free_y", nrow = 1) + labs(x = NULL, y = NULL) + theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "qc_boxplots.png"), p, width = 11, height = 3.5)
write_csv(meta |> select(-patient_id), file.path(out, "sample_qc.csv"))

# ---- PCA helper: top 2000 variable genes ----
run_pca <- function(m) {
  m <- m[order(-rowVars(m))[1:2000], ]
  pc <- prcomp(t(m))
  list(x = pc$x[, 1:5], pve = round(100 * pc$sdev[1:5]^2 / sum(pc$sdev^2), 1))
}
vmat <- assay(vsd)

# ---- Clinical features ----
na_str <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", "[Discrepancy]", "Indeterminate")
cl <- read_csv("../data_dir/brca_tcga_clinical_data.csv", show_col_types = FALSE)
cn <- read_csv("../data_dir/brca_tcga_erbb2_copy_number.csv", show_col_types = FALSE)
feat <- cl |>
  transmute(sample_id = `Sample ID`,
            Age = `Diagnosis Age`,
            Stage = str_extract(`Neoplasm Disease Stage American Joint Committee on Cancer Code`, "Stage [IV]+"),   # Stage X -> NA
            `T stage` = str_extract(`American Joint Committee on Cancer Tumor Stage Code`, "T[0-4]"),   # TX -> NA
            `N stage` = str_extract(`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`, "N[0-3]"),  # NX -> NA
            Histology = `Cancer Type Detailed`,
            ER = `ER Status By IHC`, PR = `PR status by ihc`,
            `HER2 IHC` = `IHC-HER2`,                        # keeps Equivocal as its own level
            `HER2 FISH` = `HER2 fish status`,               # tested mainly in IHC 2+ cases
            Menopause = str_remove(`Menopause Status`, " \\(.*"),
            Race = `Race Category`,
            `Fraction genome altered` = `Fraction Genome Altered`,
            `Mutation count (log)` = log1p(`Mutation Count`),
            `OS status` = `Overall Survival Status`,
            `Sample weight` = `Sample Initial Weight`) |>
  mutate(across(where(is.character), \(x) if_else(x %in% na_str, NA, x))) |>
  left_join(cn |> transmute(sample_id, `ERBB2 CN (GISTIC)` = as.character(erbb2_copy_number)), by = "sample_id")

# All samples (QC view), colored by sample type and ER status
pa <- run_pca(vmat)
grp <- meta |> left_join(feat |> select(sample_id, ER), by = "sample_id") |>
  mutate(group = case_when(sample_type != "Primary Tumor" ~ sample_type,
                           !is.na(ER) ~ paste("Tumor ER", if_else(ER == "Positive", "+", "-")),
                           TRUE ~ "Tumor ER unknown"))
p <- as_tibble(pa$x) |> bind_cols(grp) |>
  ggplot(aes(PC1, PC2, color = group)) + geom_point(size = 0.8, alpha = 0.7) +
  labs(x = paste0("PC1 (", pa$pve[1], "%)"), y = paste0("PC2 (", pa$pve[2], "%)"), color = NULL) + theme_bw()
ggsave(file.path(out, "pca_all_samples.png"), p, width = 5.5, height = 4)

# Primary tumors (clinical associations)
pt <- meta$sample_type == "Primary Tumor"
pp <- run_pca(vmat[, pt])
lg <- run_pca(log2(uq[keep, pt] + 1))                       # robustness: log2(UQ+1) PCA
robust <- apply(abs(cor(pp$x, lg$x)), 1, max)                # best |r| match per VST PC
write_csv(tibble(PC = paste0("PC", 1:5), pct_var = pp$pve, max_abs_r_log2 = round(robust, 2)),
          file.path(out, "pca_variance.csv"))

dat <- as_tibble(pp$x) |>
  mutate(sample_id = meta$sample_id[pt], Plate = meta$plate[pt], `Tissue source site` = meta$tss[pt],
         `UQ-scaled total` = meta$uq_total[pt]) |>
  left_join(feat, by = "sample_id")
write_csv(dat, file.path(out, "tumor_pcs_features.csv"))      # unlumped; used by 01.1
dat <- dat |> mutate(across(where(is.character) & !sample_id,      # rare levels (<10) -> "Other"
                            \(x) as.character(fct_lump_min(x, 10, other_level = "Other"))))

# ---- Association: adjusted R^2 of PC ~ feature (lm) ----
technical <- c("Plate", "Tissue source site", "UQ-scaled total", "Sample weight")
features <- setdiff(names(dat), c(paste0("PC", 1:5), "sample_id"))
features <- features[map_int(features, \(f) n_distinct(na.omit(dat[[f]]))) > 1]   # need >=2 levels
assoc <- expand_grid(PC = paste0("PC", 1:5), feature = features) |>
  mutate(res = map2(PC, feature, \(pc, f) {
    fit <- lm(y ~ x, tibble(y = dat[[pc]], x = dat[[f]]))
    tibble(n = nobs(fit), adj_r2 = summary(fit)$adj.r.squared, p = anova(fit)$`Pr(>F)`[1])
  })) |> unnest(res) |>
  mutate(padj = p.adjust(p, "BH"), type = if_else(feature %in% technical, "Technical", "Clinical"))
write_csv(assoc, file.path(out, "pc_clinical_assoc.csv"))

pve_lab <- setNames(paste0("PC", 1:5, "\n(", pp$pve, "%)"), paste0("PC", 1:5))
p <- assoc |> mutate(feature = fct_reorder(feature, adj_r2, max), adj_r2 = pmax(adj_r2, 0)) |>
  ggplot(aes(PC, feature, fill = adj_r2)) + geom_tile(color = "white") +
  geom_text(aes(label = round(adj_r2, 2)), size = 2.8) +
  scale_fill_gradient(low = "white", high = "firebrick", name = expression(Adj.~R^2)) +
  scale_x_discrete(labels = pve_lab) + facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) + theme_minimal()
ggsave(file.path(out, "pc_clinical_heatmap.png"), p, width = 6.5, height = 6.5)

# ---- Plate beyond biology: partial R^2 of plate given ER, HER2 IHC, histology, FGA ----
bio <- dat |> select(starts_with("PC"), Plate, ER, `HER2 IHC`, Histology, `Fraction genome altered`) |> drop_na()
plate_partial <- map_dfr(paste0("PC", 1:5), \(pc) {
  f0 <- lm(reformulate(c("ER", "`HER2 IHC`", "Histology", "`Fraction genome altered`"), pc), bio)
  f1 <- update(f0, . ~ . + Plate)
  tibble(PC = pc, n = nrow(bio), partial_r2 = 1 - deviance(f1) / deviance(f0), p = anova(f0, f1)$`Pr(>F)`[2])
})
write_csv(plate_partial, file.path(out, "plate_partial_r2.csv"))
print(plate_partial)
cat("HER2 IHC vs plate, chi-square p:", signif(chisq.test(table(dat$Plate, dat$`HER2 IHC`), simulate.p.value = TRUE)$p.value, 2), "\n")
print(assoc |> group_by(PC) |> slice_max(adj_r2, n = 3) |> select(PC, feature, n, adj_r2), n = 15)
