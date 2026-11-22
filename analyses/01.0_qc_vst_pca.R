# Step 1: QC, DESeq2 VST normalization, PCA, and clinical associations of PC1-5
# Run from repo root: Rscript analyses/01.0_qc_vst_pca.R
suppressPackageStartupMessages({library(tidyverse); library(DESeq2)})
set.seed(1)
out <- "results/01.0_qc_vst_pca"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
rsem <- read_csv("../data_dir/tcga.brca.rsem.csv", show_col_types = FALSE)
meta <- rsem |> select(1:3) |>
  mutate(sample_id = substr(bcr_patient_barcode, 1, 15),
         plate = substr(bcr_patient_barcode, 22, 25),        # sequencing plate (batch)
         tss = substr(bcr_patient_barcode, 6, 7))            # tissue source site
counts <- t(round(as.matrix(rsem[, -(1:3)])))                # RSEM expected counts -> integers, genes x samples
colnames(counts) <- meta$bcr_patient_barcode

# ---- QC ----
# Input is TCGA upper-quartile (UQ) normalized RSEM: UQ of non-zero values = 1000 in every sample.
# Read depth is therefore already scaled out; DESeq2 size factors (~1) then correct residual composition.
meta <- meta |> mutate(lib_size = colSums(counts), genes_detected = colSums(counts > 0),
                       upper_quartile = apply(counts, 2, \(x) quantile(x[x > 0], 0.75)))
keep <- rowSums(counts >= 10) >= 0.1 * ncol(counts)         # expressed in >=10% of samples
counts <- counts[keep, ]
cat("Genes kept:", sum(keep), "of", length(keep), "\n")

# ---- Normalize: DESeq2 size factors (read depth) + VST ----
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~1)
dds <- estimateSizeFactors(dds)
meta$size_factor <- sizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(out, "vsd.rds"))

qc <- meta |> select(sample_type, lib_size, upper_quartile, genes_detected, size_factor) |>
  mutate(lib_size = lib_size / 1e6) |>
  pivot_longer(-sample_type) |>
  mutate(name = recode_values(name, "lib_size" ~ "Total (M, UQ-scaled)", "upper_quartile" ~ "Upper quartile",
                              "genes_detected" ~ "Genes detected", "size_factor" ~ "DESeq2 size factor") |>
           fct_inorder())
p <- ggplot(qc, aes(sample_type, value, fill = sample_type)) + geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~name, scales = "free_y", nrow = 1) + labs(x = NULL, y = NULL) + theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "qc_boxplots.png"), p, width = 11, height = 3.5)
write_csv(meta |> select(-patient_id), file.path(out, "sample_qc.csv"))
cat("Cor(size factor, library size):", round(cor(meta$size_factor, meta$lib_size), 2), "\n")

# ---- PCA helper: top 2000 variable genes ----
run_pca <- function(m) {
  m <- m[order(-apply(m, 1, var))[1:2000], ]
  pc <- prcomp(t(m))
  list(x = pc$x[, 1:5], pve = round(100 * pc$sdev[1:5]^2 / sum(pc$sdev^2), 1))
}
vmat <- assay(vsd)

# All samples (QC view)
pa <- run_pca(vmat)
p <- as_tibble(pa$x) |> bind_cols(meta) |>
  ggplot(aes(PC1, PC2, color = sample_type)) + geom_point(size = 0.8, alpha = 0.7) +
  labs(x = paste0("PC1 (", pa$pve[1], "%)"), y = paste0("PC2 (", pa$pve[2], "%)"), color = NULL) + theme_bw()
ggsave(file.path(out, "pca_all_samples.png"), p, width = 5.5, height = 4)

# Primary tumors (clinical associations)
pt <- meta$sample_type == "Primary Tumor"
pp <- run_pca(vmat[, pt])
pcs <- as_tibble(pp$x) |> mutate(sample_id = meta$sample_id[pt], plate = meta$plate[pt],
                                 tss = meta$tss[pt], lib_size = meta$lib_size[pt])
write_csv(tibble(PC = paste0("PC", 1:5), pct_var = pp$pve), file.path(out, "pca_variance.csv"))

# ---- Clinical features ----
na_str <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", "Indeterminate", "[Discrepancy]")
cl <- read_csv("../data_dir/brca_tcga_clinical_data.csv", show_col_types = FALSE, na = c("", "NA", na_str))
cn <- read_csv("../data_dir/brca_tcga_erbb2_copy_number.csv", show_col_types = FALSE)
feat <- cl |>
  transmute(sample_id = `Sample ID`,
            Age = `Diagnosis Age`,
            Stage = str_extract(`Neoplasm Disease Stage American Joint Committee on Cancer Code`, "STAGE [IV]+"),
            `T stage` = str_extract(`American Joint Committee on Cancer Tumor Stage Code`, "T[0-4]"),
            `N stage` = str_extract(`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`, "N[0-3]"),
            Histology = `Cancer Type Detailed`,
            ER = `ER Status By IHC`, PR = `PR status by ihc`,
            `HER2 IHC` = `IHC-HER2`, `HER2 FISH` = `HER2 fish status`,
            Menopause = str_remove(`Menopause Status`, " \\(.*"),
            Race = `Race Category`,
            `Fraction genome altered` = `Fraction Genome Altered`,
            `Mutation count (log)` = log1p(`Mutation Count`),
            `OS status` = `Overall Survival Status`,
            `Sample weight` = `Sample Initial Weight`) |>
  left_join(cn |> select(sample_id, `ERBB2 copy number` = erbb2_copy_number), by = "sample_id")

dat <- pcs |> left_join(feat, by = "sample_id") |>
  mutate(`Library size` = lib_size, Plate = plate, `Tissue source site` = tss) |>
  select(-lib_size, -plate, -tss) |>
  mutate(across(where(is.character) & !sample_id, \(x) {                 # drop rare levels (<10)
    x[x %in% names(which(table(x) < 10))] <- NA; x }))

# ---- Association: R^2 of PC ~ feature (lm), BH-adjusted ANOVA p ----
features <- setdiff(names(dat), c(paste0("PC", 1:5), "sample_id"))
features <- features[map_int(features, \(f) n_distinct(na.omit(dat[[f]]))) > 1]   # need >=2 levels
assoc <- expand_grid(PC = paste0("PC", 1:5), feature = features) |>
  mutate(res = map2(PC, feature, \(pc, f) {
    d <- tibble(y = dat[[pc]], x = dat[[f]]) |> drop_na()
    fit <- lm(y ~ x, d)
    tibble(n = nrow(d), r2 = summary(fit)$r.squared, p = anova(fit)$`Pr(>F)`[1])
  })) |> unnest(res) |>
  mutate(padj = p.adjust(p, "BH"),
         type = if_else(feature %in% c("Plate", "Tissue source site", "Library size", "Sample weight"),
                        "Technical", "Clinical"))
write_csv(assoc, file.path(out, "pc_clinical_assoc.csv"))

pve_lab <- setNames(paste0(paste0("PC", 1:5), "\n(", pp$pve, "%)"), paste0("PC", 1:5))
p <- assoc |> mutate(feature = fct_reorder(feature, r2, max), sig = if_else(padj < 0.05, "*", "")) |>
  ggplot(aes(PC, feature, fill = r2)) + geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(r2, 2), sig)), size = 2.8) +
  scale_fill_gradient(low = "white", high = "firebrick", name = expression(R^2)) +
  scale_x_discrete(labels = pve_lab) + facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) + theme_minimal()
ggsave(file.path(out, "pc_clinical_heatmap.png"), p, width = 6.5, height = 6.5)

# PC1 vs PC2 colored by the strongest clinical feature of PC1/PC2
top <- assoc |> filter(type == "Clinical", PC %in% c("PC1", "PC2")) |> slice_max(r2, n = 1) |> pull(feature)
p <- ggplot(dat, aes(PC1, PC2, color = .data[[top]])) + geom_point(size = 0.8, alpha = 0.7) +
  labs(x = paste0("PC1 (", pp$pve[1], "%)"), y = paste0("PC2 (", pp$pve[2], "%)")) + theme_bw()
ggsave(file.path(out, "pca_tumor_top_feature.png"), p, width = 5.5, height = 4)

print(assoc |> group_by(PC) |> slice_max(r2, n = 4) |> select(PC, feature, n, r2, padj), n = 30)
