# Step 5.1: C5 ERBB2-low tumors vs. their matched adjacent normal tissue (11 pairs; descriptive, small n).
# Paired limma on VST (TREAT for |log2 FC| > 1) -> Hallmark GSEA (fgsea on the moderated t, checked with camera);
# MCP-counter cell scores, paired; Step 4 compact genes. Not tested here (would need the luminal pairs too): whether the
# tumor-normal change is specific to this group.
# Run from repo root: Rscript analyses/05.1_c5_low_vs_normal.R
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(limma); library(fgsea); library(msigdbr); library(tidyverse)
})
source("analyses/helpers.R")
out <- "results/05.1_c5_low_vs_normal"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Data: C5 ERBB2-low tumors with an adjacent normal from the same patient; uncorrected VST, current symbols ----
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
low <- read_csv("results/03.0_pam_clusters/clusters.csv", show_col_types = FALSE) |> filter(c5_mode == "ERBB2-low")
pairs <- as_tibble(colData(vsd)) |> mutate(col = colnames(vsd)) |>
  filter(sample_type == "Adjacent Normal" | sample_id %in% low$sample_id) |>
  group_by(patient_id) |> filter(n_distinct(sample_type) == 2) |> ungroup() |>
  mutate(tissue = factor(if_else(sample_type == "Adjacent Normal", "Normal", "Tumor"), c("Normal", "Tumor")))
stopifnot(all(table(pairs$patient_id) == 2))
pairs |> group_by(patient_id) |> summarise(same_plate = n_distinct(plate) == 1) |> count(same_plate) |> print()
write_csv(pairs |> select(patient_id, sample_id, tissue, plate, tss), file.path(out, "pairs.csv"))

v <- assay(vsd)[, pairs$col]
rownames(v) <- suppressWarnings(current_symbols(rownames(v)))

# ---- 1. Paired limma (patient as a blocking factor in the design), trend + robust; TREAT for |log2 FC| > 1 ----
X <- model.matrix(~ patient_id + tissue, pairs)
fit <- lmFit(v, X)
tr <- treat(fit, lfc = 1, trend = TRUE, robust = TRUE)
de <- eBayes(fit, trend = TRUE, robust = TRUE) |> topTable(coef = "tissueTumor", number = Inf, sort.by = "none") |>
  as_tibble(rownames = "gene") |> mutate(treat_p = tr$p.value[, "tissueTumor"], treat_padj = p.adjust(treat_p, "BH"))
write_csv(de, file.path(out, "limma.csv"))
cat("Tumor vs normal, TREAT BH p < 0.05:", sum(de$treat_padj < 0.05 & de$logFC > 0), "up,",
    sum(de$treat_padj < 0.05 & de$logFC < 0), "down\n")

# Step 4 compact-model genes: 7 vs. C1-C3 and 1 (STARD3) vs. C5 ERBB2-high; plus ERBB2
key <- read_csv("results/04.0_c5_low_enet/compact_model.csv", show_col_types = FALSE) |>
  filter(!term %in% c("(Intercept)", "Purity")) |>
  transmute(gene = symbol, comparator, key_lab = paste("Step 4 gene,", if_else(beta > 0, "higher", "lower"), "in C5 ERBB2-low vs.", comparator)) |>
  bind_rows(tibble(gene = "ERBB2")) |>     # ERBB2 shown with STARD3 (17q12 amplicon)
  mutate(key_lab = if_else(gene %in% c("ERBB2", "STARD3"), "ERBB2, and STARD3 (Step 4 gene, lower in C5 ERBB2-low vs. C5 ERBB2-high)", key_lab))
key_de <- key |> left_join(de, by = "gene"); print(key_de)
write_csv(key_de, file.path(out, "key_genes.csv"))

# Volcano: significant = TREAT BH p < 0.05; the 8 Step 4 genes and ERBB2 highlighted and labelled
vd <- de |> mutate(dir = case_when(treat_padj >= 0.05 ~ "n.s.", logFC > 0 ~ "Higher in tumor", TRUE ~ "Lower in tumor"))
vlab <- vd |> filter(dir != "n.s.", !gene %in% key$gene) |> group_by(dir) |> slice_max(abs(t), n = 8)
vkey <- vd |> inner_join(key, by = "gene")
vcount <- vd |> filter(dir != "n.s.") |> count(dir) |> mutate(x = if_else(dir == "Higher in tumor", Inf, -Inf), h = if_else(x > 0, 1.1, -0.1))
p_vol <- ggplot(vd, aes(logFC, -log10(treat_padj))) + geom_point(aes(color = dir), size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3) +
  ggrepel::geom_text_repel(data = vlab, aes(label = gene, color = dir), size = 2.5, max.overlaps = Inf, show.legend = FALSE, segment.size = 0.2) +
  geom_point(data = vkey, aes(fill = key_lab), shape = 21, size = 2.6, color = "black") +
  ggrepel::geom_label_repel(data = vkey, aes(label = gene, fill = key_lab), size = 2.8, fontface = "bold", label.size = 0.2,
                            max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, show.legend = FALSE, box.padding = 0.5) +
  geom_text(data = vcount, aes(x = x, y = -Inf, label = paste(n, "genes"), hjust = h, color = dir), vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(`Higher in tumor` = "#d6604d", `Lower in tumor` = "#4393c3", n.s. = "grey80"), name = NULL) +
  scale_fill_manual(values = c(`Step 4 gene, higher in C5 ERBB2-low vs. C1-C3` = "#fdb863",
                               `Step 4 gene, lower in C5 ERBB2-low vs. C1-C3` = "#b2abd2",
                               `ERBB2, and STARD3 (Step 4 gene, lower in C5 ERBB2-low vs. C5 ERBB2-high)` = "#66a61e"), name = NULL) +
  labs(x = "log2 fold change, tumor vs. matched normal (VST)", y = "-log10 BH-adjusted p (TREAT, |log2 FC| > 1)",
       title = sprintf("C5 ERBB2-low tumors vs. adjacent normal (%d pairs)", n_distinct(pairs$patient_id))) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)), fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out, "volcano.png"), p_vol, width = 7, height = 6.5, dpi = 200)

# ---- 2. Hallmark GSEA on the moderated t; camera check (same paired design, inter-gene correlation estimated) ----
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |> (\(m) split(m$gene_symbol, m$gs_name))()
set.seed(1)
gsea <- fgsea(hallmark, setNames(de$t, de$gene), minSize = 15, maxSize = 500, eps = 0) |> as_tibble() |>
  mutate(pathway = str_remove(pathway, "^HALLMARK_"), leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) |> arrange(-NES)
write_csv(gsea, file.path(out, "gsea_hallmark.csv"))
gsea |> filter(padj < 0.05) |> select(pathway, NES, padj) |> print(n = Inf)
cam <- map_dfr(c(Default = 0.01, Estimated = NA), .id = "correlation", \(rho)
  camera(v, ids2indices(hallmark, rownames(v)), X, contrast = "tissueTumor", inter.gene.cor = rho, trend.var = TRUE) |>
    as_tibble(rownames = "pathway") |> select(pathway, any_of("Correlation"), Direction, PValue, FDR)) |> mutate(pathway = str_remove(pathway, "^HALLMARK_"))
write_csv(cam, file.path(out, "camera_hallmark.csv"))
cat("Hallmark sets significant: fgsea", sum(gsea$padj < 0.05), "; camera default", sum(cam$FDR[cam$correlation == "Default"] < 0.05),
    "; camera estimated", sum(cam$FDR[cam$correlation == "Estimated"] < 0.05), "\n")
p_gsea <- gsea |> filter(padj < 0.05) |> ggplot(aes(NES, reorder(pathway, NES), fill = NES > 0)) + geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#d6604d", `FALSE` = "#4393c3"), guide = "none") +
  labs(x = "NES (+ = higher in tumor)", y = NULL) + theme_bw(base_size = 10)
ggsave(file.path(out, "hallmark_nes.png"), p_gsea, width = 6, height = 0.16 * sum(gsea$padj < 0.05) + 1, dpi = 200)

# ---- 3. MCP-counter cell scores (mean VST of marker genes), tumor - normal per patient; paired t-test with 95% CI, BH ----
mcp <- mcp_markers() |> filter(gene %in% rownames(v)) |>
  group_by(population) |> group_modify(\(g, k) tibble(col = colnames(v), score = colMeans(v[g$gene, , drop = FALSE]))) |>
  ungroup() |> left_join(pairs |> select(col, patient_id, tissue), by = "col")
mcp_test <- mcp |> select(-col) |> pivot_wider(names_from = tissue, values_from = score) |> group_by(population) |>
  summarise(median_normal = median(Normal), median_tumor = median(Tumor),
            tt = list(t.test(Tumor, Normal, paired = TRUE))) |>
  mutate(mean_diff = map_dbl(tt, "estimate"), lo = map_dbl(tt, \(t) t$conf.int[1]), hi = map_dbl(tt, \(t) t$conf.int[2]),
         p = map_dbl(tt, "p.value"), padj = p.adjust(p, "BH")) |> select(-tt) |> arrange(p)
write_csv(mcp_test, file.path(out, "mcp_counter_paired.csv")); print(mcp_test)

# Are the 11 paired tumors typical of the group? Purity and PAM50, paired vs. the other 72
rep_chk <- read_csv("results/02.0_her2_labels/her2_labels.csv", show_col_types = FALSE) |> select(sample_id, estimate_score) |>
  inner_join(low, by = "sample_id") |> mutate(paired = sample_id %in% pairs$sample_id)
rep_out <- rep_chk |> group_by(paired) |>
  summarise(n = n(), estimate_median = median(estimate_score), Her2 = sum(PAM50 == "Her2"), LumB = sum(PAM50 == "LumB"),
            LumA = sum(PAM50 == "LumA"), other = n - Her2 - LumB - LumA) |>
  mutate(wilcox_p_estimate = wilcox.test(estimate_score ~ paired, rep_chk)$p.value)
write_csv(rep_out, file.path(out, "paired_representativeness.csv")); print(rep_out)
pairs |> filter(tissue == "Tumor") |> count(tss) |> print()
# Matched normals available per tumor group (for the not-done group-specific comparison)
nrm <- as_tibble(colData(vsd)) |> filter(sample_type == "Adjacent Normal") |> pull(patient_id)
pairs_by_group <- read_csv("results/03.0_pam_clusters/clusters.csv", show_col_types = FALSE) |>
  mutate(group = case_when(cluster %in% c("C1", "C2", "C3") ~ "C1-C3", cluster == "C5" ~ paste("C5", c5_mode), TRUE ~ cluster)) |>
  group_by(group) |> summarise(tumors = n(), with_matched_normal = sum(substr(sample_id, 1, 12) %in% nrm))
write_csv(pairs_by_group, file.path(out, "pairs_by_group.csv")); print(pairs_by_group)

writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
