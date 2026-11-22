# Step 2.0: clinical HER2 labels (2018-style rules on reported IHC/FISH calls), HER2-low categories,
# label disagreements, and per-tumor features for 02.1 (marker genes, ESTIMATE purity score, site)
# Run from repo root: Rscript analyses/02.0_her2_labels.R
suppressPackageStartupMessages({library(tidyverse); library(tidyestimate)})
out <- "results/02.0_her2_labels"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Load: primary tumors with RNA (one per patient) ----
na_str <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", "[Discrepancy]", "Indeterminate")
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
pt <- vsd$sample_type == "Primary Tumor"
vmat <- SummarizedExperiment::assay(vsd)[, pt]; colnames(vmat) <- substr(colnames(vmat), 1, 15)
genes <- c("ERBB2", "GRB7", "STARD3", "PGAP3", "ESR1", "PGR")       # ERBB2, 17q12 amplicon, ER/PR
rna <- as_tibble(t(vmat[genes, ]), rownames = "sample_id")
# ESTIMATE score (stromal + immune ssGSEA; rank-based, so VST input is fine). Higher = lower tumor purity.
# The score-to-purity formula is calibrated for Affymetrix arrays, so the score itself is used.
est <- vmat |> filter_common_genes(id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE) |>
  estimate_score(is_affymetrix = FALSE) |> transmute(sample_id = sample, estimate_score = estimate)
feat <- read_csv("results/01.0_qc_vst_pca/tumor_pcs_features.csv", show_col_types = FALSE) |>
  select(sample_id, Era, Histology, Age, ER, PR, Race, Plate, TSS = `Tissue source site`)
cl <- read_csv("../data_dir/brca_tcga_clinical_data.csv", show_col_types = FALSE) |>
  transmute(sample_id = `Sample ID`, Sex,
            score = `HER2 ihc score`,                     # numeric IHC score 0/1/2/3
            category = `IHC-HER2`,                        # reported IHC status: Negative/Equivocal/Positive
            fish = `HER2 fish status`) |>
  mutate(across(c(category, fish), \(x) if_else(x %in% na_str, NA, x)))   # FISH Equivocal kept (does not resolve 2+)
cn <- read_csv("../data_dir/brca_tcga_erbb2_copy_number.csv", show_col_types = FALSE) |>
  transmute(sample_id, gistic = erbb2_copy_number)
d <- rna |> left_join(est, by = "sample_id") |> left_join(feat, by = "sample_id") |>
  left_join(cl, by = "sample_id") |> left_join(cn, by = "sample_id")

# ---- Labels ----
# TCGA IHC/FISH calls were made under 2007/2013 criteria; 2018-style rules are applied to those calls.
# IHC 3+ -> Positive. IHC 2+ -> FISH decides, else Unresolved. IHC 0/1+ -> Negative, but FISH+ is Unresolved
# (2018: positive only for ISH group 1, which cannot be checked here). The reported IHC status may record the
# final HER2 call rather than IHC, so reported Positive with FISH- is Unresolved.
call <- \(ihc, fish) case_when(ihc == "3+" ~ "Positive",
                               ihc == "Reported Positive" ~ if_else(fish %in% "Negative", "Unresolved", "Positive"),
                               ihc == "2+" ~ case_when(fish == "Positive" ~ "Positive", fish == "Negative" ~ "Negative",
                                                       TRUE ~ "Unresolved"),
                               ihc == "0/1+" ~ if_else(fish %in% "Positive", "Unresolved", "Negative"),
                               fish %in% c("Positive", "Negative") ~ fish)          # FISH only
d <- d |> mutate(
  ihc_s = c("0/1+", "0/1+", "2+", "3+")[score + 1],
  ihc_c = recode(category, Negative = "0/1+", Equivocal = "2+", Positive = "Reported Positive"),
  status_s = if_else(is.na(ihc_s), NA, call(ihc_s, fish)),
  status_c = if_else(is.na(ihc_c), NA, call(ihc_c, fish)),
  # Score and reported status both definitive but opposite -> Conflict (excluded); otherwise the score wins
  her2_status = case_when(status_s %in% c("Positive", "Negative") & status_c %in% c("Positive", "Negative") &
                            status_s != status_c ~ "Conflict",
                          TRUE ~ coalesce(status_s, status_c, call(NA, fish), "Missing")),
  label_source = case_when(!is.na(score) ~ "Numeric IHC score", !is.na(category) ~ "Reported IHC status only",
                           !is.na(fish) ~ "FISH only", TRUE ~ "None"),
  fish_txt = if_else(fish %in% c("Positive", "Negative"), paste0(", FISH ", fish), ""),
  rule = case_when(her2_status == "Conflict" ~ "IHC score and reported IHC status disagree",
                   !is.na(score) ~ paste0("IHC score ", ihc_s, fish_txt),
                   !is.na(category) ~ paste0("Reported IHC ", category, fish_txt),
                   !is.na(fish) ~ paste("FISH only:", fish),
                   TRUE ~ "No IHC or FISH"),
  # IHC-only target used as a sensitivity analysis in 02.1: score if present, else reported status (3+ vs 0/1+)
  ihc_binary = case_when(her2_status == "Conflict" | rule == "Reported IHC Positive, FISH Negative" ~ NA,
                         coalesce(ihc_s, ihc_c) %in% c("3+", "Reported Positive") ~ "Positive",
                         coalesce(ihc_s, ihc_c) == "0/1+" ~ "Negative"),
  # ADC-relevant HER2-low split of HER2- tumors (needs a numeric score)
  her2_low = case_when(her2_status == "Positive" ~ "HER2-positive",
                       her2_status != "Negative" ~ NA,
                       score == 0 ~ "HER2-0",
                       score == 1 ~ "HER2-low (1+)",
                       score == 2 ~ "HER2-low (2+/FISH-)",
                       TRUE ~ "Negative, no score"))

# ---- Tables ----
write_csv(d |> count(her2_status, rule, name = "n") |> arrange(her2_status, -n), file.path(out, "label_rules.csv"))

# Disagreements between sources, and the status each tumor received
ihc_lab <- \(s) paste0(s, if_else(s > 0, "+", ""))
dis <- bind_rows(
  d |> filter(status_s != status_c, !(score == 3 & category == "Positive")) |>   # 3+ reported Positive agree on IHC
    transmute(sample_id, type = paste0("IHC score ", ihc_lab(score), " vs reported IHC ", category, fish_txt)),
  d |> filter(score == 3 | (is.na(score) & category == "Positive"), fish == "Negative") |>
    transmute(sample_id, type = "IHC 3+ or reported Positive vs FISH Negative"),
  d |> filter(score <= 1 | (is.na(score) & category == "Negative"), fish == "Positive") |>
    transmute(sample_id, type = "IHC 0/1+ or reported Negative vs FISH Positive"),
  d |> filter(her2_status == "Positive", gistic < 2) |>
    transmute(sample_id, type = paste0("HER2+ (", label_source, ") vs not amplified (GISTIC < 2)")),
  d |> filter(her2_status == "Negative", gistic == 2) |>
    transmute(sample_id, type = paste0("HER2- (", label_source, ") vs amplified (GISTIC 2)"))) |>
  left_join(d |> select(sample_id, her2_status), by = "sample_id")
write_csv(dis, file.path(out, "disagreements.csv"))
dis_sum <- dis |> count(type, her2_status) |>
  pivot_wider(names_from = her2_status, values_from = n, values_fill = 0) |> mutate(Total = rowSums(pick(-type)))
write_csv(dis_sum, file.path(out, "disagreements_summary.csv"))
write_csv(d |> count(her2_status, label_source, gistic) |>
            pivot_wider(names_from = gistic, values_from = n, values_fill = 0, names_prefix = "GISTIC "),
          file.path(out, "status_by_gistic.csv"))
write_csv(d |> select(-ihc_s, -ihc_c, -status_s, -status_c, -fish_txt), file.path(out, "her2_labels.csv"))

# ---- Figure: ERBB2 mRNA by IHC score (FISH splits 2+), colored by GISTIC ----
grp <- d |> filter(!is.na(score)) |>
  mutate(group = case_when(her2_status == "Conflict" ~ "Conflict",
                           score == 2 & fish %in% c("Positive", "Negative") ~ paste0("2+ FISH", if_else(fish == "Positive", "+", "-")),
                           score == 2 ~ "2+ unresolved",
                           TRUE ~ ihc_lab(score)) |>
           factor(c("0", "1+", "2+ unresolved", "2+ FISH-", "2+ FISH+", "3+", "Conflict")),
         ERBB2 = pmax(ERBB2, 10))                                # floor at VST 10 for display (noted in legend)
cat("ERBB2 VST < 10 set to 10 (plotted):", sum(grp$ERBB2 == 10), "\n")
p <- ggplot(grp, aes(group, ERBB2)) + geom_boxplot(outlier.shape = NA, color = "grey50") +
  geom_jitter(aes(color = factor(gistic)), width = 0.2, height = 0, size = 0.8, alpha = 0.8) +
  geom_text(data = count(grp, group), aes(y = 19.1, label = paste0("n=", n)), size = 3) +
  scale_color_manual(values = c(`-2` = "#2166ac", `-1` = "#92c5de", `0` = "grey70", `1` = "#f4a582", `2` = "#b2182b"),
                     na.value = "black", name = "ERBB2 GISTIC") +
  labs(x = "HER2 IHC score", y = "ERBB2 mRNA (VST)") + theme_bw()
ggsave(file.path(out, "erbb2_by_ihc.png"), p, width = 8, height = 4)

# ---- Summary ----
cat("Primary tumors:", nrow(d), "\n"); print(count(d, her2_status)); print(count(d, her2_low))
print(count(d, ihc_binary)); print(dis_sum, n = Inf)
d |> filter(is.na(score)) |> count(category) |> print()
