# Step 2.0: clinical HER2 labels (ASCO/CAP 2018), HER2-low categories, label disagreements
# Run from repo root: Rscript analyses/02.0_her2_labels.R
suppressPackageStartupMessages(library(tidyverse))
out <- "results/02.0_her2_labels"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Load: primary tumors with RNA (one per patient) ----
na_str <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", "[Discrepancy]", "Indeterminate")
vsd <- readRDS("results/01.0_qc_vst_pca/vsd.rds")
pt <- vsd$sample_type == "Primary Tumor"
genes <- c("ERBB2", "GRB7", "STARD3", "PGAP3", "ESR1", "PGR")       # ERBB2, 17q12 amplicon, ER/PR
rna <- as_tibble(t(SummarizedExperiment::assay(vsd)[genes, pt])) |> mutate(sample_id = substr(colnames(vsd)[pt], 1, 15))
feat <- read_csv("results/01.0_qc_vst_pca/tumor_pcs_features.csv", show_col_types = FALSE) |>
  select(sample_id, Era, Histology, Age, ER, PR, Race, Plate)
cl <- read_csv("../data_dir/brca_tcga_clinical_data.csv", show_col_types = FALSE) |>
  transmute(sample_id = `Sample ID`, Sex,
            score = `HER2 ihc score`,                     # 0/1/2/3
            category = `IHC-HER2`,                        # Negative/Equivocal/Positive as reported
            fish = `HER2 fish status`) |>
  mutate(across(c(category, fish), \(x) if_else(x %in% na_str, NA, x)))   # FISH Equivocal kept (does not resolve 2+)
cn <- read_csv("../data_dir/brca_tcga_erbb2_copy_number.csv", show_col_types = FALSE) |>
  transmute(sample_id, gistic = erbb2_copy_number)
d <- rna |> left_join(feat, by = "sample_id") |> left_join(cl, by = "sample_id") |> left_join(cn, by = "sample_id")

# ---- Labels ----
# ASCO/CAP 2018 (simplified to available fields): ISH-positive or IHC 3+ -> Positive;
# IHC 0/1+ -> Negative; IHC 2+ resolved by FISH, else Unresolved. Applied separately to the
# numeric score and to the reported category; if both exist and disagree -> Conflict (excluded).
asco <- \(ihc, fish) case_when(fish == "Positive" | ihc == "3+" ~ "Positive",
                               ihc == "0/1+" ~ "Negative",
                               ihc == "2+" & fish == "Negative" ~ "Negative",
                               ihc == "2+" ~ "Unresolved",
                               fish == "Negative" ~ "Negative")              # FISH only
d <- d |> mutate(
  ihc_s = c("0/1+", "0/1+", "2+", "3+")[score + 1],
  ihc_c = recode(category, Negative = "0/1+", Equivocal = "2+", Positive = "3+"),
  status_s = if_else(is.na(ihc_s), NA, asco(ihc_s, fish)),
  status_c = if_else(is.na(ihc_c), NA, asco(ihc_c, fish)),
  status_f = asco(NA, fish),
  her2_status = case_when(!is.na(status_s) & !is.na(status_c) & status_s != status_c ~ "Conflict",
                          TRUE ~ coalesce(status_s, status_c, status_f, "Missing")),
  rule = case_when(her2_status == "Conflict" ~ "Score and category disagree",
                   !is.na(status_s) ~ paste0("IHC score ", ihc_s, if_else(fish %in% c("Positive", "Negative"), paste(", FISH", fish), "")),
                   !is.na(status_c) ~ paste0("IHC category ", category, if_else(fish %in% c("Positive", "Negative"), paste(", FISH", fish), "")),
                   !is.na(status_f) ~ paste("FISH only", fish),
                   TRUE ~ "No IHC or FISH") |> str_trim(),
  # IHC-only binary target for 02.1 (step 1): 3+ vs 0/1+, equivocal and directional conflicts excluded
  ihc_binary = case_when((ihc_s == "3+" & ihc_c == "0/1+") | (ihc_s == "0/1+" & ihc_c == "3+") ~ NA,
                         coalesce(ihc_s, ihc_c) == "3+" ~ "Positive",
                         coalesce(ihc_s, ihc_c) == "0/1+" ~ "Negative"),
  # ADC-relevant HER2-low split of clinically negative tumors (needs a numeric score)
  her2_low = case_when(her2_status == "Positive" ~ "HER2-positive",
                       her2_status != "Negative" ~ NA,
                       score == 0 ~ "HER2-0",
                       score == 1 ~ "HER2-low (1+)",
                       score == 2 ~ "HER2-low (2+/FISH-)",
                       TRUE ~ "Negative, no score"))

# ---- Tables ----
rules <- d |> count(her2_status, rule, name = "n") |> arrange(her2_status, -n)
write_csv(rules, file.path(out, "label_rules.csv"))

# Disagreements; n and how resolved. IHC score vs reported IHC status counts only when the two
# lead to different final statuses (e.g. IHC 2+ / FISH+ reported "Positive" agrees, so not counted)
dis <- bind_rows(
  d |> filter(status_s != status_c) |>
    transmute(sample_id, type = paste0("IHC score ", score, if_else(score > 0, "+", ""), " vs reported IHC status ", category), her2_status),
  d |> filter(coalesce(ihc_s, ihc_c) == "3+", fish == "Negative") |>
    transmute(sample_id, type = "IHC 3+ (score, or reported Positive) vs FISH Negative", her2_status),
  d |> filter(coalesce(ihc_s, ihc_c) == "0/1+", fish == "Positive") |>
    transmute(sample_id, type = "IHC 0/1+ (score, or reported Negative) vs FISH Positive", her2_status),
  d |> filter(her2_status == "Positive", gistic < 2) |>
    transmute(sample_id, type = "Clinical HER2+ vs ERBB2 not amplified (GISTIC < 2)", her2_status),
  d |> filter(her2_status == "Negative", gistic == 2) |>
    transmute(sample_id, type = "Clinical HER2- vs ERBB2 amplified (GISTIC 2)", her2_status))
write_csv(dis, file.path(out, "disagreements.csv"))
dis_sum <- dis |> count(type, her2_status, name = "n") |> pivot_wider(names_from = her2_status, values_from = n, values_fill = 0)
write_csv(dis_sum, file.path(out, "disagreements_summary.csv"))

conc <- d |> count(her2_status, gistic) |> pivot_wider(names_from = gistic, values_from = n, values_fill = 0, names_prefix = "GISTIC ")
write_csv(conc, file.path(out, "status_by_gistic.csv"))
write_csv(d |> select(-ihc_s, -ihc_c, -status_s, -status_c, -status_f), file.path(out, "her2_labels.csv"))

# ---- Figure: ERBB2 mRNA by IHC/FISH group, colored by GISTIC ----
grp <- d |> mutate(group = case_when(her2_status == "Conflict" ~ "Conflict",
                                     score == 2 & fish %in% c("Positive", "Negative") ~ paste0("2+ FISH", if_else(fish == "Positive", "+", "-")),
                                     !is.na(score) ~ paste0(score, "+"),
                                     TRUE ~ NA) |>
                     factor(c("0+", "1+", "2+", "2+ FISH-", "2+ FISH+", "3+", "Conflict")) |>
                     fct_recode(`0` = "0+", `2+ unresolved` = "2+"),
                   ERBB2 = pmax(ERBB2, 10)) |>                  # floor at VST 10 for display (noted in legend)
  filter(!is.na(group))                                          # no IHC score: most still have a reported IHC status (legend)
cat("ERBB2 VST < 10 set to 10 (plotted):", sum(grp$ERBB2 == 10), "\n")
d |> filter(is.na(score)) |> count(`Reported IHC status` = category) |> print()
p <- ggplot(grp, aes(group, ERBB2)) + geom_boxplot(outlier.shape = NA, color = "grey50") +
  geom_jitter(aes(color = factor(gistic)), width = 0.2, height = 0, size = 0.8, alpha = 0.8) +
  geom_text(data = count(grp, group), aes(y = 19.1, label = paste0("n=", n)), size = 3) +
  scale_color_manual(values = c(`-2` = "#2166ac", `-1` = "#92c5de", `0` = "grey70", `1` = "#f4a582", `2` = "#b2182b"),
                     na.value = "black", name = "ERBB2 GISTIC") +
  labs(x = "HER2 IHC score", y = "ERBB2 mRNA (VST)") + theme_bw()
ggsave(file.path(out, "erbb2_by_ihc.png"), p, width = 8, height = 4)

# ---- Summary ----
cat("Primary tumors:", nrow(d), "\n"); print(count(d, her2_status)); print(count(d, her2_low))
print(count(d, ihc_binary)); print(dis_sum); print(conc)