# Step 2.1: Is RNA or DNA more predictive of HER2 IHC?
# Clinical exploration, label-noise check, raw RNA vs DNA comparison, purity, cross-validated logistic models
# Run from repo root: Rscript analyses/02.1_rna_vs_dna.R
suppressPackageStartupMessages({library(tidyverse); library(pROC)})
set.seed(1)
out <- "results/02.1_rna_vs_dna"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
auc_ <- \(y, x) as.numeric(auc(y, x, direction = "<", quiet = TRUE))

# ---- Data (from 02.0) ----
d <- read_csv("results/02.0_her2_labels/her2_labels.csv", show_col_types = FALSE) |>
  mutate(y_ihc = if_else(her2_status == "Conflict" | is.na(score) | score == 2, NA, as.integer(score == 3)),  # main: IHC score 3+ vs 0/1+
         y_rep = as.integer(ihc_binary == "Positive"),            # sensitivity: reported IHC status used when no score
         y_clin = case_when(her2_status == "Positive" ~ 1L, her2_status == "Negative" ~ 0L),   # sensitivity: FISH-resolved
         y_low = case_when(her2_low %in% c("HER2-low (1+)", "HER2-low (2+/FISH-)") ~ 1L, her2_low == "HER2-0" ~ 0L),
         y_low1 = case_when(her2_low == "HER2-low (1+)" ~ 1L, her2_low == "HER2-0" ~ 0L),   # the hard clinical split: IHC 1+ vs 0
         GISTIC = factor(pmax(gistic, -1)),                       # models: categorical (-2, n = 1, merged into -1)
         Purity = estimate_score / 1000,                          # ESTIMATE score per 1000; higher = lower purity
         Age10 = Age / 10,
         Race = factor(Race, c("WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN")),   # 1 American Indian -> NA (too rare)
         Histology = if_else(Histology == "Breast Invasive Lobular Carcinoma", "Lobular", "Non-lobular") |>
           fct_relevel("Non-lobular"),
         across(c(ER, PR, Sex, Era), factor))

# ---- 1. Exploration: clinical metadata vs HER2 IHC status (main target) ----
clin <- c("Age10", "Sex", "Race", "ER", "PR", "Histology", "Era")
lrt <- \(m0, m1) anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2]
ex <- d |> filter(!is.na(y_ihc))
assoc <- map_dfr(clin, \(v) {
  a <- drop_na(ex, all_of(v))
  b <- drop_na(ex, all_of(c(v, "ERBB2", "GISTIC")))
  m_mol <- glm(y_ihc ~ ERBB2 + GISTIC, binomial, b)
  tibble(variable = v, n = nrow(a), n_pos = sum(a$y_ihc),
         p_alone = lrt(glm(y_ihc ~ 1, binomial, a), glm(reformulate(v, "y_ihc"), binomial, a)),
         p_beyond_molecular = lrt(m_mol, update(m_mol, reformulate(c(".", v)))))   # added to RNA + DNA
}) |> mutate(padj_alone = p.adjust(p_alone, "BH"), padj_beyond_molecular = p.adjust(p_beyond_molecular, "BH"))
or <- map_dfr(clin, \(v) broom::tidy(glm(reformulate(v, "y_ihc"), binomial, ex), exponentiate = TRUE, conf.int = TRUE)[-1, ] |>
                mutate(variable = v)) |>
  select(variable, term, OR = estimate, lo = conf.low, hi = conf.high, p_unadj = p.value)
write_csv(assoc, file.path(out, "clinical_assoc.csv")); write_csv(or, file.path(out, "clinical_or.csv"))
cat("Main target: males", sum(ex$Sex == "Male"), "(", sum(ex$Sex == "Male" & ex$y_ihc == 1), "positive); race missing",
    sum(is.na(ex$Race)), "\n")

# Interactions (each its own model; LRT vs main effects on the same tumors)
inter <- list(`ER x PR` = c("ER", "PR"), `Age x ER` = c("Age10", "ER"), `Race x ER` = c("Race", "ER"),
              `ERBB2 mRNA x ER` = c("ERBB2", "ER"), `ERBB2 mRNA x Race` = c("ERBB2", "Race"))
inter <- imap_dfr(inter, \(v, nm) {
  a <- drop_na(ex, all_of(v))
  m0 <- glm(reformulate(v, "y_ihc"), binomial, a)
  tibble(interaction = nm, n = nrow(a), p_interaction = lrt(m0, update(m0, reformulate(c(".", paste(v, collapse = ":"))))))
}) |> mutate(padj = p.adjust(p_interaction, "BH"))
write_csv(inter, file.path(out, "clinical_interactions.csv"))

# ---- 2. Label noise: raw ERBB2 mRNA / GISTIC vs IHC label, by label source and site (no model) ----
lab <- d |> filter(!is.na(y_rep), !is.na(gistic), label_source %in% c("Numeric IHC score", "Reported IHC status only"))
noise_auc <- \(g) map_dfr(c(`ERBB2 mRNA` = "ERBB2", `ERBB2 GISTIC` = "gistic"), \(x) {
  ci <- ci.auc(roc(g$y_rep, g[[x]], direction = "<", quiet = TRUE))
  tibble(n = nrow(g), n_pos = sum(g$y_rep), auc = ci[2], lo = ci[1], hi = ci[3])
}, .id = "predictor")
noise <- bind_rows(
  lab |> group_by(label_source) |> group_modify(\(g, ...) noise_auc(g)) |> mutate(subset = "All sites"),
  lab |> filter(!TSS %in% c("BH", "E9")) |> group_by(label_source) |> group_modify(\(g, ...) noise_auc(g)) |>
    mutate(subset = "Excluding sites BH, E9")) |> ungroup()
write_csv(noise, file.path(out, "label_noise_auc.csv"))
# Where do reported-only tumors come from? (sites with >= 3 reported-only positives)
sites <- lab |> group_by(TSS) |>
  summarise(scored = sum(label_source == "Numeric IHC score"), reported_only = sum(label_source != "Numeric IHC score"),
            reported_only_pos = sum(label_source != "Numeric IHC score" & y_rep == 1),
            reported_only_pos_amplified = sum(label_source != "Numeric IHC score" & y_rep == 1 & gistic == 2)) |>
  arrange(-reported_only_pos)
write_csv(sites, file.path(out, "label_source_by_site.csv"))
lab |> count(label_source, Era, Histology) |> group_by(label_source) |> mutate(pct = round(100 * n / sum(n))) |> print()

p_noise <- lab |> mutate(ihc_binary = paste("IHC", ihc_binary)) |>
  ggplot(aes(ihc_binary, ERBB2)) + geom_boxplot(outlier.shape = NA, color = "grey50") +
  geom_jitter(aes(color = factor(gistic)), width = 0.2, height = 0, size = 0.8, alpha = 0.8) +
  geom_text(data = \(x) count(x, label_source, ihc_binary), aes(y = 19.1, label = paste0("n=", n)), size = 3) +
  scale_color_manual(values = c(`-2` = "#2166ac", `-1` = "#92c5de", `0` = "grey70", `1` = "#f4a582", `2` = "#b2182b"),
                     name = "ERBB2 GISTIC") +
  facet_wrap(~label_source) + labs(x = NULL, y = "ERBB2 mRNA (VST)") + theme_bw()
ggsave(file.path(out, "label_noise.png"), p_noise, width = 7, height = 4)

# ---- 3. RNA vs DNA, raw values (primary test) ----
# Each is one ordered predictor, so AUC needs no fitting: compare raw-value AUCs with a paired DeLong test.
amp <- c("GRB7", "STARD3", "PGAP3")                                # 17q12 amplicon co-expressed genes
cc <- \(y, extra = NULL) d |> filter(!is.na(.data[[y]])) |>       # one complete-case set per target
  drop_na(ERBB2, gistic, ER, PR, Histology, Purity, all_of(c(amp, extra)))
targets <- c(`IHC score 3+ vs 0/1+ (main)` = "y_ihc", `IHC 3+ vs 0/1+, incl. reported status` = "y_rep",
             `Clinical status (FISH-resolved)` = "y_clin", `HER2-low vs HER2-0` = "y_low", `IHC 1+ vs 0 (HER2-)` = "y_low1")
raw_test <- \(a, y) {
  r1 <- roc(a[[y]], a$ERBB2, direction = "<", quiet = TRUE); r2 <- roc(a[[y]], a$gistic, direction = "<", quiet = TRUE)
  t <- roc.test(r1, r2, paired = TRUE, method = "delong")
  tibble(n = nrow(a), n_pos = sum(a[[y]]), auc_rna = as.numeric(r1$auc), auc_dna = as.numeric(r2$auc),
         diff = auc_rna - auc_dna, lo = t$conf.int[1], hi = t$conf.int[2], p = t$p.value)
}
raw <- bind_rows(
  imap_dfr(targets, \(y, nm) raw_test(cc(y), y) |> mutate(target = nm, subset = "All")),
  raw_test(cc("y_low") |> filter(ER == "Positive"), "y_low") |> mutate(target = "HER2-low vs HER2-0", subset = "ER+ only"),
  map_dfr(c("Plates before A213", "Plate A213 and after"), \(e)
          raw_test(cc("y_ihc") |> filter(Era == e), "y_ihc") |> mutate(target = "IHC score 3+ vs 0/1+ (main)", subset = e))) |>
  relocate(target, subset)
write_csv(raw, file.path(out, "rna_vs_dna_raw.csv"))

# ---- 4. Purity and copy number ----
# Low purity dilutes amplification signal; do non-amplified HER2+ tumors have lower purity?
# (HER2+ with a numeric IHC score, to avoid the noisy reported labels; Wilcoxon rank-sum test)
pur <- d |> filter(her2_status == "Positive", label_source == "Numeric IHC score", !is.na(gistic)) |>
  mutate(group = if_else(gistic == 2, "Amplified", "Not amplified"))
pur_test <- pur |> group_by(group) |>
  summarise(n = n(), q1 = quantile(estimate_score, 0.25), median = median(estimate_score), q3 = quantile(estimate_score, 0.75)) |>
  mutate(p_wilcoxon = wilcox.test(estimate_score ~ group, pur)$p.value,
         cohort_min = min(d$estimate_score), cohort_median = median(d$estimate_score), cohort_max = max(d$estimate_score))
write_csv(pur_test, file.path(out, "purity_her2pos_amplification.csv"))

# ---- 5. Cross-validated logistic models ----
# Repeated stratified k-fold; AUC is computed within each test fold and averaged (pooling predictions across folds
# biases AUC downward, most for tied/discrete predictors such as GISTIC). Same folds for every model in an analysis.
# Paired comparisons and CIs use the corrected repeated k-fold t-test (Nadeau & Bengio 2003; Bouckaert & Frank 2004).
make_folds <- \(y, k, reps) map(1:reps, \(r) ave(seq_along(y), y, FUN = \(i) sample(rep_len(1:k, length(i)))))
cv_model <- function(f, data, folds) {
  fits <- map(folds, \(fold) map(sort(unique(fold)), \(j) {
    test <- fold == j
    p <- predict(glm(f, binomial, data[!test, ]), data[test, ], type = "response")
    list(auc = auc_(data$y[test], p), p = set_names(p, which(test)))
  }))
  list(auc = map_dbl(list_flatten(fits), "auc"),                    # one AUC per repeat x fold
       p = map(fits, \(r) unlist(map(r, "p"))[as.character(seq_len(nrow(data)))]))   # out-of-fold predictions per repeat
}
corrected <- \(x, k) {                                               # mean, 95% CI and p (vs 0) over repeat x fold values
  se <- sqrt((1 / length(x) + 1 / (k - 1)) * var(x)); df <- length(x) - 1
  tibble(est = mean(x), lo = mean(x) - qt(0.975, df) * se, hi = mean(x) + qt(0.975, df) * se, p = 2 * pt(-abs(mean(x) / se), df))
}
models <- list(RNA = y ~ ERBB2,
               DNA = y ~ GISTIC,
               `RNA + DNA` = y ~ ERBB2 + GISTIC,
               `Molecular (RNA + DNA + ER + PR)` = y ~ ERBB2 + GISTIC + ER + PR,
               `Molecular + 17q12` = reformulate(c("ERBB2", "GISTIC", "ER", "PR", amp), "y"),
               `RNA x purity` = y ~ ERBB2 * Purity,
               `DNA x purity` = y ~ GISTIC * Purity,
               `RNA + era + histology` = y ~ ERBB2 + Era + Histology,
               `DNA + era + histology` = y ~ GISTIC + Era + Histology)
run <- function(analysis, y, mods, k = 5, reps = 10, extra = NULL) {   # exploring: 5-fold x 10 repeats
  data <- cc(y, extra) |> mutate(y = .data[[y]])
  folds <- make_folds(data$y, k, reps)
  res <- map(mods, cv_model, data = data, folds = folds)
  list(analysis = analysis, data = data, res = res, k = k,
       summary = imap_dfr(res, \(r, m) corrected(r$auc, k) |> select(-p) |>
                            mutate(analysis = analysis, model = m, n = nrow(data), n_pos = sum(data$y), k = k, .before = 1)))
}
runs <- list(
  run("IHC score 3+ vs 0/1+ (main)", "y_ihc", models),
  run("IHC 3+ vs 0/1+, incl. reported status", "y_rep", models[1:5]),
  run("Clinical status (FISH-resolved)", "y_clin", models[1:5]),
  run("HER2-low vs HER2-0", "y_low", c(models[1:3], list(`RNA + ER` = y ~ ERBB2 + ER, `DNA + ER` = y ~ GISTIC + ER))),
  run("IHC 1+ vs 0 (HER2-)", "y_low1", models[1:3]),
  run("Molecular + clinical (tumors with age and race)", "y_ihc",
      list(`Molecular (RNA + DNA + ER + PR)` = models[[4]], `Molecular + age + race` = y ~ ERBB2 + GISTIC + ER + PR + Age10 + Race),
      extra = c("Age10", "Race")),
  run("Final (10-fold)", "y_ihc", models[1:3], k = 10))
cv_auc <- map_dfr(runs, "summary")
write_csv(cv_auc, file.path(out, "cv_auc.csv"))

pairs <- list(c("RNA", "DNA"), c("RNA + DNA", "RNA"), c("Molecular (RNA + DNA + ER + PR)", "RNA + DNA"),
              c("Molecular + 17q12", "Molecular (RNA + DNA + ER + PR)"), c("RNA x purity", "RNA"), c("DNA x purity", "DNA"),
              c("RNA + ER", "DNA + ER"), c("Molecular + age + race", "Molecular (RNA + DNA + ER + PR)"),
              c("RNA + era + histology", "RNA"), c("DNA + era + histology", "DNA"))
cv_compare <- map_dfr(runs, \(a) map_dfr(pairs, \(pr) {
  if (!all(pr %in% names(a$res))) return(NULL)
  corrected(a$res[[pr[1]]]$auc - a$res[[pr[2]]]$auc, a$k) |>
    mutate(analysis = a$analysis, comparison = paste(pr, collapse = " vs "), .before = 1)
})) |>   # BH across exploration comparisons; the final 10-fold run repeats main-target tests, so it is not adjusted
  mutate(padj = if_else(analysis == "Final (10-fold)", NA, p.adjust(if_else(analysis == "Final (10-fold)", NA, p), "BH")))
write_csv(cv_compare, file.path(out, "cv_compare.csv"))

# Fitted coefficients (full data, main target and HER2-low) for the model equations in the report
coefs <- imap_dfr(list(`RNA` = list("y_ihc", models$RNA), `DNA` = list("y_ihc", models$DNA),
                       `RNA + DNA` = list("y_ihc", models$`RNA + DNA`),
                       `Molecular (RNA + DNA + ER + PR)` = list("y_ihc", models[[4]]),
                       `HER2-low: RNA` = list("y_low", models$RNA), `HER2-low: DNA` = list("y_low", models$DNA)), \(m, nm) {
  dat <- cc(m[[1]]) |> mutate(y = .data[[m[[1]]]])
  broom::tidy(glm(m[[2]], binomial, dat), conf.int = TRUE, exponentiate = FALSE) |>
    transmute(model = nm, term, beta = estimate, lo = conf.low, hi = conf.high, p = p.value, n = nrow(dat))
})
write_csv(coefs, file.path(out, "model_coefficients.csv"))

# ---- 6. Tumors where the molecular model disagrees with IHC (out-of-fold, mean over repeats) ----
main <- runs[[1]]
disc <- main$data |> mutate(p_molecular = rowMeans(do.call(cbind, main$res[["Molecular (RNA + DNA + ER + PR)"]]$p))) |>
  filter((y == 1 & p_molecular < 0.2) | (y == 0 & p_molecular > 0.8)) |>
  select(sample_id, TSS, score, category, fish, ERBB2, gistic, ER, PR, Histology, Era, estimate_score, p_molecular) |>
  arrange(-score, p_molecular)
write_csv(disc, file.path(out, "discordant_tumors.csv"))

# ---- Figures ----
fin <- runs[[length(runs)]]                                      # final 10-fold run
roc_df <- imap_dfr(fin$res, \(r, m) { o <- roc(fin$data$y, r$p[[1]], quiet = TRUE)    # repeat 1, pooled for display
  tibble(model = paste0(m, " (CV AUC ", sprintf("%.3f", mean(r$auc)), ")"), fpr = 1 - o$specificities, tpr = o$sensitivities) })
p1 <- ggplot(roc_df, aes(fpr, tpr, color = model)) + geom_path() + geom_abline(lty = 2, color = "grey60") +
  scale_color_manual(values = c("#b2182b", "#2166ac", "grey40")) + coord_equal() +
  labs(x = "False positive rate", y = "True positive rate", color = NULL) +
  theme_bw() + theme(legend.position = c(0.62, 0.15), legend.background = element_blank())
ggsave(file.path(out, "roc_final.png"), p1, width = 5, height = 5)

p2 <- cv_auc |> mutate(model = fct_rev(fct_inorder(model)), analysis = fct_inorder(analysis)) |>
  ggplot(aes(est, model)) + geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.3) +
  facet_grid(analysis ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(22)) +
  labs(x = "Cross-validated AUC (mean of test folds; corrected 95% CI)", y = NULL) +
  theme_bw() + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(out, "cv_auc_all.png"), p2, width = 8, height = 9)

# ---- Summary ----
print(assoc, width = Inf); print(inter); print(noise, width = Inf); print(head(sites, 8)); print(raw, width = Inf)
print(pur_test); print(cv_auc, n = Inf, width = Inf); print(cv_compare, n = Inf, width = Inf)
cat("Discordant tumors:", nrow(disc), "\n"); print(disc, width = Inf, n = Inf)
