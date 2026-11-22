# HER2 patient subtyping

Multiomic analysis of TCGA breast cancer (RNA-seq, ERBB2 copy number, clinical) to identify new patient populations for HER2-directed ADC treatment. Report: `report.Rmd` → `report.html` (Main results + Supplement, clickable figure/table cross-references).

Data are kept outside the repo in `../data_dir/`.

## Steps

1. **QC and normalization** — QC RNA-seq and normalize for read depth. ✅
2. **HER2+ identification** — clinical (IHC/FISH), then multimodal (RNA, copy number, other fields); compare RNA vs DNA for predicting IHC. ✅
3. **Subsets** — unsupervised learning (PAM on filtered, transformed RNA-seq).
4. **Machine learning** — logistic regression (prediction), elastic net (feature selection).
5. **Deep dive** — normal tissue comparison, pathway enrichment by tumor type vs. step 3 partitions, causal biology, pathway differences between tumor types.

## Scripts

Run from repo root (`Rscript analyses/<script>.R`); outputs go to `results/<script>/`. Render report: `Rscript -e 'rmarkdown::render("report.Rmd")'` (needs pandoc).

| Step | Script | Purpose |
|---|---|---|
| 1 | `analyses/01.0_qc_vst_pca.R` | QC, DESeq2 VST, PCA, PC1–5 vs. clinical/technical features |
| 1 | `analyses/01.1_plate_correction.R` | Plate/era assessment: composition, era balance, correction test (not applied) |
| 2 | `analyses/02.0_her2_labels.R` | HER2 labels (2018-style rules on reported IHC/FISH), HER2-low, disagreements; per-tumor features incl. ESTIMATE purity score and site |
| 2 | `analyses/02.1_rna_vs_dna.R` | Clinical exploration, label-noise check by site, RNA vs DNA (raw paired AUC), purity, CV logistic models, discordant tumors |

## Key notes for later steps

**Data facts**
- RNA input is TCGA upper-quartile (UQ) scaled RSEM, not raw counts (UQ of non-zero values = 1000 per sample). Rounded → DESeq2 size factors (~1) → blind VST. VST is an approximate stabilizer here; use **limma on VST** for differential expression, not DESeq2 count tests.
- 1,212 samples: 1,093 primary tumor, 112 adjacent normal, 7 metastasis; 16,453 genes kept. One tumor per patient; all 112 normals have a matched tumor (89 on the same plate).
- Join keys: RNA `bcr_patient_barcode` → `sample_id` = first 15 chars, matching clinical `Sample ID` and CN `sample_id`. Plate = barcode chars 22–25; tissue source site = chars 6–7.
- Clinical NA strings: `[Not Available]`, `[Not Evaluated]`, `[Unknown]`, `[Discrepancy]`, `Indeterminate`. Stage strings are `Stage IIA` style; T/N codes have substages and TX/NX.
- `results/01.0_qc_vst_pca/vsd.rds` (VST, all samples; gitignored) and `tumor_pcs_features.csv` (tumor PCs + cleaned clinical features + plate/era) are the inputs for later steps.

**Batch decision: no plate/era correction**
- Processing era: `Plates before A213` vs `Plate A213 and after` (boundary where tumor PC3/PC4 plate means step up and normals stop; 110/112 normals are early).
- Era is confounded with sample type and histology (lobular 12% vs 43%). ERBB2 CN is missing for 39% of late tumors (vs 4%). Apparent ERBB2-amp difference (13.5% vs 6.6%) vanishes within ductal tumors (14.6% vs 13.3%).
- ER/PR/HER2 IHC are balanced across eras. Plate explains PC3 0.09 / PC4 0.20 beyond biology (adjusted partial R²), over half of it era.
- Plate correction distorts small plates (single-sample plates set to cohort mean). Protecting labels inflates associations (Nygaard 2016). Would be reconsidered only with cross-plate technical replicates (or a PRPS/RUV-III approach, Molania 2023).

**Step 2 decisions and findings**
- Labels: `results/02.0_her2_labels/her2_labels.csv` (one row per primary tumor; `her2_status`, `her2_low`, `label_source`, `TSS`, `estimate_score`). TCGA IHC/FISH calls follow 2007/2013 criteria; 2018-style rules applied. IHC 0/1+ with FISH+ and reported Positive with FISH- are Unresolved.
- Reported IHC status (no numeric score) is noisy at sites BH and E9; use the **numeric IHC score** for IHC targets. Site A8 has IHC 3+/FISH+ tumors without amplification or high mRNA (annotation issue).
- RNA modestly beats DNA for IHC 3+ (AUC 0.927 vs 0.893, p = 0.045) and clearly for HER2-low vs HER2-0 (0.685 vs 0.573) and IHC 1+ vs 0 (0.652 vs 0.554). Selected model: ERBB2 mRNA alone. Clinical variables, ER/PR, 17q12 genes and purity add nothing.
- CV method: AUC per test fold then averaged (pooling biases discrete predictors); corrected repeated k-fold t-test for comparisons; GISTIC categorical in models. 5-fold x 10 repeats exploring, 10-fold x 10 repeats final.
- Use `her2_status`/`her2_low` only where `label_source == "Numeric IHC score"` when IHC precision matters; flag reported-only labels from sites BH/E9 and IHC 3+ from site A8.
- Candidate populations for HER2-directed ADCs to follow up: HER2-low (394 of 455 scored HER2- tumors); 2 IHC-negative tumors amplified with very high ERBB2 mRNA (TCGA-A2-A0YJ-01, TCGA-A8-A08J-01; re-testing candidates); TCGA "IHC 0" likely includes HER2-ultralow (not recorded).
- RNA signal for IHC 1+ vs 0 is weak (AUC 0.65): a target where multigene models (Step 4) could add value. Benchmark for Step 4: ERBB2 mRNA alone, 10-fold CV AUC 0.927 (IHC 3+ vs 0/1+).
- ESTIMATE score (`estimate_score`, higher = lower purity) is available per tumor for Steps 3–5.

**Report state (continue numbering from here)**
- Main: Sections 1–2; Figures 1–5; Tables 1–4. Next: Figure 6, Table 5.
- Supplement: Sections S1–S9; Figures S1–S5; Tables S1–S13. Next: S10, Figure S6, Table S14.
- References [1]–[17], numbered by first citation in the whole document (main then supplement). After moving text, renumber by first appearance.
- Model inline values come from `results/02.1_rna_vs_dna/model_coefficients.csv` via `b(model, term)` in the setup chunk.

**Plan for Steps 3–5**
- Step 3: PAM on uncorrected VST (refilter genes on tumors). Test clusters vs plate/era **within histology** (e.g. CMH or ductal-only). Check stability.
- Step 4: uncorrected VST; era + histology as sensitivity covariates; no 41-level plate term.
- Step 5: tumor vs normal **within patient** (paired) on uncorrected VST; never apply plate correction including normals (A089/A19E are normal-only plates).

## Workflow
- PI reviews all code; three reviewers (`.claude/agents/`: continuity, technical, reference-master) run after each step. If the agent types aren't loaded (session started outside the repo), run them as general-purpose agents told to read and adopt the agent `.md` file.
- Report rules are in `CLAUDE.md` (numbering, bold citations, S-numbered supplement, verified references with PMIDs, method traceability, symbolic equations only for specific models, rewrite over patch).
- PI preferences seen in Step 2: keep the main text short (a few sentences per result) with detail and caveats in the supplement; every result names its model/test; plan reviewer rounds until a final round finds no real issues (Step 2 took three).
- CV conventions: stratified folds, same folds for all models in an analysis, AUC per test fold then averaged, corrected repeated k-fold t-test with BH across exploration comparisons.
