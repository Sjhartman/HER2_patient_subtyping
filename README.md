# HER2 patient subtyping

Multiomic analysis of TCGA breast cancer (RNA-seq, ERBB2 copy number, clinical) to identify new patient populations for HER2-directed ADC treatment. Report: `report.Rmd` → `report.html` (Main results + Supplement, clickable figure/table cross-references).

Data are kept outside the repo in `../data_dir/`.

## Steps

1. **QC and normalization** — QC RNA-seq and normalize for read depth. ✅
2. **HER2+ identification** — clinical (IHC/FISH), then multimodal (RNA, copy number, other fields); compare RNA vs DNA for predicting IHC.
3. **Subsets** — unsupervised learning (PAM on filtered, transformed RNA-seq).
4. **Machine learning** — logistic regression (prediction), elastic net (feature selection).
5. **Deep dive** — normal tissue comparison, pathway enrichment by tumor type vs. step 3 partitions, causal biology, pathway differences between tumor types.

## Scripts

Run from repo root (`Rscript analyses/<script>.R`); outputs go to `results/<script>/`. Render report: `Rscript -e 'rmarkdown::render("report.Rmd")'` (needs pandoc).

| Step | Script | Purpose |
|---|---|---|
| 1 | `analyses/01.0_qc_vst_pca.R` | QC, DESeq2 VST, PCA, PC1–5 vs. clinical/technical features |
| 1 | `analyses/01.1_plate_correction.R` | Plate/era assessment: composition, era balance, correction test (not applied) |
| 2 | `analyses/02.0_her2_labels.R` | Clinical HER2 status (ASCO/CAP 2018), HER2-low categories, IHC-only target, label disagreements |

## Key notes for later steps (from Step 1)

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

**Plan / requirements for Step 2 (from reviewers)**
- Two separate questions: (a) clinical HER2 status (IHC + FISH, ASCO/CAP 2018), reported on its own; (b) **RNA vs DNA predicting the IHC label** — predict IHC from ERBB2 mRNA, from GISTIC CN (ordinal, 5 levels), and both.
- **No leakage:** FISH status, `HER2 cent17 ratio`, `HER2 copy number` help define the clinical label, so they are not DNA predictors of it.
- Fair comparison: same complete-case samples (CN missingness is era-dependent), cross-validated AUC, paired test (bootstrap/DeLong). Note that continuous mRNA has an advantage over 5-level GISTIC.
- IHC-equivocal: exclude from binary target or resolve by FISH. FISH covers 416 tumors, only 39% IHC-equivocal.
- **HER2-low** (ADC-relevant): `HER2 ihc score` 0 / 1+ / 2+ / 3+ = 61 / 273 / 200 / 90 (484 missing). Label conflicts exist (11 tumors scored 1+ but labelled Positive, 2 scored 3+ but labelled Negative). 246 "Negative" tumors have no score. Build an explicit label table with rules and report exclusions.
- Sensitivity covariates: era **plus histology** (era alone absorbs lobular biology); also report performance by era. Include era/purity for the CN model too.
- Lead: tumor TCGA-E9-A1NC-01 is HER2 IHC+ but ERBB2 mRNA low (RNA/IHC discordance).

**Report notes for 02.0 figure (`erbb2_by_ihc.png`) legend**
- Tumors without a numeric IHC score (n = 473) are not shown; most (290) still have a reported IHC status (`IHC-HER2`: 242 Negative, 41 Positive, 7 Equivocal), which is used in the labels and later steps.
- ERBB2 VST values below 10 are displayed at 10 (n = 2 plotted).

**PI decisions for 02.1**
- Equivocal: two steps. Main = IHC-only target (3+ vs 0/1+, equivocal excluded); sensitivity = FISH-resolved clinical status.
- Exploratory first: test clinical metadata (race, age, sex, plus ER/PR, histology, era) for association with HER2 status, including separate interaction models. Include only variables that help.
- Models: ERBB2 mRNA; GISTIC CN; both; multimodal molecular (+ ER, PR status; with vs without 17q12 amplicon genes GRB7/STARD3/PGAP3); and molecular + clinical if exploration supports it.
- CV: k = 3 folds while exploring; k = 10 for the final selected model.

**Plan for Steps 3–5**
- Step 3: PAM on uncorrected VST (refilter genes on tumors). Test clusters vs plate/era **within histology** (e.g. CMH or ductal-only). Check stability.
- Step 4: uncorrected VST; era + histology as sensitivity covariates; no 41-level plate term.
- Step 5: tumor vs normal **within patient** (paired) on uncorrected VST; never apply plate correction including normals (A089/A19E are normal-only plates).

## Workflow
- PI reviews all code; three reviewers (`.claude/agents/`: continuity, technical, reference-master) run after each step. If the agent types aren't loaded (session started outside the repo), run them as general-purpose agents told to read and adopt the agent `.md` file.
- Report rules are in `CLAUDE.md` (numbering, bold citations, S-numbered supplement, verified references with PMIDs).
