# HER2 patient subtyping

Multiomic analysis of TCGA breast cancer (RNA-seq, ERBB2 copy number, clinical) to identify new patient populations for HER2-directed ADC treatment. Report: `report.Rmd` → `report.html` (Main results + Supplement, clickable figure/table cross-references).

Data are kept outside the repo in `../data_dir/`.

## Steps

1. **QC and normalization** — QC RNA-seq and normalize for read depth. ✅
2. **HER2+ identification** — clinical (IHC/FISH), then multimodal (RNA, copy number, other fields); compare RNA vs DNA for predicting IHC. ✅
3. **Subsets** — unsupervised learning (PAM on filtered, transformed RNA-seq). ✅
4. **Machine learning** — logistic regression (prediction), elastic net (feature selection). ✅ (pending review)
5. **Deep dive** — normal tissue comparison, pathway enrichment by tumor type vs. step 3 partitions, causal biology, pathway differences between tumor types.

## Scripts

Run from repo root (`Rscript analyses/<script>.R`); outputs go to `results/<script>/`. Render report: `Rscript -e 'rmarkdown::render("report.Rmd")'` (needs pandoc).

| Step | Script | Purpose |
|---|---|---|
| 1 | `analyses/01.0_qc_vst_pca.R` | QC, DESeq2 VST, PCA, PC1–5 vs. clinical/technical features |
| 1 | `analyses/01.1_plate_correction.R` | Plate/era assessment: composition, era balance, correction test (not applied) |
| 2 | `analyses/02.0_her2_labels.R` | HER2 labels (2018-style rules on reported IHC/FISH), HER2-low, disagreements; per-tumor features incl. ESTIMATE purity score and site |
| 2 | `analyses/02.1_rna_vs_dna.R` | Clinical exploration, label-noise check by site, RNA vs DNA (raw paired AUC), purity, CV logistic models, discordant tumors |
| 3 | `analyses/03.0_pam_clusters.R` | PAM clusters (top 2000 genes, k by silhouette + consensus PAC), robustness, PAM50 (genefu), HER2/batch associations, HER2-low by cluster, C5 ERBB2 modes |
| 4 | `analyses/04.0_c5_low_enet.R` | C5 ERBB2-low vs C5 ERBB2-high and vs C1–C3 pooled: elastic net (purity unpenalized), selection stability, compact logistic model, sensitivity (no purity; + era/lobular) |
| 5 | `analyses/05.0_c5_low_gsea.R` | C5 ERBB2-low vs C1–C3 and vs C5 ERBB2-high: limma + TREAT (unadjusted and ESTIMATE-adjusted), volcano plots, GSEA (fgsea Hallmark/Reactome, camera check), MCP-counter immune/stromal scores (markers in `analyses/ref/`), ER IHC % check |
| 5 | `analyses/05.1_c5_low_vs_normal.R` | C5 ERBB2-low tumors vs matched adjacent normal (11 pairs): paired limma + TREAT, volcano with Step 4 genes, Hallmark GSEA + camera, paired MCP-counter, representativeness of the pairs |
| – | `analyses/helpers.R` | Shared helpers: current gene symbols without alias collisions; MCP-counter marker list |

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
- Restructured after Step 5: Key findings box (client-style bullets, ~180 words, no citations; CSS class .key-findings) at top; Main Sections 1–6 (5.5 = matched normal, PI: critical; 6 = Conclusions, the single implications section); Figures 1–12; Tables 1–4.
- Supplement: Sections S1–S14 (S14 = matched-normal pair details); Figures S1–S10; Tables S1–S31. Supplementary items are numbered in Supplement order (main-text items in citation order). Removed: old plate small-sample table, luminal profile table, leading-edge table, Reactome table.
- References [1]–[28], numbered by first citation in the whole document (main then supplement). After moving text, renumber by first appearance.
- Model inline values come from `results/02.1_rna_vs_dna/model_coefficients.csv` via `b(model, term)` in the setup chunk; Step 3 and 4 values are read from `results/03.0_pam_clusters/*.csv` and `results/04.0_c5_low_enet/*.csv` in the setup chunk.
- Figure colors (Step 3 on): PAM clusters C1–C5 = Dark2 (teal, orange, lavender, magenta, green); PAM50 = browns/yellow/plum/aqua; ERBB2 GISTIC = blue–red (NA white). Reuse these in later steps.

**Step 3 findings to carry forward**
- C5 ERBB2-low group (non-amplified, mostly HER2-low among scored HER2−, 78% ER+) is the least pure group: median ESTIMATE 675 vs 9 in ERBB2-high, with high stromal and immune markers (Table 6); ERBB2 mRNA diluted by infiltration (ρ = −0.30), but the purest third is still at luminal ERBB2 level, and IHC is mostly 1+/2+.
- Possible amplified-but-diluted tumors in C5 ERBB2-low: TCGA-C8-A26X-01 (GISTIC 2, IHC 1+, called HER2−: re-testing candidate), TCGA-C8-A138-01 and TCGA-C8-A1HL-01 (GISTIC 2, IHC 2+, no FISH: unresolved). In report Section 3.2 and Table S15 (all 19 amplified tumors without a HER2+ call).

**Step 4 findings to carry forward**
- C5 ERBB2-low vs C5 ERBB2-high: the 17q12 amplicon is sufficient (STARD3 top gene in 94% of training sets; 1-gene CV AUC 0.996). This is near-circular (C5 split on ERBB2 mRNA). Non-amplicon stable genes = what adds beyond the amplicon only, not all differences (ER, infiltration also differ) -> Step 5 limma.
- C5 ERBB2-low vs C1–C3 (pooled, n = 717): 43 stable genes. Lower luminal/ER genes (PGR, BCL2, MAPT, ELOVL2, GFRA1, SERPINA3, IGF1R) despite 78% ER+ IHC; higher VSTM2L, KRT86, MMP1, TMEM45A, S100A9, SOX11. Intermediate between luminal and C5-high/C4; PAM50 mixed (27 Her2, 36 LumB, 17 LumA). Not a HER2-driven group; do not call it HER2-enriched.
- Infiltration check: purity-matched C1–C3 (3 per C5-low tumor by ESTIMATE) keeps all 43 directions (differences ~96% of unmatched); within C1–C3, 7 genes reach |rho| >= 0.4 with ESTIMATE/stroma/immune (TMEM45A, SOX11, S100A9, CAVIN2, MMRN1, TMEM132C, RERG).
- 7-gene compact model: genes most often in the in-fold top 7 (VSTM2L, KRT86, MMP1, TMEM45A up; MAPT, JADE2/PHF15, BCL2 down), Firth fit; CV AUC 0.964 (procedure-level). Only VSTM2L and MAPT are consistently picked.
- Gene lists for Step 5: `results/04.0_c5_low_enet/features.csv` (stable = `sel_freq >= 0.7`; `symbol` = current HGNC). Step 5 candidate programs: ER signalling, proliferation, hypoxia (TMEM45A), immune/myeloid (S100A9).
- Packages added: `logistf` (Firth).

**Step 5 findings and conventions**
- Gene symbols: use `current_symbols()` in `analyses/helpers.R` (alias remapping alone overwrote real genes, e.g. C2orf48 -> RRM2). MCP-counter markers via `mcp_markers()` (source file lists CAVIN2 as "C").
- DE: limma-trend robust on VST, TREAT (|log2 FC| > 1) for significance; GSEA = fgsea on moderated t (ranks programs) + camera check (default 0.01 and estimated inter-gene correlation; "robust" = camera estimated FDR < 0.05). MCP-counter group tests use HC3 robust SEs.
- C5 ERBB2-low vs C1–C3: camera-confirmed in both models = low early estrogen response only (ESR1 −2.4 log2, unchanged by purity adjustment; holds within IHC ER+ tumors; staining somewhat weaker, 53% vs 72% at ≥80%, but only 1/32 ER-low <10%). Glycolysis and mTORC1 camera-confirmed only after purity adjustment (suggestive). NK MCP score uses 3 of 9 markers (low-expression filter). Immune/EMT, proliferation, hypoxia up but not robust to strict camera. Infiltrate: more T/B/monocytic/fibroblast; relative to ESTIMATE, more monocytic-lineage and fibroblast signal.
- vs C5 ERBB2-high: only the 17q12 amplicon passes TREAT; proliferation lower (fgsea, camera default only).
- Tumor vs normal (11 pairs, sites BH/E9/GI): Step 4 genes move the same way as vs luminal; not shown to be group-specific. Nice-to-have (noted in report 5.1): same comparison in the 74 luminal pairs and a difference-in-differences test.

**Plan for Steps 3–5**
- Step 3: PAM on uncorrected VST (refilter genes on tumors). Test clusters vs plate/era **within histology** (e.g. CMH or ductal-only). Check stability.
- Step 4: uncorrected VST; era + histology as sensitivity covariates; no 41-level plate term.
- Step 4 priority (PI): elastic net feature selection on **C5 (Her2) vs. the luminal clusters C1–C3** from Step 3, to find what makes C5 unique (C5 vs. all others would mostly relearn the ERBB2 amplicon/basal axis). Consider purity (ESTIMATE) as a covariate and the C5 ERBB2-low group (`c5_mode`) as a sub-analysis. Labels: `results/03.0_pam_clusters/clusters.csv`.
- Step 3 caveats to carry: clusters are regions of a continuum (k = 5 balances silhouette and PAC); residual era association among ER+ ductal tumors (later-era ER+ more often C4); HER2-low % among scored HER2− is inflated (IHC 0 likely recorded as "Negative" without score).
- Step 5: tumor vs normal **within patient** (paired) on uncorrected VST; never apply plate correction including normals (A089/A19E are normal-only plates).
- Step 5 (pinned from Step 3): differential expression (limma on VST) of C5 ERBB2-low (non-amplified, mostly HER2-low, ER+) vs. luminal clusters C1–C3 and vs. C5 ERBB2-high, to find what makes them HER2-enriched-like without amplification. Inputs: `results/03.0_pam_clusters/clusters.csv` (`cluster`, `c5_mode`).

## Workflow
- PI reviews all code; three reviewers (`.claude/agents/`: continuity, technical, reference-master) run after each step. If the agent types aren't loaded (session started outside the repo), run them as general-purpose agents told to read and adopt the agent `.md` file.
- Report rules are in `CLAUDE.md` (numbering, bold citations, S-numbered supplement, verified references with PMIDs, method traceability, symbolic equations only for specific models, rewrite over patch).
- PI preferences seen in Step 2: keep the main text short (a few sentences per result) with detail and caveats in the supplement; every result names its model/test; plan reviewer rounds until a final round finds no real issues (Step 2 took three).
- CV conventions: stratified folds, same folds for all models in an analysis, AUC per test fold then averaged, corrected repeated k-fold t-test with BH across exploration comparisons.
