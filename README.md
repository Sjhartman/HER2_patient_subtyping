# HER2 patient subtyping

Multiomic analysis of TCGA breast cancer (RNA-seq, ERBB2 copy number, clinical) to identify new patient populations for HER2-directed ADC treatment. Report: `report.Rmd` → `report.html`.

Data are kept outside the repo in `../data_dir/`.

## Steps

1. **QC and normalization** — QC RNA-seq and normalize for read depth.
2. **HER2+ identification** — clinical (IHC/FISH), then multimodal (RNA, copy number, other fields); compare RNA vs DNA for predicting IHC.
3. **Subsets** — unsupervised learning (PAM on filtered, transformed RNA-seq).
4. **Machine learning** — logistic regression (prediction), elastic net (feature selection).
5. **Deep dive** — normal tissue comparison, pathway enrichment by tumor type vs. step 3 partitions, causal biology, pathway differences between tumor types.

## Scripts

Run from repo root (`Rscript analyses/<script>.R`); outputs go to `results/<script>/`.

| Step | Script | Purpose |
|---|---|---|
| 1 | `analyses/01.0_qc_vst_pca.R` | QC, DESeq2 VST, PCA, PC1–5 vs. clinical/technical features |
| 1 | `analyses/01.1_plate_correction.R` | Plate/era assessment: composition, era balance, correction test (not applied) |
