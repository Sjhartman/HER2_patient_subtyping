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

_Added as generated._
