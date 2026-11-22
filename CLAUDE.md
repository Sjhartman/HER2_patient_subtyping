# CLAUDE.md

## Roles
- The user is the Principal Investigator (computational oncology); Claude is the co-scientist.
- Goal: a multiomic HTML report identifying new patient populations for HER2-directed ADC treatment (TCGA BRCA).
- The PI manually reviews all code and gives feedback; wait for it rather than racing ahead.

## Communication
- Be brief. Plain language, minimal jargon, no bloat in responses or code.

## Code
- R + tidyverse + ggplot2. One plain `.R` script per step, numbered by challenge step: `<step>.<sub>_<name>.R` (e.g. `01.0_qc.R`, `01.1_normalize.R`).
- Scripts in `analyses/`, outputs in `results/<script name>/`. Shared helpers unnumbered.
- Concise and legible; no excess functions or abstraction.
- Add each new script to the README under its step.
- Never add `CLAUDE.md` to `.gitignore`.

## Data (outside repo; never commit)
`../data_dir/`
- `Tempus_Pharma_RD_Coding_Challenge.v2026.1.pdf` — instructions
- `brca_tcga_clinical_data.csv` — clinical (incomplete)
- `brca_tcga_erbb2_copy_number.csv` — ERBB2 copy number (-2 to 2)
- `tcga.brca.rsem.csv` — RSEM counts; cols 1–3 metadata (`patient_id`, `bcr_patient_barcode`, `sample_type`), col 4+ genes

## Report
- `report.Rmd` → `report.html` (R Markdown).
- Opening statement: Claude Code was run in manual mode using Opus 5.5.
- Figures and tables numbered separately, with succinct legends.
- Every figure/table is cited in the text with bolded parentheses, e.g. **(Figure 1)**, **(Table 2)**.

## Outline (challenge steps 1–5)
1. QC and normalize RNA-seq (adjust for read depth).
2. Identify HER2+ patients: clinical (IHC/FISH), then multimodal (RNA, copy number, other fields); is RNA or DNA more predictive of IHC?
3. Subsets via unsupervised learning (PAM on filtered, transformed RNA-seq).
4. ML: logistic regression for prediction, elastic net for feature selection; focus on feature importance and biological validity.
5. Deep dive: normal tissue comparison, pathway enrichment per tumor type and relation to step 3 partitions, causal biology, pathway differences between tumor types.

## Reviewers (`.claude/agents/`)
- `continuity-reviewer` — logical flow, cancer biology, supporting/conflicting references, figure/table citations.
- `technical-reviewer` — methods are correct and current (well-recognized best practice).
- Run both after each step's script and report section are drafted; summarize their findings for the PI briefly.
