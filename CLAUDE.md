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
- Two parts: Main results and Supplement (supporting/technical detail, e.g. batch artifacts). Supplementary items are numbered S1, S2, … (e.g. **(Figure S1)**).
- Figures and tables numbered separately, with succinct legends.
- When revising (e.g. after review), do whatever gives the clearest result: rewrite a whole section or figure legend rather than patching sentences, which can leave text more confusing.
- Every statistical or ML result must be traceable to its method, stated before or in the same sentence. For regression/ML models that are specific (where several formulations are plausible, e.g. which terms or codings), give the model as a centered equation with symbolic coefficients (e.g. $$\text{logit}\,P(\text{IHC 3+}) = \beta_0 + \beta_1\,\text{ERBB2}$$), define each term, put fitted values in a table, and quote key values (e.g. $\beta_1$) in the text, and name the test or metric (e.g. AUC, paired DeLong test, corrected CV t-test). Do not write out self-evident equations (e.g. a univariable logistic model, the definition of AUC).
- Cite only key findings or insights (guidelines, trials, landmark biology, supporting/conflicting evidence). Do not cite software or routine methods (DESeq2, limma, pROC, DeLong, etc.).
- References section with PMIDs and DOI links (`DOI: [10.x/y](https://doi.org/10.x/y)`); cite inline as [n], numbered by first citation. Verify every reference exists. The script at the end of `report.Rmd` auto-links [n] and (Figure/Table X) citations, with back-links.
- No plate/era batch correction (era is confounded with sample type and histology); see report Step 1.3.
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
- `reference-master` — references are real, DOIs/PMIDs resolve, formatting is consistent, and each [n] supports its claim.
- Run all three after each step's script and report section are drafted; summarize their findings for the PI briefly.
