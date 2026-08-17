# ldlt-icca-survival

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

Kaplan–Meier survival analysis code accompanying our study of **living donor liver transplantation (LDLT) for intrahepatic cholangiocarcinoma (iCCA)**.

This repository contains the analysis script used to produce the survival estimates and figures reported in the manuscript, along with the resulting figures and summary tables.

## Contents

```
km_analysis.R                        Full analysis script
output/
  Figure1_Combined_KM.png            RFS and OS on shared axes, with risk table
  Figure2_RFS.png                    Recurrence-free survival
  Figure3_OS.png                     Overall survival
  survival_estimates_with_CI.csv     Time-specific survival with 95% CIs
  survival_summary_table.csv         Full summary of all reported metrics
```

## Data availability

The patient-level dataset is **not included** in this repository. The study cohort is small and consists of transplant, recurrence, and death dates, which cannot be shared publicly without compromising patient confidentiality. Requests for access may be directed to the corresponding author, subject to institutional review board approval and a data use agreement.

To run the script on your own cohort, provide a `KM_data.csv` in the working directory with six columns in this order:

| Position | Meaning | Format |
|---|---|---|
| 1 | Patient identifier | any |
| 2 | Recurrence indicator | `1` = recurred, `0` = censored |
| 3 | Death indicator | `1` = died, `0` = alive at last contact |
| 4 | Date of transplant | `MM/DD/YY` |
| 5 | Date of recurrence, or of last imaging without recurrence | `MM/DD/YY` |
| 6 | Date of death, or of last contact if alive | `MM/DD/YY` |

Column *names* are reassigned positionally by the script, so headers may differ; column *order* must match. The file is read with `fileEncoding = "UTF-8-BOM"`.

## Analysis

Time to event is measured from the date of transplant.

- **Recurrence-free survival (RFS)** — from transplant to radiographic recurrence, censored at last imaging without recurrence.
- **Overall survival (OS)** — from transplant to death, censored at last known contact.

Survival is estimated by the Kaplan–Meier product-limit method. Confidence intervals use the **complementary log-log transformation**, which is preferred over the plain-scale interval in small cohorts because it constrains bounds to the interval [0, 1] and behaves better in the tails. Median survival and its 95% CI are obtained from `survminer::surv_median()`. Follow-up duration is reported as a median with interquartile range over observed follow-up times.

The script also reports time-specific survival at 0.5, 1, 2, 3, 4, and 5 years, per-patient time to recurrence, and the number of patients surviving beyond each annual landmark. Time-specific estimates are computed with `extend = FALSE`, so no estimate is reported at a timepoint beyond the last observed follow-up.

Given the cohort size, confidence intervals are wide and these estimates should be read as descriptive rather than inferential.

## Requirements

R (≥ 4.0) with:

```r
install.packages(c("survival", "survminer", "dplyr", "ggplot2"))
```

## Usage

```bash
Rscript km_analysis.R
```

Figures are drawn to the active graphics device and the summary table is printed to the console. The `write.csv()` call that saves the summary table is commented out at the end of the script; uncomment it to write `survival_summary_table.csv`.

## Citation

Manuscript under review. Citation details will be added upon publication.

## License

MIT — see [LICENSE](LICENSE).
