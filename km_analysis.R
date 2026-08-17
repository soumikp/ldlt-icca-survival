## ============================================================
## KAPLAN-MEIER SURVIVAL ANALYSIS
## Living donor liver transplantation for intrahepatic
## cholangiocarcinoma: recurrence-free and overall survival
## ============================================================
##
## Expects "KM_data.csv" in the working directory. The patient-level
## data are not distributed with this repository; see README.md for
## the required column schema.
##
## Usage: Rscript km_analysis.R

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

## =========================
## 1. Read and prepare data
## =========================

df <- read.csv("KM_data.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

colnames(df) <- c(
  "Pt_no",
  "Recurred",
  "Dead",
  "Date_Tx",
  "Date_LastScan_or_Recurrence",
  "Date_Death_or_LastAlive"
)

cat("=== INPUT DATA ===\n")
print(df)
cat("\nTotal recipients analyzed:", nrow(df), "\n\n")

## Convert dates
df$Date_Tx <- as.Date(df$Date_Tx, format = "%m/%d/%y")
df$Date_LastScan_or_Recurrence <- as.Date(df$Date_LastScan_or_Recurrence, format = "%m/%d/%y")
df$Date_Death_or_LastAlive <- as.Date(df$Date_Death_or_LastAlive, format = "%m/%d/%y")

## Event variables
df$RFS_event <- as.numeric(df$Recurred)
df$OS_event  <- as.numeric(df$Dead)

## Time variables
df$RFS_days <- as.numeric(df$Date_LastScan_or_Recurrence - df$Date_Tx)
df$OS_days  <- as.numeric(df$Date_Death_or_LastAlive - df$Date_Tx)

df$RFS_years <- df$RFS_days / 365.25
df$OS_years  <- df$OS_days / 365.25

df$RFS_months <- df$RFS_years * 12
df$OS_months  <- df$OS_years * 12

## =========================
## 2. Data validation
## =========================

cat("=== DATA VALIDATION ===\n")
cat("Missing OS times:", sum(is.na(df$OS_years)), "\n")
cat("Negative OS times:", sum(df$OS_years < 0, na.rm = TRUE), "\n")
cat("Invalid OS event values:", sum(!df$OS_event %in% c(0, 1)), "\n\n")

## =========================
## 3. Event summary
## =========================

cat("=== EVENT SUMMARY ===\n")
cat("Recurrence events:", sum(df$RFS_event, na.rm = TRUE), "of", nrow(df), "\n")
cat("Death events:", sum(df$OS_event, na.rm = TRUE), "of", nrow(df), "\n\n")

## =========================
## 4. Kaplan-Meier fits (log-log CI for stats, not plotted)
## =========================

fit_rfs <- survfit(Surv(RFS_years, RFS_event) ~ 1, data = df,
                   conf.int = 0.95, conf.type = "log-log")

fit_os <- survfit(Surv(OS_years, OS_event) ~ 1, data = df,
                  conf.int = 0.95, conf.type = "log-log")

## =========================
## 5. Combined KM plot (NO CI shading)
## =========================

max_followup <- max(c(df$RFS_years, df$OS_years), na.rm = TRUE)

p_combined <- ggsurvplot_combine(
  fit = list(RFS = fit_rfs, OS = fit_os),
  data = df,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.28,
  risk.table.fontsize = 4,
  palette = c("#0072B2", "#D55E00"),
  xlab = "Years after transplant",
  ylab = "Survival probability",
  surv.scale = "percent",
  break.time.by = 1,
  xlim = c(0, ceiling(max_followup)),
  ggtheme = theme_classic(base_size = 12),
  legend.title = "",
  legend.labs = c("Recurrence-Free Survival", "Overall Survival"),
  censor.shape = 124,
  censor.size = 3
)

print(p_combined)

## =========================
## 6. Individual KM plots (NO CI shading)
## =========================

p_rfs <- ggsurvplot(
  fit_rfs, data = df,
  conf.int = FALSE,
  risk.table = TRUE, risk.table.height = 0.25,
  palette = "#0072B2",
  xlab = "Years after transplant", ylab = "Recurrence-free survival",
  surv.scale = "percent", break.time.by = 1,
  xlim = c(0, ceiling(max_followup)),
  ggtheme = theme_classic(base_size = 12),
  title = "Recurrence-Free Survival",
  legend = "none", censor.shape = 124, censor.size = 3
)

p_os <- ggsurvplot(
  fit_os, data = df,
  conf.int = FALSE,
  risk.table = TRUE, risk.table.height = 0.25,
  palette = "#D55E00",
  xlab = "Years after transplant", ylab = "Overall survival",
  surv.scale = "percent", break.time.by = 1,
  xlim = c(0, ceiling(max_followup)),
  ggtheme = theme_classic(base_size = 12),
  title = "Overall Survival",
  legend = "none", censor.shape = 124, censor.size = 3
)

print(p_rfs)
print(p_os)

## =========================
## 7. Median survival with 95% CI (reported in text)
## =========================

med_rfs <- surv_median(fit_rfs)
med_os  <- surv_median(fit_os)

cat("=== MEDIAN SURVIVAL ===\n")

if (is.na(med_rfs$median)) {
  cat("Median RFS: Not reached (NR)\n")
} else {
  cat("Median RFS:", round(med_rfs$median, 2), "years\n")
  cat("95% CI: [", round(med_rfs$lower, 2), " - ", round(med_rfs$upper, 2), "] years\n", sep = "")
}

if (is.na(med_os$median)) {
  cat("Median OS: Not reached (NR)\n")
} else {
  cat("Median OS:", round(med_os$median, 2), "years\n")
  cat("95% CI: [", round(med_os$lower, 2), " - ", round(med_os$upper, 2), "] years\n", sep = "")
}
cat("\n")

## =========================
## 8. Time-specific survival with 95% CI
## =========================

time_points <- c(0.5, 1, 2, 3, 4, 5)

rfs_summary <- summary(fit_rfs, times = time_points, extend = FALSE)
os_summary  <- summary(fit_os,  times = time_points, extend = FALSE)

cat("=== TIME-SPECIFIC SURVIVAL RATES (95% CI) ===\n\n")

cat("Recurrence-Free Survival:\n")
for (i in seq_along(rfs_summary$time)) {
  cat(sprintf("  %.1f-year: %.1f%% (95%% CI: %.1f%% - %.1f%%), n at risk = %d\n",
              rfs_summary$time[i],
              rfs_summary$surv[i] * 100,
              rfs_summary$lower[i] * 100,
              rfs_summary$upper[i] * 100,
              rfs_summary$n.risk[i]))
}

cat("\nOverall Survival:\n")
for (i in seq_along(os_summary$time)) {
  cat(sprintf("  %.1f-year: %.1f%% (95%% CI: %.1f%% - %.1f%%), n at risk = %d\n",
              os_summary$time[i],
              os_summary$surv[i] * 100,
              os_summary$lower[i] * 100,
              os_summary$upper[i] * 100,
              os_summary$n.risk[i]))
}
cat("\n")

## =========================
## 9. MEDIAN FOLLOW-UP as median (IQR)
## Based on observed follow-up times
## =========================

median_fu <- median(df$OS_months, na.rm = TRUE)
q1_fu     <- quantile(df$OS_months, 0.25, na.rm = TRUE)
q3_fu     <- quantile(df$OS_months, 0.75, na.rm = TRUE)
min_fu    <- min(df$OS_months, na.rm = TRUE)
max_fu    <- max(df$OS_months, na.rm = TRUE)

cat("=== MEDIAN FOLLOW-UP (observed) ===\n")
cat(sprintf("Median follow-up: %.1f months (IQR: %.1f - %.1f)\n",
            median_fu, q1_fu, q3_fu))
cat(sprintf("Range: %.1f - %.1f months\n", min_fu, max_fu))
cat(sprintf("In years: %.2f (IQR: %.2f - %.2f)\n\n",
            median_fu/12, q1_fu/12, q3_fu/12))

## =========================
## 10. Time to recurrence
## =========================

df_rec <- subset(df, RFS_event == 1)

if (nrow(df_rec) > 0) {
  cat("=== TIME TO RECURRENCE (n = ", nrow(df_rec), ") ===\n", sep = "")
  for (i in seq_len(nrow(df_rec))) {
    cat(sprintf("  Pt %s: %.1f months (%.2f years) post-transplant\n",
                df_rec$Pt_no[i], df_rec$RFS_months[i], df_rec$RFS_years[i]))
  }
}
cat("\n")

## =========================
## 10B. NUMBER OF PATIENTS SURVIVING BEYOND EACH TIMEPOINT
## =========================

cat("=== PATIENTS SURVIVING BEYOND EACH TIMEPOINT ===\n\n")

year_marks  <- c(12, 24, 36, 48, 60)   # months
year_labels <- c("1 year", "2 years", "3 years", "4 years", "5 years")

## Patients who actually survived beyond each timepoint
## (either alive with follow-up past it, or died after it)
cat("Patients who survived beyond each timepoint:\n")
for (i in seq_along(year_marks)) {
  n_survived <- sum(df$OS_months >= year_marks[i], na.rm = TRUE)
  cat(sprintf("  Survived >= %s: %d of %d patients\n",
              year_labels[i], n_survived, nrow(df)))
}
cat("\n")

## Detailed per-patient status (sorted by follow-up length)
cat("Per-patient follow-up status (sorted):\n")
df_sorted <- df[order(df$OS_months), ]
for (i in seq_len(nrow(df_sorted))) {
  status <- ifelse(df_sorted$OS_event[i] == 1, "DIED", "alive")
  cat(sprintf("  Pt %s: %.1f months (%.2f years) - %s\n",
              df_sorted$Pt_no[i],
              df_sorted$OS_months[i],
              df_sorted$OS_years[i],
              status))
}
cat("\n")

## =========================
## 11. Helper functions for CI extraction
## =========================

get_survival <- function(fit, time) {
  x <- summary(fit, times = time, extend = FALSE)
  if (length(x$surv) == 0) return(NA_real_)
  return(x$surv[1] * 100)
}

get_survival_lower <- function(fit, time) {
  x <- summary(fit, times = time, extend = FALSE)
  if (length(x$lower) == 0) return(NA_real_)
  return(x$lower[1] * 100)
}

get_survival_upper <- function(fit, time) {
  x <- summary(fit, times = time, extend = FALSE)
  if (length(x$upper) == 0) return(NA_real_)
  return(x$upper[1] * 100)
}

format_ci <- function(fit, time) {
  est <- get_survival(fit, time)
  lo  <- get_survival_lower(fit, time)
  up  <- get_survival_upper(fit, time)
  if (is.na(est)) return("NR")
  sprintf("%.1f (95%% CI: %.1f-%.1f)", est, lo, up)
}

## =========================
## 12. Summary table with CIs and survivor counts
## =========================

median_rfs_value <- ifelse(is.na(med_rfs$median), "NR",
  sprintf("%.2f (95%% CI: %.2f-%.2f)", med_rfs$median, med_rfs$lower, med_rfs$upper))

median_os_value <- ifelse(is.na(med_os$median), "NR",
  sprintf("%.2f (95%% CI: %.2f-%.2f)", med_os$median, med_os$lower, med_os$upper))

## Survivor counts for the table
n_surv_1yr <- sum(df$OS_months >= 12, na.rm = TRUE)
n_surv_2yr <- sum(df$OS_months >= 24, na.rm = TRUE)
n_surv_3yr <- sum(df$OS_months >= 36, na.rm = TRUE)
n_surv_4yr <- sum(df$OS_months >= 48, na.rm = TRUE)
n_surv_5yr <- sum(df$OS_months >= 60, na.rm = TRUE)

summary_table <- data.frame(
  Metric = c(
    "Total patients (n)",
    "Recurrence events",
    "Death events",
    "Median follow-up (months)",
    "Follow-up Q1 (months)",
    "Follow-up Q3 (months)",
    "Follow-up min (months)",
    "Follow-up max (months)",
    "Patients survived >= 1 year",
    "Patients survived >= 2 years",
    "Patients survived >= 3 years",
    "Patients survived >= 4 years",
    "Patients survived >= 5 years",
    "Median RFS (years)",
    "Median OS (years)",
    "6-month RFS (%)",
    "1-year RFS (%)",
    "2-year RFS (%)",
    "3-year RFS (%)",
    "4-year RFS (%)",
    "5-year RFS (%)",
    "6-month OS (%)",
    "1-year OS (%)",
    "2-year OS (%)",
    "3-year OS (%)",
    "4-year OS (%)",
    "5-year OS (%)"
  ),
  Value = c(
    nrow(df),
    sum(df$RFS_event, na.rm = TRUE),
    sum(df$OS_event, na.rm = TRUE),
    round(median_fu, 1),
    round(q1_fu, 1),
    round(q3_fu, 1),
    round(min_fu, 1),
    round(max_fu, 1),
    n_surv_1yr,
    n_surv_2yr,
    n_surv_3yr,
    n_surv_4yr,
    n_surv_5yr,
    median_rfs_value,
    median_os_value,
    format_ci(fit_rfs, 0.5),
    format_ci(fit_rfs, 1),
    format_ci(fit_rfs, 2),
    format_ci(fit_rfs, 3),
    format_ci(fit_rfs, 4),
    format_ci(fit_rfs, 5),
    format_ci(fit_os, 0.5),
    format_ci(fit_os, 1),
    format_ci(fit_os, 2),
    format_ci(fit_os, 3),
    format_ci(fit_os, 4),
    format_ci(fit_os, 5)
  ),
  stringsAsFactors = FALSE
)

cat("=== SUMMARY TABLE ===\n")
print(summary_table, row.names = FALSE)

## =========================
## 13. Save summary
## =========================

write.csv(summary_table, "survival_summary_table.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")