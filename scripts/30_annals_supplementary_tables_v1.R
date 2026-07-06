#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

outdir <- "outputs/annals_supplement_v1"
tabdir <- file.path(outdir, "tables")
srcdir <- file.path(outdir, "source_data")
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

db_levels <- c("MIMIC-IV", "eICU", "SICdb")

fmt_n <- function(x) {
  ifelse(is.na(x), "NA", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
}
fmt_pct <- function(num, den, digits = 1) {
  if (is.na(num) || is.na(den) || den == 0) return("NA")
  paste0(sprintf(paste0("%.", digits, "f"), 100 * num / den), "%")
}
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}
fmt_or <- function(or, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}
fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
fmt_med_iqr <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  paste0(
    sprintf(paste0("%.", digits, "f"), q[2]),
    " [",
    sprintf(paste0("%.", digits, "f"), q[1]),
    "-",
    sprintf(paste0("%.", digits, "f"), q[3]),
    "]"
  )
}

write_md_table <- function(df, path, title = NULL, note = NULL) {
  lines <- character()
  if (!is.null(title)) lines <- c(lines, paste0("# ", title), "")
  lines <- c(lines, paste(names(df), collapse = " | "))
  lines <- c(lines, paste(rep("---", ncol(df)), collapse = " | "))
  for (i in seq_len(nrow(df))) {
    vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
    vals <- gsub("\\|", "/", vals)
    lines <- c(lines, paste(vals, collapse = " | "))
  }
  if (!is.null(note)) lines <- c(lines, "", note)
  writeLines(lines, path)
}

complete_n <- function(dat, vars) {
  vars <- vars[vars %in% names(dat)]
  if (!length(vars)) return(NA_integer_)
  sum(complete.cases(dat[, vars, drop = FALSE]))
}

patient_level_path <- Sys.getenv(
  "HRC_PATIENT_LEVEL_ANALYSIS_DATA",
  "outputs_patient_level/analysis_dataset_core_hrc.csv"
)
if (!file.exists(patient_level_path)) {
  stop(
    "Patient-level analysis file not found. Generate it locally and set ",
    "HRC_PATIENT_LEVEL_ANALYSIS_DATA, but do not commit it to the public repository."
  )
}

analysis_df <- read_csv(patient_level_path, show_col_types = FALSE) %>%
  mutate(database = factor(database, levels = db_levels))

primary_n <- analysis_df %>%
  count(database, name = "primary_analytic_n")

# Supplementary Table S1 ------------------------------------------------------
strobe <- tribble(
  ~Item, ~Reporting_element, ~Manuscript_or_supplement_location, ~Status,
  "1", "Title and abstract identify the design, data sources, and major finding.", "Title; Abstract", "Drafted; line numbers pending final formatting",
  "2", "Scientific background and rationale.", "Introduction, paragraphs 1-3", "Drafted",
  "3", "Specific objectives and hypothesis.", "End of Introduction", "Drafted",
  "4", "Study design presented early.", "Methods: Study design and data sources; Figure 1", "Drafted",
  "5", "Setting, databases, eligibility period, and database roles.", "Methods: Study design and data sources; Supplementary Table S2", "Needs final database version/date details",
  "6", "Participants, eligibility criteria, exclusions, and cohort flow.", "Methods: Study population; Figure 1; Supplementary Table S3", "Drafted",
  "7", "Variables and definitions.", "Methods: Variables; Supplementary Table S2", "Drafted",
  "8", "Data sources and measurement rules.", "Methods: Variable harmonization; Supplementary Tables S2 and S4", "Drafted",
  "9", "Bias and missingness considerations.", "Methods: Sensitivity analyses; Supplementary Tables S4-S6", "Drafted",
  "10", "Study size.", "Results: Cohort construction; Figure 1; Supplementary Table S3", "Drafted",
  "11", "Quantitative variable handling.", "Methods: HRC construction; Supplementary Table S5", "Drafted",
  "12", "Statistical methods, confounding adjustment, subgroup and sensitivity analyses.", "Methods: Statistical analysis; Supplementary Tables S5-S8", "Drafted",
  "13", "Participant flow and missing data at each stage.", "Figure 1; Supplementary Tables S3-S4", "Drafted",
  "14", "Descriptive data and baseline characteristics.", "Table 1; Supplementary Table S4", "Drafted",
  "15", "Outcome data.", "Results; Supplementary Tables S7-S8", "Drafted",
  "16", "Main results with estimates and precision.", "Results; Figure 2; Supplementary Table S7", "Drafted",
  "17", "Other analyses, sensitivity analyses, and robustness checks.", "Results; Figure 4; Supplementary Tables S6-S8", "Drafted",
  "18", "Key results summarized with study objectives.", "Discussion, opening paragraph", "Drafted",
  "19", "Limitations.", "Discussion: Strengths and limitations", "Drafted",
  "20", "Interpretation calibrated to observational design.", "Discussion", "Drafted",
  "21", "Generalisability.", "Discussion", "Drafted",
  "22", "Funding and role of funders.", "Declarations", "Needs author input"
)
write_csv(strobe, file.path(tabdir, "supp_table_s1_strobe_mapping.csv"))
write_md_table(
  strobe,
  file.path(tabdir, "supp_table_s1_strobe_mapping.md"),
  "Supplementary Table S1. STROBE reporting map",
  "Note: This is a manuscript-specific STROBE map. Final page and line numbers should be inserted after journal formatting."
)

# Supplementary Table S2 ------------------------------------------------------
harmonization <- tribble(
  ~Domain, ~Role_in_analysis, ~`MIMIC-IV`, ~eICU, ~SICdb, ~Harmonized_rule_or_direction,
  "Unit of analysis", "Cohort unit", "Adult ICU stay; first eligible support stay per patient", "Adult ICU stay; first eligible support stay per patient", "Adult ICU case; first eligible support case per patient", "One index support episode per patient/case was retained for the primary analysis.",
  "Index time", "Time zero", "Initiation of early vasoactive/hemodynamic support", "Initiation of early vasoactive support", "Initiation of early vasoactive support", "All windows were anchored to the first qualifying early support time.",
  "Baseline window", "Expected-recovery conditioning", "0-6 h after index support", "0-6 h after index support", "0-6 h after index support", "Baseline physiology was summarized before the response window.",
  "Response window", "Observed recovery", "6-24 h after index support", "6-24 h after index support", "6-24 h after index support", "Observed recovery was measured before the 24 h landmark.",
  "Landmark", "Temporal ordering", "24 h after index support", "24 h after index support", "24 h after index support", "Primary outcomes were evaluated after the landmark.",
  "MAP", "Core HRC circulatory domain", "Mean arterial pressure, mmHg", "Mean arterial pressure, mmHg", "Mean arterial pressure, mmHg", "Higher response minus baseline indicates better pressure recovery.",
  "Vasoactive support", "Core HRC circulatory-support domain", "Norepinephrine-equivalent average dose where available", "Harmonized vasoactive burden from medication records", "Harmonized vasoactive burden / agent burden from treatment records", "Lower response burden relative to baseline indicates better support de-escalation; database-specific raw units were not directly compared.",
  "Urine output", "Core HRC renal-response domain", "Urine output normalized to mL/kg/h", "Urine output normalized to mL/kg/h", "Urine output normalized to mL/kg/h", "Higher response relative to baseline indicates better renal/perfusion recovery.",
  "Creatinine", "Secondary organ-recovery signal", "Baseline/response and 24-72 h creatinine summaries", "Baseline/response and 24-72 h creatinine summaries", "Baseline/response and 24-72 h creatinine summaries", "Creatinine worsening was used for secondary organ nonrecovery analyses, not the core HRC construct.",
  "Lactate", "Secondary metabolic signal", "Available subset", "Available subset", "Available subset", "Lactate nonclearance and hyperlactatemia were secondary endpoints due variable availability.",
  "Severity", "Expected-recovery and outcome adjustment", "First-day SOFA and SOFA domains", "Acute Physiology Score and APACHE score", "SAPS3", "Severity was harmonized by role rather than identical raw scale.",
  "Demographics", "Adjustment", "Age, sex, weight", "Age, sex, weight", "Age, sex, weight", "Used for expected-recovery and outcome adjustment where available.",
  "Hospital mortality", "Primary outcome", "Post-landmark hospital death", "Post-landmark hospital death", "Post-landmark hospital death", "Deaths before 24 h were handled in landmark sensitivity analyses.",
  "Organ nonrecovery", "Key secondary outcome", "Persistent vasoactive support, oliguria, creatinine worsening, lactate nonclearance/hyperlactatemia, RRT/ventilation modules where available", "Persistent vasoactive support, oliguria, creatinine worsening, lactate nonclearance/hyperlactatemia", "Persistent vasoactive support proxy, oliguria, creatinine worsening, lactate nonclearance/hyperlactatemia", "Domain definitions preserved physiologic directionality across databases.",
  "Fluid exposure", "Sensitivity covariate", "0-24 h net fluid balance where available", "0-24 h fluid balance where available", "0-24 h fluid exposure rate; net balance not consistently available", "Fluid adjustment tested whether HRC associations were explained by fluid exposure.",
  "CVP/MPP", "MIMIC-IV mechanistic module", "CVP and mean perfusion pressure available in a subset", "Not used for primary cross-database validation", "Not used for primary cross-database validation", "Used only for MIMIC-IV mechanistic validation because cross-database availability was insufficient.",
  "Mechanical circulatory support", "MIMIC-IV mechanistic module", "IABP, Impella, ECMO from procedure/device records", "Not primary due inconsistent device capture", "Not primary due inconsistent device capture", "MCS analyses were mechanistic, not part of the cross-database HRC definition.",
  "Clustering", "Robustness", "ICU stay/patient structure", "Hospital/ward structure", "Hospital-unit/case structure", "Cluster-aware and database-stratified analyses were used as robustness checks."
)
write_csv(harmonization, file.path(tabdir, "supp_table_s2_variable_harmonization.csv"))
write_md_table(
  harmonization,
  file.path(tabdir, "supp_table_s2_variable_harmonization.md"),
  "Supplementary Table S2. Cross-database variable harmonization dictionary",
  "Note: Harmonization prioritized physiologic directionality and timing rather than identical raw measurement scales."
)

# Supplementary Table S3 ------------------------------------------------------
flow_order <- c(
  "Adult ICU stays/cases",
  "Early vasoactive support",
  "24-h landmark eligible",
  "First support stay/case",
  "Core HRC analytic cohort"
)
flow <- read_csv("outputs/annals_main_v5/source_data/figure1_flow_counts.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    step = factor(step, levels = flow_order)
  ) %>%
  arrange(database, step) %>%
  group_by(database) %>%
  mutate(
    previous_n = lag(n),
    excluded_from_previous = if_else(is.na(previous_n), NA_real_, previous_n - n),
    retained_from_previous = if_else(is.na(previous_n), NA_real_, n / previous_n)
  ) %>%
  ungroup()
s3 <- flow %>%
  transmute(
    Database = as.character(database),
    Step = as.character(step),
    `Source metric` = metric,
    `n retained` = fmt_n(n),
    `Excluded from previous step` = ifelse(is.na(excluded_from_previous), "NA", fmt_n(excluded_from_previous)),
    `Retained from previous step` = ifelse(is.na(retained_from_previous), "NA", paste0(sprintf("%.1f", 100 * retained_from_previous), "%"))
  )
write_csv(flow, file.path(srcdir, "supp_table_s3_cohort_flow_source.csv"))
write_csv(s3, file.path(tabdir, "supp_table_s3_cohort_flow.csv"))
write_md_table(
  s3,
  file.path(tabdir, "supp_table_s3_cohort_flow.md"),
  "Supplementary Table S3. Cohort construction and exclusions by database",
  "Note: The final row is the primary analytic HRC cohort used in the manuscript. The SICdb formal HRC scoring set contained 8,275 rows; the primary merged analytic set used for final outcome/sensitivity analyses retained 8,274 rows."
)

# Supplementary Table S4 ------------------------------------------------------
formal_paths <- c(
  "MIMIC-IV" = "outputs/mimic_formal/mimic_formal_cohort.csv",
  "eICU" = "outputs/eicu_formal/eicu_formal_cohort.csv",
  "SICdb" = "outputs/sicdb_formal/sicdb_formal_cohort.csv"
)
formal <- lapply(formal_paths, read_csv, show_col_types = FALSE)
formal_den <- tibble(
  database = names(formal),
  formal_n = vapply(formal, nrow, integer(1))
)

availability_rows <- list()
for (db in names(formal)) {
  dat <- formal[[db]]
  vaso_vars <- if (db == "MIMIC-IV") {
    c("baseline_neq_avg", "response_neq_avg")
  } else {
    c("baseline_vaso_burden", "response_vaso_burden")
  }
  core_vars <- if (db == "MIMIC-IV") {
    c("delta_map", "log_neq_reduction", "log_uo_recovery", "hospital_mortality_after_landmark",
      "baseline_map_mean", "baseline_neq_avg", "baseline_uo_ml_kg_h")
  } else {
    c("delta_map", "log_vaso_burden_reduction", "log_uo_recovery", "hospital_mortality_after_landmark",
      "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h")
  }
  den <- nrow(dat)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Baseline and response MAP", complete_n = complete_n(dat, c("baseline_map_mean", "response_map_mean")), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Baseline and response vasoactive support", complete_n = complete_n(dat, vaso_vars), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Baseline and response urine output, weight-normalized", complete_n = complete_n(dat, c("baseline_uo_ml_kg_h", "response_uo_ml_kg_h")), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "All core observed HRC requirements", complete_n = complete_n(dat, core_vars), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Baseline and response creatinine", complete_n = complete_n(dat, c("baseline_creatinine_mean", "response_creatinine_mean")), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Baseline and response lactate", complete_n = complete_n(dat, c("baseline_lactate_mean", "response_lactate_mean")), denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Formal early-support cohort", domain = "Post-landmark hospital mortality", complete_n = complete_n(dat, "hospital_mortality_after_landmark"), denominator_n = den)
}

post_avail <- analysis_df %>%
  group_by(database) %>%
  summarise(
    primary_n = n(),
    organ_nonrecovery = sum(!is.na(organ_nonrecovery_24_72)),
    fluid_covariate = sum(!is.na(fluid_covariate)),
    cvp_baseline_response = sum(!is.na(baseline_cvp_mean) & !is.na(response_cvp_mean)),
    lactate_24_72 = sum(!is.na(post24_lactate_mean)),
    creatinine_24_72 = sum(!is.na(post24_creatinine_mean)),
    .groups = "drop"
  )
for (i in seq_len(nrow(post_avail))) {
  db <- as.character(post_avail$database[i])
  den <- post_avail$primary_n[i]
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Primary HRC analytic cohort", domain = "24-72 h organ nonrecovery evaluable", complete_n = post_avail$organ_nonrecovery[i], denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Primary HRC analytic cohort", domain = "0-24 h fluid covariate", complete_n = post_avail$fluid_covariate[i], denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Primary HRC analytic cohort", domain = "24-72 h lactate summary", complete_n = post_avail$lactate_24_72[i], denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Primary HRC analytic cohort", domain = "24-72 h creatinine summary", complete_n = post_avail$creatinine_24_72[i], denominator_n = den)
  availability_rows[[length(availability_rows) + 1]] <- tibble(database = db, assessment_set = "Primary HRC analytic cohort", domain = "Baseline and response CVP", complete_n = ifelse(db == "MIMIC-IV", post_avail$cvp_baseline_response[i], NA_integer_), denominator_n = ifelse(db == "MIMIC-IV", den, NA_integer_))
}
availability <- bind_rows(availability_rows) %>%
  mutate(
    database = factor(database, levels = db_levels),
    availability_percent = ifelse(is.na(complete_n) | is.na(denominator_n), NA_real_, complete_n / denominator_n)
  ) %>%
  arrange(database, assessment_set, domain)
s4 <- availability %>%
  transmute(
    Database = as.character(database),
    `Assessment set` = assessment_set,
    Domain = domain,
    `Complete n` = fmt_n(complete_n),
    `Denominator n` = fmt_n(denominator_n),
    `Complete %` = ifelse(is.na(availability_percent), "NA", paste0(sprintf("%.1f", 100 * availability_percent), "%"))
  )
write_csv(availability, file.path(srcdir, "supp_table_s4_variable_availability_source.csv"))
write_csv(s4, file.path(tabdir, "supp_table_s4_variable_availability.csv"))
write_md_table(
  s4,
  file.path(tabdir, "supp_table_s4_variable_availability.md"),
  "Supplementary Table S4. Variable availability and missingness by database",
  "Note: CVP was used only for the MIMIC-IV mechanistic module. SICdb fluid adjustment used a fluid-exposure signal rather than consistently measured net fluid balance."
)

# Supplementary Table S5 ------------------------------------------------------
model_specs <- tribble(
  ~Database, ~Component, ~Observed_response_variable, ~Expected_recovery_model, ~Imputation, ~Cross_fitting, ~HRC_contribution,
  "MIMIC-IV", "MAP recovery", "delta_map", "Linear regression with natural splines for age, weight, first-day SOFA, baseline MAP, baseline norepinephrine-equivalent burden, and baseline urine output; sex and SOFA respiration/cardiovascular/renal domains entered linearly.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "MIMIC-IV", "Vasopressor reduction", "log_neq_reduction", "Same MIMIC-IV expected-recovery model structure, outcome changed to log-scale norepinephrine-equivalent reduction.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "MIMIC-IV", "Urine output recovery", "log_uo_recovery", "Same MIMIC-IV expected-recovery model structure, outcome changed to log-scale urine-output recovery.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "eICU", "MAP recovery", "delta_map", "Linear regression with natural splines for age, weight, Acute Physiology Score, APACHE score, baseline MAP, baseline vasoactive burden, and baseline urine output; sex entered linearly.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "eICU", "Vasopressor-burden reduction", "log_vaso_burden_reduction", "Same eICU expected-recovery model structure, outcome changed to log-scale vasoactive-burden reduction.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "eICU", "Urine output recovery", "log_uo_recovery", "Same eICU expected-recovery model structure, outcome changed to log-scale urine-output recovery.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "SICdb", "MAP recovery", "delta_map", "Linear regression with natural splines for age, weight, SAPS3, baseline MAP, baseline vasoactive burden, and baseline urine output; sex and heart-surgery indicator entered linearly.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "SICdb", "Vasopressor-burden reduction", "log_vaso_burden_reduction", "Same SICdb expected-recovery model structure, outcome changed to log-scale vasoactive-burden reduction.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC.",
  "SICdb", "Urine output recovery", "log_uo_recovery", "Same SICdb expected-recovery model structure, outcome changed to log-scale urine-output recovery.", "Median imputation for expected-model covariates within the database.", "5-fold out-of-fold prediction.", "Standardized residual enters equal-weighted core HRC."
)
write_csv(model_specs, file.path(tabdir, "supp_table_s5_expected_recovery_model_specification.csv"))
write_md_table(
  model_specs,
  file.path(tabdir, "supp_table_s5_expected_recovery_model_specification.md"),
  "Supplementary Table S5. Expected-recovery model specification for HRC construction",
  "Note: The HRC framework used the same domain structure across databases, with database-specific recalibration of expected-recovery coefficients."
)

# Supplementary Table S6 ------------------------------------------------------
component_paths <- c(
  "MIMIC-IV" = "outputs/mimic_formal/mimic_hrc_formal_component_model_metrics.csv",
  "eICU" = "outputs/eicu_formal/eicu_hrc_formal_component_model_metrics.csv",
  "SICdb" = "outputs/sicdb_formal/sicdb_hrc_formal_component_model_metrics.csv"
)
component_metrics <- bind_rows(lapply(names(component_paths), function(db) {
  read_csv(component_paths[[db]], show_col_types = FALSE) %>% mutate(database = db)
})) %>%
  left_join(primary_n %>% mutate(database = as.character(database)), by = "database") %>%
  select(database, component, observed_variable, scoring_n = n, primary_analytic_n, oof_rmse, oof_r2)
s6a <- component_metrics %>%
  transmute(
    Database = database,
    Component = component,
    `Observed response variable` = observed_variable,
    `HRC scoring n` = fmt_n(scoring_n),
    `Primary analytic n` = fmt_n(primary_analytic_n),
    `Out-of-fold RMSE` = fmt_num(oof_rmse, 3),
    `Out-of-fold R2` = fmt_num(oof_r2, 3)
  )

score_paths <- c(
  "MIMIC-IV" = "outputs/mimic_formal/mimic_hrc_formal_scores.csv",
  "eICU" = "outputs/eicu_formal/eicu_hrc_formal_scores.csv",
  "SICdb" = "outputs/sicdb_formal/sicdb_hrc_formal_scores.csv"
)
resid_diag <- bind_rows(lapply(names(score_paths), function(db) {
  x <- read_csv(score_paths[[db]], show_col_types = FALSE)
  cor_mat <- suppressWarnings(cor(
    x[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")],
    method = "spearman",
    use = "pairwise.complete.obs"
  ))
  tibble(
    database = db,
    scoring_n = nrow(x),
    hrc_mean = mean(x$hrc_core_z, na.rm = TRUE),
    hrc_sd = sd(x$hrc_core_z, na.rm = TRUE),
    hrc_median_iqr = fmt_med_iqr(x$hrc_core_z, 2),
    map_vaso_spearman = cor_mat["hrc_map_residual_z", "hrc_vaso_residual_z"],
    map_uo_spearman = cor_mat["hrc_map_residual_z", "hrc_uo_residual_z"],
    vaso_uo_spearman = cor_mat["hrc_vaso_residual_z", "hrc_uo_residual_z"]
  )
})) %>%
  left_join(primary_n %>% mutate(database = as.character(database)), by = "database") %>%
  left_join(read_csv("outputs/final_methodology_extensions/simplified_hrc_correlation.csv", show_col_types = FALSE) %>%
              transmute(database, simplified_hrc_spearman = spearman_rho), by = "database")
s6b <- resid_diag %>%
  transmute(
    Database = database,
    `HRC scoring n` = fmt_n(scoring_n),
    `Primary analytic n` = fmt_n(primary_analytic_n),
    `HRC mean` = fmt_num(hrc_mean, 3),
    `HRC SD` = fmt_num(hrc_sd, 3),
    `HRC median [IQR]` = hrc_median_iqr,
    `MAP-vasopressor residual rho` = fmt_num(map_vaso_spearman, 3),
    `MAP-urine residual rho` = fmt_num(map_uo_spearman, 3),
    `Vasopressor-urine residual rho` = fmt_num(vaso_uo_spearman, 3),
    `Simplified HRC correlation rho` = fmt_num(simplified_hrc_spearman, 3)
  )
write_csv(component_metrics, file.path(srcdir, "supp_table_s6a_component_model_performance_source.csv"))
write_csv(s6a, file.path(tabdir, "supp_table_s6a_component_model_performance.csv"))
write_md_table(
  s6a,
  file.path(tabdir, "supp_table_s6a_component_model_performance.md"),
  "Supplementary Table S6A. Expected-recovery component model performance",
  "Note: Out-of-fold metrics were estimated from 5-fold cross-fitted expected-recovery models."
)
write_csv(resid_diag, file.path(srcdir, "supp_table_s6b_hrc_residual_diagnostics_source.csv"))
write_csv(s6b, file.path(tabdir, "supp_table_s6b_hrc_residual_diagnostics.csv"))
write_md_table(
  s6b,
  file.path(tabdir, "supp_table_s6b_hrc_residual_diagnostics.md"),
  "Supplementary Table S6B. HRC residual and domain-correlation diagnostics",
  "Note: Residual-domain correlations are Spearman correlations among standardized out-of-fold residuals. Simplified HRC correlation compares the primary residualized HRC with the simplified raw-domain construction."
)

# Supplementary Table S7 ------------------------------------------------------
linear_paths <- c(
  "MIMIC-IV" = "outputs/mimic_formal/mimic_hrc_formal_linear_or.csv",
  "eICU" = "outputs/eicu_formal/eicu_hrc_formal_linear_or.csv",
  "SICdb" = "outputs/sicdb_formal/sicdb_hrc_formal_linear_or.csv"
)
full_mortality <- bind_rows(lapply(names(linear_paths), function(db) {
  read_csv(linear_paths[[db]], show_col_types = FALSE) %>% mutate(database = db)
})) %>%
  select(database, everything())
write_csv(full_mortality, file.path(srcdir, "supp_table_s7_primary_mortality_full_coefficients_source.csv"))
mortality_summary_raw <- read_csv("outputs/cross_database_validation/mimic_eicu_sicdb_hrc_meta_summary.csv", show_col_types = FALSE) %>%
  filter(database %in% c(db_levels, "Fixed-effect summary")) %>%
  mutate(weight_percent = ifelse(database %in% db_levels, 100 * weight / sum(weight[database %in% db_levels], na.rm = TRUE), NA_real_))
mortality_summary <- mortality_summary_raw %>%
  transmute(
    Database = database,
    Exposure = exposure,
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value),
    `Fixed-effect weight` = ifelse(is.na(weight_percent), "NA", paste0(fmt_num(weight_percent, 1), "%"))
  )
mortality_random <- read_csv("outputs/cross_database_validation/mimic_eicu_sicdb_hrc_meta_model_stats.csv", show_col_types = FALSE) %>%
  filter(model == "DerSimonian-Laird random effects") %>%
  transmute(
    Database = "Random-effects summary",
    Exposure = "HRC per 1 SD increase",
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value),
    `Fixed-effect weight` = "NA"
  )
mortality_het <- read_csv("outputs/cross_database_validation/mimic_eicu_sicdb_hrc_meta_model_stats.csv", show_col_types = FALSE) %>%
  filter(model == "Heterogeneity") %>%
  transmute(
    Database = "Between-database heterogeneity",
    Exposure = paste0("Q=", fmt_num(q_stat, 2), "; df=", q_df, "; I2=", fmt_num(i2_percent, 1), "%"),
    `OR (95% CI)` = "NA",
    `P value` = fmt_p(p_value),
    `Fixed-effect weight` = "NA"
  )
s7 <- bind_rows(mortality_summary, mortality_random, mortality_het)
write_csv(s7, file.path(tabdir, "supp_table_s7_primary_mortality_models.csv"))
write_md_table(
  s7,
  file.path(tabdir, "supp_table_s7_primary_mortality_models.md"),
  "Supplementary Table S7. Primary HRC-mortality models and meta-analysis",
  "Note: Database-specific estimates are adjusted odds ratios per 1-SD higher HRC. Full model coefficients are provided in the source CSV."
)

# Supplementary Table S8 ------------------------------------------------------
organ_db <- read_csv("outputs/mechanism_clinical_utility/post24_endpoint_hrc_or.csv", show_col_types = FALSE) %>%
  filter(endpoint == "organ_nonrecovery_24_72") %>%
  transmute(
    Database = database,
    Endpoint = "24-72 h organ nonrecovery",
    n = fmt_n(n),
    Events = fmt_n(events),
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value)
  )
organ_meta <- read_csv("outputs/mechanism_clinical_utility/post24_endpoint_hrc_meta.csv", show_col_types = FALSE) %>%
  filter(analysis == "organ_nonrecovery_24_72") %>%
  transmute(
    Database = "Random-effects summary",
    Endpoint = "24-72 h organ nonrecovery",
    n = "NA",
    Events = "NA",
    `OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high),
    `P value` = "NA"
  )
organ_het <- read_csv("outputs/mechanism_clinical_utility/post24_endpoint_hrc_meta.csv", show_col_types = FALSE) %>%
  filter(analysis == "organ_nonrecovery_24_72") %>%
  transmute(
    Database = "Between-database heterogeneity",
    Endpoint = paste0("Q p=", fmt_p(q_p), "; I2=", fmt_num(i2_percent, 1), "%"),
    n = "NA",
    Events = "NA",
    `OR (95% CI)` = "NA",
    `P value` = "NA"
  )
s8 <- bind_rows(organ_db, organ_meta, organ_het)
write_csv(s8, file.path(tabdir, "supp_table_s8_organ_nonrecovery_models.csv"))
write_md_table(
  s8,
  file.path(tabdir, "supp_table_s8_organ_nonrecovery_models.md"),
  "Supplementary Table S8. HRC and post-landmark organ nonrecovery",
  "Note: Database-specific estimates are adjusted odds ratios per 1-SD higher HRC."
)

# Index -----------------------------------------------------------------------
supp_index <- tribble(
  ~Item, ~Title, ~Files, ~Purpose,
  "Supplementary Table S1", "STROBE reporting map", "supp_table_s1_strobe_mapping.md/csv", "Reporting-compliance map; final line numbers pending.",
  "Supplementary Table S2", "Cross-database variable harmonization dictionary", "supp_table_s2_variable_harmonization.md/csv", "Defines how physiologic variables were mapped across MIMIC-IV, eICU, and SICdb.",
  "Supplementary Table S3", "Cohort construction and exclusions by database", "supp_table_s3_cohort_flow.md/csv", "Makes Figure 1 cohort counts auditable.",
  "Supplementary Table S4", "Variable availability and missingness by database", "supp_table_s4_variable_availability.md/csv", "Shows feasibility and missingness for core and secondary domains.",
  "Supplementary Table S5", "Expected-recovery model specification", "supp_table_s5_expected_recovery_model_specification.md/csv", "Locks HRC construction model inputs and cross-fitting strategy.",
  "Supplementary Table S6A", "Expected-recovery component model performance", "supp_table_s6a_component_model_performance.md/csv", "Reports OOF RMSE and R2 for each HRC component.",
  "Supplementary Table S6B", "HRC residual and domain-correlation diagnostics", "supp_table_s6b_hrc_residual_diagnostics.md/csv", "Shows standardized HRC behavior and domain coherence.",
  "Supplementary Table S7", "Primary HRC-mortality models and meta-analysis", "supp_table_s7_primary_mortality_models.md/csv", "Supports Figure 2 and primary inference.",
  "Supplementary Table S8", "HRC and post-landmark organ nonrecovery", "supp_table_s8_organ_nonrecovery_models.md/csv", "Supports Figure 3 and the key secondary outcome."
)
write_csv(supp_index, file.path(outdir, "supplementary_tables_index_v1.csv"))
write_md_table(
  supp_index,
  file.path(outdir, "supplementary_tables_index_v1.md"),
  "Annals supplementary table package v1"
)

cat("Supplementary table package written to ", outdir, "\n", sep = "")
