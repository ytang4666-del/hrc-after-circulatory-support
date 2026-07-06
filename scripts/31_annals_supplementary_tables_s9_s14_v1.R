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
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}
fmt_pct_value <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(sprintf(paste0("%.", digits, "f"), x), "%"))
}
fmt_or <- function(or, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}
fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
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

or_table <- function(df, database_col = "database", analysis_col = "analysis") {
  df %>%
    transmute(
      Database = .data[[database_col]],
      Analysis = .data[[analysis_col]],
      n = fmt_n(n),
      Events = fmt_n(ifelse("events" %in% names(df), events, deaths)),
      `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
      `P value` = fmt_p(p_value)
    )
}

# Supplementary Table S9 ------------------------------------------------------
audit <- read_csv("outputs/methodology_stress_tests/issue2_landmark_selection_audit.csv", show_col_types = FALSE)
profile <- read_csv("outputs/landmark_sensitivity/landmark_early_event_profile.csv", show_col_types = FALSE)
pred <- read_csv("outputs/landmark_sensitivity/early_death_discharge_prediction_models.csv", show_col_types = FALSE)

profile_wide <- profile %>%
  select(database, group, n, percent_of_early_support, mortality_anytime, mean_age, mean_severity, median_index_hour) %>%
  pivot_wider(
    names_from = group,
    values_from = c(n, percent_of_early_support, mortality_anytime, mean_age, mean_severity, median_index_hour),
    names_sep = "__"
  )
pred_wide <- pred %>%
  select(database, event, model_auc) %>%
  pivot_wider(names_from = event, values_from = model_auc, names_prefix = "auc__")

s9_source <- audit %>%
  left_join(profile_wide, by = "database") %>%
  left_join(pred_wide, by = "database") %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  arrange(database)

s9a <- s9_source %>%
  transmute(
    Database = as.character(database),
    `Early-support cohort n` = fmt_n(early_support_n),
    `24-h landmark eligible n (%)` = paste0(fmt_n(landmark_eligible_n), " (", fmt_pct_value(100 * landmark_eligible_n / early_support_n), ")"),
    `Excluded before 24 h n (%)` = paste0(fmt_n(excluded_before_landmark_n), " (", fmt_pct_value(excluded_before_landmark_percent), ")"),
    `Early death n (%)` = paste0(fmt_n(n__early_death), " (", fmt_pct_value(percent_of_early_support__early_death), ")"),
    `Alive ICU discharge before 24 h n (%)` = paste0(fmt_n(n__early_discharge_alive), " (", fmt_pct_value(percent_of_early_support__early_discharge_alive), ")"),
    `Landmark eligibility AUC` = fmt_num(auc__landmark_eligible, 3),
    `Early-death AUC` = fmt_num(auc__early_death_before_landmark, 3),
    `Early-discharge AUC` = fmt_num(auc__early_icu_discharge_alive_before_landmark, 3)
  )
s9b <- profile %>%
  mutate(
    database = factor(database, levels = db_levels),
    group = factor(group, levels = c("landmark_eligible", "early_death", "early_discharge_alive"))
  ) %>%
  arrange(database, group) %>%
  transmute(
    Database = as.character(database),
    Group = recode(as.character(group),
      "landmark_eligible" = "24-h landmark eligible",
      "early_death" = "Died before 24 h",
      "early_discharge_alive" = "Alive ICU discharge before 24 h"
    ),
    n = fmt_n(n),
    `% of early-support cohort` = fmt_pct_value(percent_of_early_support),
    `Any-time mortality` = fmt_pct_value(100 * mortality_anytime),
    `Mean age` = fmt_num(mean_age, 1),
    `Mean severity score` = fmt_num(mean_severity, 1),
    `Median index hour` = fmt_num(median_index_hour, 2)
  )

write_csv(s9_source, file.path(srcdir, "supp_table_s9a_landmark_eligibility_audit_source.csv"))
write_csv(s9a, file.path(tabdir, "supp_table_s9a_landmark_eligibility_audit.csv"))
write_md_table(
  s9a,
  file.path(tabdir, "supp_table_s9a_landmark_eligibility_audit.md"),
  "Supplementary Table S9A. Landmark eligibility audit",
  "Note: HRC required a 0-24 h response window; patients who died or left the ICU before 24 h were profiled and addressed in sensitivity analyses."
)
write_csv(s9b, file.path(tabdir, "supp_table_s9b_early_event_profile.csv"))
write_md_table(
  s9b,
  file.path(tabdir, "supp_table_s9b_early_event_profile.md"),
  "Supplementary Table S9B. Early event profile before the 24 h landmark",
  "Note: Severity scores are database-specific and should be interpreted within database."
)

# Supplementary Table S10 -----------------------------------------------------
ipcw <- read_csv("outputs/landmark_sensitivity/landmark_ipcw_hrc_or.csv", show_col_types = FALSE)
early_composite <- read_csv("outputs/landmark_sensitivity/early_event_composite_worst_response_or.csv", show_col_types = FALSE)
s10_source <- bind_rows(ipcw, early_composite) %>%
  mutate(
    database = factor(database, levels = db_levels),
    analysis = factor(analysis, levels = c(
      "Landmark IPCW HRC",
      "Early deaths assigned worst observed HRC",
      "Early deaths worst + early discharges best HRC"
    ))
  ) %>%
  arrange(database, analysis)
s10 <- s10_source %>%
  transmute(
    Database = as.character(database),
    Analysis = as.character(analysis),
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    `OR per 1-SD higher HRC (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value)
  )
write_csv(s10_source, file.path(srcdir, "supp_table_s10_landmark_sensitivity_source.csv"))
write_csv(s10, file.path(tabdir, "supp_table_s10_landmark_sensitivity.csv"))
write_md_table(
  s10,
  file.path(tabdir, "supp_table_s10_landmark_sensitivity.md"),
  "Supplementary Table S10. IPCW and early-event composite sensitivity analyses",
  "Note: The composite analyses assigned early deaths the worst observed HRC value; the second composite also assigned alive early ICU discharges the best observed HRC value. These are stress tests of the landmark design, not primary estimands."
)

# Supplementary Table S11 -----------------------------------------------------
tw_or <- read_csv("outputs/time_window_sensitivity/time_window_hrc_or.csv", show_col_types = FALSE)
tw_meta <- read_csv("outputs/time_window_sensitivity/time_window_hrc_meta.csv", show_col_types = FALSE)
window_labels <- c(
  "primary_0_6_6_24" = "Primary: baseline 0-6 h, response 6-24 h",
  "early_0_6_6_12" = "Early: baseline 0-6 h, response 6-12 h",
  "late_0_12_12_24" = "Late: baseline 0-12 h, response 12-24 h"
)
s11_db <- tw_or %>%
  mutate(
    database = factor(database, levels = db_levels),
    window_label = recode(window_strategy, !!!window_labels)
  ) %>%
  arrange(window_strategy, database) %>%
  transmute(
    Window = window_label,
    Database = as.character(database),
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value),
    `Random-effects OR (95% CI)` = "NA",
    `I2` = "NA",
    row_order = match(as.character(database), db_levels)
  )
s11_meta <- tw_meta %>%
  mutate(window_label = recode(analysis, !!!window_labels)) %>%
  transmute(
    Window = window_label,
    Database = "Random-effects summary",
    n = "NA",
    Deaths = "NA",
    `OR (95% CI)` = "NA",
    `P value` = fmt_p(random_p),
    `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high),
    `I2` = paste0(fmt_num(i2_percent, 1), "%"),
    row_order = length(db_levels) + 1
  )
s11 <- bind_rows(s11_db, s11_meta) %>%
  mutate(Window = factor(Window, levels = window_labels)) %>%
  arrange(Window, row_order) %>%
  select(-row_order)
write_csv(bind_rows(
  tw_or %>% mutate(source = "database_specific"),
  tw_meta %>% mutate(database = NA_character_, window_strategy = analysis, n = NA_real_, deaths = NA_real_, beta = NA_real_, se = NA_real_, odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high, p_value = random_p, source = "meta") %>%
    select(names(tw_or), source)
), file.path(srcdir, "supp_table_s11_time_window_sensitivity_source.csv"))
write_csv(s11, file.path(tabdir, "supp_table_s11_time_window_sensitivity.csv"))
write_md_table(
  s11,
  file.path(tabdir, "supp_table_s11_time_window_sensitivity.md"),
  "Supplementary Table S11. Alternative response-window analyses",
  "Note: HRC was rebuilt separately for each response-window strategy; estimates are not simple refits of the primary HRC score."
)

# Supplementary Table S12 -----------------------------------------------------
ipw_cc <- read_csv("outputs/methodology_stress_tests/issue5_ipw_complete_case_or.csv", show_col_types = FALSE)
mice <- read_csv("outputs/high_priority_methodology/issue6_mice_covariate_imputation_or.csv", show_col_types = FALSE)
smd <- read_csv("outputs/methodology_stress_tests/issue5_complete_case_smd.csv", show_col_types = FALSE) %>%
  group_by(database) %>%
  summarise(
    max_abs_smd = max(abs(smd_complete_vs_incomplete), na.rm = TRUE),
    variable_with_max_smd = variable[which.max(abs(smd_complete_vs_incomplete))],
    .groups = "drop"
  )
s12_source <- bind_rows(ipw_cc, mice) %>%
  left_join(smd, by = "database") %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  arrange(database, analysis)
s12 <- s12_source %>%
  transmute(
    Database = as.character(database),
    Analysis = analysis,
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value),
    `Max complete-vs-incomplete abs SMD` = fmt_num(max_abs_smd, 3),
    `Variable with max SMD` = variable_with_max_smd
  )
write_csv(s12_source, file.path(srcdir, "supp_table_s12_missing_data_sensitivity_source.csv"))
write_csv(s12, file.path(tabdir, "supp_table_s12_missing_data_sensitivity.csv"))
write_md_table(
  s12,
  file.path(tabdir, "supp_table_s12_missing_data_sensitivity.md"),
  "Supplementary Table S12. Missing-data and complete-case sensitivity analyses",
  "Note: The max standardized mean difference compares HRC-complete versus HRC-incomplete patients in each formal cohort."
)

# Supplementary Table S13 -----------------------------------------------------
no_leak <- read_csv("outputs/high_priority_methodology/issue1_no_leakage_hrc_or.csv", show_col_types = FALSE)
gam <- read_csv("outputs/high_priority_methodology/issue4_gam_expected_recovery_or.csv", show_col_types = FALSE)
index_meta <- read_csv("outputs/high_priority_methodology/high_priority_extension_meta.csv", show_col_types = FALSE)
measurement <- read_csv("outputs/extra_methodology_tests/measurement_frequency_adjusted_hrc_or.csv", show_col_types = FALSE) %>%
  filter(analysis == "Adjusted for measurement frequency")
fluid <- read_csv("outputs/final_methodology_extensions/fluid_adjusted_hrc_or.csv", show_col_types = FALSE) %>%
  filter(endpoint == "outcome") %>%
  transmute(database, analysis, n, deaths = events, beta, se, odds_ratio, ci_low, ci_high, p_value)
simplified <- read_csv("outputs/final_methodology_extensions/simplified_hrc_or.csv", show_col_types = FALSE) %>%
  filter(endpoint == "outcome") %>%
  transmute(database, analysis, n, deaths = events, beta, se, odds_ratio, ci_low, ci_high, p_value)
pca <- read_csv("outputs/extra_methodology_tests/pca_weighted_hrc_or.csv", show_col_types = FALSE) %>%
  select(database, analysis, n, deaths, beta, se, odds_ratio, ci_low, ci_high, p_value, pca_variance_explained, pca_loadings)
index_adj <- read_csv("outputs/high_priority_methodology/high_priority_extension_meta.csv", show_col_types = FALSE) %>%
  filter(analysis == "Primary HRC with index-hour adjustment")

index_hour_source <- read_csv("outputs/high_priority_methodology/issue5_index_hour_adjusted_or.csv", show_col_types = FALSE)

s13a_source <- bind_rows(
  no_leak %>% mutate(variant_note = NA_character_),
  gam %>% mutate(variant_note = NA_character_),
  index_hour_source %>% mutate(variant_note = NA_character_),
  measurement %>% mutate(variant_note = NA_character_),
  fluid %>% mutate(variant_note = NA_character_),
  simplified %>% mutate(variant_note = NA_character_),
  pca %>% mutate(variant_note = paste0("PC1 variance ", fmt_num(100 * pca_variance_explained, 1), "%; ", pca_loadings)) %>%
    select(database, analysis, n, deaths, beta, se, odds_ratio, ci_low, ci_high, p_value, variant_note)
) %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  arrange(analysis, database)
s13a <- s13a_source %>%
  transmute(
    Database = as.character(database),
    Variant = analysis,
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    `OR (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value),
    Note = ifelse(is.na(variant_note), "", variant_note)
  )

alt_meta <- bind_rows(
  index_meta %>% filter(analysis %in% c("No-leakage baseline HRC", "GAM expected-recovery HRC", "Primary HRC with index-hour adjustment")) %>%
    transmute(Analysis = analysis, `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high), I2 = paste0(fmt_num(i2_percent, 1), "%")),
  read_csv("outputs/final_methodology_extensions/fluid_adjusted_hrc_meta.csv", show_col_types = FALSE) %>%
    filter(analysis == "Fluid-adjusted outcome") %>%
    transmute(Analysis = analysis, `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high), I2 = paste0(fmt_num(i2_percent, 1), "%")),
  read_csv("outputs/final_methodology_extensions/simplified_hrc_meta.csv", show_col_types = FALSE) %>%
    filter(analysis == "Simplified HRC outcome") %>%
    transmute(Analysis = analysis, `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high), I2 = paste0(fmt_num(i2_percent, 1), "%"))
)
domain_meta <- read_csv("outputs/methodology_stress_tests/issue3_domain_ablation_meta.csv", show_col_types = FALSE) %>%
  transmute(
    Analysis = paste0("Domain ablation: ", analysis),
    `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high),
    I2 = paste0(fmt_num(i2_percent, 1), "%")
  )
s13b <- bind_rows(alt_meta, domain_meta)

cutoff_meta <- read_csv("outputs/high_priority_methodology/issue7_low_hrc_cutoff_meta.csv", show_col_types = FALSE)
s13c <- cutoff_meta %>%
  transmute(
    `Low-HRC definition` = recode(analysis,
      "Low-HRC cutoff p10" = "Lowest 10%",
      "Low-HRC cutoff p20" = "Lowest 20%",
      "Low-HRC cutoff p25" = "Lowest 25% (primary clinical category)",
      "Low-HRC cutoff p33" = "Lowest 33%"
    ),
    Exposure = exposure,
    `Random-effects OR (95% CI)` = fmt_or(random_or, random_ci_low, random_ci_high),
    `I2` = paste0(fmt_num(i2_percent, 1), "%"),
    `Q-test P value` = fmt_p(q_p)
  )

write_csv(s13a_source, file.path(srcdir, "supp_table_s13a_alternative_hrc_constructions_source.csv"))
write_csv(s13a, file.path(tabdir, "supp_table_s13a_alternative_hrc_constructions.csv"))
write_md_table(
  s13a,
  file.path(tabdir, "supp_table_s13a_alternative_hrc_constructions.md"),
  "Supplementary Table S13A. Alternative HRC construction and adjustment analyses",
  "Note: Estimates are database-specific mortality associations. Sensitivity-analysis sample sizes may differ from the final primary analytic cohort because each score was rebuilt under its own data requirements."
)
write_csv(s13b, file.path(tabdir, "supp_table_s13b_alternative_hrc_meta_and_domain_ablation.csv"))
write_md_table(
  s13b,
  file.path(tabdir, "supp_table_s13b_alternative_hrc_meta_and_domain_ablation.md"),
  "Supplementary Table S13B. Meta-analytic sensitivity and domain-ablation summaries",
  "Note: Domain-ablation analyses test whether the HRC association was driven by one component alone."
)
write_csv(s13c, file.path(tabdir, "supp_table_s13c_low_hrc_cutoff_sensitivity.csv"))
write_md_table(
  s13c,
  file.path(tabdir, "supp_table_s13c_low_hrc_cutoff_sensitivity.md"),
  "Supplementary Table S13C. Low-HRC cutoff sensitivity",
  "Note: Continuous HRC remained the primary inferential exposure; categorical low-HRC definitions were secondary clinical translations."
)

# Supplementary Table S14 -----------------------------------------------------
neg <- read_csv("outputs/final_methodology_extensions/permuted_hrc_negative_control.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    endpoint = factor(endpoint, levels = c("outcome", "organ_nonrecovery_24_72"))
  ) %>%
  arrange(endpoint, database)
s14 <- neg %>%
  transmute(
    Database = as.character(database),
    Endpoint = recode(endpoint,
      "outcome" = "Post-landmark hospital mortality",
      "organ_nonrecovery_24_72" = "24-72 h organ nonrecovery"
    ),
    n = fmt_n(n),
    Events = fmt_n(events),
    `Observed OR` = fmt_num(actual_or, 2),
    `Permuted-null median OR` = fmt_num(null_median_or, 2),
    `Permuted-null 95% range` = paste0(fmt_num(null_ci_low, 2), "-", fmt_num(null_ci_high, 2)),
    `Empirical P value` = fmt_p(empirical_p),
    Permutations = fmt_n(permutations)
  )
write_csv(neg, file.path(srcdir, "supp_table_s14_negative_control_source.csv"))
write_csv(s14, file.path(tabdir, "supp_table_s14_negative_control.csv"))
write_md_table(
  s14,
  file.path(tabdir, "supp_table_s14_negative_control.md"),
  "Supplementary Table S14. Permuted-HRC negative-control analysis",
  "Note: HRC labels were randomly permuted within each database 200 times. The observed association was compared with the database-specific permuted-null distribution."
)

# Update supplementary table index -------------------------------------------
index_path <- file.path(outdir, "supplementary_tables_index_v1.csv")
base_index <- if (file.exists(index_path)) {
  read_csv(index_path, show_col_types = FALSE)
} else {
  tibble(Item = character(), Title = character(), Files = character(), Purpose = character())
}
new_index <- tribble(
  ~Item, ~Title, ~Files, ~Purpose,
  "Supplementary Table S9A", "Landmark eligibility audit", "supp_table_s9a_landmark_eligibility_audit.md/csv", "Quantifies landmark eligibility, early death, early discharge, and selection-model AUCs.",
  "Supplementary Table S9B", "Early event profile before the 24 h landmark", "supp_table_s9b_early_event_profile.md/csv", "Profiles patients excluded by early death or early ICU discharge.",
  "Supplementary Table S10", "IPCW and early-event composite sensitivity analyses", "supp_table_s10_landmark_sensitivity.md/csv", "Tests whether the HRC association depends on the 24 h landmark design.",
  "Supplementary Table S11", "Alternative response-window analyses", "supp_table_s11_time_window_sensitivity.md/csv", "Tests whether the HRC association depends on the 6-24 h response window.",
  "Supplementary Table S12", "Missing-data and complete-case sensitivity analyses", "supp_table_s12_missing_data_sensitivity.md/csv", "Reports IPW complete-case and MICE covariate-imputation sensitivity checks.",
  "Supplementary Table S13A", "Alternative HRC construction and adjustment analyses", "supp_table_s13a_alternative_hrc_constructions.md/csv", "Reports no-leakage, GAM, index-hour, measurement-frequency, fluid-adjusted, simplified, and PCA-weighted sensitivity analyses.",
  "Supplementary Table S13B", "Meta-analytic sensitivity and domain-ablation summaries", "supp_table_s13b_alternative_hrc_meta_and_domain_ablation.md/csv", "Summarizes random-effects estimates and component-ablation checks.",
  "Supplementary Table S13C", "Low-HRC cutoff sensitivity", "supp_table_s13c_low_hrc_cutoff_sensitivity.md/csv", "Shows that categorical low-HRC results are not dependent on a single percentile threshold.",
  "Supplementary Table S14", "Permuted-HRC negative-control analysis", "supp_table_s14_negative_control.md/csv", "Tests whether random assignment of HRC values reproduces the observed associations."
)
index_all <- bind_rows(base_index, new_index) %>%
  distinct(Item, .keep_all = TRUE)
write_csv(index_all, index_path)
write_md_table(index_all, file.path(outdir, "supplementary_tables_index_v1.md"), "Annals supplementary table package v1")

cat("Supplementary Tables S9-S14 written to ", outdir, "\n", sep = "")
