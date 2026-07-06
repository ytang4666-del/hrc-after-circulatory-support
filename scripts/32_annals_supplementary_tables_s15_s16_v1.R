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

fmt_n <- function(x) {
  ifelse(is.na(x), "NA", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
}
fmt_pct_value <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(sprintf(paste0("%.", digits, "f"), x), "%"))
}
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
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

endpoint_label <- function(x) {
  recode(
    x,
    "outcome" = "Post-landmark hospital mortality",
    "organ_nonrecovery_24_72" = "24-72 h organ nonrecovery",
    "support_or_organ_nonrecovery_24_72" = "24-72 h support or organ nonrecovery",
    .default = x
  )
}

# Supplementary Table S15: MIMIC-IV CVP/MPP mechanistic module ----------------
cvp_order <- c(
  "nonlow_HRC_no_congestion_lowMPP",
  "congestion_lowMPP_only",
  "low_HRC_only",
  "low_HRC_plus_congestion_lowMPP"
)
cvp_labels <- c(
  "nonlow_HRC_no_congestion_lowMPP" = "Non-low HRC without venous congestion/low MPP",
  "congestion_lowMPP_only" = "Venous congestion/low MPP only",
  "low_HRC_only" = "Low HRC only",
  "low_HRC_plus_congestion_lowMPP" = "Low HRC + venous congestion/low MPP"
)

s15a_source <- read_csv(
  "outputs/mechanism_clinical_utility/mimic_cvp_mpp_joint_phenotype.csv",
  show_col_types = FALSE
) %>%
  mutate(
    group = factor(group, levels = cvp_order),
    phenotype = recode(as.character(group), !!!cvp_labels)
  ) %>%
  arrange(group)

s15a <- s15a_source %>%
  transmute(
    Phenotype = phenotype,
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    Mortality = fmt_pct_value(100 * mortality),
    `24-72 h organ nonrecovery` = fmt_pct_value(100 * post24_organ_nonrecovery),
    `Median response CVP, mmHg` = fmt_num(median_response_cvp, 1),
    `Median response MPP, mmHg` = fmt_num(median_response_mpp, 1),
    `Median HRC` = fmt_num(median_hrc, 2)
  )

s15b_source <- read_csv(
  "outputs/mechanism_clinical_utility/mimic_cvp_mpp_hrc_or.csv",
  show_col_types = FALSE
) %>%
  mutate(endpoint = factor(endpoint, levels = c("outcome", "organ_nonrecovery_24_72"))) %>%
  arrange(endpoint)

s15b <- s15b_source %>%
  transmute(
    Database = database,
    Endpoint = endpoint_label(as.character(endpoint)),
    n = fmt_n(n),
    Events = fmt_n(events),
    `OR per 1-SD higher HRC (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value)
  )

write_csv(s15a_source, file.path(srcdir, "supp_table_s15a_cvp_mpp_joint_phenotype_source.csv"))
write_csv(s15b_source, file.path(srcdir, "supp_table_s15b_cvp_mpp_hrc_associations_source.csv"))
write_csv(s15a, file.path(tabdir, "supp_table_s15a_cvp_mpp_joint_phenotype.csv"))
write_md_table(
  s15a,
  file.path(tabdir, "supp_table_s15a_cvp_mpp_joint_phenotype.md"),
  "Supplementary Table S15A. MIMIC-IV CVP/MPP joint mechanistic phenotype",
  "Note: Venous congestion/low MPP was evaluated only in the MIMIC-IV subset with paired CVP and MPP data. This table is mechanistic support, not part of the cross-database HRC definition."
)
write_csv(s15b, file.path(tabdir, "supp_table_s15b_cvp_mpp_hrc_associations.csv"))
write_md_table(
  s15b,
  file.path(tabdir, "supp_table_s15b_cvp_mpp_hrc_associations.md"),
  "Supplementary Table S15B. HRC associations in the MIMIC-IV CVP/MPP subcohort",
  "Note: Estimates are adjusted associations within the CVP/MPP-evaluable subcohort; ORs below 1 indicate lower risk with higher HRC."
)

# Supplementary Table S16: MIMIC-IV CS/MCS mechanistic module -----------------
subgroup_order <- c(
  "Overall",
  "No CS/MCS",
  "Cardiogenic shock ICD",
  "CS or MCS",
  "Any MCS during ICU",
  "MCS during 0-24h",
  "IABP",
  "Impella",
  "ECMO"
)
subgroup_labels <- c(
  "Overall" = "Overall MIMIC-IV formal cohort",
  "No CS/MCS" = "No cardiogenic shock/MCS",
  "Cardiogenic shock ICD" = "Cardiogenic shock ICD",
  "CS or MCS" = "Cardiogenic shock or any MCS",
  "Any MCS during ICU" = "Any MCS during ICU stay",
  "MCS during 0-24h" = "MCS initiated/recorded during 0-24 h",
  "IABP" = "IABP",
  "Impella" = "Impella",
  "ECMO" = "ECMO"
)

s16a_source <- read_csv(
  "outputs/mimic_formal/mimic_hrc_formal_subgroup_summary.csv",
  show_col_types = FALSE
) %>%
  mutate(
    subgroup = factor(subgroup, levels = subgroup_order),
    subgroup_label = recode(as.character(subgroup), !!!subgroup_labels)
  ) %>%
  arrange(subgroup)

s16a <- s16a_source %>%
  transmute(
    Subgroup = subgroup_label,
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    Mortality = fmt_pct_value(100 * mortality_rate),
    `Low-HRC n` = fmt_n(low_hrc_n),
    `Low-HRC mortality` = fmt_pct_value(100 * low_hrc_mortality),
    `Non-low-HRC mortality` = fmt_pct_value(100 * nonlow_hrc_mortality),
    `Median HRC` = fmt_num(median_hrc, 2)
  )

s16b_source <- read_csv(
  "outputs/mechanism_clinical_utility/mimic_mcs_mechanism_hrc_or.csv",
  show_col_types = FALSE
) %>%
  mutate(endpoint = factor(endpoint, levels = c(
    "outcome",
    "organ_nonrecovery_24_72",
    "support_or_organ_nonrecovery_24_72"
  ))) %>%
  arrange(endpoint)

s16b <- s16b_source %>%
  transmute(
    Database = database,
    Endpoint = endpoint_label(as.character(endpoint)),
    Analysis = analysis,
    n = fmt_n(n),
    Events = fmt_n(events),
    `OR per 1-SD higher HRC (95% CI)` = fmt_or(odds_ratio, ci_low, ci_high),
    `P value` = fmt_p(p_value)
  )

discordance_order <- c(
  "MAP_restored_organ_recovered",
  "MAP_restored_organ_not_recovered",
  "MAP_not_restored"
)
discordance_labels <- c(
  "MAP_restored_organ_recovered" = "MAP restored and organ recovery",
  "MAP_restored_organ_not_recovered" = "MAP restored without organ recovery",
  "MAP_not_restored" = "MAP not restored"
)

s16c_source <- read_csv(
  "outputs/mechanism_clinical_utility/mimic_mcs_discordant_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(
    discordant_group = factor(discordant_group, levels = discordance_order),
    discordance_label = recode(as.character(discordant_group), !!!discordance_labels)
  ) %>%
  arrange(discordant_group)

s16c <- s16c_source %>%
  transmute(
    Group = discordance_label,
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    Mortality = fmt_pct_value(100 * mortality),
    `24-72 h organ nonrecovery events` = fmt_n(post24_organ_nonrecovery_events),
    `24-72 h organ nonrecovery` = fmt_pct_value(100 * post24_organ_nonrecovery_rate),
    `Median HRC` = fmt_num(median_hrc, 2)
  )

s16d_source <- read_csv(
  "outputs/methodology_stress_tests/issue8_mimic_mcs_bootstrap_or.csv",
  show_col_types = FALSE
) %>%
  mutate(subgroup = factor(subgroup, levels = subgroup_order)) %>%
  arrange(subgroup)

s16d <- s16d_source %>%
  transmute(
    Subgroup = recode(as.character(subgroup), !!!subgroup_labels),
    n = fmt_n(n),
    Deaths = fmt_n(deaths),
    `OR per 1-SD higher HRC (bootstrap 95% CI)` = fmt_or(odds_ratio, bootstrap_ci_low, bootstrap_ci_high),
    `Bootstrap replicates` = fmt_n(bootstrap_success),
    Note = note
  )

write_csv(s16a_source, file.path(srcdir, "supp_table_s16a_mimic_cs_mcs_subgroups_source.csv"))
write_csv(s16b_source, file.path(srcdir, "supp_table_s16b_mcs_endpoint_associations_source.csv"))
write_csv(s16c_source, file.path(srcdir, "supp_table_s16c_mcs_pressure_organ_discordance_source.csv"))
write_csv(s16d_source, file.path(srcdir, "supp_table_s16d_mcs_bootstrap_sensitivity_source.csv"))

write_csv(s16a, file.path(tabdir, "supp_table_s16a_mimic_cs_mcs_subgroups.csv"))
write_md_table(
  s16a,
  file.path(tabdir, "supp_table_s16a_mimic_cs_mcs_subgroups.md"),
  "Supplementary Table S16A. MIMIC-IV cardiogenic shock and MCS subgroup profile",
  "Note: MCS categories are not mutually exclusive at the device level; these analyses are mechanistic subgroups rather than the primary cross-database validation cohort."
)
write_csv(s16b, file.path(tabdir, "supp_table_s16b_mcs_endpoint_associations.csv"))
write_md_table(
  s16b,
  file.path(tabdir, "supp_table_s16b_mcs_endpoint_associations.md"),
  "Supplementary Table S16B. MCS 0-24 h endpoint associations",
  "Note: The MCS 0-24 h module is exploratory and underpowered for mortality alone; the support/organ nonrecovery endpoint captures persistent refractory physiology."
)
write_csv(s16c, file.path(tabdir, "supp_table_s16c_mcs_pressure_organ_discordance.csv"))
write_md_table(
  s16c,
  file.path(tabdir, "supp_table_s16c_mcs_pressure_organ_discordance.md"),
  "Supplementary Table S16C. Pressure-organ recovery discordance among MCS patients",
  "Note: This table separates pressure restoration from downstream organ recovery, supporting the manuscript's claim that HRC is not merely a MAP-restoration measure."
)
write_csv(s16d, file.path(tabdir, "supp_table_s16d_mcs_bootstrap_sensitivity.csv"))
write_md_table(
  s16d,
  file.path(tabdir, "supp_table_s16d_mcs_bootstrap_sensitivity.md"),
  "Supplementary Table S16D. MCS/CS subgroup bootstrap sensitivity",
  "Note: Bootstrap confidence intervals used a reduced adjustment model because device-specific MCS strata were small."
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
  "Supplementary Table S15A", "MIMIC-IV CVP/MPP joint mechanistic phenotype", "supp_table_s15a_cvp_mpp_joint_phenotype.md/csv", "Shows mortality and organ nonrecovery across low-HRC and venous congestion/low MPP phenotypes.",
  "Supplementary Table S15B", "HRC associations in the MIMIC-IV CVP/MPP subcohort", "supp_table_s15b_cvp_mpp_hrc_associations.md/csv", "Reports adjusted HRC associations in the CVP/MPP-evaluable mechanistic subcohort.",
  "Supplementary Table S16A", "MIMIC-IV cardiogenic shock and MCS subgroup profile", "supp_table_s16a_mimic_cs_mcs_subgroups.md/csv", "Describes CS/MCS/device strata and low-HRC mortality gradients.",
  "Supplementary Table S16B", "MCS 0-24 h endpoint associations", "supp_table_s16b_mcs_endpoint_associations.md/csv", "Reports HRC associations among patients receiving MCS during the first 24 h.",
  "Supplementary Table S16C", "Pressure-organ recovery discordance among MCS patients", "supp_table_s16c_mcs_pressure_organ_discordance.md/csv", "Shows that organ recovery failure persists even when MAP is restored in MCS patients.",
  "Supplementary Table S16D", "MCS/CS subgroup bootstrap sensitivity", "supp_table_s16d_mcs_bootstrap_sensitivity.md/csv", "Reports bootstrap confidence intervals for small CS/MCS/device subgroups."
)

index_all <- base_index %>%
  filter(!Item %in% new_index$Item) %>%
  bind_rows(new_index)

write_csv(index_all, index_path)
write_md_table(index_all, file.path(outdir, "supplementary_tables_index_v1.md"), "Annals supplementary table package v1")

cat("Supplementary Tables S15-S16 written to ", outdir, "\n", sep = "")
