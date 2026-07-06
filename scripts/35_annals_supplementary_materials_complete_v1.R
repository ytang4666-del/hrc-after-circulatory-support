suppressPackageStartupMessages({
  library(stringr)
})

project_root <- getwd()

out_path <- "outputs/annals_manuscript/supplementary_materials_complete_v1.md"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

sap_path <- "outputs/annals_manuscript/supplementary_statistical_analysis_plan_v1.md"
table_index_path <- "outputs/annals_supplement_v1/supplementary_tables_index_v1.md"
figure_index_path <- "outputs/annals_supplement_figures_v1/supplementary_figure_index_v1.md"
figure_legend_path <- "outputs/annals_manuscript/supplementary_figure_legends_v1.md"

table_files <- file.path(
  "outputs/annals_supplement_v1/tables",
  c(
    "supp_table_s1_strobe_mapping.md",
    "supp_table_s2_variable_harmonization.md",
    "supp_table_s3_cohort_flow.md",
    "supp_table_s4_variable_availability.md",
    "supp_table_s5_expected_recovery_model_specification.md",
    "supp_table_s6a_component_model_performance.md",
    "supp_table_s6b_hrc_residual_diagnostics.md",
    "supp_table_s7_primary_mortality_models.md",
    "supp_table_s8_organ_nonrecovery_models.md",
    "supp_table_s9a_landmark_eligibility_audit.md",
    "supp_table_s9b_early_event_profile.md",
    "supp_table_s10_landmark_sensitivity.md",
    "supp_table_s11_time_window_sensitivity.md",
    "supp_table_s12_missing_data_sensitivity.md",
    "supp_table_s13a_alternative_hrc_constructions.md",
    "supp_table_s13b_alternative_hrc_meta_and_domain_ablation.md",
    "supp_table_s13c_low_hrc_cutoff_sensitivity.md",
    "supp_table_s14_negative_control.md",
    "supp_table_s15a_cvp_mpp_joint_phenotype.md",
    "supp_table_s15b_cvp_mpp_hrc_associations.md",
    "supp_table_s16a_mimic_cs_mcs_subgroups.md",
    "supp_table_s16b_mcs_endpoint_associations.md",
    "supp_table_s16c_mcs_pressure_organ_discordance.md",
    "supp_table_s16d_mcs_bootstrap_sensitivity.md",
    "supp_table_s17_decision_curve_summary.md"
  )
)

required_files <- c(sap_path, table_index_path, table_files, figure_index_path, figure_legend_path)
absent <- required_files[!file.exists(required_files)]
if (length(absent) > 0) {
  stop("Required source files not found:\n", paste(absent, collapse = "\n"))
}

read_md <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

section_break <- function(title = NULL) {
  if (is.null(title)) {
    "\n\n---\n\n"
  } else {
    paste0("\n\n---\n\n", title, "\n\n")
  }
}

figure_legend_text <- read_md(figure_legend_path)
figure_legend_count <- str_count(figure_legend_text, "(?m)^## Supplementary Figure S[0-9]+\\.")

cover <- paste(
  "# Supplementary Materials v1",
  "",
  "Project: Hemodynamic Recovery Capacity After Early Circulatory Support Across Three Critical Care Databases",
  "",
  "Target journal: Annals of Intensive Care",
  "",
  "Status: Complete drafting package for author review and journal formatting. This file assembles the statistical protocol, supplementary tables, supplementary figure index, and supplementary figure legends.",
  "",
  "## One-sentence supplement argument",
  "",
  "These supplementary materials document cohort construction, variable harmonization, HRC construction, primary and secondary models, landmark and missing-data sensitivity analyses, mechanistic MIMIC-IV modules, and negative-control analyses supporting the multicohort HRC study.",
  "",
  "## How to read this supplement",
  "",
  "1. Supplementary Statistical Analysis Plan: locked analytic design, HRC construction, modeling principles, and sensitivity hierarchy.",
  "2. Supplementary Table Index: complete map of S1-S17 and their purpose.",
  "3. Supplementary Tables S1-S17: auditable cohort, variable, model, sensitivity, and mechanistic results.",
  "4. Supplementary Figure Index: complete map of S1-S10 figure files and source data.",
  "5. Supplementary Figure Legends S1-S10: journal-ready legends with interpretation guardrails.",
  "6. Data and Code Availability Notes: submission-facing notes to finalize with repository and database-access details.",
  "",
  "## Abbreviations",
  "",
  "CS, cardiogenic shock; CVP, central venous pressure; ECMO, extracorporeal membrane oxygenation; HRC, Hemodynamic Recovery Capacity; IABP, intra-aortic balloon pump; ICU, intensive care unit; IPCW, inverse probability weighting for 24 h landmark eligibility; MAP, mean arterial pressure; MCS, mechanical circulatory support; MPP, mean perfusion pressure; OR, odds ratio.",
  sep = "\n"
)

data_code_notes <- paste(
  "## Data and Code Availability Notes",
  "",
  "The present supplement records the analysis outputs and source data files generated for manuscript drafting. Final public code repository, database credential requirements, and accession or version details should be completed before submission.",
  "",
  "MIMIC-IV, eICU, and SICdb analyses require appropriate database access, credentialing, and compliance with each database's data-use terms. Reproducible source data for the displayed tables and figures are stored in the project output folders generated by the analysis scripts.",
  "",
  "Primary manuscript display items are stored under `outputs/annals_main_v5/`. Supplementary table source files are stored under `outputs/annals_supplement_v1/`. Supplementary figure files and source data are stored under `outputs/annals_supplement_figures_v1/`.",
  "",
  "Before submission, the final author team should add the public code repository URL, software environment file, database version numbers, and any required institutional or data-use statements.",
  sep = "\n"
)

assembly_qa <- paste(
  "## Assembly QA",
  "",
  paste0("- Supplementary table markdown files included: ", length(table_files), " files covering numbered table groups S1-S17."),
  paste0("- Supplementary figure legends included: ", figure_legend_count, " legends covering S1-S10."),
  "- Source files checked before assembly: all required markdown inputs were present.",
  paste0("- Assembly script: `scripts/35_annals_supplementary_materials_complete_v1.R`."),
  paste0("- Generated file: `", out_path, "`."),
  sep = "\n"
)

parts <- c(
  cover,
  section_break("# Supplementary Statistical Analysis Plan"),
  read_md(sap_path),
  section_break("# Supplementary Table Index"),
  read_md(table_index_path)
)

for (tf in table_files) {
  parts <- c(parts, section_break(), read_md(tf))
}

parts <- c(
  parts,
  section_break("# Supplementary Figure Index"),
  read_md(figure_index_path),
  section_break("# Supplementary Figure Legends"),
  figure_legend_text,
  section_break("# Data and Code Availability"),
  data_code_notes,
  section_break("# Assembly Quality Control"),
  assembly_qa
)

out_text <- paste(parts, collapse = "\n")
writeLines(out_text, out_path, useBytes = TRUE)

message("Generated: ", out_path)
message("Included table files: ", length(table_files))
message("Included supplementary figure legends: ", figure_legend_count)
