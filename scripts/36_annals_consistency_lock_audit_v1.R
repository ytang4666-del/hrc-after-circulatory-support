fmt_bool <- function(x) if (isTRUE(x)) "PASS" else "FAIL"

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

main_path <- "outputs/annals_manuscript/manuscript_complete_v1_annals.md"
figure_legend_path <- "outputs/annals_manuscript/figure_legends_v1.md"
supp_path <- "outputs/annals_manuscript/supplementary_materials_complete_v1.md"
index_path <- "outputs/annals_manuscript/supplementary_tables_figures_index_v1.md"
out_path <- "outputs/annals_manuscript/full_consistency_lock_audit_v1.md"

scan_files <- c(main_path, figure_legend_path, supp_path, index_path)
missing_inputs <- scan_files[!file.exists(scan_files)]
if (length(missing_inputs) > 0) {
  stop("Missing required input files:\n", paste(missing_inputs, collapse = "\n"))
}

main <- read_text(main_path)
all_text <- paste(vapply(scan_files, read_text, character(1)), collapse = "\n")

anchor_checks <- c(
  "Primary MIMIC-IV HRC cohort n = 16,757" = "16,757 patients in MIMIC-IV",
  "Primary eICU HRC cohort n = 9,275" = "9,275 in eICU",
  "Primary SICdb HRC cohort n = 8,274" = "8,274 in SICdb",
  "Early-support MIMIC-IV denominator n = 24,176" = "24,176 adult ICU stays/cases in MIMIC-IV",
  "Early-support eICU denominator n = 22,427" = "22,427 in eICU",
  "Early-support SICdb denominator n = 11,328" = "11,328 in SICdb",
  "Mortality OR MIMIC-IV = 0.76 (0.73 to 0.80)" = "OR 0.76, 95% CI 0.73 to 0.80",
  "Mortality OR eICU = 0.75 (0.71 to 0.80)" = "OR 0.75, 95% CI 0.71 to 0.80",
  "Mortality OR SICdb = 0.67 (0.61 to 0.72)" = "OR 0.67, 95% CI 0.61 to 0.72",
  "Mortality random-effects OR = 0.73 (0.68 to 0.78)" = "0.73 (95% CI 0.68 to 0.78)",
  "Organ nonrecovery pooled OR = 0.61 (0.57 to 0.66)" = "0.61 (95% CI 0.57 to 0.66)",
  "Landmark MIMIC-IV exclusion = 3,340 of 24,176" = "3,340 of 24,176",
  "Landmark eICU exclusion = 4,174 of 22,427" = "4,174 of 22,427",
  "Landmark SICdb exclusion = 1,807 of 11,328" = "1,807 of 11,328",
  "MCS mortality OR = 0.85 (0.69 to 1.06)" = "OR 0.85, 95% CI 0.69 to 1.06",
  "MCS support/organ nonrecovery OR = 0.70 (0.52 to 0.95)" = "OR 0.70, 95% CI 0.52 to 0.95"
)

anchor_table <- data.frame(
  Check = names(anchor_checks),
  Required_text = unname(anchor_checks),
  Status = vapply(unname(anchor_checks), function(x) fmt_bool(grepl(x, main, fixed = TRUE)), character(1)),
  stringsAsFactors = FALSE
)

forbidden <- c(
  "annals_main_v4",
  "_v4",
  "center clustering",
  "formal early-support cohorts included",
  "organ-not-recovered phenotype also",
  "highest-risk physiology",
  "8,275 in SICdb",
  "HRC was constructed from harmonized circulatory and renal response domains"
)
forbidden_table <- data.frame(
  Pattern = forbidden,
  Present = vapply(forbidden, function(x) grepl(x, all_text, fixed = TRUE), logical(1)),
  stringsAsFactors = FALSE
)
forbidden_table$Status <- ifelse(forbidden_table$Present, "FAIL", "PASS")

main_lines <- readLines(main_path, warn = FALSE, encoding = "UTF-8")
placeholder_pat <- "\\[refs\\]|\\[To be completed|\\[AUTHOR TO CONFIRM|\\[initials\\]|will be provided|should be stated here|pending confirmation|Reference list to be assembled|approximately 35"
placeholder_lines <- grep(placeholder_pat, main_lines, value = TRUE)

path_tokens <- unique(unlist(regmatches(all_text, gregexpr("`outputs/[^`]+`", all_text))))
path_tokens <- gsub("`", "", path_tokens)
path_missing <- path_tokens[!file.exists(path_tokens)]

term_ledger <- data.frame(
  Canonical_term = c(
    "Hemodynamic Recovery Capacity (HRC)",
    "early circulatory support",
    "early vasoactive support",
    "core HRC cohort",
    "post-landmark hospital mortality",
    "post-landmark organ nonrecovery",
    "low-HRC phenotype",
    "mechanistic module",
    "MIMIC-IV",
    "eICU",
    "SICdb"
  ),
  Locked_usage = c(
    "Model-based residualized physiologic recovery measure, standardized within database.",
    "General clinical framing; primary cohort operationalized through early vasoactive support.",
    "Operational denominator for primary cohort construction.",
    "Patients with sufficient 0-24 h data to construct HRC and assess post-landmark outcomes.",
    "Primary outcome; avoid unqualified mortality in model descriptions.",
    "Key secondary outcome; operational database endpoint, not prospectively adjudicated organ failure.",
    "Secondary clinical translation only; continuous HRC remains primary exposure.",
    "MIMIC-IV-only support/physiology analyses; not cross-database primary inference.",
    "Use exact capitalization and hyphen.",
    "Use exact capitalization.",
    "Use exact capitalization."
  ),
  stringsAsFactors = FALSE
)

md_table <- function(df) {
  if (nrow(df) == 0) return("")
  cols <- names(df)
  rows <- apply(df, 1, function(x) paste(x, collapse = " | "))
  paste(
    paste(cols, collapse = " | "),
    paste(rep("---", length(cols)), collapse = " | "),
    paste(rows, collapse = "\n"),
    sep = "\n"
  )
}

overall_pass <- all(anchor_table$Status == "PASS") &&
  all(forbidden_table$Status == "PASS") &&
  length(path_missing) == 0

report <- c(
  "# Full Consistency Lock Audit v1",
  "",
  "Project: Hemodynamic Recovery Capacity After Early Circulatory Support Across Three Critical Care Databases",
  "",
  paste0("Generated from script: `scripts/36_annals_consistency_lock_audit_v1.R`"),
  "",
  paste0("Overall machine-check status: ", if (overall_pass) "PASS for numeric/path/version consistency" else "REVIEW REQUIRED"),
  "",
  "## One-sentence argument",
  "",
  "In ICU patients receiving early vasoactive circulatory support, HRC captured residual physiologic recovery after baseline severity and support intensity were accounted for, and lower HRC was reproducibly associated with post-landmark mortality and organ nonrecovery across MIMIC-IV, eICU, and SICdb, with MIMIC-IV mechanism modules kept secondary.",
  "",
  "## Numeric Anchor Checks",
  "",
  md_table(anchor_table),
  "",
  "## Forbidden Residual Phrases",
  "",
  md_table(forbidden_table[, c("Pattern", "Status")]),
  "",
  "## Path Check",
  "",
  paste0("- Referenced `outputs/...` paths scanned: ", length(path_tokens)),
  paste0("- Missing referenced `outputs/...` paths: ", length(path_missing)),
  if (length(path_missing) > 0) paste(path_missing, collapse = "\n") else "- No missing referenced output paths detected.",
  "",
  "## Terminology Ledger",
  "",
  md_table(term_ledger),
  "",
  "## Remaining Author/Input Items",
  "",
  if (length(placeholder_lines) == 0) {
    "- No manuscript placeholders detected by the audit pattern."
  } else {
    c(
      paste0("- Placeholder/input lines detected in main manuscript: ", length(placeholder_lines)),
      paste0("- ", placeholder_lines)
    )
  },
  "",
  "## Manual Consistency Decisions Locked in This Pass",
  "",
  "- Cohort denominators are now separated into early-support denominator, 24 h landmark/first-stay formal modeling cohort, and core HRC analytic cohort.",
  "- MIMIC-IV CVP/MPP and MCS analyses are consistently framed as mechanistic modules, not primary cross-database validation.",
  "- Continuous HRC remains the primary exposure; low-HRC remains secondary clinical translation.",
  "- HRC wording remains associative and phenotype-oriented; causal and treatment-guidance claims remain bounded.",
  "- The unverified `center clustering` robustness claim was removed from the main manuscript and component discussion draft.",
  "",
  "## Next Work Unit",
  "",
  "Resolve author-confirmation fields, create the public code repository/release, and convert the manuscript into journal submission format."
)

writeLines(report, out_path, useBytes = TRUE)
message("Generated: ", out_path)
message("Overall status: ", if (overall_pass) "PASS" else "REVIEW REQUIRED")
