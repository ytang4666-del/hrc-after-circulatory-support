#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}
fmt_range <- function(lo, hi, digits = 2) {
  paste0(sprintf(paste0("%.", digits, "f"), lo), "-", sprintf(paste0("%.", digits, "f"), hi))
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

summary_source <- read_csv(
  "outputs/mechanism_clinical_utility/decision_curve_summary.csv",
  show_col_types = FALSE
) %>%
  mutate(
    database = factor(database, levels = db_levels),
    outcome = factor(outcome, levels = c("outcome", "organ_nonrecovery_24_72"))
  ) %>%
  arrange(outcome, database)

threshold_source <- read_csv(
  "outputs/mechanism_clinical_utility/decision_curve_net_benefit.csv",
  show_col_types = FALSE
) %>%
  mutate(
    database = factor(database, levels = db_levels),
    outcome = factor(outcome, levels = c("outcome", "organ_nonrecovery_24_72"))
  ) %>%
  arrange(outcome, database, threshold)

s17 <- summary_source %>%
  transmute(
    Database = as.character(database),
    Outcome = label,
    n = fmt_n(n),
    Events = fmt_n(events),
    `Threshold range` = fmt_range(threshold_min, threshold_max, 2),
    `Mean delta net benefit` = fmt_num(mean_delta_net_benefit, 4),
    `Positive thresholds / total` = paste0(positive_delta_thresholds, "/", total_thresholds)
  )

write_csv(summary_source, file.path(srcdir, "supp_table_s17_decision_curve_summary_source.csv"))
write_csv(threshold_source, file.path(srcdir, "supp_table_s17_decision_curve_threshold_source.csv"))
write_csv(s17, file.path(tabdir, "supp_table_s17_decision_curve_summary.csv"))
write_md_table(
  s17,
  file.path(tabdir, "supp_table_s17_decision_curve_summary.md"),
  "Supplementary Table S17. Optional decision-curve summary",
  "Note: Delta net benefit compares an HRC-augmented model with the baseline clinical model across prespecified threshold ranges. This table is optional clinical-utility support and is not part of the primary inference."
)

index_path <- file.path(outdir, "supplementary_tables_index_v1.csv")
base_index <- if (file.exists(index_path)) {
  read_csv(index_path, show_col_types = FALSE)
} else {
  tibble(Item = character(), Title = character(), Files = character(), Purpose = character())
}

s17_index <- tibble(
  Item = "Supplementary Table S17",
  Title = "Optional decision-curve summary",
  Files = "supp_table_s17_decision_curve_summary.md/csv",
  Purpose = "Summarizes whether HRC-augmented models improve net benefit across clinically plausible thresholds."
)

index_all <- base_index %>%
  filter(Item != "Supplementary Table S17") %>%
  bind_rows(s17_index)

write_csv(index_all, index_path)
write_md_table(index_all, file.path(outdir, "supplementary_tables_index_v1.md"), "Annals supplementary table package v1")

cat("Supplementary Table S17 written to ", outdir, "\n", sep = "")
