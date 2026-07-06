#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(svglite)
  library(ragg)
  library(grid)
})

outdir <- "outputs/annals_main_v5"
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
srcdir <- file.path(outdir, "source_data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

pal <- list(
  ink = "#1b2430",
  muted = "#667085",
  rule = "#d8dee6",
  grid = "#edf0f2",
  blue = "#2f6f9f",
  blue2 = "#8ab6d6",
  blue_light = "#eaf3fb",
  green = "#3f7f5f",
  green2 = "#97c7aa",
  green_light = "#eaf5ee",
  purple = "#6e5a8a",
  purple2 = "#b2a7c8",
  purple_light = "#f0ebf6",
  red = "#ba4a48",
  red_light = "#f7e5e3",
  gold = "#b07d2b",
  gold_light = "#f5ead2"
)

db_levels <- c("MIMIC-IV", "eICU", "SICdb")
db_cols <- c("MIMIC-IV" = pal$blue, "eICU" = pal$green, "SICdb" = pal$purple)
db_fills <- c("MIMIC-IV" = pal$blue_light, "eICU" = pal$green_light, "SICdb" = pal$purple_light)

theme_annals <- function(base_size = 7.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal$ink),
      axis.ticks = element_line(linewidth = 0.30, colour = pal$ink),
      axis.text = element_text(size = base_size - 0.8, colour = pal$ink),
      axis.title = element_text(size = base_size, colour = pal$ink),
      plot.title = element_text(size = base_size + 0.9, face = "bold", colour = pal$ink),
      plot.subtitle = element_text(size = base_size - 0.3, colour = pal$muted),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.8, colour = pal$ink),
      legend.key.size = unit(3.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.35, face = "bold", colour = pal$ink),
      panel.grid.major = element_line(linewidth = 0.16, colour = pal$grid),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_annals())

fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x, accuracy = 1) paste0(sprintf(paste0("%.", accuracy, "f"), 100 * x), "%")
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
fmt_med_iqr <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  qs <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  paste0(
    sprintf(paste0("%.", digits, "f"), qs[2]),
    " [",
    sprintf(paste0("%.", digits, "f"), qs[1]),
    "-",
    sprintf(paste0("%.", digits, "f"), qs[3]),
    "]"
  )
}

save_pub <- function(plot, filename, width_mm = 183, height_mm = 125, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 220)
  print(plot)
  dev.off()
}

write_md_table <- function(df, path, title = NULL) {
  lines <- character()
  if (!is.null(title)) lines <- c(lines, paste0("# ", title), "")
  lines <- c(lines, paste(names(df), collapse = " | "))
  lines <- c(lines, paste(rep("---", ncol(df)), collapse = " | "))
  for (i in seq_len(nrow(df))) {
    vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
    vals <- gsub("\\|", "/", vals)
    lines <- c(lines, paste(vals, collapse = " | "))
  }
  writeLines(lines, path)
}

read_summary <- function(path, database, metric_map) {
  x <- read_csv(path, show_col_types = FALSE)
  tibble(
    database = database,
    step = names(metric_map),
    metric = unname(metric_map),
    n = x$value[match(unname(metric_map), x$metric)]
  )
}

flow <- bind_rows(
  read_summary(
    "outputs/mimic_formal/mimic_formal_cohort_summary.csv",
    "MIMIC-IV",
    c(
      "Adult ICU stays/cases" = "adult_valid_icu_stays",
      "Early vasoactive support" = "early_vaso_support_stays",
      "24-h landmark eligible" = "landmark_eligible_stays",
      "First support stay/case" = "formal_first_patient_support_stays",
      "Core HRC complete" = "core_hrc_complete_rows"
    )
  ),
  read_summary(
    "outputs/eicu_formal/eicu_formal_cohort_summary.csv",
    "eICU",
    c(
      "Adult ICU stays/cases" = "adult_valid_icu_stays",
      "Early vasoactive support" = "early_vasoactive_support_stays",
      "24-h landmark eligible" = "landmark_eligible_stays",
      "First support stay/case" = "formal_first_patient_support_stays",
      "Core HRC complete" = "core_hrc_complete_rows"
    )
  ),
  read_summary(
    "outputs/sicdb_formal/sicdb_formal_cohort_summary.csv",
    "SICdb",
    c(
      "Adult ICU stays/cases" = "adult_valid_icu_cases",
      "Early vasoactive support" = "early_vasoactive_support_cases",
      "24-h landmark eligible" = "landmark_eligible_cases",
      "First support stay/case" = "formal_first_patient_support_cases",
      "Core HRC complete" = "core_hrc_complete_rows"
    )
  )
) %>%
  mutate(database = factor(database, levels = db_levels))

analysis_df <- read_csv("outputs/final_methodology_extensions/final_extension_analysis_dataset.csv", show_col_types = FALSE) %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  filter(!is.na(database), !is.na(hrc_core_z))
if (identical(Sys.getenv("HRC_EXPORT_PATIENT_LEVEL"), "1")) {
  dir.create("outputs_patient_level", showWarnings = FALSE, recursive = TRUE)
  write_csv(analysis_df, file.path("outputs_patient_level", "analysis_dataset_core_hrc.csv"))
} else {
  message("Skipping export of patient-level analysis_dataset_core_hrc.csv for public-release safety.")
}

analytic_n <- analysis_df %>%
  count(database, name = "analytic_n")
flow <- flow %>%
  left_join(analytic_n, by = "database") %>%
  mutate(
    n = if_else(metric == "core_hrc_complete_rows", as.double(analytic_n), n),
    step = if_else(metric == "core_hrc_complete_rows", "Core HRC analytic cohort", step)
  ) %>%
  select(-analytic_n)
write_csv(flow, file.path(srcdir, "figure1_flow_counts.csv"))

# Table 1 ---------------------------------------------------------------------
table1_vars <- c(
  "Core HRC cohort, n", "Age, years", "Male sex", "Body weight, kg",
  "Primary severity score", "Index time after ICU admission, h",
  "Baseline MAP, mmHg", "Response MAP, mmHg",
  "Baseline vasopressor burden", "Response vasopressor burden",
  "Baseline urine output, mL/kg/h", "Response urine output, mL/kg/h",
  "HRC, z score", "Low-HRC phenotype", "Post-landmark hospital death",
  "24-72h organ nonrecovery"
)
table1_value <- function(df, var) {
  switch(
    var,
    "Core HRC cohort, n" = fmt_n(nrow(df)),
    "Age, years" = fmt_med_iqr(df$age, 1),
    "Male sex" = paste0(fmt_n(sum(df$male == 1, na.rm = TRUE)), " (", fmt_pct(mean(df$male == 1, na.rm = TRUE)), ")"),
    "Body weight, kg" = fmt_med_iqr(df$weight_kg, 1),
    "Primary severity score" = fmt_med_iqr(df$severity_primary, 1),
    "Index time after ICU admission, h" = fmt_med_iqr(df$index_hour_from_icu, 1),
    "Baseline MAP, mmHg" = fmt_med_iqr(df$baseline_map_mean, 1),
    "Response MAP, mmHg" = fmt_med_iqr(df$response_map_mean, 1),
    "Baseline vasopressor burden" = fmt_med_iqr(df$baseline_vaso_burden, 2),
    "Response vasopressor burden" = fmt_med_iqr(df$response_vaso_burden, 2),
    "Baseline urine output, mL/kg/h" = fmt_med_iqr(df$baseline_uo_ml_kg_h, 2),
    "Response urine output, mL/kg/h" = fmt_med_iqr(df$response_uo_ml_kg_h, 2),
    "HRC, z score" = fmt_med_iqr(df$hrc_core_z, 2),
    "Low-HRC phenotype" = paste0(fmt_n(sum(df$low_hrc_q1 == 1, na.rm = TRUE)), " (", fmt_pct(mean(df$low_hrc_q1 == 1, na.rm = TRUE)), ")"),
    "Post-landmark hospital death" = paste0(fmt_n(sum(df$outcome == 1, na.rm = TRUE)), " (", fmt_pct(mean(df$outcome == 1, na.rm = TRUE)), ")"),
    "24-72h organ nonrecovery" = paste0(fmt_n(sum(df$organ_nonrecovery_24_72 == 1, na.rm = TRUE)), " (", fmt_pct(mean(df$organ_nonrecovery_24_72 == 1, na.rm = TRUE)), ")"),
    NA_character_
  )
}
table1 <- tibble(Characteristic = table1_vars)
for (db in db_levels) {
  db_df <- analysis_df %>% filter(database == db)
  table1[[db]] <- vapply(table1_vars, function(v) table1_value(db_df, v), character(1))
}
write_csv(table1, file.path(tabdir, "table1_baseline_characteristics.csv"))
write_md_table(table1, file.path(tabdir, "table1_baseline_characteristics.md"), "Table 1. Baseline characteristics of the core HRC cohorts")
cat(
  c(
    "",
    "Values are median [IQR] or n (%). HRC denotes hemodynamic recovery capacity.",
    "Primary severity score and vasopressor burden are database-specific harmonized variables and should be interpreted within, rather than directly compared between, databases.",
    "Percentages for 24-72h organ nonrecovery use patients with available post-landmark organ-recovery assessment as the denominator."
  ),
  file = file.path(tabdir, "table1_baseline_characteristics.md"),
  sep = "\n",
  append = TRUE
)

# Figure 1 --------------------------------------------------------------------
card <- function(xmin, xmax, ymin, ymax, label, fill = "white", colour = pal$rule, fontface = "plain", size = 2.5) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = colour, linewidth = 0.45),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, size = size, lineheight = 0.93,
             colour = pal$ink, fontface = fontface)
  )
}
arrow_seg <- function(x, y, xend, yend) {
  annotate("segment", x = x, y = y, xend = xend, yend = yend,
           linewidth = 0.35, colour = pal$muted,
           arrow = arrow(type = "closed", length = unit(1.5, "mm")))
}

flow_steps <- c(
  "Adult ICU stays/cases",
  "Early vasoactive support",
  "24-h landmark eligible",
  "First support stay/case",
  "Core HRC analytic cohort"
)
flow_stage_labels <- c(
  "Adult ICU",
  "Early vasoactive\nsupport",
  "24-h landmark\neligible",
  "First support\ncase",
  "Core HRC\nanalytic"
)
flow_stage <- flow %>%
  mutate(
    stage = match(step, flow_steps),
    x = c(12, 31.5, 51, 70.5, 90)[stage],
    y = c("MIMIC-IV" = 17.8, "eICU" = 11.0, "SICdb" = 4.2)[as.character(database)],
    stage_label = flow_stage_labels[stage],
    label = paste0(stage_label, "\n", fmt_n(n)),
    fill = if_else(stage == length(flow_steps), unname(db_fills[as.character(database)]), "white")
  )
flow_arrows <- flow_stage %>%
  arrange(database, stage) %>%
  group_by(database) %>%
  mutate(xend = lead(x), yend = lead(y)) %>%
  ungroup() %>%
  filter(!is.na(xend))

figure1 <- ggplot() +
  annotate("text", x = 3, y = 91, label = "A  Study time structure",
           hjust = 0, size = 3.15, fontface = "bold", colour = pal$ink) +
  annotate("rect", xmin = 6, xmax = 25, ymin = 69, ymax = 81, fill = "#f5f7fa", colour = pal$rule, linewidth = 0.35) +
  annotate("rect", xmin = 25, xmax = 56, ymin = 69, ymax = 81, fill = pal$blue_light, colour = pal$blue, linewidth = 0.40) +
  annotate("rect", xmin = 56, xmax = 91, ymin = 69, ymax = 81, fill = "#f7f7f7", colour = pal$rule, linewidth = 0.35) +
  annotate("segment", x = 6, xend = 94, y = 66, yend = 66,
           linewidth = 0.36, colour = pal$ink,
           arrow = arrow(type = "closed", length = unit(1.35, "mm"))) +
  annotate("segment", x = c(6, 25, 56), xend = c(6, 25, 56), y = 63.8, yend = 68.2,
           linewidth = 0.3, colour = pal$ink) +
  annotate("text", x = c(15.5, 40.5, 73.5), y = 75,
           label = c("Baseline\n0-6 h", "Response\n6-24 h", "Post-landmark outcomes\nfrom 24 h"),
           size = 2.40, lineheight = 0.92, colour = pal$ink) +
  annotate("text", x = c(6, 25, 56, 94), y = 60.5,
           label = c("Index support", "6 h", "24 h landmark", "Hospital discharge"),
           size = 2.22, colour = pal$ink) +
  annotate("text", x = 3, y = 51.5, label = "B  HRC construction",
           hjust = 0, size = 3.15, fontface = "bold", colour = pal$ink) +
  annotate("rect", xmin = 8, xmax = 31, ymin = 33, ymax = 45, fill = "white", colour = pal$rule, linewidth = 0.38) +
  annotate("rect", xmin = 39, xmax = 62, ymin = 33, ymax = 45, fill = "white", colour = pal$rule, linewidth = 0.38) +
  annotate("rect", xmin = 72, xmax = 92, ymin = 33, ymax = 45, fill = pal$blue_light, colour = pal$blue, linewidth = 0.42) +
  annotate("text", x = 19.5, y = 39, label = "Observed recovery\nMAP, vasopressor,\nurine output",
           size = 2.30, lineheight = 0.92, colour = pal$ink) +
  annotate("text", x = 50.5, y = 39, label = "Expected recovery\nbaseline state,\nseverity, support",
           size = 2.30, lineheight = 0.92, colour = pal$ink) +
  annotate("text", x = 82, y = 39, label = "HRC\nresidualized\nz score",
           size = 2.48, lineheight = 0.92, fontface = "bold", colour = pal$ink) +
  annotate("text", x = 35, y = 39, label = "-", size = 4.0, colour = pal$muted) +
  annotate("text", x = 67, y = 39, label = "=", size = 3.3, colour = pal$muted) +
  annotate("text", x = 3, y = 24, label = "C  Cohort flow by database",
           hjust = 0, size = 3.15, fontface = "bold", colour = pal$ink) +
  geom_segment(
    data = flow_arrows,
    aes(x = x + 6.0, xend = xend - 6.0, y = y, yend = yend),
    linewidth = 0.30, colour = pal$muted,
    arrow = arrow(type = "closed", length = unit(1.05, "mm"))
  ) +
  geom_label(
    data = flow_stage,
    aes(x = x, y = y, label = label, fill = fill),
    size = 2.22, lineheight = 0.86, label.size = 0.28,
    label.padding = unit(1.05, "mm"), family = "Arial",
    colour = pal$ink,
    show.legend = FALSE
  ) +
  geom_text(
    data = distinct(flow_stage, database, y),
    aes(x = 0.7, y = y, label = as.character(database)),
    hjust = 0, size = 2.55, fontface = "bold", colour = pal$ink
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 94), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(4, 7, 6, 7)
  )
save_pub(figure1, file.path(figdir, "figure1_annals_hrc_concept_flow_v5"), width_mm = 183, height_mm = 112)

# Figure 2 --------------------------------------------------------------------
quartile_files <- tibble(
  database = db_levels,
  path = c(
    "outputs/mimic_formal/mimic_hrc_formal_quartile_mortality.csv",
    "outputs/eicu_formal/eicu_hrc_formal_quartile_mortality.csv",
    "outputs/sicdb_formal/sicdb_hrc_formal_quartile_mortality.csv"
  )
)
quart <- bind_rows(lapply(seq_len(nrow(quartile_files)), function(i) {
  read_csv(quartile_files$path[i], show_col_types = FALSE) %>%
    mutate(database = quartile_files$database[i])
})) %>%
  mutate(
    database = factor(database, levels = db_levels),
    hrc_quartile = factor(hrc_quartile, levels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
                          labels = c("Q1\nlowest", "Q2", "Q3", "Q4\nhighest"))
  )
write_csv(quart, file.path(srcdir, "figure2_quartile_mortality.csv"))

spline_files <- tibble(
  database = db_levels,
  path = c(
    "outputs/mimic_formal/mimic_hrc_formal_spline_predictions.csv",
    "outputs/eicu_formal/eicu_hrc_formal_spline_predictions.csv",
    "outputs/sicdb_formal/sicdb_hrc_formal_spline_predictions.csv"
  )
)
spline <- bind_rows(lapply(seq_len(nrow(spline_files)), function(i) {
  read_csv(spline_files$path[i], show_col_types = FALSE) %>% mutate(database = spline_files$database[i])
})) %>%
  mutate(database = factor(database, levels = db_levels))
write_csv(spline, file.path(srcdir, "figure2_spline_source.csv"))

forest_db <- read_csv("outputs/cross_database_validation/mimic_eicu_sicdb_hrc_meta_summary.csv", show_col_types = FALSE) %>%
  filter(database %in% db_levels)
forest_re <- read_csv("outputs/cross_database_validation/mimic_eicu_sicdb_hrc_meta_model_stats.csv", show_col_types = FALSE) %>%
  filter(model == "DerSimonian-Laird random effects") %>%
  transmute(database = "Random-effects summary", exposure = "HRC per 1 SD increase", beta, se, odds_ratio, ci_low, ci_high, p_value, weight = NA_real_)
mortality_forest <- bind_rows(forest_db, forest_re) %>%
  mutate(
    database = factor(database, levels = rev(c(db_levels, "Random-effects summary"))),
    label = fmt_or(odds_ratio, ci_low, ci_high)
  )
write_csv(mortality_forest, file.path(srcdir, "figure2_mortality_forest_source.csv"))

p_q <- ggplot(quart, aes(x = hrc_quartile, y = mortality_rate, group = database, colour = database)) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = db_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.06, 0.18))) +
  labs(title = "A  Mortality gradient by HRC quartile", x = "HRC quartile", y = "Post-landmark death") +
  theme_annals() +
  theme(legend.position = "top")

p_spline <- ggplot(spline, aes(x = hrc_core_z, y = mortality_probability, colour = database, fill = database)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.65) +
  scale_colour_manual(values = db_cols) +
  scale_fill_manual(values = db_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.06, 0.18))) +
  labs(title = "B  Adjusted dose-response", x = "HRC, z score", y = "Post-landmark death") +
  theme_annals() +
  theme(legend.position = "none")

p_forest <- ggplot(mortality_forest, aes(y = database, x = odds_ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.42, colour = pal$ink) +
  geom_point(aes(fill = database == "Random-effects summary"), shape = 21, size = 2.0, colour = pal$ink) +
  geom_text(aes(x = 1.02, label = label), hjust = 0, size = 2.25, colour = pal$ink) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = pal$gold_light), guide = "none") +
  scale_x_log10(limits = c(0.54, 1.36), breaks = c(0.6, 0.75, 1, 1.25)) +
  labs(title = "C  Mortality association", x = "OR per 1-SD higher HRC", y = NULL) +
  theme_annals() +
  theme(panel.grid.major.y = element_blank())

figure2 <- (p_q | p_spline) / p_forest +
  plot_layout(heights = c(1, 0.86))
save_pub(figure2, file.path(figdir, "figure2_annals_hrc_mortality_v5"), width_mm = 183, height_mm = 116)

# Figure 3 --------------------------------------------------------------------
organ_db <- read_csv("outputs/mechanism_clinical_utility/post24_endpoint_hrc_or.csv", show_col_types = FALSE) %>%
  filter(endpoint == "organ_nonrecovery_24_72") %>%
  transmute(database, odds_ratio, ci_low, ci_high, n, events)
organ_meta <- read_csv("outputs/mechanism_clinical_utility/post24_endpoint_hrc_meta.csv", show_col_types = FALSE) %>%
  filter(analysis == "organ_nonrecovery_24_72") %>%
  transmute(database = "Random-effects summary", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high, n = NA_real_, events = NA_real_)
organ_forest <- bind_rows(organ_db, organ_meta) %>%
  mutate(database = factor(database, levels = rev(c(db_levels, "Random-effects summary"))), label = fmt_or(odds_ratio, ci_low, ci_high))
write_csv(organ_forest, file.path(srcdir, "figure3_organ_forest_source.csv"))

discord <- read_csv("outputs/mechanism_clinical_utility/discordant_pressure_organ_recovery_rates.csv", show_col_types = FALSE) %>%
  filter(discordant_group %in% c("MAP_restored_organ_recovered", "MAP_restored_organ_not_recovered")) %>%
  mutate(
    database = factor(database, levels = db_levels),
    phenotype = recode(discordant_group,
      "MAP_restored_organ_recovered" = "Organ recovered",
      "MAP_restored_organ_not_recovered" = "Organ not recovered"
    ),
    phenotype = factor(phenotype, levels = c("Organ recovered", "Organ not recovered"))
  )
discord_long <- discord %>%
  select(database, phenotype, n, mortality, post24_organ_nonrecovery_rate) %>%
  pivot_longer(c(mortality, post24_organ_nonrecovery_rate), names_to = "endpoint", values_to = "rate") %>%
  mutate(endpoint = recode(endpoint,
    "mortality" = "Hospital death",
    "post24_organ_nonrecovery_rate" = "24-72 h organ nonrecovery"
  ))
write_csv(discord_long, file.path(srcdir, "figure3_pressure_organ_discordance_source.csv"))

p_oforest <- ggplot(organ_forest, aes(y = database, x = odds_ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.42, colour = pal$ink) +
  geom_point(aes(fill = database == "Random-effects summary"), shape = 21, size = 2.0, colour = pal$ink) +
  geom_text(aes(x = 1.02, label = label), hjust = 0, size = 2.25, colour = pal$ink) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = pal$gold_light), guide = "none") +
  scale_x_log10(limits = c(0.42, 1.35), breaks = c(0.5, 0.65, 0.8, 1.0, 1.25)) +
  labs(title = "A  Association with 24-72 h organ nonrecovery", x = "OR per 1-SD higher HRC", y = NULL) +
  theme_annals() +
  theme(panel.grid.major.y = element_blank())

p_discord <- ggplot(discord_long, aes(x = database, y = rate, colour = phenotype, shape = phenotype)) +
  geom_point(position = position_dodge(width = 0.52), size = 2.15, stroke = 0.35) +
  geom_text(
    aes(label = percent(rate, accuracy = 0.1)),
    position = position_dodge(width = 0.52),
    vjust = -0.9, size = 2.18, colour = pal$ink, show.legend = FALSE
  ) +
  facet_wrap(~endpoint, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("Organ recovered" = pal$green, "Organ not recovered" = pal$red)) +
  scale_shape_manual(values = c("Organ recovered" = 16, "Organ not recovered" = 17)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.06, 0.18))) +
  labs(title = "B  Risk among patients with restored MAP", x = NULL, y = NULL) +
  theme_annals() +
  theme(legend.position = "top", axis.text.x = element_text(size = 6.8), panel.grid.major.x = element_blank())

figure3 <- p_oforest / p_discord +
  plot_layout(heights = c(0.82, 1.12))
save_pub(figure3, file.path(figdir, "figure3_annals_organ_nonrecovery_v5"), width_mm = 183, height_mm = 112)

# Figure 4 --------------------------------------------------------------------
cvp <- read_csv("outputs/mechanism_clinical_utility/mimic_cvp_mpp_joint_phenotype.csv", show_col_types = FALSE) %>%
  mutate(
    label = recode(group,
      "nonlow_HRC_no_congestion_lowMPP" = "Reference",
      "congestion_lowMPP_only" = "Congestion/\nlow MPP",
      "low_HRC_only" = "Low HRC",
      "low_HRC_plus_congestion_lowMPP" = "Low HRC +\ncongestion/low MPP"
    ),
    phenotype = factor(label, levels = c("Reference", "Congestion/\nlow MPP", "Low HRC", "Low HRC +\ncongestion/low MPP"))
  )
mcs <- read_csv("outputs/mechanism_clinical_utility/mimic_mcs_mechanism_hrc_or.csv", show_col_types = FALSE) %>%
  filter(endpoint %in% c("outcome", "organ_nonrecovery_24_72", "support_or_organ_nonrecovery_24_72")) %>%
  mutate(
    endpoint_label = recode(endpoint,
      "outcome" = "Mortality",
      "organ_nonrecovery_24_72" = "Organ nonrecovery",
      "support_or_organ_nonrecovery_24_72" = "Support/organ nonrecovery"
    ),
    endpoint_label = factor(endpoint_label, levels = rev(c("Mortality", "Organ nonrecovery", "Support/organ nonrecovery"))),
    label = fmt_or(odds_ratio, ci_low, ci_high)
  )

hk <- read_csv("outputs/high_priority_methodology/issue3_hartung_knapp_prediction_meta.csv", show_col_types = FALSE)
time_meta <- read_csv("outputs/time_window_sensitivity/time_window_hrc_meta.csv", show_col_types = FALSE)
fluid_meta <- read_csv("outputs/final_methodology_extensions/fluid_adjusted_hrc_meta.csv", show_col_types = FALSE)
simple_meta <- read_csv("outputs/final_methodology_extensions/simplified_hrc_meta.csv", show_col_types = FALSE)
hp_meta <- read_csv("outputs/high_priority_methodology/high_priority_extension_meta.csv", show_col_types = FALSE)
robust <- bind_rows(
  hk %>% filter(analysis == "Primary main HRC") %>% transmute(analysis = "Primary HRC\nHartung-Knapp", odds_ratio = random_or, ci_low = hk_ci_low, ci_high = hk_ci_high),
  time_meta %>% filter(analysis == "early_0_6_6_12") %>% transmute(analysis = "Early response\n6-12 h", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high),
  time_meta %>% filter(analysis == "late_0_12_12_24") %>% transmute(analysis = "Late response\n12-24 h", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high),
  fluid_meta %>% filter(analysis == "Fluid-adjusted outcome") %>% transmute(analysis = "Fluid-adjusted\nmortality", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high),
  simple_meta %>% filter(analysis == "Simplified HRC outcome") %>% transmute(analysis = "Simplified HRC\nmortality", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high),
  hp_meta %>% filter(analysis == "No-leakage baseline HRC") %>% transmute(analysis = "No-leakage\nbaseline HRC", odds_ratio = random_or, ci_low = random_ci_low, ci_high = random_ci_high)
) %>% mutate(analysis = factor(analysis, levels = rev(analysis)), label = fmt_or(odds_ratio, ci_low, ci_high))
neg <- read_csv("outputs/final_methodology_extensions/permuted_hrc_negative_control.csv", show_col_types = FALSE) %>%
  filter(endpoint == "outcome") %>%
  mutate(database = factor(database, levels = db_levels))
write_csv(cvp, file.path(srcdir, "figure4_cvp_mpp_source.csv"))
write_csv(mcs, file.path(srcdir, "figure4_mcs_source.csv"))
write_csv(robust, file.path(srcdir, "figure4_robustness_source.csv"))
write_csv(neg, file.path(srcdir, "figure4_negative_control_source.csv"))

p_cvp <- ggplot(cvp, aes(x = post24_organ_nonrecovery, y = mortality, size = n, fill = phenotype)) +
  geom_point(shape = 21, alpha = 0.88, colour = pal$ink, stroke = 0.25) +
  geom_text(aes(label = label), size = 2.22, vjust = -0.9, colour = pal$ink, lineheight = 0.88, show.legend = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.22, 0.62)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.035, 0.22)) +
  scale_size_continuous(range = c(2.1, 5.0), guide = "none") +
  scale_fill_manual(values = c(pal$blue_light, pal$gold_light, pal$purple_light, pal$red_light), guide = "none") +
  labs(title = "A  CVP/MPP phenotype in MIMIC-IV", x = "24-72 h organ nonrecovery", y = "Mortality") +
  coord_cartesian(clip = "off") +
  theme_annals() +
  theme(plot.margin = margin(5, 12, 4, 7))

p_mcs <- ggplot(mcs, aes(y = endpoint_label, x = odds_ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.42, colour = pal$ink) +
  geom_point(shape = 21, size = 2.1, fill = "white", colour = pal$ink) +
  geom_text(aes(x = 1.08, label = label), hjust = 0, size = 2.25, colour = pal$ink) +
  scale_x_log10(limits = c(0.43, 1.70), breaks = c(0.5, 0.7, 1.0, 1.4)) +
  labs(title = "B  MCS 0-24 h module", x = "OR per 1-SD higher HRC", y = NULL) +
  theme_annals() +
  theme(panel.grid.major.y = element_blank())

p_robust <- ggplot(robust, aes(y = analysis, x = odds_ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.13, linewidth = 0.40, colour = pal$ink) +
  geom_point(shape = 21, size = 1.95, fill = "white", colour = pal$ink) +
  geom_text(aes(x = 1.14, label = label), hjust = 0, size = 2.25, colour = pal$ink) +
  scale_x_log10(limits = c(0.45, 1.78), breaks = c(0.5, 0.75, 1, 1.5)) +
  labs(title = "C  Sensitivity analyses", x = "Pooled OR", y = NULL) +
  theme_annals() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.3))

p_neg <- ggplot(neg, aes(y = database)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = pal$muted) +
  geom_errorbarh(aes(xmin = null_ci_low, xmax = null_ci_high), height = 0.20, linewidth = 0.55, colour = pal$muted) +
  geom_point(aes(x = null_median_or), shape = 21, size = 1.95, fill = "white", colour = pal$muted) +
  geom_point(aes(x = actual_or), shape = 21, size = 2.1, fill = pal$red_light, colour = pal$red) +
  annotate("text", x = 0.66, y = 3.38, label = "actual", size = 2.25, colour = pal$red) +
  annotate("text", x = 1.03, y = 3.38, label = "permuted null", size = 2.25, colour = pal$muted) +
  scale_x_log10(limits = c(0.55, 1.25), breaks = c(0.6, 0.8, 1.0, 1.2)) +
  labs(title = "D  Permuted-HRC falsification", x = "OR", y = NULL) +
  theme_annals() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = 4), colour = pal$ink)
  )

left_col <- p_cvp / p_robust + plot_layout(heights = c(1, 1))
right_col <- p_mcs / p_neg + plot_layout(heights = c(1, 1))
figure4 <- left_col | right_col
save_pub(figure4, file.path(figdir, "figure4_annals_mechanism_robustness_v5"), width_mm = 183, height_mm = 120)

contracts <- tibble(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
  core_conclusion = c(
    "HRC is a residualized recovery construct derived from a common post-support time structure.",
    "Lower HRC identifies higher post-landmark hospital mortality across databases.",
    "Lower HRC identifies post-support organ nonrecovery and pressure-organ discordance.",
    "Mechanistic and robustness analyses support HRC as a physiologic recovery signal rather than a model artifact."
  ),
  display_role = c(
    "Conceptual framework and cohort availability",
    "Mortality gradient, dose-response, and forest plot",
    "Organ nonrecovery forest plot and MAP-restored discordance",
    "CVP/MPP module, MCS module, sensitivity analyses, and negative control"
  )
)
write_csv(contracts, file.path(outdir, "annals_v5_figure_contracts.csv"))
write_md_table(contracts, file.path(outdir, "annals_v5_figure_contracts.md"), "Annals v5 figure contracts")

cat("Annals v5 display items written to ", outdir, "\n", sep = "")
