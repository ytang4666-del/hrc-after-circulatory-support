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

outdir <- "outputs/annals_supplement_figures_v1"
figdir <- file.path(outdir, "figures")
srcdir <- file.path(outdir, "source_data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
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
  gold_light = "#f5ead2",
  grey_light = "#f6f7f9"
)

db_levels <- c("MIMIC-IV", "eICU", "SICdb")
db_cols <- c("MIMIC-IV" = pal$blue, "eICU" = pal$green, "SICdb" = pal$purple)
db_fills <- c("MIMIC-IV" = pal$blue_light, "eICU" = pal$green_light, "SICdb" = pal$purple_light)

theme_annals <- function(base_size = 7.8) {
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
      legend.text = element_text(size = base_size - 0.9, colour = pal$ink),
      legend.key.size = unit(3.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.3, face = "bold", colour = pal$ink),
      panel.grid.major = element_line(linewidth = 0.16, colour = pal$grid),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_annals())

fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x, accuracy = 1) paste0(sprintf(paste0("%.", accuracy, "f"), 100 * x), "%")
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)

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

lab_window <- function(x) {
  recode(
    x,
    "primary_0_6_6_24" = "0-6 / 6-24 h",
    "early_0_6_6_12" = "0-6 / 6-12 h",
    "late_0_12_12_24" = "0-12 / 12-24 h",
    .default = x
  )
}

clean_analysis <- function(x) {
  recode(
    x,
    "Landmark IPCW HRC" = "24-h landmark\nIPCW",
    "Early deaths assigned worst observed HRC" = "Early deaths\nworst HRC",
    "Early deaths worst + early discharges best HRC" = "Early deaths worst;\nearly discharges best",
    "IPW complete-case HRC" = "IPW complete\ncase",
    "IPW complete-case sensitivity" = "IPW complete\ncase",
    "MICE covariate imputation HRC" = "MICE covariate\nimputation",
    "MICE covariate-imputation sensitivity" = "MICE covariate\nimputation",
    "Adjusted for measurement frequency" = "Measurement-frequency\nadjusted",
    "GAM expected-recovery HRC" = "GAM expected\nrecovery",
    "HRC adjusted for fluid balance/exposure" = "Fluid-adjusted\nHRC",
    "No-leakage baseline HRC" = "No-leakage\nbaseline HRC",
    "PCA-weighted HRC" = "PCA-weighted\nHRC",
    "Primary HRC with index-hour adjustment" = "Index-hour\nadjusted",
    "Simplified raw-domain HRC" = "Simplified\nraw-domain HRC",
    .default = x
  )
}

col_or_na <- function(dat, nm) {
  if (nm %in% names(dat)) dat[[nm]] else rep(NA_real_, nrow(dat))
}

read_score <- function(path, database) {
  dat <- read_csv(path, show_col_types = FALSE)
  tibble(
    database = database,
    hrc_core_z = dat$hrc_core_z,
    hrc_map_residual_z = dat$hrc_map_residual_z,
    hrc_vaso_residual_z = dat$hrc_vaso_residual_z,
    hrc_uo_residual_z = dat$hrc_uo_residual_z,
    delta_map = dat$delta_map,
    expected_delta_map = dat$expected_delta_map,
    log_vaso_reduction = if ("log_neq_reduction" %in% names(dat)) dat$log_neq_reduction else dat$log_vaso_burden_reduction,
    expected_log_vaso_reduction = if ("expected_log_neq_reduction" %in% names(dat)) dat$expected_log_neq_reduction else dat$expected_log_vaso_burden_reduction,
    log_uo_recovery = dat$log_uo_recovery,
    expected_log_uo_recovery = dat$expected_log_uo_recovery
  )
}

scores <- bind_rows(
  read_score("outputs/mimic_formal/mimic_hrc_formal_scores.csv", "MIMIC-IV"),
  read_score("outputs/eicu_formal/eicu_hrc_formal_scores.csv", "eICU"),
  read_score("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", "SICdb")
) %>%
  mutate(database = factor(database, levels = db_levels))

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

# Supplementary Figure S1 -----------------------------------------------------
flow <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s3_cohort_flow_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    step = factor(step, levels = c(
      "Adult ICU stays/cases",
      "Early vasoactive support",
      "24-h landmark eligible",
      "First support stay/case",
      "Core HRC analytic cohort"
    )),
    step_short = recode(
      as.character(step),
      "Adult ICU stays/cases" = "Adult ICU\nstays/cases",
      "Early vasoactive support" = "Early\nsupport",
      "24-h landmark eligible" = "24-h\neligible",
      "First support stay/case" = "First support\ncase",
      "Core HRC analytic cohort" = "Core HRC\nanalytic"
    )
  ) %>%
  mutate(step_short = factor(step_short, levels = c(
    "Adult ICU\nstays/cases",
    "Early\nsupport",
    "24-h\neligible",
    "First support\ncase",
    "Core HRC\nanalytic"
  ))) %>%
  group_by(database) %>%
  arrange(step, .by_group = TRUE) %>%
  mutate(
    previous_n = lag(n),
    excluded_n = if_else(is.na(previous_n), NA_real_, previous_n - n),
    retained_percent = if_else(is.na(previous_n), NA_real_, n / previous_n)
  ) %>%
  ungroup()
write_csv(flow, file.path(srcdir, "supp_figure_s1_detailed_cohort_flow_source.csv"))

p_s1a <- ggplot(flow, aes(n, step_short, group = 1, colour = database)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.9) +
  geom_text(aes(label = fmt_n(n)), hjust = -0.12, size = 2.20, show.legend = FALSE) +
  facet_wrap(~database, nrow = 1) +
  scale_colour_manual(values = db_cols) +
  scale_x_log10(
    labels = label_comma(),
    breaks = c(8000, 15000, 30000, 60000, 120000, 220000),
    limits = c(7000, 320000),
    expand = expansion(mult = c(0.02, 0.16))
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Retained cohort size", x = "Patients/stays, log scale", y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.6),
    legend.position = "none",
    plot.margin = margin(4, 20, 4, 5)
  )

retention <- flow %>%
  filter(!is.na(retained_percent)) %>%
  mutate(retained_label = paste0(fmt_n(n), "\n", fmt_pct(retained_percent)))
p_s1b <- ggplot(retention, aes(step_short, database, fill = retained_percent)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = retained_label), size = 2.15, lineheight = 0.9, colour = pal$ink) +
  scale_fill_gradient(low = pal$red_light, high = pal$blue, labels = percent_format(accuracy = 1), limits = c(0.35, 1)) +
  labs(title = "Retention from previous step", x = NULL, y = NULL, fill = "Retained") +
  theme(axis.text.x = element_text(size = 6.7), legend.position = "right")

fig_s1 <- (p_s1a / p_s1b) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = pal$ink))
save_pub(fig_s1, file.path(figdir, "supp_figure_s1_detailed_cohort_flow_v1"), 183, 135)

# Supplementary Figure S2 -----------------------------------------------------
component_long <- scores %>%
  select(database, hrc_map_residual_z, hrc_vaso_residual_z, hrc_uo_residual_z) %>%
  pivot_longer(
    cols = starts_with("hrc_"),
    names_to = "component",
    values_to = "residual_z"
  ) %>%
  mutate(component = recode(component,
    "hrc_map_residual_z" = "MAP recovery",
    "hrc_vaso_residual_z" = "Vasopressor de-escalation",
    "hrc_uo_residual_z" = "Urine output recovery"
  )) %>%
  mutate(component = factor(component, levels = c("MAP recovery", "Vasopressor de-escalation", "Urine output recovery")))
write_csv(component_long, file.path(srcdir, "supp_figure_s2_component_residuals_source.csv"))

cor_diag <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s6b_hrc_residual_diagnostics_source.csv", show_col_types = FALSE) %>%
  mutate(database = factor(database, levels = db_levels)) %>%
  select(database, map_vaso_spearman, map_uo_spearman, vaso_uo_spearman) %>%
  pivot_longer(-database, names_to = "pair", values_to = "rho") %>%
  mutate(pair = recode(pair,
    "map_vaso_spearman" = "MAP vs vaso",
    "map_uo_spearman" = "MAP vs urine",
    "vaso_uo_spearman" = "Vaso vs urine"
  )) %>%
  mutate(pair = factor(pair, levels = c("MAP vs vaso", "MAP vs urine", "Vaso vs urine")))
write_csv(cor_diag, file.path(srcdir, "supp_figure_s2_component_correlations_source.csv"))

p_s2a <- ggplot(component_long, aes(residual_z, fill = database, colour = database)) +
  geom_density(alpha = 0.22, linewidth = 0.45, adjust = 1.1) +
  facet_wrap(~component, nrow = 1) +
  scale_fill_manual(values = db_fills) +
  scale_colour_manual(values = db_cols) +
  coord_cartesian(xlim = c(-3.2, 3.2)) +
  labs(title = "Standardized residual components", x = "Component residual, z score", y = "Density") +
  theme(legend.position = "top")

p_s2b <- ggplot(cor_diag, aes(pair, database, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5, colour = pal$ink) +
  scale_fill_gradient2(low = pal$red_light, mid = "white", high = pal$blue, midpoint = 0, limits = c(-0.15, 0.45)) +
  labs(title = "Residual-domain coherence", x = NULL, y = NULL, fill = "Spearman rho") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

fig_s2 <- (p_s2a / p_s2b) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = pal$ink))
save_pub(fig_s2, file.path(figdir, "supp_figure_s2_hrc_component_distributions_v1"), 183, 130)

# Supplementary Figure S3 -----------------------------------------------------
hrc_summary <- scores %>%
  group_by(database) %>%
  summarise(
    n = n(),
    mean = mean(hrc_core_z, na.rm = TRUE),
    sd = sd(hrc_core_z, na.rm = TRUE),
    q25 = quantile(hrc_core_z, 0.25, na.rm = TRUE),
    median = median(hrc_core_z, na.rm = TRUE),
    q75 = quantile(hrc_core_z, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(hrc_summary, file.path(srcdir, "supp_figure_s3_hrc_distribution_summary_source.csv"))

p_s3a <- ggplot(scores, aes(hrc_core_z, fill = database, colour = database)) +
  geom_density(alpha = 0.22, linewidth = 0.45, adjust = 1.1) +
  scale_fill_manual(values = db_fills) +
  scale_colour_manual(values = db_cols) +
  coord_cartesian(xlim = c(-3.2, 3.2)) +
  labs(title = "Within-database HRC distributions", x = "HRC, z score", y = "Density") +
  theme(legend.position = "top")

p_s3b <- ggplot(scores, aes(database, hrc_core_z, fill = database)) +
  geom_violin(width = 0.7, linewidth = 0.25, colour = pal$rule, alpha = 0.8, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = pal$muted) +
  scale_fill_manual(values = db_fills) +
  coord_cartesian(ylim = c(-3.2, 3.2)) +
  labs(title = "Median and spread after standardization", x = NULL, y = "HRC, z score") +
  theme(legend.position = "none")

fig_s3 <- (p_s3a | p_s3b) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = pal$ink))
save_pub(fig_s3, file.path(figdir, "supp_figure_s3_hrc_distribution_v1"), 183, 105)

# Supplementary Figure S4 -----------------------------------------------------
scores_with_row <- scores %>%
  group_by(database) %>%
  mutate(source_row = row_number()) %>%
  ungroup()

expected_long <- bind_rows(
  scores_with_row %>%
    transmute(database, source_row, component = "MAP recovery", observed = delta_map, expected = expected_delta_map),
  scores_with_row %>%
    transmute(database, source_row, component = "Vasopressor de-escalation", observed = log_vaso_reduction, expected = expected_log_vaso_reduction),
  scores_with_row %>%
    transmute(database, source_row, component = "Urine output recovery", observed = log_uo_recovery, expected = expected_log_uo_recovery)
) %>%
  mutate(component = factor(component, levels = c("MAP recovery", "Vasopressor de-escalation", "Urine output recovery"))) %>%
  filter(is.finite(observed), is.finite(expected)) %>%
  group_by(database, component) %>%
  mutate(
    observed_z = as.numeric(scale(observed)),
    expected_z = as.numeric(scale(expected))
  ) %>%
  ungroup() %>%
  filter(is.finite(observed_z), is.finite(expected_z))

set.seed(20260705)
expected_sample <- expected_long %>%
  group_by(database, component) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), 1200))) %>%
  ungroup()
write_csv(expected_sample, file.path(srcdir, "supp_figure_s4_expected_observed_recovery_sample_source.csv"))

p_s4 <- ggplot(expected_sample, aes(expected_z, observed_z, colour = database)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", colour = pal$muted) +
  geom_point(alpha = 0.20, size = 0.42) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.42) +
  facet_grid(component ~ database, scales = "free") +
  scale_colour_manual(values = db_cols) +
  coord_cartesian(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
  labs(title = "Observed recovery versus expected recovery", x = "Expected recovery, within-domain z score", y = "Observed recovery, within-domain z score") +
  theme(legend.position = "none")
save_pub(p_s4, file.path(figdir, "supp_figure_s4_expected_observed_recovery_v1"), 183, 145)

# Supplementary Figure S5 -----------------------------------------------------
spline <- read_csv("outputs/annals_main_v5/source_data/figure2_spline_source.csv", show_col_types = FALSE) %>%
  mutate(database = factor(database, levels = db_levels))
write_csv(spline, file.path(srcdir, "supp_figure_s5_full_mortality_splines_source.csv"))

p_s5 <- ggplot(spline, aes(hrc_core_z, mortality_probability, colour = database, fill = database)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.22, linewidth = 0) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~database, nrow = 1) +
  scale_colour_manual(values = db_cols) +
  scale_fill_manual(values = db_fills) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(-2.8, 2.8)) +
  labs(title = "Full HRC-mortality spline curves", x = "HRC, z score", y = "Adjusted mortality probability") +
  theme(legend.position = "none")
save_pub(p_s5, file.path(figdir, "supp_figure_s5_full_mortality_splines_v1"), 183, 85)

# Supplementary Figure S6 -----------------------------------------------------
organ_component_long <- analysis_df %>%
  group_by(database) %>%
  mutate(hrc_quartile = ntile(hrc_core_z, 4)) %>%
  ungroup() %>%
  select(
    database, hrc_quartile,
    persistent_vaso_24_72,
    oliguria_24_72,
    creatinine_worsening_24_72,
    lactate_nonclearance_24_72,
    organ_nonrecovery_24_72
  ) %>%
  pivot_longer(
    cols = -c(database, hrc_quartile),
    names_to = "component",
    values_to = "event"
  ) %>%
  filter(!is.na(event)) %>%
  group_by(database, hrc_quartile, component) %>%
  summarise(n = n(), events = sum(event == 1, na.rm = TRUE), rate = mean(event == 1, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    hrc_quartile = factor(hrc_quartile, levels = 1:4, labels = c("Q1\nlowest", "Q2", "Q3", "Q4\nhighest")),
    component = recode(component,
      "persistent_vaso_24_72" = "Persistent vasoactive support",
      "oliguria_24_72" = "Oliguria",
      "creatinine_worsening_24_72" = "Creatinine worsening",
      "lactate_nonclearance_24_72" = "Lactate nonclearance",
      "organ_nonrecovery_24_72" = "Composite organ nonrecovery"
    )
  )
write_csv(organ_component_long, file.path(srcdir, "supp_figure_s6_organ_component_rates_source.csv"))

p_s6 <- ggplot(organ_component_long, aes(hrc_quartile, component, fill = rate)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = paste0(fmt_pct(rate), "\n", events, "/", n)), size = 2.20, lineheight = 0.87, colour = pal$ink) +
  facet_wrap(~database, nrow = 1) +
  scale_fill_gradient(low = "white", high = pal$red, labels = percent_format(accuracy = 1)) +
  labs(title = "Organ nonrecovery components by HRC quartile", x = "Within-database HRC quartile", y = NULL, fill = "Rate") +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "right")
save_pub(p_s6, file.path(figdir, "supp_figure_s6_organ_nonrecovery_components_v1"), 183, 105)

# Supplementary Figure S7 -----------------------------------------------------
tw <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s11_time_window_sensitivity_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = if_else(is.na(database), "Random-effects", database),
    database = factor(database, levels = c("MIMIC-IV", "eICU", "SICdb", "Random-effects")),
    window = factor(lab_window(window_strategy), levels = c("0-6 / 6-24 h", "0-6 / 6-12 h", "0-12 / 12-24 h"))
  )
write_csv(tw, file.path(srcdir, "supp_figure_s7_time_window_sensitivity_source.csv"))

p_s7 <- ggplot(tw, aes(odds_ratio, database, colour = database)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, linewidth = 0.45) +
  geom_point(aes(shape = database == "Random-effects"), size = 2.0, fill = "white") +
  facet_wrap(~window, nrow = 1) +
  scale_colour_manual(values = c(db_cols, "Random-effects" = pal$ink)) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18)) +
  scale_x_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1.0), labels = number_format(accuracy = 0.01), limits = c(0.58, 1.05)) +
  labs(title = "Alternative response-window sensitivity", x = "OR per 1-SD higher HRC", y = NULL) +
  theme(legend.position = "none")
save_pub(p_s7, file.path(figdir, "supp_figure_s7_time_window_sensitivity_v1"), 183, 82)

# Supplementary Figure S8 -----------------------------------------------------
early_profile <- read_csv("outputs/landmark_sensitivity/landmark_early_event_profile.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    group = factor(group, levels = c("landmark_eligible", "early_death", "early_discharge_alive")),
    group_label = recode(as.character(group),
      "landmark_eligible" = "24-h landmark eligible",
      "early_death" = "Died before 24 h",
      "early_discharge_alive" = "Alive ICU discharge before 24 h"
    )
  )
landmark_sens <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s10_landmark_sensitivity_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    analysis_label = factor(clean_analysis(analysis), levels = c(
      "24-h landmark\nIPCW",
      "Early deaths\nworst HRC",
      "Early deaths worst;\nearly discharges best"
    ))
  )
write_csv(early_profile, file.path(srcdir, "supp_figure_s8_early_event_profile_source.csv"))
write_csv(landmark_sens, file.path(srcdir, "supp_figure_s8_landmark_sensitivity_source.csv"))

p_s8a <- ggplot(early_profile, aes(database, percent_of_early_support / 100, fill = group_label)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(percent_of_early_support >= 5, paste0(fmt_pct(percent_of_early_support / 100), "\n", fmt_n(n)), "")),
            position = position_stack(vjust = 0.5), size = 2.20, lineheight = 0.88, colour = pal$ink) +
  scale_fill_manual(values = c(
    "24-h landmark eligible" = pal$blue_light,
    "Died before 24 h" = pal$red_light,
    "Alive ICU discharge before 24 h" = pal$gold_light
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Pre-landmark disposition", x = NULL, y = "% of early-support cohort") +
  theme(legend.position = "top")

p_s8b <- ggplot(landmark_sens, aes(odds_ratio, analysis_label, colour = database)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, position = position_dodge(width = 0.55), linewidth = 0.42) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_colour_manual(values = db_cols) +
  scale_x_log10(breaks = c(0.35, 0.45, 0.6, 0.8, 1), labels = number_format(accuracy = 0.01), limits = c(0.33, 1.02)) +
  labs(title = "Landmark sensitivity", x = "OR per 1-SD higher HRC", y = NULL) +
  theme(legend.position = "top")

fig_s8 <- (p_s8a | p_s8b) +
  plot_layout(widths = c(0.95, 1.35)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = pal$ink))
save_pub(fig_s8, file.path(figdir, "supp_figure_s8_landmark_early_event_sensitivity_v1"), 183, 110)

# Supplementary Figure S9 -----------------------------------------------------
missing_sens <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s12_missing_data_sensitivity_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    analysis_label = factor(clean_analysis(analysis), levels = c("IPW complete\ncase", "MICE covariate\nimputation"))
  )
alt_sens <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s13a_alternative_hrc_constructions_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    analysis_label = factor(clean_analysis(analysis), levels = rev(c(
      "No-leakage\nbaseline HRC",
      "GAM expected\nrecovery",
      "Index-hour\nadjusted",
      "Measurement-frequency\nadjusted",
      "Fluid-adjusted\nHRC",
      "Simplified\nraw-domain HRC",
      "PCA-weighted\nHRC"
    )))
  )
write_csv(missing_sens, file.path(srcdir, "supp_figure_s9_missing_data_sensitivity_source.csv"))
write_csv(alt_sens, file.path(srcdir, "supp_figure_s9_alternative_hrc_sensitivity_source.csv"))

p_s9a <- ggplot(missing_sens, aes(odds_ratio, analysis_label, colour = database)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14, position = position_dodge(width = 0.52), linewidth = 0.42) +
  geom_point(position = position_dodge(width = 0.52), size = 1.8) +
  scale_colour_manual(values = db_cols) +
  scale_x_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1), labels = number_format(accuracy = 0.01), limits = c(0.58, 1.02)) +
  labs(title = "Missing-data sensitivity", x = "OR per 1-SD higher HRC", y = NULL) +
  theme(legend.position = "top")

p_s9b <- ggplot(alt_sens, aes(odds_ratio, analysis_label, colour = database)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, position = position_dodge(width = 0.55), linewidth = 0.38) +
  geom_point(position = position_dodge(width = 0.55), size = 1.55) +
  scale_colour_manual(values = db_cols) +
  scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = number_format(accuracy = 0.01), limits = c(0.48, 1.02)) +
  labs(title = "Alternative HRC construction and adjustment", x = "OR per 1-SD higher HRC", y = NULL) +
  theme(legend.position = "top", axis.text.y = element_text(size = 6.3))

fig_s9 <- (p_s9a / p_s9b) +
  plot_layout(heights = c(0.65, 1.45)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = pal$ink))
save_pub(fig_s9, file.path(figdir, "supp_figure_s9_missing_alternative_hrc_sensitivity_v1"), 183, 150)

# Supplementary Figure S10 ----------------------------------------------------
neg <- read_csv("outputs/annals_supplement_v1/source_data/supp_table_s14_negative_control_source.csv", show_col_types = FALSE) %>%
  mutate(
    database = factor(database, levels = db_levels),
    endpoint_label = factor(recode(endpoint,
      "outcome" = "Hospital mortality",
      "organ_nonrecovery_24_72" = "24-72 h organ nonrecovery"
    ), levels = c("Hospital mortality", "24-72 h organ nonrecovery"))
  )
write_csv(neg, file.path(srcdir, "supp_figure_s10_negative_control_source.csv"))

p_s10 <- ggplot(neg, aes(database, null_median_or)) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", colour = pal$muted) +
  geom_errorbar(aes(ymin = null_ci_low, ymax = null_ci_high), width = 0.18, linewidth = 0.45, colour = pal$muted) +
  geom_point(size = 2.1, colour = pal$muted, fill = "white", shape = 21) +
  geom_point(aes(y = actual_or, colour = database), size = 2.3) +
  geom_text(aes(y = actual_or, label = sprintf("%.2f", actual_or)), vjust = 1.65, size = 2.15, colour = pal$ink) +
  facet_wrap(~endpoint_label, nrow = 1) +
  scale_colour_manual(values = db_cols) +
  scale_y_continuous(limits = c(0.48, 1.08), breaks = c(0.5, 0.7, 0.9, 1.0)) +
  labs(title = "Permuted-HRC negative-control analysis", x = NULL, y = "OR; grey interval is permuted-null 95% range") +
  theme(legend.position = "none")
save_pub(p_s10, file.path(figdir, "supp_figure_s10_negative_control_v1"), 183, 90)

# Index and QA notes ----------------------------------------------------------
fig_index <- tribble(
  ~Item, ~Title, ~Files, ~Purpose,
  "Supplementary Figure S1", "Detailed cohort flow by database", "supp_figure_s1_detailed_cohort_flow_v1.png/svg/pdf/tiff", "Audits cohort attrition and 24-h landmark retention.",
  "Supplementary Figure S2", "HRC component distributions and residual-domain coherence", "supp_figure_s2_hrc_component_distributions_v1.png/svg/pdf/tiff", "Shows standardized residual components and correlations between recovery domains.",
  "Supplementary Figure S3", "HRC distribution across databases", "supp_figure_s3_hrc_distribution_v1.png/svg/pdf/tiff", "Shows within-database standardization of the final HRC score.",
  "Supplementary Figure S4", "Observed versus expected recovery", "supp_figure_s4_expected_observed_recovery_v1.png/svg/pdf/tiff", "Visualizes residualized HRC construction across component domains.",
  "Supplementary Figure S5", "Full mortality spline curves by database", "supp_figure_s5_full_mortality_splines_v1.png/svg/pdf/tiff", "Expands the HRC-mortality dose-response display.",
  "Supplementary Figure S6", "Organ nonrecovery components by HRC quartile", "supp_figure_s6_organ_nonrecovery_components_v1.png/svg/pdf/tiff", "Shows that lower HRC tracks multiple downstream organ nonrecovery components.",
  "Supplementary Figure S7", "Alternative response-window sensitivity", "supp_figure_s7_time_window_sensitivity_v1.png/svg/pdf/tiff", "Tests whether the HRC association depends on the primary 6-24 h window.",
  "Supplementary Figure S8", "Landmark and early-event sensitivity", "supp_figure_s8_landmark_early_event_sensitivity_v1.png/svg/pdf/tiff", "Defends the 24-h landmark design with early-event profiling and IPCW/composite sensitivity.",
  "Supplementary Figure S9", "Missing-data and alternative HRC construction sensitivity", "supp_figure_s9_missing_alternative_hrc_sensitivity_v1.png/svg/pdf/tiff", "Shows robustness to missing-data handling and alternative score construction.",
  "Supplementary Figure S10", "Permuted-HRC negative-control analysis", "supp_figure_s10_negative_control_v1.png/svg/pdf/tiff", "Shows observed associations are not reproduced by random HRC assignment."
)
write_csv(fig_index, file.path(outdir, "supplementary_figure_index_v1.csv"))
write_md_table(fig_index, file.path(outdir, "supplementary_figure_index_v1.md"), "Annals supplementary figure package v1")

qa <- c(
  "# Supplementary figure QA notes v1",
  "",
  "Backend: R/ggplot2/patchwork only.",
  "Archetype: quantitative-grid supplementary defense figures.",
  "Core conclusion: the HRC association is supported by auditable cohort construction, residualized component behavior, consistent dose-response, and sensitivity analyses.",
  "",
  "Figure-level claims:",
  "- S1: cohort attrition and landmark retention are transparent by database.",
  "- S2: HRC residual components are standardized and show modest, interpretable domain coherence.",
  "- S3: final HRC is centered and scaled within each database.",
  "- S4: observed recovery is evaluated against expected recovery, not raw delta alone.",
  "- S5: mortality risk falls with higher HRC across all three databases.",
  "- S6: lower HRC tracks higher downstream organ nonrecovery component rates.",
  "- S7: the association is not dependent on a single response window.",
  "- S8: the 24-h landmark is supported by IPCW and early-event composite sensitivity analyses.",
  "- S9: findings persist across missing-data and alternative HRC-construction strategies.",
  "- S10: permuted HRC does not reproduce the observed associations.",
  "",
  "Exports: SVG/PDF/TIFF/PNG produced for every supplementary figure. Source data CSVs were written for every panel."
)
writeLines(qa, file.path(outdir, "supplementary_figure_qa_notes_v1.md"))

cat("Supplementary Figures S1-S10 written to ", outdir, "\n", sep = "")
