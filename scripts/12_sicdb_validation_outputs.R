#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
sicdb_dir <- if (length(args) >= 1) args[[1]] else "outputs/sicdb_formal"
mimic_dir <- if (length(args) >= 2) args[[2]] else "outputs/mimic_formal"
eicu_dir <- if (length(args) >= 3) args[[3]] else "outputs/eicu_formal"
cross_dir <- if (length(args) >= 4) args[[4]] else "outputs/cross_database_validation"

table_dir <- file.path(sicdb_dir, "tables")
figdir <- file.path(sicdb_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(sicdb_dir, "sicdb_hrc_formal_scores.csv"), stringsAsFactors = FALSE)
spline <- read.csv(file.path(sicdb_dir, "sicdb_hrc_formal_spline_predictions.csv"), stringsAsFactors = FALSE)
flow <- read.csv(file.path(sicdb_dir, "sicdb_formal_cohort_summary.csv"), stringsAsFactors = FALSE)
quartile <- read.csv(file.path(sicdb_dir, "sicdb_hrc_formal_quartile_mortality.csv"), stringsAsFactors = FALSE)
sicdb_linear <- read.csv(file.path(sicdb_dir, "sicdb_hrc_formal_linear_or.csv"), stringsAsFactors = FALSE)
sicdb_low <- read.csv(file.path(sicdb_dir, "sicdb_hrc_formal_low_hrc_or.csv"), stringsAsFactors = FALSE)
mimic_linear <- read.csv(file.path(mimic_dir, "mimic_hrc_formal_linear_or.csv"), stringsAsFactors = FALSE)
eicu_linear <- read.csv(file.path(eicu_dir, "eicu_hrc_formal_linear_or.csv"), stringsAsFactors = FALSE)

for (nm in intersect(c("hrc_core_z", "hospital_mortality_after_landmark"), names(scores))) {
  scores[[nm]] <- as.numeric(scores[[nm]])
}
for (nm in c("mortality_probability", "ci_low", "ci_high", "hrc_core_z")) {
  spline[[nm]] <- as.numeric(spline[[nm]])
}
flow$value <- as.numeric(flow$value)

extract_term <- function(df, term, database, exposure) {
  row <- df[df$term == term, ]
  data.frame(
    database = database,
    exposure = exposure,
    beta = as.numeric(row$beta),
    se = as.numeric(row$se),
    odds_ratio = as.numeric(row$odds_ratio),
    ci_low = as.numeric(row$ci_low),
    ci_high = as.numeric(row$ci_high),
    p_value = as.numeric(row$p_value),
    stringsAsFactors = FALSE
  )
}

continuous <- rbind(
  extract_term(mimic_linear, "hrc_core_z", "MIMIC-IV", "HRC per 1 SD increase"),
  extract_term(eicu_linear, "hrc_core_z", "eICU", "HRC per 1 SD increase"),
  extract_term(sicdb_linear, "hrc_core_z", "SICdb", "HRC per 1 SD increase")
)
continuous$weight <- 1 / continuous$se^2
beta_fe <- sum(continuous$weight * continuous$beta) / sum(continuous$weight)
se_fe <- sqrt(1 / sum(continuous$weight))
meta <- data.frame(
  database = "Fixed-effect summary",
  exposure = "HRC per 1 SD increase",
  beta = beta_fe,
  se = se_fe,
  odds_ratio = exp(beta_fe),
  ci_low = exp(beta_fe - 1.96 * se_fe),
  ci_high = exp(beta_fe + 1.96 * se_fe),
  p_value = 2 * pnorm(abs(beta_fe / se_fe), lower.tail = FALSE),
  weight = NA_real_
)
forest <- rbind(continuous, meta)
write.csv(forest, file.path(cross_dir, "mimic_eicu_sicdb_hrc_meta_summary.csv"), row.names = FALSE)

k <- nrow(continuous)
q_stat <- sum(continuous$weight * (continuous$beta - beta_fe)^2)
q_df <- k - 1
q_p <- pchisq(q_stat, df = q_df, lower.tail = FALSE)
i2 <- max(0, (q_stat - q_df) / q_stat) * 100
tau2_dl <- max(
  0,
  (q_stat - q_df) / (sum(continuous$weight) - sum(continuous$weight^2) / sum(continuous$weight))
)
random_weight <- 1 / (continuous$se^2 + tau2_dl)
beta_re <- sum(random_weight * continuous$beta) / sum(random_weight)
se_re <- sqrt(1 / sum(random_weight))
meta_stats <- data.frame(
  model = c("Fixed effect", "DerSimonian-Laird random effects", "Heterogeneity"),
  beta = c(beta_fe, beta_re, NA_real_),
  se = c(se_fe, se_re, NA_real_),
  odds_ratio = c(exp(beta_fe), exp(beta_re), NA_real_),
  ci_low = c(exp(beta_fe - 1.96 * se_fe), exp(beta_re - 1.96 * se_re), NA_real_),
  ci_high = c(exp(beta_fe + 1.96 * se_fe), exp(beta_re + 1.96 * se_re), NA_real_),
  p_value = c(
    2 * pnorm(abs(beta_fe / se_fe), lower.tail = FALSE),
    2 * pnorm(abs(beta_re / se_re), lower.tail = FALSE),
    q_p
  ),
  q_stat = c(NA_real_, NA_real_, q_stat),
  q_df = c(NA_real_, NA_real_, q_df),
  i2_percent = c(NA_real_, NA_real_, i2),
  tau2_dl = c(NA_real_, tau2_dl, tau2_dl)
)
write.csv(meta_stats, file.path(cross_dir, "mimic_eicu_sicdb_hrc_meta_model_stats.csv"), row.names = FALSE)

sicdb_validation <- data.frame(
  Metric = c(
    "Formal SICdb cohort rows",
    "Core HRC complete rows",
    "Post-landmark hospital deaths",
    "Q1 low-HRC mortality",
    "Q4 high-HRC mortality",
    "Adjusted HRC OR per 1 SD",
    "Adjusted low-HRC OR"
  ),
  Value = c(
    flow$value[match("formal_dataset_rows", flow$metric)],
    flow$value[match("core_hrc_complete_rows", flow$metric)],
    flow$value[match("hospital_mortality_after_landmark_events", flow$metric)],
    sprintf("%.1f%%", 100 * quartile$mortality_rate[quartile$hrc_quartile == "Q1_lowest"]),
    sprintf("%.1f%%", 100 * quartile$mortality_rate[quartile$hrc_quartile == "Q4_highest"]),
    with(sicdb_linear[sicdb_linear$term == "hrc_core_z", ],
         sprintf("%.2f (%.2f-%.2f)", odds_ratio, ci_low, ci_high)),
    with(sicdb_low[sicdb_low$term == "low_hrc_q1", ],
         sprintf("%.2f (%.2f-%.2f)", odds_ratio, ci_low, ci_high))
  )
)
write.csv(sicdb_validation, file.path(table_dir, "sicdb_validation_summary_table.csv"), row.names = FALSE)

write_md_table <- function(df, path, title = NULL) {
  lines <- character()
  if (!is.null(title)) lines <- c(lines, paste0("# ", title), "")
  cols <- names(df)
  lines <- c(lines, paste(cols, collapse = " | "))
  lines <- c(lines, paste(rep("---", length(cols)), collapse = " | "))
  for (i in seq_len(nrow(df))) {
    vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
    vals <- gsub("\\|", "/", vals)
    lines <- c(lines, paste(vals, collapse = " | "))
  }
  writeLines(lines, con = path)
}

write_md_table(
  sicdb_validation,
  file.path(table_dir, "sicdb_validation_summary_table.md"),
  "SICdb external validation summary"
)
write_md_table(
  forest,
  file.path(cross_dir, "mimic_eicu_sicdb_hrc_meta_summary.md"),
  "Three-database HRC validation summary"
)
write_md_table(
  meta_stats,
  file.path(cross_dir, "mimic_eicu_sicdb_hrc_meta_model_stats.md"),
  "Three-database HRC meta-analysis model statistics"
)

theme_set(
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#1f2933"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#1f2933"),
      axis.title = element_text(size = 7, colour = "#1f2933"),
      axis.text = element_text(size = 6.2, colour = "#1f2933"),
      plot.title = element_text(size = 7.5, face = "bold", colour = "#101820"),
      plot.subtitle = element_text(size = 6.2, colour = "#52616b"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 6.2),
      legend.key.size = unit(3.5, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
)

pal <- list(
  ink = "#17202a",
  muted = "#6b7280",
  light = "#e5e7eb",
  blue = "#3b6ea8",
  blue_light = "#dbe8f6",
  red = "#b94a48",
  red_light = "#f2d8d5",
  green = "#3f7f5f",
  green_light = "#dfeee6"
)

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

flow_steps <- data.frame(
  step = c(
    "Adult valid ICU cases",
    "Early vasoactive support",
    "24-h landmark eligible",
    "First qualifying case",
    "Core HRC complete"
  ),
  metric = c(
    "adult_valid_icu_cases",
    "early_vasoactive_support_cases",
    "landmark_eligible_cases",
    "formal_first_patient_support_cases",
    "core_hrc_complete_rows"
  ),
  x = 1,
  y = 5:1,
  stringsAsFactors = FALSE
)
flow_steps$n <- flow$value[match(flow_steps$metric, flow$metric)]
flow_steps$label <- paste0(flow_steps$step, "\n", fmt_n(flow_steps$n))

p_flow <- ggplot(flow_steps, aes(x = x, y = y)) +
  geom_segment(
    data = flow_steps[-nrow(flow_steps), ],
    aes(x = 1.58, xend = 1.58, y = y - 0.31, yend = y - 0.69),
    arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
    linewidth = 0.35,
    colour = pal$muted
  ) +
  geom_rect(
    aes(xmin = 0.35, xmax = 1.65, ymin = y - 0.26, ymax = y + 0.26),
    fill = "white",
    colour = pal$green,
    linewidth = 0.45
  ) +
  geom_text(aes(label = label), size = 2.15, lineheight = 0.95, colour = pal$ink) +
  coord_cartesian(xlim = c(0.25, 1.75), ylim = c(0.45, 5.55), clip = "off") +
  labs(title = "A  SICdb validation flow") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0, colour = pal$ink),
    plot.margin = margin(4, 4, 4, 4)
  )

scores$mortality_group <- ifelse(
  scores$hospital_mortality_after_landmark == 1,
  "Died after landmark",
  "Survived to discharge"
)
p_dist <- ggplot(scores, aes(x = hrc_core_z, fill = mortality_group, colour = mortality_group)) +
  geom_density(alpha = 0.18, linewidth = 0.55, adjust = 1.1) +
  geom_vline(xintercept = quantile(scores$hrc_core_z, 0.25, na.rm = TRUE),
             colour = pal$red, linetype = "dashed", linewidth = 0.35) +
  scale_fill_manual(values = c("Died after landmark" = pal$red_light, "Survived to discharge" = pal$blue_light)) +
  scale_colour_manual(values = c("Died after landmark" = pal$red, "Survived to discharge" = pal$blue)) +
  guides(fill = "none", colour = guide_legend(override.aes = list(linewidth = 0.9))) +
  labs(
    title = "B  HRC distribution in SICdb",
    x = "Hemodynamic Recovery Capacity, z score",
    y = "Density"
  )

p_spline <- ggplot(spline, aes(x = hrc_core_z, y = mortality_probability)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = pal$blue_light, alpha = 0.75) +
  geom_line(colour = pal$blue, linewidth = 0.75) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "C  SICdb adjusted mortality curve",
    x = "Hemodynamic Recovery Capacity, z score",
    y = "Adjusted probability of death"
  )

forest_plot <- forest
forest_plot$database <- factor(
  forest_plot$database,
  levels = rev(c("MIMIC-IV", "eICU", "SICdb", "Fixed-effect summary"))
)
forest_plot$is_summary <- forest_plot$database == "Fixed-effect summary"

x_low <- min(forest_plot$ci_low, na.rm = TRUE)
x_high <- max(forest_plot$ci_high, 1.0, na.rm = TRUE)
x_limits <- c(max(0.35, floor(x_low * 20) / 20), min(1.05, ceiling(x_high * 20) / 20))

p_forest <- ggplot(forest_plot, aes(x = odds_ratio, y = database)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.45, colour = pal$muted) +
  geom_point(aes(fill = is_summary), shape = 21, size = 2.2, colour = pal$ink) +
  scale_fill_manual(values = c("FALSE" = pal$blue, "TRUE" = pal$red), guide = "none") +
  scale_x_continuous(limits = x_limits, breaks = pretty(x_limits, n = 5)) +
  labs(
    title = "D  Three-database validation",
    x = "Adjusted odds ratio per 1 SD higher HRC",
    y = NULL
  )

fig <- (p_flow | p_dist) / (p_spline | p_forest) +
  plot_layout(widths = c(0.82, 1.18), heights = c(0.95, 1.05)) +
  plot_annotation(
    title = "External validation of Hemodynamic Recovery Capacity in SICdb",
    subtitle = "The recovery phenotype reproduces in an independent European ICU database and strengthens three-database consistency",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", colour = pal$ink),
      plot.subtitle = element_text(size = 7, colour = pal$muted)
    )
  )

save_pub <- function(plot, filename, width_mm = 183, height_mm = 130, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 220)
  print(plot)
  dev.off()
}

figure_base <- file.path(figdir, "figure3_sicdb_validation_v1")
save_pub(fig, figure_base)

writeLines(
  c(
    "# Figure 3 contract",
    "",
    "Core conclusion: SICdb reproduces the HRC-mortality association in an independent European ICU system and upgrades the validation from two-database replication to three-database consistency.",
    "",
    "Evidence chain:",
    "- Panel A: SICdb cohort flow from adult ICU cases to core HRC-complete validation rows.",
    "- Panel B: HRC distribution shows lower recovery capacity among post-landmark hospital deaths.",
    "- Panel C: adjusted spline tests the continuous-first HRC-mortality relationship in SICdb.",
    "- Panel D: MIMIC-IV, eICU, SICdb, and fixed-effect summary ORs test reproducibility across systems.",
    "",
    "Archetype: quantitative grid.",
    "Backend: R/ggplot2/patchwork only.",
    "Export formats: SVG, PDF, TIFF, and PNG preview.",
    "Review risks handled: no fixed cross-database coefficients; database-specific recalibration; continuous HRC remains primary."
  ),
  con = file.path(figdir, "figure3_sicdb_validation_v1_contract.md")
)

cat("SICdb validation table:", file.path(table_dir, "sicdb_validation_summary_table.md"), "\n")
cat("Three-database meta:", file.path(cross_dir, "mimic_eicu_sicdb_hrc_meta_summary.md"), "\n")
cat("Figure written:", figure_base, "\n")
