#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "outputs/mimic_formal"
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(outdir, "mimic_hrc_formal_scores.csv"), stringsAsFactors = FALSE)
spline <- read.csv(file.path(outdir, "mimic_hrc_formal_spline_predictions.csv"), stringsAsFactors = FALSE)
subgroups <- read.csv(file.path(outdir, "mimic_hrc_formal_subgroup_summary.csv"), stringsAsFactors = FALSE)
flow <- read.csv(file.path(outdir, "mimic_formal_cohort_summary.csv"), stringsAsFactors = FALSE)

num_cols <- c(
  "hrc_core_z", "hospital_mortality_after_landmark", "low_hrc_q1",
  "mortality_probability", "ci_low", "ci_high", "n", "deaths",
  "mortality_rate", "low_hrc_mortality", "nonlow_hrc_mortality"
)
for (nm in intersect(num_cols, names(scores))) scores[[nm]] <- as.numeric(scores[[nm]])
for (nm in intersect(num_cols, names(spline))) spline[[nm]] <- as.numeric(spline[[nm]])
for (nm in intersect(num_cols, names(subgroups))) subgroups[[nm]] <- as.numeric(subgroups[[nm]])
flow$value <- as.numeric(flow$value)

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
      legend.spacing.x = unit(2.2, "mm"),
      strip.text = element_text(size = 6.5, face = "bold"),
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
  green = "#3f7f5f"
)

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * x), "%")

flow_steps <- data.frame(
  step = c(
    "Adult valid ICU stays",
    "Early vasoactive support",
    "24-h landmark eligible",
    "First qualifying stay",
    "Core HRC complete"
  ),
  metric = c(
    "adult_valid_icu_stays",
    "early_vaso_support_stays",
    "landmark_eligible_stays",
    "formal_first_patient_support_stays",
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
    colour = pal$blue,
    linewidth = 0.45
  ) +
  geom_text(aes(label = label), size = 2.15, lineheight = 0.95, colour = pal$ink) +
  coord_cartesian(xlim = c(0.25, 1.75), ylim = c(0.45, 5.55), clip = "off") +
  labs(title = "A  MIMIC cohort flow") +
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
  guides(
    fill = "none",
    colour = guide_legend(override.aes = list(linewidth = 0.9))
  ) +
  labs(
    title = "B  HRC distribution",
    x = "Hemodynamic Recovery Capacity, z score",
    y = "Density"
  ) +
  theme(legend.position = "top")

p_spline <- ggplot(spline, aes(x = hrc_core_z, y = mortality_probability)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = pal$blue_light, alpha = 0.75) +
  geom_line(colour = pal$blue, linewidth = 0.75) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "C  Adjusted mortality curve",
    x = "Hemodynamic Recovery Capacity, z score",
    y = "Adjusted probability of death"
  )

subgroup_order <- c(
  "Overall", "No CS/MCS", "Cardiogenic shock ICD",
  "Any MCS during ICU", "IABP", "Impella", "ECMO"
)
sg <- subgroups[subgroups$subgroup %in% subgroup_order, ]
sg$subgroup <- factor(sg$subgroup, levels = rev(subgroup_order))
sg_long <- rbind(
  data.frame(
    subgroup = sg$subgroup,
    hrc_status = "Non-low HRC",
    mortality = sg$nonlow_hrc_mortality,
    n = sg$n
  ),
  data.frame(
    subgroup = sg$subgroup,
    hrc_status = "Low HRC",
    mortality = sg$low_hrc_mortality,
    n = sg$n
  )
)
sg_long$hrc_status <- factor(sg_long$hrc_status, levels = c("Non-low HRC", "Low HRC"))

p_subgroup <- ggplot(sg_long, aes(x = mortality, y = subgroup, colour = hrc_status)) +
  geom_line(aes(group = subgroup), colour = pal$light, linewidth = 1.4) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Non-low HRC" = pal$blue, "Low HRC" = pal$red)) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "D  Low-HRC phenotype across subgroups",
    x = "Hospital mortality after landmark",
    y = NULL
  ) +
  theme(legend.position = "top")

fig <- (p_flow | p_dist) / (p_spline | p_subgroup) +
  plot_layout(widths = c(0.82, 1.18), heights = c(0.95, 1.05)) +
  plot_annotation(
    title = "Hemodynamic Recovery Capacity in MIMIC-IV",
    subtitle = "A residualized recovery phenotype after early hemodynamic support",
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

figure_base <- file.path(figdir, "figure1_mimic_hrc_v1")
save_pub(fig, figure_base)

writeLines(
  c(
    "# Figure 1 contract",
    "",
    "Core conclusion: HRC captures a residualized physiologic recovery gradient after early hemodynamic support in MIMIC-IV.",
    "",
    "Panel A: cohort flow from adult ICU stays to the core HRC complete analysis set.",
    "Panel B: HRC distribution by post-landmark hospital survival status.",
    "Panel C: adjusted spline curve linking HRC to post-landmark hospital mortality.",
    "Panel D: low-HRC versus non-low-HRC mortality across cardiogenic shock and MCS subgroups.",
    "",
    "Export formats: SVG, PDF, TIFF, and PNG preview."
  ),
  con = file.path(figdir, "figure1_mimic_hrc_v1_contract.md")
)

cat("Figure written:", figure_base, "\n")
