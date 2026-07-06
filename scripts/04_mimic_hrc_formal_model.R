#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "outputs/mimic_formal/mimic_formal_cohort.csv"
outdir <- if (length(args) >= 2) args[[2]] else "outputs/mimic_formal"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260704)

d <- read.csv(input_path, stringsAsFactors = FALSE)

numeric_cols <- c(
  "subject_id", "hadm_id", "stay_id", "anchor_age", "weight_kg",
  "first_day_sofa", "sofa_respiration", "sofa_coagulation", "sofa_liver",
  "sofa_cardiovascular", "sofa_cns", "sofa_renal", "sofa2_total",
  "hospital_expire_flag", "hospital_mortality_after_landmark",
  "cardiogenic_shock_icd", "mcs_any_icu", "mcs_any_0_24h", "mcs_pre_index",
  "iabp_any_icu", "impella_any_icu", "ecmo_any_icu",
  "iabp_0_24h", "impella_0_24h", "ecmo_0_24h", "cs_or_mcs_subgroup",
  "baseline_map_mean", "response_map_mean", "delta_map",
  "baseline_neq_avg", "response_neq_avg", "delta_neq_reduction",
  "log_neq_reduction", "baseline_uo_ml", "response_uo_ml",
  "baseline_uo_ml_kg_h", "response_uo_ml_kg_h", "delta_uo_ml_kg_h",
  "log_uo_recovery", "baseline_creatinine_mean", "response_creatinine_mean",
  "delta_creatinine_reduction", "baseline_lactate_mean",
  "response_lactate_mean", "delta_lactate_reduction"
)
for (col in intersect(numeric_cols, names(d))) {
  d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
}

d$gender_male <- ifelse(d$gender == "M", 1, 0)
d$baseline_neq_log <- log1p(pmax(d$baseline_neq_avg, 0))
d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))

required <- c(
  "delta_map", "log_neq_reduction", "log_uo_recovery",
  "hospital_mortality_after_landmark",
  "baseline_map_mean", "baseline_neq_avg", "baseline_uo_ml_kg_h"
)
core <- d[complete.cases(d[, required]), ]
core <- core[is.finite(core$delta_map) &
               is.finite(core$log_neq_reduction) &
               is.finite(core$log_uo_recovery), ]

impute_median <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

zscore <- function(x) {
  as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

for (col in c(
  "anchor_age", "gender_male", "weight_kg", "first_day_sofa",
  "sofa_respiration", "sofa_cardiovascular", "sofa_renal",
  "baseline_map_mean", "baseline_neq_log", "baseline_uo_log",
  "baseline_neq_avg", "baseline_uo_ml_kg_h"
)) {
  core[[col]] <- impute_median(core[[col]])
}

folds <- sample(rep(1:5, length.out = nrow(core)))

expected_formula <- function(outcome) {
  as.formula(paste0(
    outcome,
    " ~ ns(anchor_age, df = 3) + gender_male + ns(weight_kg, df = 3) + ",
    "ns(first_day_sofa, df = 4) + sofa_respiration + ",
    "sofa_cardiovascular + sofa_renal + ",
    "ns(baseline_map_mean, df = 4) + ns(baseline_neq_log, df = 4) + ",
    "ns(baseline_uo_log, df = 4)"
  ))
}

crossfit_component <- function(dat, outcome, folds) {
  expected <- rep(NA_real_, nrow(dat))
  fold_rmse <- rep(NA_real_, max(folds))
  for (fold in sort(unique(folds))) {
    train <- dat[folds != fold, ]
    test <- dat[folds == fold, ]
    fit <- lm(expected_formula(outcome), data = train)
    expected[folds == fold] <- as.numeric(predict(fit, newdata = test))
    fold_rmse[fold] <- sqrt(mean((test[[outcome]] - expected[folds == fold])^2, na.rm = TRUE))
  }
  residual <- dat[[outcome]] - expected
  list(
    expected = expected,
    residual = residual,
    residual_z = zscore(residual),
    rmse = mean(fold_rmse, na.rm = TRUE),
    r2_oof = 1 - sum(residual^2, na.rm = TRUE) /
      sum((dat[[outcome]] - mean(dat[[outcome]], na.rm = TRUE))^2, na.rm = TRUE)
  )
}

map_cf <- crossfit_component(core, "delta_map", folds)
vaso_cf <- crossfit_component(core, "log_neq_reduction", folds)
uo_cf <- crossfit_component(core, "log_uo_recovery", folds)

core$expected_delta_map <- map_cf$expected
core$expected_log_neq_reduction <- vaso_cf$expected
core$expected_log_uo_recovery <- uo_cf$expected
core$hrc_map_residual_z <- map_cf$residual_z
core$hrc_vaso_residual_z <- vaso_cf$residual_z
core$hrc_uo_residual_z <- uo_cf$residual_z
core$hrc_core_raw <- rowMeans(
  core[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")],
  na.rm = FALSE
)
core$hrc_core_z <- zscore(core$hrc_core_raw)

quartile_cut <- quantile(core$hrc_core_z, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
core$hrc_quartile <- cut(
  core$hrc_core_z,
  breaks = c(-Inf, quartile_cut, Inf),
  labels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
  include.lowest = TRUE
)
core$low_hrc_q1 <- ifelse(core$hrc_quartile == "Q1_lowest", 1, 0)

score_cols <- c(
  "subject_id", "hadm_id", "stay_id", "index_time", "landmark_time",
  "first_careunit", "last_careunit", "anchor_age", "gender", "weight_kg",
  "first_day_sofa", "sofa_respiration", "sofa_cardiovascular", "sofa_renal",
  "hospital_mortality_after_landmark",
  "cardiogenic_shock_icd", "mcs_any_icu", "mcs_any_0_24h",
  "iabp_any_icu", "impella_any_icu", "ecmo_any_icu", "cs_or_mcs_subgroup",
  "hrc_core_z", "hrc_core_raw", "hrc_quartile", "low_hrc_q1",
  "hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z",
  "delta_map", "expected_delta_map",
  "log_neq_reduction", "expected_log_neq_reduction",
  "log_uo_recovery", "expected_log_uo_recovery",
  "baseline_map_mean", "baseline_neq_avg", "baseline_uo_ml_kg_h",
  "response_map_mean", "response_neq_avg", "response_uo_ml_kg_h"
)
write.csv(
  core[, score_cols],
  file = file.path(outdir, "mimic_hrc_formal_scores.csv"),
  row.names = FALSE
)

component_metrics <- data.frame(
  component = c("MAP recovery", "Vasopressor reduction", "Urine output recovery"),
  observed_variable = c("delta_map", "log_neq_reduction", "log_uo_recovery"),
  n = nrow(core),
  oof_rmse = c(map_cf$rmse, vaso_cf$rmse, uo_cf$rmse),
  oof_r2 = c(map_cf$r2_oof, vaso_cf$r2_oof, uo_cf$r2_oof)
)
write.csv(
  component_metrics,
  file = file.path(outdir, "mimic_hrc_formal_component_model_metrics.csv"),
  row.names = FALSE
)

quartile_summary <- do.call(rbind, lapply(levels(core$hrc_quartile), function(q) {
  x <- core[core$hrc_quartile == q, ]
  data.frame(
    hrc_quartile = q,
    n = nrow(x),
    deaths = sum(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
    mortality_rate = mean(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
    median_hrc = median(x$hrc_core_z, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  quartile_summary,
  file = file.path(outdir, "mimic_hrc_formal_quartile_mortality.csv"),
  row.names = FALSE
)

adjust_covars <- c(
  "anchor_age", "first_day_sofa", "baseline_neq_log",
  "baseline_map_mean", "baseline_uo_log"
)
linear_fit <- glm(
  hospital_mortality_after_landmark ~ hrc_core_z + ns(anchor_age, df = 3) +
    ns(first_day_sofa, df = 4) + ns(baseline_neq_log, df = 4) +
    ns(baseline_map_mean, df = 4) + ns(baseline_uo_log, df = 4),
  data = core,
  family = binomial()
)
coef_table <- summary(linear_fit)$coefficients
ci_table <- confint.default(linear_fit)
linear_rows <- data.frame(
  term = rownames(coef_table),
  beta = coef_table[, "Estimate"],
  se = coef_table[, "Std. Error"],
  odds_ratio = exp(coef_table[, "Estimate"]),
  ci_low = exp(ci_table[, 1]),
  ci_high = exp(ci_table[, 2]),
  p_value = coef_table[, "Pr(>|z|)"],
  row.names = NULL
)
write.csv(
  linear_rows,
  file = file.path(outdir, "mimic_hrc_formal_linear_or.csv"),
  row.names = FALSE
)

spline_fit <- glm(
  hospital_mortality_after_landmark ~ ns(hrc_core_z, df = 4) +
    ns(anchor_age, df = 3) + ns(first_day_sofa, df = 4) +
    ns(baseline_neq_log, df = 4) + ns(baseline_map_mean, df = 4) +
    ns(baseline_uo_log, df = 4),
  data = core,
  family = binomial()
)
hrc_grid <- seq(
  quantile(core$hrc_core_z, 0.01, na.rm = TRUE),
  quantile(core$hrc_core_z, 0.99, na.rm = TRUE),
  length.out = 200
)
newdata <- data.frame(
  hrc_core_z = hrc_grid,
  anchor_age = median(core$anchor_age, na.rm = TRUE),
  first_day_sofa = median(core$first_day_sofa, na.rm = TRUE),
  baseline_neq_log = median(core$baseline_neq_log, na.rm = TRUE),
  baseline_map_mean = median(core$baseline_map_mean, na.rm = TRUE),
  baseline_uo_log = median(core$baseline_uo_log, na.rm = TRUE)
)
pred <- predict(spline_fit, newdata = newdata, type = "link", se.fit = TRUE)
spline_pred <- data.frame(
  hrc_core_z = hrc_grid,
  mortality_probability = plogis(pred$fit),
  ci_low = plogis(pred$fit - 1.96 * pred$se.fit),
  ci_high = plogis(pred$fit + 1.96 * pred$se.fit)
)
write.csv(
  spline_pred,
  file = file.path(outdir, "mimic_hrc_formal_spline_predictions.csv"),
  row.names = FALSE
)

summarize_group <- function(label, x) {
  data.frame(
    subgroup = label,
    n = nrow(x),
    deaths = sum(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
    mortality_rate = mean(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
    low_hrc_n = sum(x$low_hrc_q1 == 1, na.rm = TRUE),
    low_hrc_deaths = sum(x$hospital_mortality_after_landmark == 1 & x$low_hrc_q1 == 1, na.rm = TRUE),
    low_hrc_mortality = mean(x$hospital_mortality_after_landmark[x$low_hrc_q1 == 1] == 1, na.rm = TRUE),
    nonlow_hrc_mortality = mean(x$hospital_mortality_after_landmark[x$low_hrc_q1 == 0] == 1, na.rm = TRUE),
    median_hrc = median(x$hrc_core_z, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
subgroups <- list(
  "Overall" = core,
  "No CS/MCS" = core[core$cs_or_mcs_subgroup == 0, ],
  "Cardiogenic shock ICD" = core[core$cardiogenic_shock_icd == 1, ],
  "Any MCS during ICU" = core[core$mcs_any_icu == 1, ],
  "MCS during 0-24h" = core[core$mcs_any_0_24h == 1, ],
  "IABP" = core[core$iabp_any_icu == 1, ],
  "Impella" = core[core$impella_any_icu == 1, ],
  "ECMO" = core[core$ecmo_any_icu == 1, ],
  "CS or MCS" = core[core$cs_or_mcs_subgroup == 1, ]
)
subgroup_summary <- do.call(rbind, Map(summarize_group, names(subgroups), subgroups))
write.csv(
  subgroup_summary,
  file = file.path(outdir, "mimic_hrc_formal_subgroup_summary.csv"),
  row.names = FALSE
)

flow_counts <- read.csv(file.path(outdir, "mimic_formal_cohort_summary.csv"), stringsAsFactors = FALSE)
write.csv(
  flow_counts,
  file = file.path(outdir, "mimic_hrc_formal_flow_counts.csv"),
  row.names = FALSE
)

summary_md <- c(
  "# MIMIC formal HRC model",
  "",
  "## Construct",
  "",
  "- Expected recovery was estimated within MIMIC using 5-fold cross-fitted natural-spline regression models.",
  "- Core recovery components were MAP recovery, log-scale vasopressor reduction, and log-scale urine-output recovery.",
  "- HRC was calculated as the equal-weighted mean of standardized out-of-fold residuals, then standardized again.",
  "- Higher HRC indicates better-than-expected physiologic recovery after early hemodynamic support.",
  "",
  "## Counts",
  "",
  paste0("- formal cohort rows: ", nrow(d)),
  paste0("- core HRC rows: ", nrow(core)),
  paste0("- CS or MCS subgroup rows in scored cohort: ", nrow(core[core$cs_or_mcs_subgroup == 1, ])),
  "",
  "## Component model performance",
  "",
  paste(
    apply(component_metrics, 1, function(row) {
      paste0(
        "- ", row[["component"]], ": OOF R2=", sprintf("%.3f", as.numeric(row[["oof_r2"]])),
        ", OOF RMSE=", sprintf("%.3f", as.numeric(row[["oof_rmse"]]))
      )
    }),
    collapse = "\n"
  ),
  "",
  "## Mortality by HRC quartile",
  "",
  paste(
    apply(quartile_summary, 1, function(row) {
      paste0(
        "- ", row[["hrc_quartile"]], ": n=", row[["n"]],
        ", deaths=", row[["deaths"]],
        ", mortality=", sprintf("%.3f", as.numeric(row[["mortality_rate"]]))
      )
    }),
    collapse = "\n"
  ),
  "",
  "## Linear adjusted association",
  "",
  paste(
    apply(linear_rows[linear_rows$term == "hrc_core_z", ], 1, function(row) {
      paste0(
        "- HRC per 1 SD: OR=", sprintf("%.3f", as.numeric(row[["odds_ratio"]])),
        " (", sprintf("%.3f", as.numeric(row[["ci_low"]])),
        "-", sprintf("%.3f", as.numeric(row[["ci_high"]])),
        "), p=", signif(as.numeric(row[["p_value"]]), 3)
      )
    }),
    collapse = "\n"
  )
)
writeLines(summary_md, con = file.path(outdir, "mimic_hrc_formal_model_summary.md"))

cat("Formal cohort rows:", nrow(d), "\n")
cat("Core HRC rows:", nrow(core), "\n")
cat("Scores:", file.path(outdir, "mimic_hrc_formal_scores.csv"), "\n")
cat("Model summary:", file.path(outdir, "mimic_hrc_formal_model_summary.md"), "\n")
