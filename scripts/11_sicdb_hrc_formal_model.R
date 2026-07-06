#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "outputs/sicdb_formal/sicdb_formal_cohort.csv"
outdir <- if (length(args) >= 2) args[[2]] else "outputs/sicdb_formal"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260704)

d <- read.csv(input_path, stringsAsFactors = FALSE)

numeric_cols <- c(
  "caseid", "patientid", "icu_offset_sec", "time_of_stay_sec",
  "icu_out_offset_sec", "index_offset_sec", "landmark_offset_sec",
  "age", "male", "weight_kg", "saps3", "heart_surgery_cpb_min",
  "heart_surgery_end_sec", "hospital_mortality_after_landmark",
  "baseline_map_mean", "response_map_mean", "delta_map",
  "baseline_vaso_burden", "response_vaso_burden",
  "delta_vaso_burden_reduction", "log_vaso_burden_reduction",
  "baseline_vaso_agents_n", "response_vaso_agents_n",
  "baseline_uo_ml", "response_uo_ml", "baseline_uo_ml_kg_h",
  "response_uo_ml_kg_h", "delta_uo_ml_kg_h", "log_uo_recovery",
  "baseline_creatinine_mean", "response_creatinine_mean",
  "delta_creatinine_reduction", "baseline_lactate_mean",
  "response_lactate_mean", "delta_lactate_reduction"
)
for (col in intersect(numeric_cols, names(d))) {
  d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
}

d$heart_surgery <- ifelse(d$heartsurgeryadditionaldata == "740", 1, 0)
d$weight_kg[d$weight_kg <= 0 | d$weight_kg > 300] <- NA_real_
d$baseline_vaso_log <- log1p(pmax(d$baseline_vaso_burden, 0))
d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))

required <- c(
  "delta_map", "log_vaso_burden_reduction", "log_uo_recovery",
  "hospital_mortality_after_landmark",
  "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h"
)
core <- d[complete.cases(d[, required]), ]
core <- core[is.finite(core$delta_map) &
               is.finite(core$log_vaso_burden_reduction) &
               is.finite(core$log_uo_recovery), ]

impute_median <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

for (col in c(
  "age", "male", "weight_kg", "saps3", "heart_surgery",
  "baseline_map_mean", "baseline_vaso_log", "baseline_uo_log",
  "baseline_vaso_burden", "baseline_uo_ml_kg_h"
)) {
  core[[col]] <- impute_median(core[[col]])
}

folds <- sample(rep(1:5, length.out = nrow(core)))

expected_formula <- function(outcome) {
  as.formula(paste0(
    outcome,
    " ~ ns(age, df = 3) + male + ns(weight_kg, df = 3) + ",
    "ns(saps3, df = 4) + heart_surgery + ",
    "ns(baseline_map_mean, df = 4) + ns(baseline_vaso_log, df = 4) + ",
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
vaso_cf <- crossfit_component(core, "log_vaso_burden_reduction", folds)
uo_cf <- crossfit_component(core, "log_uo_recovery", folds)

core$expected_delta_map <- map_cf$expected
core$expected_log_vaso_burden_reduction <- vaso_cf$expected
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
  "caseid", "patientid", "hospitalunit", "hospitalunit_label",
  "age", "male", "weight_kg", "saps3", "heart_surgery",
  "index_offset_sec", "landmark_offset_sec",
  "hospital_mortality_after_landmark",
  "hrc_core_z", "hrc_core_raw", "hrc_quartile", "low_hrc_q1",
  "hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z",
  "delta_map", "expected_delta_map",
  "log_vaso_burden_reduction", "expected_log_vaso_burden_reduction",
  "log_uo_recovery", "expected_log_uo_recovery",
  "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h",
  "response_map_mean", "response_vaso_burden", "response_uo_ml_kg_h"
)
write.csv(
  core[, score_cols],
  file = file.path(outdir, "sicdb_hrc_formal_scores.csv"),
  row.names = FALSE
)

component_metrics <- data.frame(
  component = c("MAP recovery", "Vasopressor-burden reduction", "Urine output recovery"),
  observed_variable = c("delta_map", "log_vaso_burden_reduction", "log_uo_recovery"),
  n = nrow(core),
  oof_rmse = c(map_cf$rmse, vaso_cf$rmse, uo_cf$rmse),
  oof_r2 = c(map_cf$r2_oof, vaso_cf$r2_oof, uo_cf$r2_oof)
)
write.csv(
  component_metrics,
  file = file.path(outdir, "sicdb_hrc_formal_component_model_metrics.csv"),
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
  file = file.path(outdir, "sicdb_hrc_formal_quartile_mortality.csv"),
  row.names = FALSE
)

linear_fit <- glm(
  hospital_mortality_after_landmark ~ hrc_core_z + ns(age, df = 3) +
    male + ns(saps3, df = 4) + heart_surgery +
    ns(baseline_vaso_log, df = 4) + ns(baseline_map_mean, df = 4) +
    ns(baseline_uo_log, df = 4),
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
  file = file.path(outdir, "sicdb_hrc_formal_linear_or.csv"),
  row.names = FALSE
)

low_fit <- glm(
  hospital_mortality_after_landmark ~ low_hrc_q1 + ns(age, df = 3) +
    male + ns(saps3, df = 4) + heart_surgery +
    ns(baseline_vaso_log, df = 4) + ns(baseline_map_mean, df = 4) +
    ns(baseline_uo_log, df = 4),
  data = core,
  family = binomial()
)
low_coef <- summary(low_fit)$coefficients
low_ci <- confint.default(low_fit)
low_rows <- data.frame(
  term = rownames(low_coef),
  beta = low_coef[, "Estimate"],
  se = low_coef[, "Std. Error"],
  odds_ratio = exp(low_coef[, "Estimate"]),
  ci_low = exp(low_ci[, 1]),
  ci_high = exp(low_ci[, 2]),
  p_value = low_coef[, "Pr(>|z|)"],
  row.names = NULL
)
write.csv(
  low_rows,
  file = file.path(outdir, "sicdb_hrc_formal_low_hrc_or.csv"),
  row.names = FALSE
)

spline_fit <- glm(
  hospital_mortality_after_landmark ~ ns(hrc_core_z, df = 4) + ns(age, df = 3) +
    male + ns(saps3, df = 4) + heart_surgery +
    ns(baseline_vaso_log, df = 4) + ns(baseline_map_mean, df = 4) +
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
  age = median(core$age, na.rm = TRUE),
  male = median(core$male, na.rm = TRUE),
  saps3 = median(core$saps3, na.rm = TRUE),
  heart_surgery = median(core$heart_surgery, na.rm = TRUE),
  baseline_vaso_log = median(core$baseline_vaso_log, na.rm = TRUE),
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
  file = file.path(outdir, "sicdb_hrc_formal_spline_predictions.csv"),
  row.names = FALSE
)

cat("SICdb core HRC rows:", nrow(core), "\n")
cat("SICdb HRC OR per SD:", linear_rows$odds_ratio[linear_rows$term == "hrc_core_z"], "\n")
cat("SICdb low-HRC OR:", low_rows$odds_ratio[low_rows$term == "low_hrc_q1"], "\n")
