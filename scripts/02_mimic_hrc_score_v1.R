#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "outputs/mimic_feasibility/mimic_hrc_feasibility_dataset.csv"
outdir <- if (length(args) >= 2) args[[2]] else "outputs/mimic_feasibility"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

analysis <- read.csv(input_path, stringsAsFactors = FALSE)

numeric_cols <- c(
  "anchor_age", "weight_kg", "first_day_sofa", "sofa_respiration",
  "sofa_cardiovascular", "sofa_renal", "sofa2_total",
  "hospital_expire_flag", "hospital_mortality_after_landmark",
  "baseline_map_mean", "response_map_mean", "delta_map",
  "baseline_neq_avg", "response_neq_avg", "delta_neq_reduction",
  "baseline_uo_ml", "response_uo_ml", "baseline_uo_ml_kg_h",
  "response_uo_ml_kg_h", "delta_uo_ml_kg_h",
  "baseline_creatinine_mean", "response_creatinine_mean",
  "delta_creatinine_reduction",
  "baseline_lactate_mean", "response_lactate_mean",
  "delta_lactate_reduction"
)
for (col in intersect(numeric_cols, names(analysis))) {
  analysis[[col]] <- suppressWarnings(as.numeric(analysis[[col]]))
}

analysis$gender_male <- ifelse(analysis$gender == "M", 1, 0)

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

fit_component <- function(dat, outcome, covars) {
  model_formula <- as.formula(paste(outcome, "~", paste(covars, collapse = " + ")))
  fit <- lm(model_formula, data = dat)
  observed <- dat[[outcome]]
  expected <- as.numeric(predict(fit, newdata = dat))
  residual <- observed - expected
  list(
    fit = fit,
    expected = expected,
    residual = residual,
    residual_z = zscore(residual)
  )
}

core_observed <- c("delta_map", "delta_neq_reduction", "delta_uo_ml_kg_h")
core <- analysis[complete.cases(analysis[, c(core_observed, "hospital_mortality_after_landmark")]), ]

covars <- c(
  "anchor_age", "gender_male", "weight_kg", "first_day_sofa",
  "sofa_cardiovascular", "sofa_renal",
  "baseline_map_mean", "baseline_neq_avg", "baseline_uo_ml_kg_h"
)
for (col in covars) {
  core[[col]] <- impute_median(core[[col]])
}

map_fit <- fit_component(core, "delta_map", covars)
vaso_fit <- fit_component(core, "delta_neq_reduction", covars)
uo_fit <- fit_component(core, "delta_uo_ml_kg_h", covars)

core$hrc_map_residual_z <- map_fit$residual_z
core$hrc_vaso_residual_z <- vaso_fit$residual_z
core$hrc_uo_residual_z <- uo_fit$residual_z
core$hrc_core_raw <- rowMeans(
  core[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")],
  na.rm = FALSE
)
core$hrc_core_z <- zscore(core$hrc_core_raw)

quartile_cut <- quantile(core$hrc_core_z, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
core$hrc_quartile <- cut(
  core$hrc_core_z,
  breaks = c(-Inf, quartile_cut, Inf),
  labels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
  include.lowest = TRUE
)
core$low_hrc_q1 <- ifelse(core$hrc_quartile == "Q1_lowest", 1, 0)

score_cols <- c(
  "subject_id", "hadm_id", "stay_id", "index_time", "landmark_time",
  "hospital_mortality_after_landmark",
  "hrc_core_z", "hrc_core_raw", "hrc_quartile", "low_hrc_q1",
  "hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z",
  "delta_map", "delta_neq_reduction", "delta_uo_ml_kg_h",
  "baseline_map_mean", "baseline_neq_avg", "baseline_uo_ml_kg_h",
  "first_day_sofa", "sofa_cardiovascular", "sofa_renal"
)
write.csv(
  core[, score_cols],
  file = file.path(outdir, "mimic_hrc_v1_scores.csv"),
  row.names = FALSE
)

quartile_rows <- do.call(
  rbind,
  lapply(levels(core$hrc_quartile), function(q) {
    x <- core[core$hrc_quartile == q, ]
    data.frame(
      hrc_quartile = q,
      n = nrow(x),
      deaths = sum(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
      mortality_rate = mean(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
      median_hrc = median(x$hrc_core_z, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  quartile_rows,
  file = file.path(outdir, "mimic_hrc_v1_quartile_mortality.csv"),
  row.names = FALSE
)

glm_covars <- c(
  "hrc_core_z", "anchor_age", "first_day_sofa",
  "baseline_neq_avg", "baseline_map_mean", "baseline_uo_ml_kg_h"
)
for (col in glm_covars) {
  core[[col]] <- impute_median(core[[col]])
}
glm_fit <- glm(
  hospital_mortality_after_landmark ~ hrc_core_z + anchor_age + first_day_sofa +
    baseline_neq_avg + baseline_map_mean + baseline_uo_ml_kg_h,
  data = core,
  family = binomial()
)
coef_table <- summary(glm_fit)$coefficients
ci_table <- confint.default(glm_fit)
or_rows <- data.frame(
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
  or_rows,
  file = file.path(outdir, "mimic_hrc_v1_logistic_or.csv"),
  row.names = FALSE
)

creatinine <- core[complete.cases(core[, c("delta_creatinine_reduction", "baseline_creatinine_mean")]), ]
if (nrow(creatinine) > 100) {
  creatinine_covars <- c(covars, "baseline_creatinine_mean")
  for (col in creatinine_covars) {
    creatinine[[col]] <- impute_median(creatinine[[col]])
  }
  creat_fit <- fit_component(creatinine, "delta_creatinine_reduction", creatinine_covars)
  creatinine$hrc_creatinine_residual_z <- creat_fit$residual_z
  creatinine$hrc_core_plus_creatinine_raw <- rowMeans(
    cbind(
      creatinine$hrc_map_residual_z,
      creatinine$hrc_vaso_residual_z,
      creatinine$hrc_uo_residual_z,
      creatinine$hrc_creatinine_residual_z
    ),
    na.rm = FALSE
  )
  creatinine$hrc_core_plus_creatinine_z <- zscore(creatinine$hrc_core_plus_creatinine_raw)
  write.csv(
    creatinine[, c(
      "subject_id", "hadm_id", "stay_id", "hrc_core_z",
      "hrc_core_plus_creatinine_z", "hrc_creatinine_residual_z",
      "delta_creatinine_reduction", "baseline_creatinine_mean",
      "hospital_mortality_after_landmark"
    )],
    file = file.path(outdir, "mimic_hrc_v1_creatinine_sensitivity_scores.csv"),
    row.names = FALSE
  )
}

summary_lines <- c(
  "# MIMIC HRC v1.0 pilot scoring",
  "",
  "This is a MIMIC-only pilot implementation of the residualized HRC construct.",
  "",
  "## Core definition",
  "",
  "- Core observed recovery variables: `delta_map`, `delta_neq_reduction`, `delta_uo_ml_kg_h`.",
  "- Expected recovery model: linear model within MIMIC, conditional on age, sex, weight, baseline SOFA domains, baseline MAP, baseline NE-equivalent, and baseline urine output.",
  "- Component HRC: observed minus expected recovery, standardized to z-scores.",
  "- Core HRC: equal-weighted mean of the three component residual z-scores, then standardized.",
  "",
  "## Counts",
  "",
  paste0("- input rows: ", nrow(analysis)),
  paste0("- core HRC rows: ", nrow(core)),
  paste0("- creatinine sensitivity rows: ", ifelse(exists("creatinine"), nrow(creatinine), 0)),
  "",
  "## Mortality by HRC quartile",
  "",
  paste(
    apply(quartile_rows, 1, function(row) {
      paste0(
        "- ", row[["hrc_quartile"]], ": n=", row[["n"]],
        ", deaths=", row[["deaths"]],
        ", mortality=", sprintf("%.3f", as.numeric(row[["mortality_rate"]]))
      )
    }),
    collapse = "\n"
  ),
  "",
  "## Adjusted logistic association",
  "",
  paste(
    apply(or_rows, 1, function(row) {
      paste0(
        "- ", row[["term"]], ": OR=", sprintf("%.3f", as.numeric(row[["odds_ratio"]])),
        " (", sprintf("%.3f", as.numeric(row[["ci_low"]])),
        "-", sprintf("%.3f", as.numeric(row[["ci_high"]])),
        "), p=", signif(as.numeric(row[["p_value"]]), 3)
      )
    }),
    collapse = "\n"
  )
)
writeLines(summary_lines, con = file.path(outdir, "mimic_hrc_v1_model_summary.md"))

cat("Input rows:", nrow(analysis), "\n")
cat("Core HRC rows:", nrow(core), "\n")
cat("Scores:", file.path(outdir, "mimic_hrc_v1_scores.csv"), "\n")
cat("Quartile mortality:", file.path(outdir, "mimic_hrc_v1_quartile_mortality.csv"), "\n")
cat("Logistic OR:", file.path(outdir, "mimic_hrc_v1_logistic_or.csv"), "\n")
cat("Summary:", file.path(outdir, "mimic_hrc_v1_model_summary.md"), "\n")
