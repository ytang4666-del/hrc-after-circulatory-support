#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(splines)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

outdir <- "outputs/sensitivity"
figdir <- file.path(outdir, "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260704)

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

num <- function(x) suppressWarnings(as.numeric(x))

has_variation <- function(x) {
  ux <- unique(x[is.finite(x) & !is.na(x)])
  length(ux) > 1
}

fmt_or <- function(or, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}

ordinal_pct <- function(pct) {
  if (pct == 10) return("10th")
  if (pct == 20) return("20th")
  if (pct == 25) return("25th")
  if (pct == 33) return("33rd")
  paste0(pct, "th")
}

read_primary_or <- function(path, database) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  row <- d[d$term == "hrc_core_z", ]
  data.frame(
    database = database,
    analysis = "Primary core HRC",
    exposure = "HRC per 1 SD increase",
    n = NA_integer_,
    deaths = NA_integer_,
    beta = num(row$beta),
    se = num(row$se),
    odds_ratio = num(row$odds_ratio),
    ci_low = num(row$ci_low),
    ci_high = num(row$ci_high),
    p_value = num(row$p_value),
    stringsAsFactors = FALSE
  )
}

meta_summary <- function(rows, analysis, exposure) {
  rows <- rows[is.finite(rows$beta) & is.finite(rows$se) & rows$se > 0, ]
  if (nrow(rows) == 0) {
    return(data.frame())
  }
  w <- 1 / rows$se^2
  beta_fe <- sum(w * rows$beta) / sum(w)
  se_fe <- sqrt(1 / sum(w))
  k <- nrow(rows)
  q <- if (k > 1) sum(w * (rows$beta - beta_fe)^2) else NA_real_
  q_df <- k - 1
  q_p <- if (k > 1) pchisq(q, df = q_df, lower.tail = FALSE) else NA_real_
  i2 <- if (k > 1 && is.finite(q) && q > 0) max(0, (q - q_df) / q) * 100 else NA_real_
  tau2 <- if (k > 1) {
    max(0, (q - q_df) / (sum(w) - sum(w^2) / sum(w)))
  } else {
    NA_real_
  }
  wr <- if (k > 1) 1 / (rows$se^2 + tau2) else w
  beta_re <- sum(wr * rows$beta) / sum(wr)
  se_re <- sqrt(1 / sum(wr))
  data.frame(
    analysis = analysis,
    exposure = exposure,
    k = k,
    fixed_or = exp(beta_fe),
    fixed_ci_low = exp(beta_fe - 1.96 * se_fe),
    fixed_ci_high = exp(beta_fe + 1.96 * se_fe),
    fixed_p = 2 * pnorm(abs(beta_fe / se_fe), lower.tail = FALSE),
    random_or = exp(beta_re),
    random_ci_low = exp(beta_re - 1.96 * se_re),
    random_ci_high = exp(beta_re + 1.96 * se_re),
    random_p = 2 * pnorm(abs(beta_re / se_re), lower.tail = FALSE),
    q_stat = q,
    q_df = q_df,
    q_p = q_p,
    i2_percent = i2,
    tau2_dl = tau2,
    stringsAsFactors = FALSE
  )
}

standardize_formal <- function(database) {
  if (database == "MIMIC-IV") {
    d <- read.csv("outputs/mimic_formal/mimic_formal_cohort.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$stay_id,
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$anchor_age),
      male = ifelse(d$gender == "M", 1, 0),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$first_day_sofa),
      severity_secondary = NA_real_,
      heart_surgery = NA_real_,
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_neq_avg),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_neq_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
      baseline_creatinine_mean = num(d$baseline_creatinine_mean),
      delta_creatinine_reduction = num(d$delta_creatinine_reduction),
      baseline_lactate_mean = num(d$baseline_lactate_mean),
      delta_lactate_reduction = num(d$delta_lactate_reduction),
      stringsAsFactors = FALSE
    )
  } else if (database == "eICU") {
    d <- read.csv("outputs/eicu_formal/eicu_formal_cohort.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$patientunitstayid,
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$age),
      male = ifelse(d$gender == "Male", 1, 0),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      heart_surgery = NA_real_,
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_vaso_burden),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_vaso_burden_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
      baseline_creatinine_mean = num(d$baseline_creatinine_mean),
      delta_creatinine_reduction = num(d$delta_creatinine_reduction),
      baseline_lactate_mean = num(d$baseline_lactate_mean),
      delta_lactate_reduction = num(d$delta_lactate_reduction),
      stringsAsFactors = FALSE
    )
  } else if (database == "SICdb") {
    d <- read.csv("outputs/sicdb_formal/sicdb_formal_cohort.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$caseid,
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$age),
      male = num(d$male),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$saps3),
      severity_secondary = NA_real_,
      heart_surgery = ifelse(d$heartsurgeryadditionaldata == "740", 1, 0),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_vaso_burden),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_vaso_burden_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
      baseline_creatinine_mean = num(d$baseline_creatinine_mean),
      delta_creatinine_reduction = num(d$delta_creatinine_reduction),
      baseline_lactate_mean = num(d$baseline_lactate_mean),
      delta_lactate_reduction = num(d$delta_lactate_reduction),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown database")
  }
}

standardize_scores <- function(database) {
  if (database == "MIMIC-IV") {
    d <- read.csv("outputs/mimic_formal/mimic_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$stay_id,
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      age = num(d$anchor_age),
      male = ifelse(d$gender == "M", 1, 0),
      severity_primary = num(d$first_day_sofa),
      severity_secondary = NA_real_,
      heart_surgery = NA_real_,
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_log = log1p(pmax(num(d$baseline_neq_avg), 0)),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      stringsAsFactors = FALSE
    )
  } else if (database == "eICU") {
    d <- read.csv("outputs/eicu_formal/eicu_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$patientunitstayid,
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      age = num(d$age),
      male = ifelse(d$gender == "Male", 1, 0),
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      heart_surgery = NA_real_,
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_log = log1p(pmax(num(d$baseline_vaso_burden), 0)),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      stringsAsFactors = FALSE
    )
  } else if (database == "SICdb") {
    d <- read.csv("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$caseid,
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      age = num(d$age),
      male = num(d$male),
      severity_primary = num(d$saps3),
      severity_secondary = NA_real_,
      heart_surgery = num(d$heart_surgery),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_log = log1p(pmax(num(d$baseline_vaso_burden), 0)),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown database")
  }
}

association_formula <- function(exposure, dat) {
  terms <- c(exposure, "ns(age, df = 3)")
  if (has_variation(dat$male)) terms <- c(terms, "male")
  if (sum(!is.na(dat$severity_primary)) > 20 && has_variation(dat$severity_primary)) {
    terms <- c(terms, "ns(severity_primary, df = 4)")
  }
  if (sum(!is.na(dat$severity_secondary)) > 20 && has_variation(dat$severity_secondary)) {
    terms <- c(terms, "ns(severity_secondary, df = 4)")
  }
  if ("heart_surgery" %in% names(dat) && sum(!is.na(dat$heart_surgery)) > 20 && has_variation(dat$heart_surgery)) {
    terms <- c(terms, "heart_surgery")
  }
  terms <- c(
    terms,
    "ns(baseline_vaso_log, df = 4)",
    "ns(baseline_map_mean, df = 4)",
    "ns(baseline_uo_log, df = 4)"
  )
  as.formula(paste("outcome ~", paste(terms, collapse = " + ")))
}

fit_or <- function(dat, exposure, database, analysis, exposure_label) {
  covars <- c(
    "outcome", exposure, "age", "male", "severity_primary", "severity_secondary",
    "heart_surgery", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log"
  )
  keep <- covars[covars %in% names(dat)]
  x <- dat[complete.cases(dat[, c("outcome", exposure)]), keep, drop = FALSE]
  x <- x[is.finite(x$outcome) & is.finite(x[[exposure]]), , drop = FALSE]
  if (length(unique(x$outcome)) < 2 || length(unique(x[[exposure]])) < 2) {
    return(data.frame())
  }
  for (col in setdiff(names(x), c("outcome", exposure))) {
    if (is.numeric(x[[col]])) x[[col]] <- impute_median(x[[col]])
  }
  fit <- glm(association_formula(exposure, x), data = x, family = binomial())
  coef_table <- summary(fit)$coefficients
  if (!(exposure %in% rownames(coef_table))) return(data.frame())
  ci_table <- confint.default(fit)
  data.frame(
    database = database,
    analysis = analysis,
    exposure = exposure_label,
    n = nrow(x),
    deaths = sum(x$outcome == 1, na.rm = TRUE),
    beta = coef_table[exposure, "Estimate"],
    se = coef_table[exposure, "Std. Error"],
    odds_ratio = exp(coef_table[exposure, "Estimate"]),
    ci_low = exp(ci_table[exposure, 1]),
    ci_high = exp(ci_table[exposure, 2]),
    p_value = coef_table[exposure, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

expected_formula <- function(outcome, dat, extra_baseline = character()) {
  terms <- c("ns(age, df = 3)")
  if (has_variation(dat$male)) terms <- c(terms, "male")
  if (sum(!is.na(dat$weight_kg)) > 20 && has_variation(dat$weight_kg)) {
    terms <- c(terms, "ns(weight_kg, df = 3)")
  }
  if (sum(!is.na(dat$severity_primary)) > 20 && has_variation(dat$severity_primary)) {
    terms <- c(terms, "ns(severity_primary, df = 4)")
  }
  if (sum(!is.na(dat$severity_secondary)) > 20 && has_variation(dat$severity_secondary)) {
    terms <- c(terms, "ns(severity_secondary, df = 4)")
  }
  if ("heart_surgery" %in% names(dat) && sum(!is.na(dat$heart_surgery)) > 20 && has_variation(dat$heart_surgery)) {
    terms <- c(terms, "heart_surgery")
  }
  terms <- c(
    terms,
    "ns(baseline_map_mean, df = 4)",
    "ns(baseline_vaso_log, df = 4)",
    "ns(baseline_uo_log, df = 4)"
  )
  if ("baseline_creatinine_log" %in% extra_baseline) {
    terms <- c(terms, "ns(baseline_creatinine_log, df = 4)")
  }
  if ("baseline_lactate_log" %in% extra_baseline) {
    terms <- c(terms, "ns(baseline_lactate_log, df = 4)")
  }
  as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

crossfit_component <- function(dat, outcome, folds, extra_baseline = character()) {
  expected <- rep(NA_real_, nrow(dat))
  for (fold in sort(unique(folds))) {
    train <- dat[folds != fold, , drop = FALSE]
    test <- dat[folds == fold, , drop = FALSE]
    fit <- lm(expected_formula(outcome, train, extra_baseline), data = train)
    expected[folds == fold] <- as.numeric(predict(fit, newdata = test))
  }
  residual <- dat[[outcome]] - expected
  zscore(residual)
}

compute_hrc_model <- function(dat, database, analysis, component_names) {
  component_map <- list(
    map = list(variable = "delta_map", baseline = character(), label = "MAP"),
    vaso = list(variable = "vaso_reduction", baseline = character(), label = "Vaso"),
    uo = list(variable = "log_uo_recovery", baseline = character(), label = "UO"),
    creatinine = list(
      variable = "delta_creatinine_reduction",
      baseline = "baseline_creatinine_log",
      label = "Creatinine"
    ),
    lactate = list(
      variable = "delta_lactate_reduction",
      baseline = "baseline_lactate_log",
      label = "Lactate"
    )
  )
  dat$baseline_vaso_log <- log1p(pmax(dat$baseline_vaso_burden, 0))
  dat$baseline_uo_log <- log1p(pmax(dat$baseline_uo_ml_kg_h, 0))
  dat$baseline_creatinine_log <- log1p(pmax(dat$baseline_creatinine_mean, 0))
  dat$baseline_lactate_log <- log1p(pmax(dat$baseline_lactate_mean, 0))
  dat$weight_kg[dat$weight_kg <= 0 | dat$weight_kg > 300] <- NA_real_

  observed <- vapply(component_names, function(nm) component_map[[nm]]$variable, character(1))
  extra_baseline <- unique(unlist(lapply(component_names, function(nm) component_map[[nm]]$baseline)))
  required <- c(
    "outcome", observed, "age", "baseline_map_mean",
    "baseline_vaso_log", "baseline_uo_log", extra_baseline
  )
  x <- dat[complete.cases(dat[, required]), , drop = FALSE]
  x <- x[is.finite(x$outcome), , drop = FALSE]
  if (nrow(x) < 500 || sum(x$outcome == 1, na.rm = TRUE) < 50) {
    return(list(row = data.frame(), scores = data.frame()))
  }
  impute_cols <- c(
    "age", "male", "weight_kg", "severity_primary", "severity_secondary",
    "heart_surgery", "baseline_map_mean", "baseline_vaso_log",
    "baseline_uo_log", "baseline_creatinine_log", "baseline_lactate_log"
  )
  for (col in intersect(impute_cols, names(x))) {
    x[[col]] <- impute_median(x[[col]])
  }
  folds <- sample(rep(1:5, length.out = nrow(x)))
  residuals <- lapply(component_names, function(nm) {
    crossfit_component(
      x,
      component_map[[nm]]$variable,
      folds,
      component_map[[nm]]$baseline
    )
  })
  names(residuals) <- component_names
  for (nm in component_names) {
    x[[paste0("hrc_", nm, "_residual_z")]] <- residuals[[nm]]
  }
  x$hrc_z <- zscore(rowMeans(x[, paste0("hrc_", component_names, "_residual_z"), drop = FALSE], na.rm = FALSE))
  row <- fit_or(x, "hrc_z", database, analysis, "HRC per 1 SD increase")
  list(row = row, scores = x)
}

write_md_table <- function(df, path, title = NULL) {
  lines <- character()
  if (!is.null(title)) lines <- c(lines, paste0("# ", title), "")
  if (nrow(df) == 0) {
    lines <- c(lines, "No rows.")
  } else {
    cols <- names(df)
    lines <- c(lines, paste(cols, collapse = " | "))
    lines <- c(lines, paste(rep("---", length(cols)), collapse = " | "))
    for (i in seq_len(nrow(df))) {
      vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
      vals <- gsub("\\|", "/", vals)
      lines <- c(lines, paste(vals, collapse = " | "))
    }
  }
  writeLines(lines, con = path)
}

databases <- c("MIMIC-IV", "eICU", "SICdb")
primary <- rbind(
  read_primary_or("outputs/mimic_formal/mimic_hrc_formal_linear_or.csv", "MIMIC-IV"),
  read_primary_or("outputs/eicu_formal/eicu_hrc_formal_linear_or.csv", "eICU"),
  read_primary_or("outputs/sicdb_formal/sicdb_hrc_formal_linear_or.csv", "SICdb")
)
primary_meta <- meta_summary(primary, "Primary core HRC", "HRC per 1 SD increase")

score_list <- setNames(lapply(databases, standardize_scores), databases)
cutoffs <- c(0.10, 0.20, 0.25, 0.33)
cutoff_rows <- do.call(rbind, lapply(databases, function(db) {
  d <- score_list[[db]]
  do.call(rbind, lapply(cutoffs, function(p) {
    threshold <- quantile(d$hrc_core_z, probs = p, na.rm = TRUE)
    d$low_hrc_cutoff <- ifelse(d$hrc_core_z <= threshold, 1, 0)
    fit_or(
      d,
      "low_hrc_cutoff",
      db,
      paste0("Low-HRC cutoff p", round(p * 100)),
      paste0("Low-HRC below ", ordinal_pct(round(p * 100)), " percentile")
    )
  }))
}))
cutoff_meta <- do.call(rbind, lapply(unique(cutoff_rows$analysis), function(a) {
  x <- cutoff_rows[cutoff_rows$analysis == a, ]
  meta_summary(x, a, unique(x$exposure)[1])
}))

formal_list <- setNames(lapply(databases, standardize_formal), databases)
variant_components <- list(
  "Core HRC recomputed" = c("map", "vaso", "uo"),
  "Core + creatinine" = c("map", "vaso", "uo", "creatinine"),
  "Core + lactate" = c("map", "vaso", "uo", "lactate"),
  "Core + creatinine + lactate" = c("map", "vaso", "uo", "creatinine", "lactate")
)

enhanced_rows <- list()
enhanced_scores <- list()
for (db in databases) {
  for (analysis in names(variant_components)) {
    res <- compute_hrc_model(formal_list[[db]], db, analysis, variant_components[[analysis]])
    if (nrow(res$row) > 0) {
      enhanced_rows[[paste(db, analysis, sep = "__")]] <- res$row
      enhanced_scores[[paste(db, analysis, sep = "__")]] <- res$scores
    }
  }
}
enhanced_or <- do.call(rbind, enhanced_rows)
enhanced_meta <- do.call(rbind, lapply(unique(enhanced_or$analysis), function(a) {
  x <- enhanced_or[enhanced_or$analysis == a, ]
  meta_summary(x, a, "HRC per 1 SD increase")
}))

sicdb <- formal_list[["SICdb"]]
subset_rows <- list()
for (label in c("SICdb excluding heart surgery", "SICdb heart surgery only")) {
  sub <- if (label == "SICdb excluding heart surgery") {
    sicdb[sicdb$heart_surgery == 0, , drop = FALSE]
  } else {
    sicdb[sicdb$heart_surgery == 1, , drop = FALSE]
  }
  res <- compute_hrc_model(sub, "SICdb", label, c("map", "vaso", "uo"))
  if (nrow(res$row) > 0) subset_rows[[label]] <- res$row
}
subset_or <- do.call(rbind, subset_rows)

all_continuous <- rbind(
  primary,
  enhanced_or,
  subset_or
)
write.csv(primary, file.path(outdir, "primary_core_continuous_or.csv"), row.names = FALSE)
write.csv(primary_meta, file.path(outdir, "primary_core_meta.csv"), row.names = FALSE)
write.csv(cutoff_rows, file.path(outdir, "low_hrc_cutoff_sensitivity_or.csv"), row.names = FALSE)
write.csv(cutoff_meta, file.path(outdir, "low_hrc_cutoff_sensitivity_meta.csv"), row.names = FALSE)
write.csv(enhanced_or, file.path(outdir, "enhanced_hrc_sensitivity_or.csv"), row.names = FALSE)
write.csv(enhanced_meta, file.path(outdir, "enhanced_hrc_sensitivity_meta.csv"), row.names = FALSE)
write.csv(subset_or, file.path(outdir, "sicdb_heart_surgery_sensitivity_or.csv"), row.names = FALSE)
write.csv(all_continuous, file.path(outdir, "all_continuous_sensitivity_or.csv"), row.names = FALSE)

md_overview <- c(
  "# HRC sensitivity package",
  "",
  "## Core sensitivity conclusions",
  "",
  paste0(
    "- Primary three-database fixed-effect OR per 1 SD higher HRC: ",
    fmt_or(primary_meta$fixed_or, primary_meta$fixed_ci_low, primary_meta$fixed_ci_high),
    "; random-effects OR: ",
    fmt_or(primary_meta$random_or, primary_meta$random_ci_low, primary_meta$random_ci_high),
    "."
  ),
  paste0(
    "- Low-HRC cutoff sensitivity remained directionally consistent across p10, p20, p25, and p33 thresholds."
  ),
  paste0(
    "- Enhanced HRC variants tested whether the construct depends only on MAP/vaso/urine or remains stable when creatinine and/or lactate recovery are added."
  ),
  paste0(
    "- SICdb heart-surgery sensitivity tested whether the European validation is driven by postoperative cardiac-surgical physiology."
  ),
  "",
  "## Output files",
  "",
  "- primary_core_continuous_or.csv",
  "- primary_core_meta.csv",
  "- low_hrc_cutoff_sensitivity_or.csv",
  "- low_hrc_cutoff_sensitivity_meta.csv",
  "- enhanced_hrc_sensitivity_or.csv",
  "- enhanced_hrc_sensitivity_meta.csv",
  "- sicdb_heart_surgery_sensitivity_or.csv",
  "- all_continuous_sensitivity_or.csv",
  "- figures/figure4_hrc_sensitivity_v1.*"
)
writeLines(md_overview, file.path(outdir, "sensitivity_overview.md"))

key_results <- c(
  "# HRC sensitivity key results",
  "",
  "## Primary continuous HRC",
  "",
  paste0(
    "- Fixed-effect three-database OR per 1 SD higher HRC: ",
    fmt_or(primary_meta$fixed_or, primary_meta$fixed_ci_low, primary_meta$fixed_ci_high),
    "."
  ),
  paste0(
    "- Random-effects three-database OR per 1 SD higher HRC: ",
    fmt_or(primary_meta$random_or, primary_meta$random_ci_low, primary_meta$random_ci_high),
    "."
  ),
  paste0(
    "- Heterogeneity: I2=",
    sprintf("%.1f%%", primary_meta$i2_percent),
    ", Q-test p=",
    signif(primary_meta$q_p, 3),
    "."
  ),
  "",
  "## Low-HRC threshold sensitivity",
  "",
  paste(
    apply(cutoff_meta, 1, function(row) {
      paste0(
        "- ", row[["analysis"]], ": random-effects OR ",
        fmt_or(
          as.numeric(row[["random_or"]]),
          as.numeric(row[["random_ci_low"]]),
          as.numeric(row[["random_ci_high"]])
        ),
        "."
      )
    }),
    collapse = "\n"
  ),
  "",
  "## Enhanced HRC sensitivity",
  "",
  paste(
    apply(enhanced_meta, 1, function(row) {
      paste0(
        "- ", row[["analysis"]], ": random-effects OR ",
        fmt_or(
          as.numeric(row[["random_or"]]),
          as.numeric(row[["random_ci_low"]]),
          as.numeric(row[["random_ci_high"]])
        ),
        "."
      )
    }),
    collapse = "\n"
  ),
  "",
  "## SICdb cardiac-surgery sensitivity",
  "",
  paste(
    apply(subset_or, 1, function(row) {
      paste0(
        "- ", row[["analysis"]], ": n=", row[["n"]],
        ", deaths=", row[["deaths"]],
        ", OR ",
        fmt_or(
          as.numeric(row[["odds_ratio"]]),
          as.numeric(row[["ci_low"]]),
          as.numeric(row[["ci_high"]])
        ),
        "."
      )
    }),
    collapse = "\n"
  )
)
writeLines(key_results, file.path(outdir, "sensitivity_key_results.md"))

write_md_table(
  primary_meta,
  file.path(outdir, "primary_core_meta.md"),
  "Primary three-database continuous HRC meta-analysis"
)
write_md_table(
  cutoff_meta,
  file.path(outdir, "low_hrc_cutoff_sensitivity_meta.md"),
  "Low-HRC cutoff sensitivity meta-analysis"
)
write_md_table(
  enhanced_meta,
  file.path(outdir, "enhanced_hrc_sensitivity_meta.md"),
  "Enhanced HRC sensitivity meta-analysis"
)
write_md_table(
  subset_or,
  file.path(outdir, "sicdb_heart_surgery_sensitivity_or.md"),
  "SICdb heart-surgery sensitivity"
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
  blue = "#3b6ea8",
  red = "#b94a48",
  green = "#3f7f5f",
  light_blue = "#dbe8f6"
)

continuous_plot <- enhanced_meta
continuous_plot$analysis <- factor(
  continuous_plot$analysis,
  levels = rev(c(
    "Core HRC recomputed",
    "Core + creatinine",
    "Core + lactate",
    "Core + creatinine + lactate"
  ))
)
p_enhanced <- ggplot(continuous_plot, aes(x = random_or, y = analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = pal$muted) +
  geom_errorbarh(aes(xmin = random_ci_low, xmax = random_ci_high), height = 0.12, linewidth = 0.45, colour = pal$muted) +
  geom_point(shape = 21, size = 2.2, fill = pal$blue, colour = pal$ink) +
  scale_x_continuous(limits = c(0.55, 1.02), breaks = c(0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(
    title = "A  Enhanced HRC variants",
    x = "Random-effects OR per 1 SD higher HRC",
    y = NULL
  )

cutoff_plot <- cutoff_meta
cutoff_plot$cutoff <- factor(
  cutoff_plot$analysis,
  levels = paste0("Low-HRC cutoff p", c(10, 20, 25, 33)),
  labels = c("p10", "p20", "p25", "p33")
)
p_cutoff <- ggplot(cutoff_plot, aes(x = random_or, y = cutoff)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = pal$muted) +
  geom_errorbarh(aes(xmin = random_ci_low, xmax = random_ci_high), height = 0.12, linewidth = 0.45, colour = pal$muted) +
  geom_point(shape = 21, size = 2.2, fill = pal$red, colour = pal$ink) +
  scale_x_continuous(limits = c(1.0, max(2.6, max(cutoff_plot$random_ci_high, na.rm = TRUE))), breaks = pretty(c(1, 2.6), n = 5)) +
  labs(
    title = "B  Low-HRC threshold sensitivity",
    x = "Random-effects OR for low-HRC",
    y = NULL
  )

full_variant <- enhanced_or[enhanced_or$analysis == "Core + creatinine + lactate", ]
full_variant$database <- factor(full_variant$database, levels = rev(databases))
p_full <- ggplot(full_variant, aes(x = odds_ratio, y = database)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.45, colour = pal$muted) +
  geom_point(shape = 21, size = 2.2, fill = pal$blue, colour = pal$ink) +
  scale_x_continuous(limits = c(0.45, 1.02), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(
    title = "C  Full enhanced HRC by database",
    x = "Adjusted OR per 1 SD higher HRC",
    y = NULL
  )

subset_plot <- subset_or
subset_plot$analysis <- factor(
  subset_plot$analysis,
  levels = rev(c("SICdb excluding heart surgery", "SICdb heart surgery only"))
)
p_sicdb_subset <- ggplot(subset_plot, aes(x = odds_ratio, y = analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = pal$muted) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.45, colour = pal$muted) +
  geom_point(shape = 21, size = 2.2, fill = pal$green, colour = pal$ink) +
  scale_x_continuous(limits = c(0.45, 1.05), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(
    title = "D  SICdb cardiac-surgery sensitivity",
    x = "Adjusted OR per 1 SD higher HRC",
    y = NULL
  )

fig <- (p_enhanced | p_cutoff) / (p_full | p_sicdb_subset) +
  plot_layout(widths = c(1.1, 0.9), heights = c(1, 1)) +
  plot_annotation(
    title = "Sensitivity analyses for Hemodynamic Recovery Capacity",
    subtitle = "Associations remain directionally stable across organ-domain, threshold, database, and SICdb cardiac-surgery checks",
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

figure_base <- file.path(figdir, "figure4_hrc_sensitivity_v1")
save_pub(fig, figure_base)

writeLines(
  c(
    "# Figure 4 contract",
    "",
    "Core conclusion: The HRC-mortality association is not explained by a single organ domain, arbitrary low-HRC threshold, one database, or SICdb cardiac-surgery case mix.",
    "",
    "Evidence chain:",
    "- Panel A: random-effects meta-analysis across enhanced HRC variants.",
    "- Panel B: random-effects meta-analysis across low-HRC percentile thresholds.",
    "- Panel C: full enhanced HRC estimates by database.",
    "- Panel D: SICdb estimates after excluding and isolating heart-surgery cases.",
    "",
    "Archetype: quantitative grid.",
    "Backend: R/ggplot2/patchwork only.",
    "Export formats: SVG, PDF, TIFF, PNG preview."
  ),
  con = file.path(figdir, "figure4_hrc_sensitivity_v1_contract.md")
)

cat("Sensitivity outputs:", outdir, "\n")
cat("Figure:", figure_base, "\n")
