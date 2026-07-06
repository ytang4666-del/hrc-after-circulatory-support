#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

outdir <- "outputs/methodology_stress_tests"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260704)

num <- function(x) suppressWarnings(as.numeric(x))

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

impute_median <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

has_variation <- function(x) {
  ux <- unique(x[is.finite(x) & !is.na(x)])
  length(ux) > 1
}

fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)

auc_rank <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calibration_slope <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  lp <- qlogis(p)
  fit <- try(glm(y ~ lp, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  unname(coef(fit)["lp"])
}

model_metrics <- function(y, p, fit = NULL) {
  data.frame(
    auc = auc_rank(y, p),
    brier = mean((y - p)^2, na.rm = TRUE),
    calibration_slope = calibration_slope(y, p),
    log_likelihood = if (is.null(fit)) NA_real_ else as.numeric(logLik(fit)),
    stringsAsFactors = FALSE
  )
}

e_value_from_estimate <- function(x) {
  if (!is.finite(x) || x < 1) return(NA_real_)
  x + sqrt(x * (x - 1))
}

e_value_or <- function(or, lo, hi) {
  if (or < 1) {
    est <- 1 / or
    bound <- if (hi < 1) 1 / hi else 1
  } else {
    est <- or
    bound <- if (lo > 1) lo else 1
  }
  c(e_value = e_value_from_estimate(est), e_value_ci = e_value_from_estimate(bound))
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
  writeLines(lines, path)
}

association_terms <- function(dat, exposure = NULL, include_hrc = TRUE) {
  terms <- character()
  if (include_hrc && !is.null(exposure)) terms <- c(terms, exposure)
  terms <- c(terms, "ns(age, df = 3)")
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  if ("severity_primary" %in% names(dat) && has_variation(dat$severity_primary)) {
    terms <- c(terms, "ns(severity_primary, df = 4)")
  }
  if ("severity_secondary" %in% names(dat) && has_variation(dat$severity_secondary)) {
    terms <- c(terms, "ns(severity_secondary, df = 4)")
  }
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) {
    terms <- c(terms, "heart_surgery")
  }
  terms <- c(terms, "ns(baseline_vaso_log, df = 4)", "ns(baseline_map_mean, df = 4)", "ns(baseline_uo_log, df = 4)")
  terms
}

association_formula <- function(dat, exposure = "hrc_core_z", include_hrc = TRUE) {
  as.formula(paste("outcome ~", paste(association_terms(dat, exposure, include_hrc), collapse = " + ")))
}

prepare_model_data <- function(dat, exposure = "hrc_core_z") {
  cols <- unique(c(
    "outcome", exposure, "age", "male", "severity_primary", "severity_secondary",
    "heart_surgery", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log"
  ))
  cols <- cols[cols %in% names(dat)]
  x <- dat[complete.cases(dat[, c("outcome", exposure)]), cols, drop = FALSE]
  x <- x[is.finite(x$outcome) & is.finite(x[[exposure]]), , drop = FALSE]
  for (col in setdiff(names(x), c("outcome", exposure))) {
    if (is.numeric(x[[col]])) x[[col]] <- impute_median(x[[col]])
  }
  x
}

fit_or <- function(dat, exposure, database, analysis, exposure_label, weights = NULL) {
  x <- prepare_model_data(dat, exposure)
  if (nrow(x) == 0 || length(unique(x$outcome)) < 2 || length(unique(x[[exposure]])) < 2) {
    return(data.frame())
  }
  if (!is.null(weights)) {
    x$ipw <- weights[match(rownames(x), names(weights))]
    x$ipw <- impute_median(x$ipw)
    fit <- suppressWarnings(glm(association_formula(x, exposure, TRUE), data = x, weights = ipw, family = quasibinomial()))
  } else {
    fit <- glm(association_formula(x, exposure, TRUE), data = x, family = binomial())
  }
  tab <- summary(fit)$coefficients
  if (!(exposure %in% rownames(tab))) return(data.frame())
  ci <- confint.default(fit)
  p_col <- if ("Pr(>|t|)" %in% colnames(tab)) "Pr(>|t|)" else "Pr(>|z|)"
  data.frame(
    database = database,
    analysis = analysis,
    exposure = exposure_label,
    n = nrow(x),
    deaths = sum(x$outcome == 1),
    beta = tab[exposure, "Estimate"],
    se = tab[exposure, "Std. Error"],
    odds_ratio = exp(tab[exposure, "Estimate"]),
    ci_low = exp(ci[exposure, 1]),
    ci_high = exp(ci[exposure, 2]),
    p_value = tab[exposure, p_col],
    stringsAsFactors = FALSE
  )
}

standardize_scores <- function(database) {
  if (database == "MIMIC-IV") {
    d <- read.csv("outputs/mimic_formal/mimic_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$stay_id,
      cluster = "MIMIC",
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      hrc_map_residual_z = num(d$hrc_map_residual_z),
      hrc_vaso_residual_z = num(d$hrc_vaso_residual_z),
      hrc_uo_residual_z = num(d$hrc_uo_residual_z),
      age = num(d$anchor_age),
      male = ifelse(d$gender == "M", 1, 0),
      severity_primary = num(d$first_day_sofa),
      severity_secondary = NA_real_,
      heart_surgery = NA_real_,
      baseline_vaso_log = log1p(pmax(num(d$baseline_neq_avg), 0)),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      cs_or_mcs_subgroup = num(d$cs_or_mcs_subgroup),
      mcs_any_icu = num(d$mcs_any_icu),
      mcs_any_0_24h = num(d$mcs_any_0_24h),
      iabp_any_icu = num(d$iabp_any_icu),
      impella_any_icu = num(d$impella_any_icu),
      ecmo_any_icu = num(d$ecmo_any_icu),
      stringsAsFactors = FALSE
    )
  } else if (database == "eICU") {
    d <- read.csv("outputs/eicu_formal/eicu_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$patientunitstayid,
      cluster = d$hospitalid,
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      hrc_map_residual_z = num(d$hrc_map_residual_z),
      hrc_vaso_residual_z = num(d$hrc_vaso_residual_z),
      hrc_uo_residual_z = num(d$hrc_uo_residual_z),
      age = num(d$age),
      male = ifelse(d$gender == "Male", 1, 0),
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      heart_surgery = NA_real_,
      baseline_vaso_log = log1p(pmax(num(d$baseline_vaso_burden), 0)),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      stringsAsFactors = FALSE
    )
  } else if (database == "SICdb") {
    d <- read.csv("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    data.frame(
      database = database,
      id = d$caseid,
      cluster = d$hospitalunit,
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      hrc_map_residual_z = num(d$hrc_map_residual_z),
      hrc_vaso_residual_z = num(d$hrc_vaso_residual_z),
      hrc_uo_residual_z = num(d$hrc_uo_residual_z),
      age = num(d$age),
      male = num(d$male),
      severity_primary = num(d$saps3),
      severity_secondary = NA_real_,
      heart_surgery = num(d$heart_surgery),
      baseline_vaso_log = log1p(pmax(num(d$baseline_vaso_burden), 0)),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown database")
  }
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
      severity_primary = num(d$first_day_sofa),
      severity_secondary = NA_real_,
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_neq_avg),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_neq_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
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
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_vaso_burden),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_vaso_burden_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
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
      severity_primary = num(d$saps3),
      severity_secondary = NA_real_,
      heart_surgery = ifelse(d$heartsurgeryadditionaldata == "740", 1, 0),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_vaso_burden = num(d$baseline_vaso_burden),
      baseline_uo_ml_kg_h = num(d$baseline_uo_ml_kg_h),
      delta_map = num(d$delta_map),
      vaso_reduction = num(d$log_vaso_burden_reduction),
      log_uo_recovery = num(d$log_uo_recovery),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown database")
  }
}

cluster_vcov_glm <- function(fit, cluster) {
  x <- model.matrix(fit)
  y <- fit$y
  mu <- fitted(fit)
  prior <- fit$prior.weights
  if (is.null(prior)) prior <- rep(1, length(y))
  ok <- !is.na(cluster)
  x <- x[ok, , drop = FALSE]
  y <- y[ok]
  mu <- mu[ok]
  prior <- prior[ok]
  cluster <- as.factor(cluster[ok])
  w <- as.vector(prior * mu * (1 - mu))
  bread <- try(solve(crossprod(x, x * w)), silent = TRUE)
  if (inherits(bread, "try-error")) return(NULL)
  score <- x * as.vector(prior * (y - mu))
  u <- rowsum(score, cluster)
  meat <- crossprod(u)
  g <- nrow(u)
  n <- nrow(x)
  k <- ncol(x)
  correction <- if (g > 1 && n > k) (g / (g - 1)) * ((n - 1) / (n - k)) else 1
  correction * bread %*% meat %*% bread
}

fit_cluster_or <- function(dat, database, cluster_var = "cluster") {
  x <- prepare_model_data(dat, "hrc_core_z")
  cluster <- dat[[cluster_var]][as.integer(rownames(x))]
  fit <- glm(association_formula(x, "hrc_core_z", TRUE), data = x, family = binomial())
  vc <- cluster_vcov_glm(fit, cluster)
  tab <- summary(fit)$coefficients
  beta <- tab["hrc_core_z", "Estimate"]
  naive_se <- tab["hrc_core_z", "Std. Error"]
  robust_se <- if (is.null(vc)) NA_real_ else sqrt(diag(vc))["hrc_core_z"]
  data.frame(
    database = database,
    n = nrow(x),
    deaths = sum(x$outcome == 1),
    clusters = length(unique(cluster[!is.na(cluster)])),
    beta = beta,
    naive_or = exp(beta),
    naive_ci_low = exp(beta - 1.96 * naive_se),
    naive_ci_high = exp(beta + 1.96 * naive_se),
    cluster_robust_se = robust_se,
    cluster_robust_or = exp(beta),
    cluster_robust_ci_low = exp(beta - 1.96 * robust_se),
    cluster_robust_ci_high = exp(beta + 1.96 * robust_se),
    cluster_robust_p = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE),
    note = ifelse(length(unique(cluster[!is.na(cluster)])) < 20, "Few clusters; interpret cluster-robust SE cautiously.", ""),
    stringsAsFactors = FALSE
  )
}

cv_incremental_metrics <- function(dat, database, folds_n = 5) {
  x <- prepare_model_data(dat, "hrc_core_z")
  folds <- sample(rep(seq_len(folds_n), length.out = nrow(x)))
  p_base <- rep(NA_real_, nrow(x))
  p_hrc <- rep(NA_real_, nrow(x))
  for (fold in seq_len(folds_n)) {
    train <- x[folds != fold, , drop = FALSE]
    test <- x[folds == fold, , drop = FALSE]
    base_fit <- glm(association_formula(train, "hrc_core_z", FALSE), data = train, family = binomial())
    hrc_fit <- glm(association_formula(train, "hrc_core_z", TRUE), data = train, family = binomial())
    p_base[folds == fold] <- predict(base_fit, newdata = test, type = "response")
    p_hrc[folds == fold] <- predict(hrc_fit, newdata = test, type = "response")
  }
  base_full <- glm(association_formula(x, "hrc_core_z", FALSE), data = x, family = binomial())
  hrc_full <- glm(association_formula(x, "hrc_core_z", TRUE), data = x, family = binomial())
  lrt <- anova(base_full, hrc_full, test = "Chisq")
  base_m <- model_metrics(x$outcome, p_base, base_full)
  hrc_m <- model_metrics(x$outcome, p_hrc, hrc_full)
  data.frame(
    database = database,
    n = nrow(x),
    deaths = sum(x$outcome == 1),
    baseline_auc = base_m$auc,
    hrc_auc = hrc_m$auc,
    delta_auc = hrc_m$auc - base_m$auc,
    baseline_brier = base_m$brier,
    hrc_brier = hrc_m$brier,
    delta_brier = hrc_m$brier - base_m$brier,
    baseline_calibration_slope = base_m$calibration_slope,
    hrc_calibration_slope = hrc_m$calibration_slope,
    lrt_chisq = lrt$Deviance[2],
    lrt_p = lrt$`Pr(>Chi)`[2],
    stringsAsFactors = FALSE
  )
}

smd_one <- function(x, group) {
  x <- num(x)
  group <- as.integer(group)
  a <- x[group == 1]
  b <- x[group == 0]
  m1 <- mean(a, na.rm = TRUE)
  m0 <- mean(b, na.rm = TRUE)
  v1 <- var(a, na.rm = TRUE)
  v0 <- var(b, na.rm = TRUE)
  denom <- sqrt((v1 + v0) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (m1 - m0) / denom
}

complete_case_diagnostics <- function(formal, database) {
  formal$complete_core <- complete.cases(formal[, c("delta_map", "vaso_reduction", "log_uo_recovery", "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h")])
  vars <- c("outcome", "age", "male", "severity_primary", "severity_secondary", "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h")
  rows <- do.call(rbind, lapply(vars[vars %in% names(formal)], function(v) {
    data.frame(
      database = database,
      variable = v,
      complete_mean = mean(num(formal[[v]])[formal$complete_core], na.rm = TRUE),
      incomplete_mean = mean(num(formal[[v]])[!formal$complete_core], na.rm = TRUE),
      smd_complete_vs_incomplete = smd_one(formal[[v]], formal$complete_core),
      stringsAsFactors = FALSE
    )
  }))
  rows
}

ipw_complete_case <- function(formal, scores, database) {
  formal$complete_core <- as.integer(complete.cases(formal[, c("delta_map", "vaso_reduction", "log_uo_recovery", "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h")]))
  formal$baseline_vaso_log <- log1p(pmax(formal$baseline_vaso_burden, 0))
  formal$baseline_uo_log <- log1p(pmax(formal$baseline_uo_ml_kg_h, 0))
  pred_cols <- c("age", "male", "severity_primary", "severity_secondary", "baseline_map_mean", "baseline_vaso_log", "baseline_uo_log")
  pred_cols <- pred_cols[pred_cols %in% names(formal)]
  for (v in pred_cols) {
    formal[[paste0(v, "_missing")]] <- as.integer(is.na(formal[[v]]))
    formal[[v]] <- impute_median(formal[[v]])
  }
  missing_cols <- paste0(pred_cols, "_missing")
  miss_terms <- missing_cols[vapply(formal[, missing_cols, drop = FALSE], has_variation, logical(1))]
  value_terms <- pred_cols[vapply(formal[, pred_cols, drop = FALSE], has_variation, logical(1))]
  form <- as.formula(paste("complete_core ~", paste(c(value_terms, miss_terms), collapse = " + ")))
  fit <- glm(form, data = formal, family = binomial())
  p_complete <- pmin(pmax(predict(fit, type = "response"), 0.02), 0.98)
  names(p_complete) <- formal$id
  scores$baseline_vaso_log <- scores$baseline_vaso_log
  scores$baseline_uo_log <- scores$baseline_uo_log
  rownames(scores) <- as.character(scores$id)
  raw_w <- mean(formal$complete_core == 1) / p_complete[as.character(scores$id)]
  lo <- quantile(raw_w, 0.01, na.rm = TRUE)
  hi <- quantile(raw_w, 0.99, na.rm = TRUE)
  w <- pmin(pmax(raw_w, lo), hi)
  names(w) <- rownames(scores)
  fit_or(scores, "hrc_core_z", database, "IPW complete-case HRC", "HRC per 1 SD increase", weights = w)
}

domain_ablation <- function(scores, database) {
  variants <- list(
    "MAP residual only" = "hrc_map_residual_z",
    "Vaso residual only" = "hrc_vaso_residual_z",
    "Urine residual only" = "hrc_uo_residual_z",
    "MAP + Vaso" = c("hrc_map_residual_z", "hrc_vaso_residual_z"),
    "MAP + Urine" = c("hrc_map_residual_z", "hrc_uo_residual_z"),
    "Vaso + Urine" = c("hrc_vaso_residual_z", "hrc_uo_residual_z"),
    "Core all domains" = c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")
  )
  do.call(rbind, lapply(names(variants), function(label) {
    x <- scores
    x$hrc_variant_z <- zscore(rowMeans(x[, variants[[label]], drop = FALSE], na.rm = FALSE))
    fit_or(x, "hrc_variant_z", database, label, "HRC variant per 1 SD increase")
  }))
}

meta_summary <- function(rows, analysis, exposure) {
  rows <- rows[is.finite(rows$beta) & is.finite(rows$se) & rows$se > 0, ]
  if (nrow(rows) == 0) return(data.frame())
  w <- 1 / rows$se^2
  beta_fe <- sum(w * rows$beta) / sum(w)
  se_fe <- sqrt(1 / sum(w))
  q <- if (nrow(rows) > 1) sum(w * (rows$beta - beta_fe)^2) else NA_real_
  q_df <- nrow(rows) - 1
  q_p <- if (nrow(rows) > 1) pchisq(q, df = q_df, lower.tail = FALSE) else NA_real_
  i2 <- if (nrow(rows) > 1 && is.finite(q) && q > 0) max(0, (q - q_df) / q) * 100 else NA_real_
  tau2 <- if (nrow(rows) > 1) max(0, (q - q_df) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (rows$se^2 + tau2)
  beta_re <- sum(wr * rows$beta) / sum(wr)
  se_re <- sqrt(1 / sum(wr))
  data.frame(
    analysis = analysis,
    exposure = exposure,
    k = nrow(rows),
    random_or = exp(beta_re),
    random_ci_low = exp(beta_re - 1.96 * se_re),
    random_ci_high = exp(beta_re + 1.96 * se_re),
    random_p = 2 * pnorm(abs(beta_re / se_re), lower.tail = FALSE),
    i2_percent = i2,
    q_p = q_p,
    stringsAsFactors = FALSE
  )
}

sicdb_interaction <- function(scores) {
  x <- prepare_model_data(scores, "hrc_core_z")
  x$heart_surgery <- scores$heart_surgery[as.integer(rownames(x))]
  x$heart_surgery <- impute_median(x$heart_surgery)
  form <- update(association_formula(x, "hrc_core_z", TRUE), . ~ . + hrc_core_z:heart_surgery)
  fit <- glm(form, data = x, family = binomial())
  tab <- summary(fit)$coefficients
  term <- "hrc_core_z:heart_surgery"
  data.frame(
    database = "SICdb",
    analysis = "HRC by heart-surgery interaction",
    n = nrow(x),
    deaths = sum(x$outcome == 1),
    beta_interaction = tab[term, "Estimate"],
    interaction_or = exp(tab[term, "Estimate"]),
    p_interaction = tab[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

mcs_bootstrap <- function(scores, subgroup_col, label, b = 500) {
  x <- scores[scores[[subgroup_col]] == 1, , drop = FALSE]
  x <- prepare_model_data(x, "hrc_core_z")
  if (nrow(x) < 30 || sum(x$outcome == 1) < 10) return(data.frame())
  simple_form <- outcome ~ hrc_core_z + age + severity_primary + baseline_vaso_log + baseline_map_mean + baseline_uo_log
  fit <- glm(simple_form, data = x, family = binomial())
  beta <- coef(fit)["hrc_core_z"]
  boot_beta <- rep(NA_real_, b)
  for (i in seq_len(b)) {
    idx <- sample(seq_len(nrow(x)), replace = TRUE)
    boot_fit <- try(glm(simple_form, data = x[idx, , drop = FALSE], family = binomial()), silent = TRUE)
    if (!inherits(boot_fit, "try-error")) boot_beta[i] <- coef(boot_fit)["hrc_core_z"]
  }
  boot_beta <- boot_beta[is.finite(boot_beta)]
  data.frame(
    subgroup = label,
    n = nrow(x),
    deaths = sum(x$outcome == 1),
    odds_ratio = exp(beta),
    bootstrap_ci_low = exp(quantile(boot_beta, 0.025, na.rm = TRUE)),
    bootstrap_ci_high = exp(quantile(boot_beta, 0.975, na.rm = TRUE)),
    bootstrap_success = length(boot_beta),
    note = "Bootstrap CI from reduced adjustment model due small MCS subgroups.",
    stringsAsFactors = FALSE
  )
}

scores <- setNames(lapply(c("MIMIC-IV", "eICU", "SICdb"), standardize_scores), c("MIMIC-IV", "eICU", "SICdb"))
formal <- setNames(lapply(c("MIMIC-IV", "eICU", "SICdb"), standardize_formal), c("MIMIC-IV", "eICU", "SICdb"))

primary_sources <- data.frame(
  database = c("MIMIC-IV", "eICU", "SICdb", "Random-effects summary"),
  odds_ratio = c(0.761880166729433, 0.754964987982257, 0.665950585639046, 0.731850877662616),
  ci_low = c(0.729533643763184, 0.713853331758085, 0.613626695925942, 0.68265216342514),
  ci_high = c(0.79566089023866, 0.798444313029004, 0.722736128427034, 0.784595341276285)
)
ev <- t(mapply(e_value_or, primary_sources$odds_ratio, primary_sources$ci_low, primary_sources$ci_high))
evalue <- cbind(primary_sources, as.data.frame(ev))
write.csv(evalue, file.path(outdir, "issue1_unmeasured_confounding_evalues.csv"), row.names = FALSE)
write_md_table(evalue, file.path(outdir, "issue1_unmeasured_confounding_evalues.md"), "Issue 1: E-values for unmeasured confounding")

summary_files <- list(
  "MIMIC-IV" = "outputs/mimic_formal/mimic_formal_cohort_summary.csv",
  "eICU" = "outputs/eicu_formal/eicu_formal_cohort_summary.csv",
  "SICdb" = "outputs/sicdb_formal/sicdb_formal_cohort_summary.csv"
)
landmark <- do.call(rbind, lapply(names(summary_files), function(db) {
  s <- read.csv(summary_files[[db]], stringsAsFactors = FALSE)
  metric <- setNames(num(s$value), s$metric)
  early <- if (db == "MIMIC-IV") {
    metric["early_vaso_support_stays"]
  } else if (db == "eICU") {
    metric["early_vasoactive_support_stays"]
  } else {
    metric["early_vasoactive_support_cases"]
  }
  landmark_n <- if (db == "MIMIC-IV") metric["landmark_eligible_stays"] else if (db == "eICU") metric["landmark_eligible_stays"] else metric["landmark_eligible_cases"]
  formal_n <- if (db == "SICdb") metric["formal_first_patient_support_cases"] else metric["formal_first_patient_support_stays"]
  data.frame(
    database = db,
    early_support_n = early,
    landmark_eligible_n = landmark_n,
    excluded_before_landmark_n = early - landmark_n,
    excluded_before_landmark_percent = 100 * (early - landmark_n) / early,
    first_patient_support_n = formal_n,
    stringsAsFactors = FALSE
  )
}))
write.csv(landmark, file.path(outdir, "issue2_landmark_selection_audit.csv"), row.names = FALSE)
write_md_table(landmark, file.path(outdir, "issue2_landmark_selection_audit.md"), "Issue 2: Landmark selection audit")

incremental <- do.call(rbind, lapply(names(scores), function(db) cv_incremental_metrics(scores[[db]], db)))
write.csv(incremental, file.path(outdir, "issue4_incremental_model_performance.csv"), row.names = FALSE)
write_md_table(incremental, file.path(outdir, "issue4_incremental_model_performance.md"), "Issue 4: Baseline model versus baseline plus HRC")

cc_smd <- do.call(rbind, lapply(names(formal), function(db) complete_case_diagnostics(formal[[db]], db)))
write.csv(cc_smd, file.path(outdir, "issue5_complete_case_smd.csv"), row.names = FALSE)

ipw <- do.call(rbind, lapply(names(scores), function(db) ipw_complete_case(formal[[db]], scores[[db]], db)))
write.csv(ipw, file.path(outdir, "issue5_ipw_complete_case_or.csv"), row.names = FALSE)
write_md_table(ipw, file.path(outdir, "issue5_ipw_complete_case_or.md"), "Issue 5: IPW complete-case sensitivity")

cluster <- rbind(
  fit_cluster_or(scores[["eICU"]], "eICU"),
  fit_cluster_or(scores[["SICdb"]], "SICdb")
)
write.csv(cluster, file.path(outdir, "issue6_cluster_robust_or.csv"), row.names = FALSE)
write_md_table(cluster, file.path(outdir, "issue6_cluster_robust_or.md"), "Issue 6: Center/cluster-robust sensitivity")

domain <- do.call(rbind, lapply(names(scores), function(db) domain_ablation(scores[[db]], db)))
domain_meta <- do.call(rbind, lapply(unique(domain$analysis), function(a) {
  x <- domain[domain$analysis == a, ]
  meta_summary(x, a, "HRC variant per 1 SD increase")
}))
write.csv(domain, file.path(outdir, "issue3_domain_ablation_or.csv"), row.names = FALSE)
write.csv(domain_meta, file.path(outdir, "issue3_domain_ablation_meta.csv"), row.names = FALSE)
write_md_table(domain_meta, file.path(outdir, "issue3_domain_ablation_meta.md"), "Issue 3: HRC domain and vasopressor-harmonization sensitivity")

interaction <- sicdb_interaction(scores[["SICdb"]])
write.csv(interaction, file.path(outdir, "issue7_sicdb_heart_surgery_interaction.csv"), row.names = FALSE)
write_md_table(interaction, file.path(outdir, "issue7_sicdb_heart_surgery_interaction.md"), "Issue 7: SICdb cardiac-surgery interaction")

mcs_rows <- do.call(rbind, list(
  mcs_bootstrap(scores[["MIMIC-IV"]], "cs_or_mcs_subgroup", "CS or MCS"),
  mcs_bootstrap(scores[["MIMIC-IV"]], "mcs_any_icu", "Any MCS during ICU"),
  mcs_bootstrap(scores[["MIMIC-IV"]], "mcs_any_0_24h", "MCS during 0-24h"),
  mcs_bootstrap(scores[["MIMIC-IV"]], "iabp_any_icu", "IABP"),
  mcs_bootstrap(scores[["MIMIC-IV"]], "impella_any_icu", "Impella"),
  mcs_bootstrap(scores[["MIMIC-IV"]], "ecmo_any_icu", "ECMO")
))
write.csv(mcs_rows, file.path(outdir, "issue8_mimic_mcs_bootstrap_or.csv"), row.names = FALSE)
write_md_table(mcs_rows, file.path(outdir, "issue8_mimic_mcs_bootstrap_or.md"), "Issue 8: MIMIC MCS bootstrap sensitivity")

key <- c(
  "# Methodology stress-test key results",
  "",
  "## Issue 1: Unmeasured confounding",
  paste0("- Random-effects summary E-value: ", sprintf("%.2f", evalue$e_value[evalue$database == "Random-effects summary"]), "; CI-bound E-value: ", sprintf("%.2f", evalue$e_value_ci[evalue$database == "Random-effects summary"]), "."),
  "",
  "## Issue 2: Landmark selection",
  paste(apply(landmark, 1, function(r) paste0("- ", r[["database"]], ": excluded before landmark ", r[["excluded_before_landmark_n"]], " / ", r[["early_support_n"]], " (", sprintf("%.1f%%", as.numeric(r[["excluded_before_landmark_percent"]])), ").")), collapse = "\n"),
  "",
  "## Issue 3: Domain/vasopressor harmonization",
  paste(apply(domain_meta[domain_meta$analysis %in% c("Vaso residual only", "Core all domains", "Vaso + Urine"), ], 1, function(r) paste0("- ", r[["analysis"]], ": random-effects OR ", fmt_or(as.numeric(r[["random_or"]]), as.numeric(r[["random_ci_low"]]), as.numeric(r[["random_ci_high"]])), ".")), collapse = "\n"),
  "",
  "## Issue 4: Dynamic severity score concern",
  paste(apply(incremental, 1, function(r) paste0("- ", r[["database"]], ": delta AUC ", sprintf("%.3f", as.numeric(r[["delta_auc"]])), ", delta Brier ", sprintf("%.4f", as.numeric(r[["delta_brier"]])), ", LRT p=", signif(as.numeric(r[["lrt_p"]]), 3), ".")), collapse = "\n"),
  "",
  "## Issue 5: Complete-case/IPW",
  paste(apply(ipw, 1, function(r) paste0("- ", r[["database"]], ": IPW OR ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["ci_low"]]), as.numeric(r[["ci_high"]])), ".")), collapse = "\n"),
  "",
  "## Issue 6: Center clustering",
  paste(apply(cluster, 1, function(r) paste0(
    "- ", r[["database"]], ": cluster-robust OR ",
    fmt_or(
      as.numeric(r[["cluster_robust_or"]]),
      as.numeric(r[["cluster_robust_ci_low"]]),
      as.numeric(r[["cluster_robust_ci_high"]])
    ),
    ", clusters=", r[["clusters"]], "."
  )), collapse = "\n"),
  "",
  "## Issue 7: SICdb cardiac surgery",
  paste0("- HRC x heart-surgery interaction OR ", sprintf("%.2f", interaction$interaction_or), ", p=", signif(interaction$p_interaction, 3), "."),
  "",
  "## Issue 8: MCS small samples",
  paste(apply(mcs_rows, 1, function(r) paste0("- ", r[["subgroup"]], ": n=", r[["n"]], ", deaths=", r[["deaths"]], ", bootstrap OR ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["bootstrap_ci_low"]]), as.numeric(r[["bootstrap_ci_high"]])), ".")), collapse = "\n")
)
writeLines(key, file.path(outdir, "methodology_stress_test_key_results.md"))

cat("Methodology stress-test outputs:", outdir, "\n")
cat("Key results:", file.path(outdir, "methodology_stress_test_key_results.md"), "\n")
