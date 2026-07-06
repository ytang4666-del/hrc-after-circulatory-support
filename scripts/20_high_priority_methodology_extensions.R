#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(splines)
  library(survival)
  library(cmprsk)
  library(mgcv)
  library(mice)
})

outdir <- "outputs/high_priority_methodology"
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

fmt_p <- function(p) {
  ifelse(
    !is.finite(p),
    "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3g", p))
  )
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

term_for_numeric <- function(dat, var, df = 4) {
  if (!(var %in% names(dat))) return(character())
  x <- dat[[var]]
  ux <- unique(x[is.finite(x) & !is.na(x)])
  if (length(ux) <= 1) return(character())
  if (length(ux) <= df) return(var)
  qs <- unique(as.numeric(quantile(x, probs = seq(0, 1, length.out = df + 2), na.rm = TRUE)))
  if (length(qs) < df + 2) return(var)
  paste0("ns(", var, ", df = ", df, ")")
}

gam_term_for_numeric <- function(dat, var, k = 5) {
  if (!(var %in% names(dat))) return(character())
  x <- dat[[var]]
  ux <- unique(x[is.finite(x) & !is.na(x)])
  if (length(ux) <= 1) return(character())
  if (length(ux) <= 4) return(var)
  paste0("s(", var, ", k = ", min(k, length(ux) - 1), ")")
}

meta_summary <- function(rows, analysis, exposure = "HRC per 1 SD increase") {
  rows <- rows[is.finite(rows$beta) & is.finite(rows$se) & rows$se > 0, ]
  if (nrow(rows) == 0) return(data.frame())
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
    random_or = exp(beta_re),
    random_ci_low = exp(beta_re - 1.96 * se_re),
    random_ci_high = exp(beta_re + 1.96 * se_re),
    q_stat = q,
    q_p = q_p,
    i2_percent = i2,
    tau2_dl = tau2,
    stringsAsFactors = FALSE
  )
}

hartung_knapp_summary <- function(rows, analysis) {
  rows <- rows[is.finite(rows$beta) & is.finite(rows$se) & rows$se > 0, ]
  if (nrow(rows) < 2) return(data.frame())
  base <- meta_summary(rows, analysis)
  tau2 <- base$tau2_dl[1]
  wr <- 1 / (rows$se^2 + tau2)
  beta_re <- sum(wr * rows$beta) / sum(wr)
  k <- nrow(rows)
  hk_scale <- sum(wr * (rows$beta - beta_re)^2) / (k - 1)
  hk_se <- sqrt(max(hk_scale, 1) / sum(wr))
  hk_t <- qt(0.975, df = k - 1)
  pred_t <- if (k > 2) qt(0.975, df = k - 2) else qt(0.975, df = k - 1)
  pred_se <- sqrt(tau2 + hk_se^2)
  data.frame(
    analysis = analysis,
    k = k,
    random_or = exp(beta_re),
    hk_ci_low = exp(beta_re - hk_t * hk_se),
    hk_ci_high = exp(beta_re + hk_t * hk_se),
    prediction_low = exp(beta_re - pred_t * pred_se),
    prediction_high = exp(beta_re + pred_t * pred_se),
    i2_percent = base$i2_percent[1],
    tau2_dl = tau2,
    stringsAsFactors = FALSE
  )
}

read_formal <- function(database) {
  if (database == "MIMIC-IV") {
    d <- read.csv("outputs/mimic_formal/mimic_formal_cohort.csv", stringsAsFactors = FALSE)
    t_index <- as.POSIXct(d$index_time, tz = "UTC")
    t_icu <- as.POSIXct(d$icu_intime, tz = "UTC")
    out <- data.frame(
      database = database,
      id = as.character(d$stay_id),
      cluster_id = d$first_careunit,
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$anchor_age),
      male = ifelse(d$gender == "M", 1, 0),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$first_day_sofa),
      severity_secondary = num(d$sofa2_total),
      heart_surgery = NA_real_,
      index_hour_from_icu = as.numeric(difftime(t_index, t_icu, units = "hours")),
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
    out <- data.frame(
      database = database,
      id = as.character(d$patientunitstayid),
      cluster_id = as.character(d$hospitalid),
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$age),
      male = ifelse(d$gender == "Male", 1, 0),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      heart_surgery = NA_real_,
      index_hour_from_icu = num(d$index_offset) / 60,
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
    out <- data.frame(
      database = database,
      id = as.character(d$caseid),
      cluster_id = as.character(d$hospitalunit),
      outcome = num(d$hospital_mortality_after_landmark),
      age = num(d$age),
      male = num(d$male),
      weight_kg = num(d$weight_kg),
      severity_primary = num(d$saps3),
      severity_secondary = NA_real_,
      heart_surgery = ifelse(d$heartsurgeryadditionaldata == "740", 1, 0),
      index_hour_from_icu = (num(d$index_offset_sec) - num(d$icu_offset_sec)) / 3600,
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
  out$weight_kg[out$weight_kg <= 0 | out$weight_kg > 300] <- NA_real_
  out$baseline_vaso_log <- log1p(pmax(out$baseline_vaso_burden, 0))
  out$baseline_uo_log <- log1p(pmax(out$baseline_uo_ml_kg_h, 0))
  out
}

read_scores <- function(database) {
  if (database == "MIMIC-IV") {
    s <- read.csv("outputs/mimic_formal/mimic_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = database,
      id = as.character(s$stay_id),
      hrc_core_z = num(s$hrc_core_z),
      outcome = num(s$hospital_mortality_after_landmark),
      stringsAsFactors = FALSE
    )
  } else if (database == "eICU") {
    s <- read.csv("outputs/eicu_formal/eicu_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = database,
      id = as.character(s$patientunitstayid),
      hrc_core_z = num(s$hrc_core_z),
      outcome = num(s$hospital_mortality_after_landmark),
      stringsAsFactors = FALSE
    )
  } else {
    s <- read.csv("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = database,
      id = as.character(s$caseid),
      hrc_core_z = num(s$hrc_core_z),
      outcome = num(s$hospital_mortality_after_landmark),
      stringsAsFactors = FALSE
    )
  }
  merge(read_formal(database), out[, c("id", "hrc_core_z")], by = "id", all.x = FALSE, all.y = FALSE)
}

base_expected_terms <- function(dat, include_severity = FALSE, gam = FALSE) {
  f <- if (gam) gam_term_for_numeric else term_for_numeric
  terms <- c(f(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, f(dat, "weight_kg", 3))
  if (include_severity) {
    terms <- c(terms, f(dat, "severity_primary", 4))
    terms <- c(terms, f(dat, "severity_secondary", 4))
  }
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) terms <- c(terms, "heart_surgery")
  terms <- c(terms, f(dat, "index_hour_from_icu", 4))
  terms <- c(terms, f(dat, "baseline_map_mean", 4))
  terms <- c(terms, f(dat, "baseline_vaso_log", 4))
  terms <- c(terms, f(dat, "baseline_uo_log", 4))
  unique(terms[nzchar(terms)])
}

mortality_terms <- function(dat, include_severity = TRUE, include_index = TRUE) {
  terms <- c(term_for_numeric(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  if (include_severity) {
    terms <- c(terms, term_for_numeric(dat, "severity_primary", 4))
    terms <- c(terms, term_for_numeric(dat, "severity_secondary", 4))
  }
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) terms <- c(terms, "heart_surgery")
  if (include_index) terms <- c(terms, term_for_numeric(dat, "index_hour_from_icu", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  unique(terms[nzchar(terms)])
}

crossfit_component <- function(dat, outcome, folds, terms, model = c("lm", "gam")) {
  model <- match.arg(model)
  expected <- rep(NA_real_, nrow(dat))
  fold_rmse <- rep(NA_real_, max(folds))
  form <- as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
  for (fold in sort(unique(folds))) {
    train <- dat[folds != fold, , drop = FALSE]
    test <- dat[folds == fold, , drop = FALSE]
    fit <- if (model == "gam") {
      suppressWarnings(mgcv::gam(form, data = train, method = "REML"))
    } else {
      lm(form, data = train)
    }
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

build_hrc <- function(dat, database, analysis, include_severity = FALSE, model = c("lm", "gam")) {
  model <- match.arg(model)
  required <- c(
    "delta_map", "vaso_reduction", "log_uo_recovery", "outcome",
    "baseline_map_mean", "baseline_vaso_burden", "baseline_uo_ml_kg_h"
  )
  x <- dat[complete.cases(dat[, required]), , drop = FALSE]
  x <- x[
    is.finite(x$delta_map) & is.finite(x$vaso_reduction) &
      is.finite(x$log_uo_recovery) & is.finite(x$outcome),
    ,
    drop = FALSE
  ]
  covars <- c(
    "age", "male", "weight_kg", "severity_primary", "severity_secondary",
    "heart_surgery", "index_hour_from_icu", "baseline_map_mean",
    "baseline_vaso_log", "baseline_uo_log", "baseline_vaso_burden",
    "baseline_uo_ml_kg_h"
  )
  for (col in intersect(covars, names(x))) x[[col]] <- impute_median(x[[col]])
  folds <- sample(rep(1:5, length.out = nrow(x)))
  terms <- base_expected_terms(x, include_severity = include_severity, gam = model == "gam")
  map_cf <- crossfit_component(x, "delta_map", folds, terms, model = model)
  vaso_cf <- crossfit_component(x, "vaso_reduction", folds, terms, model = model)
  uo_cf <- crossfit_component(x, "log_uo_recovery", folds, terms, model = model)
  x$hrc_map_residual_z <- map_cf$residual_z
  x$hrc_vaso_residual_z <- vaso_cf$residual_z
  x$hrc_uo_residual_z <- uo_cf$residual_z
  x$hrc_core_raw <- rowMeans(
    x[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")],
    na.rm = FALSE
  )
  x$hrc_core_z <- zscore(x$hrc_core_raw)
  list(
    data = x,
    metrics = data.frame(
      database = database,
      analysis = analysis,
      component = c("MAP recovery", "Vasopressor-burden reduction", "Urine output recovery"),
      observed_variable = c("delta_map", "vaso_reduction", "log_uo_recovery"),
      n = nrow(x),
      oof_rmse = c(map_cf$rmse, vaso_cf$rmse, uo_cf$rmse),
      oof_r2 = c(map_cf$r2_oof, vaso_cf$r2_oof, uo_cf$r2_oof),
      stringsAsFactors = FALSE
    )
  )
}

fit_logit_or <- function(dat, exposure = "hrc_core_z", analysis = "", database = "", include_severity = TRUE, include_index = TRUE) {
  terms <- c(exposure, mortality_terms(dat, include_severity = include_severity, include_index = include_index))
  terms <- unique(terms[nzchar(terms)])
  cols <- unique(c("outcome", all.vars(as.formula(paste("~", paste(terms, collapse = "+"))))))
  x <- dat[complete.cases(dat[, c("outcome", exposure)]), cols, drop = FALSE]
  x <- x[is.finite(x$outcome) & is.finite(x[[exposure]]), , drop = FALSE]
  for (v in setdiff(names(x), c("outcome", exposure))) if (is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  fit <- glm(as.formula(paste("outcome ~", paste(terms, collapse = " + "))), data = x, family = binomial())
  tab <- summary(fit)$coefficients
  ci <- confint.default(fit)
  data.frame(
    database = database,
    analysis = analysis,
    n = nrow(x),
    deaths = sum(x$outcome == 1, na.rm = TRUE),
    beta = tab[exposure, "Estimate"],
    se = tab[exposure, "Std. Error"],
    odds_ratio = exp(tab[exposure, "Estimate"]),
    ci_low = exp(ci[exposure, 1]),
    ci_high = exp(ci[exposure, 2]),
    p_value = tab[exposure, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

fit_mice_or <- function(dat, database) {
  vars <- c(
    "outcome", "hrc_core_z", "age", "male", "severity_primary", "severity_secondary",
    "heart_surgery", "index_hour_from_icu", "baseline_vaso_log",
    "baseline_map_mean", "baseline_uo_log"
  )
  x <- dat[, intersect(vars, names(dat)), drop = FALSE]
  x <- x[complete.cases(x[, c("outcome", "hrc_core_z")]), , drop = FALSE]
  keep <- vapply(x, function(col) !(all(is.na(col)) || !has_variation(num(col))), logical(1))
  keep[c("outcome", "hrc_core_z")] <- TRUE
  x <- x[, keep, drop = FALSE]
  for (v in names(x)) x[[v]] <- num(x[[v]])
  methods <- make.method(x)
  methods[c("outcome", "hrc_core_z")] <- ""
  pred <- make.predictorMatrix(x)
  pred[, c("outcome", "hrc_core_z")] <- pred[, c("outcome", "hrc_core_z")]
  pred[c("outcome", "hrc_core_z"), ] <- 0
  imp <- mice(x, m = 5, maxit = 5, method = methods, predictorMatrix = pred, printFlag = FALSE, seed = 20260704)
  rhs <- paste(setdiff(names(x), "outcome"), collapse = " + ")
  pooled <- pool(with(imp, glm(as.formula(paste("outcome ~", rhs)), family = binomial())))
  ptab <- summary(pooled)
  row <- ptab[ptab$term == "hrc_core_z", ]
  data.frame(
    database = database,
    analysis = "MICE covariate-imputation sensitivity",
    n = nrow(x),
    deaths = sum(x$outcome == 1, na.rm = TRUE),
    beta = row$estimate,
    se = row$std.error,
    odds_ratio = exp(row$estimate),
    ci_low = exp(row$estimate - 1.96 * row$std.error),
    ci_high = exp(row$estimate + 1.96 * row$std.error),
    p_value = row$p.value,
    stringsAsFactors = FALSE
  )
}

fit_low_cutoffs <- function(dat, database) {
  cutoffs <- c(0.10, 0.20, 0.25, 0.33)
  do.call(rbind, lapply(cutoffs, function(p) {
    x <- dat
    threshold <- as.numeric(quantile(x$hrc_core_z, probs = p, na.rm = TRUE))
    x$low_hrc <- ifelse(x$hrc_core_z <= threshold, 1, 0)
    row <- fit_logit_or(
      x,
      exposure = "low_hrc",
      analysis = paste0("Low-HRC cutoff p", sprintf("%02.0f", p * 100)),
      database = database,
      include_severity = TRUE,
      include_index = TRUE
    )
    row$cutoff_percentile <- p
    row$threshold <- threshold
    row
  }))
}

fit_index_strata <- function(dat, database) {
  breaks <- c(0, 6, 12, 24)
  labels <- c("index_0_6h", "index_6_12h", "index_12_24h")
  dat$index_stratum <- cut(dat$index_hour_from_icu, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)
  rows <- list()
  for (lab in labels) {
    x <- dat[dat$index_stratum == lab, , drop = FALSE]
    if (nrow(x) < 300 || sum(x$outcome == 1, na.rm = TRUE) < 30) next
    rows[[lab]] <- fit_logit_or(
      x,
      exposure = "hrc_core_z",
      analysis = paste0("Index-delay stratum: ", lab),
      database = database,
      include_severity = TRUE,
      include_index = FALSE
    )
  }
  do.call(rbind, rows)
}

read_survival <- function(database) {
  file <- switch(
    database,
    "MIMIC-IV" = file.path(outdir, "mimic_survival_source.csv"),
    "eICU" = file.path(outdir, "eicu_survival_source.csv"),
    "SICdb" = file.path(outdir, "sicdb_survival_source.csv")
  )
  d <- read.csv(file, stringsAsFactors = FALSE)
  d$id <- as.character(d$stay_id)
  d$time_to_hospital_exit_days <- num(d$time_to_hospital_exit_days)
  d$death_event <- num(d$death_event)
  d$alive_discharge_event <- num(d$alive_discharge_event)
  d
}

fit_survival_models <- function(dat, database) {
  surv <- read_survival(database)
  x <- merge(dat, surv[, c("id", "time_to_hospital_exit_days", "death_event", "alive_discharge_event")], by = "id")
  x <- x[is.finite(x$time_to_hospital_exit_days) & x$time_to_hospital_exit_days > 0, , drop = FALSE]
  x$fstatus <- ifelse(x$death_event == 1, 1, ifelse(x$alive_discharge_event == 1, 2, 0))
  terms <- c("hrc_core_z", mortality_terms(x, include_severity = TRUE, include_index = TRUE))
  terms <- unique(terms[nzchar(terms)])
  covars <- setdiff(all.vars(as.formula(paste("~", paste(terms, collapse = "+")))), character())
  keep <- unique(c("time_to_hospital_exit_days", "death_event", "fstatus", covars))
  y <- x[, keep, drop = FALSE]
  y <- y[complete.cases(y[, c("time_to_hospital_exit_days", "death_event", "fstatus", "hrc_core_z")]), , drop = FALSE]
  for (v in setdiff(names(y), c("time_to_hospital_exit_days", "death_event", "fstatus", "hrc_core_z"))) {
    if (is.numeric(y[[v]])) y[[v]] <- impute_median(y[[v]])
  }
  cox_fit <- coxph(
    as.formula(paste("Surv(time_to_hospital_exit_days, death_event == 1) ~", paste(terms, collapse = " + "))),
    data = y
  )
  cox_tab <- summary(cox_fit)$coefficients
  cox_ci <- confint.default(cox_fit)
  mm <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = y)[, -1, drop = FALSE]
  fg_fit <- cmprsk::crr(
    ftime = y$time_to_hospital_exit_days,
    fstatus = y$fstatus,
    cov1 = mm,
    failcode = 1,
    cencode = 0
  )
  fg_tab <- summary(fg_fit)$coef
  fg_row <- "hrc_core_z"
  if (!(fg_row %in% rownames(fg_tab))) {
    fg_row <- grep("^hrc_core_z$", rownames(fg_tab), value = TRUE)[1]
  }
  rbind(
    data.frame(
      database = database,
      analysis = "Cause-specific Cox death model",
      n = nrow(y),
      deaths = sum(y$death_event == 1, na.rm = TRUE),
      beta = cox_tab["hrc_core_z", "coef"],
      se = cox_tab["hrc_core_z", "se(coef)"],
      hazard_ratio = exp(cox_tab["hrc_core_z", "coef"]),
      ci_low = exp(cox_ci["hrc_core_z", 1]),
      ci_high = exp(cox_ci["hrc_core_z", 2]),
      p_value = cox_tab["hrc_core_z", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    ),
    data.frame(
      database = database,
      analysis = "Fine-Gray death subdistribution model",
      n = nrow(y),
      deaths = sum(y$death_event == 1, na.rm = TRUE),
      beta = fg_tab[fg_row, "coef"],
      se = fg_tab[fg_row, "se(coef)"],
      hazard_ratio = exp(fg_tab[fg_row, "coef"]),
      ci_low = exp(fg_tab[fg_row, "coef"] - 1.96 * fg_tab[fg_row, "se(coef)"]),
      ci_high = exp(fg_tab[fg_row, "coef"] + 1.96 * fg_tab[fg_row, "se(coef)"]),
      p_value = fg_tab[fg_row, "p-value"],
      stringsAsFactors = FALSE
    )
  )
}

fit_pooled_interaction <- function(score_list) {
  x <- do.call(rbind, score_list)
  keep <- c(
    "database", "outcome", "hrc_core_z", "age", "male", "severity_primary",
    "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log", "index_hour_from_icu"
  )
  x <- x[, intersect(keep, names(x)), drop = FALSE]
  x <- x[complete.cases(x[, c("outcome", "hrc_core_z", "database")]), , drop = FALSE]
  for (v in setdiff(names(x), c("database", "outcome", "hrc_core_z"))) if (is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  x$database <- factor(x$database, levels = c("MIMIC-IV", "eICU", "SICdb"))
  base_terms <- c(
    "hrc_core_z", "database", "ns(age, df = 3)", "male",
    "ns(severity_primary, df = 4)", "ns(baseline_vaso_log, df = 4)",
    "ns(baseline_map_mean, df = 4)", "ns(baseline_uo_log, df = 4)",
    "ns(index_hour_from_icu, df = 4)"
  )
  no_int <- glm(as.formula(paste("outcome ~", paste(base_terms, collapse = " + "))), data = x, family = binomial())
  with_int <- glm(
    as.formula(paste("outcome ~", paste(c(base_terms, "hrc_core_z:database"), collapse = " + "))),
    data = x,
    family = binomial()
  )
  lrt <- anova(no_int, with_int, test = "LRT")
  data.frame(
    analysis = "Pooled HRC-by-database interaction",
    n = nrow(x),
    deaths = sum(x$outcome == 1, na.rm = TRUE),
    lrt_chisq = lrt$Deviance[2],
    df = lrt$Df[2],
    p_value = lrt$`Pr(>Chi)`[2],
    stringsAsFactors = FALSE
  )
}

databases <- c("MIMIC-IV", "eICU", "SICdb")

message("Reading primary scores...")
score_data <- setNames(lapply(databases, read_scores), databases)
formal_data <- setNames(lapply(databases, read_formal), databases)

message("Issue 1: no-leakage HRC rebuild...")
no_leakage <- lapply(databases, function(db) {
  built <- build_hrc(formal_data[[db]], db, "No-leakage baseline HRC", include_severity = FALSE, model = "lm")
  list(
    or = fit_logit_or(built$data, "hrc_core_z", "No-leakage baseline HRC", db, include_severity = FALSE, include_index = TRUE),
    metrics = built$metrics
  )
})
no_leakage_or <- do.call(rbind, lapply(no_leakage, `[[`, "or"))
no_leakage_metrics <- do.call(rbind, lapply(no_leakage, `[[`, "metrics"))

message("Issue 4: flexible GAM expected-recovery HRC...")
gam_sens <- lapply(databases, function(db) {
  built <- build_hrc(formal_data[[db]], db, "GAM expected-recovery HRC", include_severity = FALSE, model = "gam")
  list(
    or = fit_logit_or(built$data, "hrc_core_z", "GAM expected-recovery HRC", db, include_severity = FALSE, include_index = TRUE),
    metrics = built$metrics
  )
})
gam_or <- do.call(rbind, lapply(gam_sens, `[[`, "or"))
gam_metrics <- do.call(rbind, lapply(gam_sens, `[[`, "metrics"))

message("Issue 5: index-time adjustment and strata...")
index_adjusted_or <- do.call(rbind, lapply(databases, function(db) {
  fit_logit_or(score_data[[db]], "hrc_core_z", "Primary HRC with index-hour adjustment", db, include_severity = TRUE, include_index = TRUE)
}))
index_strata_or <- do.call(rbind, lapply(databases, function(db) fit_index_strata(score_data[[db]], db)))

message("Issue 6: MICE covariate-imputation sensitivity...")
mice_or <- do.call(rbind, lapply(databases, function(db) fit_mice_or(score_data[[db]], db)))

message("Issue 7: low-HRC cutoff sensitivity...")
low_cutoff_or <- do.call(rbind, lapply(databases, function(db) fit_low_cutoffs(score_data[[db]], db)))
low_cutoff_meta <- do.call(rbind, lapply(split(low_cutoff_or, low_cutoff_or$analysis), function(x) {
  meta_summary(x, unique(x$analysis), "low-HRC versus non-low HRC")
}))

message("Issue 2: survival and competing-risk sensitivity...")
survival_rows <- do.call(rbind, lapply(databases, function(db) fit_survival_models(score_data[[db]], db)))

message("Issue 3: Hartung-Knapp meta-analysis and prediction intervals...")
main_or <- do.call(rbind, lapply(databases, function(db) {
  file <- switch(
    db,
    "MIMIC-IV" = "outputs/mimic_formal/mimic_hrc_formal_linear_or.csv",
    "eICU" = "outputs/eicu_formal/eicu_hrc_formal_linear_or.csv",
    "SICdb" = "outputs/sicdb_formal/sicdb_hrc_formal_linear_or.csv"
  )
  x <- read.csv(file, stringsAsFactors = FALSE)
  row <- x[x$term == "hrc_core_z", ]
  data.frame(database = db, analysis = "Primary main HRC", beta = num(row$beta), se = num(row$se), stringsAsFactors = FALSE)
}))
tw <- read.csv("outputs/time_window_sensitivity/time_window_hrc_or.csv", stringsAsFactors = FALSE)
tw_meta_input <- data.frame(
  database = tw$database,
  analysis = paste0("Time-window: ", tw$window_strategy),
  beta = num(tw$beta),
  se = num(tw$se),
  stringsAsFactors = FALSE
)
hk_meta <- do.call(rbind, lapply(split(rbind(main_or, tw_meta_input), rbind(main_or, tw_meta_input)$analysis), function(x) {
  hartung_knapp_summary(x, unique(x$analysis))
}))
loo_meta <- do.call(rbind, lapply(split(main_or, main_or$analysis), function(x) {
  do.call(rbind, lapply(x$database, function(drop_db) {
    y <- x[x$database != drop_db, ]
    m <- meta_summary(y, paste0(unique(x$analysis), " leave-one-out excluding ", drop_db))
    m$excluded_database <- drop_db
    m
  }))
}))

message("Issue 8: recalibrated construct-replication audit...")
pooled_interaction <- fit_pooled_interaction(score_data)

all_or <- rbind(no_leakage_or, gam_or, index_adjusted_or, mice_or)
all_meta <- do.call(rbind, lapply(split(all_or, all_or$analysis), function(x) {
  meta_summary(x, unique(x$analysis), "HRC per 1 SD increase")
}))

write.csv(no_leakage_or, file.path(outdir, "issue1_no_leakage_hrc_or.csv"), row.names = FALSE)
write.csv(no_leakage_metrics, file.path(outdir, "issue1_no_leakage_component_metrics.csv"), row.names = FALSE)
write.csv(survival_rows, file.path(outdir, "issue2_survival_competing_risk.csv"), row.names = FALSE)
write.csv(hk_meta, file.path(outdir, "issue3_hartung_knapp_prediction_meta.csv"), row.names = FALSE)
write.csv(loo_meta, file.path(outdir, "issue3_leave_one_database_out_meta.csv"), row.names = FALSE)
write.csv(gam_or, file.path(outdir, "issue4_gam_expected_recovery_or.csv"), row.names = FALSE)
write.csv(gam_metrics, file.path(outdir, "issue4_gam_component_metrics.csv"), row.names = FALSE)
write.csv(index_adjusted_or, file.path(outdir, "issue5_index_hour_adjusted_or.csv"), row.names = FALSE)
write.csv(index_strata_or, file.path(outdir, "issue5_index_delay_stratified_or.csv"), row.names = FALSE)
write.csv(mice_or, file.path(outdir, "issue6_mice_covariate_imputation_or.csv"), row.names = FALSE)
write.csv(low_cutoff_or, file.path(outdir, "issue7_low_hrc_cutoff_or.csv"), row.names = FALSE)
write.csv(low_cutoff_meta, file.path(outdir, "issue7_low_hrc_cutoff_meta.csv"), row.names = FALSE)
write.csv(pooled_interaction, file.path(outdir, "issue8_construct_replication_interaction.csv"), row.names = FALSE)
write.csv(all_meta, file.path(outdir, "high_priority_extension_meta.csv"), row.names = FALSE)

or_md <- all_or
or_md$OR_95CI <- fmt_or(or_md$odds_ratio, or_md$ci_low, or_md$ci_high)
or_md$p_value <- fmt_p(or_md$p_value)
write_md_table(
  or_md[, c("database", "analysis", "n", "deaths", "OR_95CI", "p_value")],
  file.path(outdir, "high_priority_extension_or.md"),
  "High-priority HRC methodology extension ORs"
)

surv_md <- survival_rows
surv_md$HR_95CI <- fmt_or(surv_md$hazard_ratio, surv_md$ci_low, surv_md$ci_high)
surv_md$p_value <- fmt_p(surv_md$p_value)
write_md_table(
  surv_md[, c("database", "analysis", "n", "deaths", "HR_95CI", "p_value")],
  file.path(outdir, "issue2_survival_competing_risk.md"),
  "Issue 2: survival and competing-risk sensitivity"
)

hk_md <- hk_meta
hk_md$HK_OR_95CI <- fmt_or(hk_md$random_or, hk_md$hk_ci_low, hk_md$hk_ci_high)
hk_md$prediction_interval <- fmt_or(hk_md$random_or, hk_md$prediction_low, hk_md$prediction_high)
hk_md$i2_percent <- sprintf("%.1f", hk_md$i2_percent)
write_md_table(
  hk_md[, c("analysis", "k", "HK_OR_95CI", "prediction_interval", "i2_percent")],
  file.path(outdir, "issue3_hartung_knapp_prediction_meta.md"),
  "Issue 3: Hartung-Knapp meta-analysis and prediction intervals"
)

cut_md <- low_cutoff_or
cut_md$OR_95CI <- fmt_or(cut_md$odds_ratio, cut_md$ci_low, cut_md$ci_high)
cut_md$p_value <- fmt_p(cut_md$p_value)
write_md_table(
  cut_md[, c("database", "analysis", "cutoff_percentile", "threshold", "n", "deaths", "OR_95CI", "p_value")],
  file.path(outdir, "issue7_low_hrc_cutoff_or.md"),
  "Issue 7: low-HRC secondary cutoff sensitivity"
)

interaction_md <- pooled_interaction
interaction_md$p_value <- fmt_p(interaction_md$p_value)
write_md_table(
  interaction_md,
  file.path(outdir, "issue8_construct_replication_interaction.md"),
  "Issue 8: construct replication interaction audit"
)

key_lines <- c(
  "# High-priority methodology extensions",
  "",
  "## Implemented fixes",
  "",
  "1. No-leakage baseline HRC: rebuilt HRC without first-day SOFA/APACHE/SAPS-style severity in expected recovery models.",
  "2. Survival/competing-risk sensitivity: post-landmark cause-specific Cox and Fine-Gray models using hospital death versus alive discharge.",
  "3. Meta-analysis heterogeneity: Hartung-Knapp random-effects confidence intervals, prediction intervals, and leave-one-database-out meta-analysis.",
  "4. Flexible expected-recovery model: GAM-based cross-fitted HRC sensitivity.",
  "5. Index-time delay: adjusted and stratified models by support-start delay after ICU admission.",
  "6. Covariate missingness: MICE covariate-imputation sensitivity.",
  "7. Low-HRC cutoff: treated as secondary and tested 10th, 20th, 25th, and 33rd percentile cutoffs.",
  "8. External validation wording: added statistical construct-replication audit with HRC-by-database interaction.",
  "",
  "## Key continuous-HRC sensitivity ORs",
  "",
  paste(apply(or_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["analysis"]], ": OR=", row[["OR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n"),
  "",
  "## Survival/competing-risk results",
  "",
  paste(apply(surv_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["analysis"]], ": HR=", row[["HR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n"),
  "",
  "## Construct-replication interaction",
  "",
  paste0(
    "- Pooled HRC-by-database interaction LRT p=",
    fmt_p(pooled_interaction$p_value),
    ". This supports wording the analysis as recalibrated construct replication rather than fixed-coefficient score transport."
  )
)
writeLines(key_lines, file.path(outdir, "high_priority_methodology_key_results.md"))

cat("Wrote high-priority methodology extension outputs to", outdir, "\n")
