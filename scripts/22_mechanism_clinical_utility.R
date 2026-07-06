#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

outdir <- "outputs/mechanism_clinical_utility"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260704)

num <- function(x) suppressWarnings(as.numeric(x))

impute_median <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

has_variation <- function(x) {
  ux <- unique(x[is.finite(x) & !is.na(x)])
  length(ux) > 1
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

fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)

fmt_p <- function(p) {
  ifelse(!is.finite(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3g", p)))
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

bind_rows <- function(x) {
  x <- x[!vapply(x, function(y) is.null(y) || nrow(y) == 0, logical(1))]
  if (length(x) == 0) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- lapply(x, function(y) {
    missing <- setdiff(cols, names(y))
    for (m in missing) y[[m]] <- NA
    y[, cols, drop = FALSE]
  })
  do.call(rbind, out)
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
  tau2 <- if (k > 1) max(0, (q - q_df) / (sum(w) - sum(w^2) / sum(w))) else NA_real_
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
    q_p = q_p,
    i2_percent = i2,
    tau2_dl = tau2,
    stringsAsFactors = FALSE
  )
}

read_dataset <- function(database) {
  if (database == "MIMIC-IV") {
    f <- read.csv("outputs/mimic_formal/mimic_formal_cohort.csv", stringsAsFactors = FALSE)
    s <- read.csv("outputs/mimic_formal/mimic_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    p <- read.csv(file.path(outdir, "mimic_post24_mechanism.csv"), stringsAsFactors = FALSE)
    t_index <- as.POSIXct(f$index_time, tz = "UTC")
    t_icu <- as.POSIXct(f$icu_intime, tz = "UTC")
    d <- data.frame(
      database = database,
      id = as.character(f$stay_id),
      patient_id = as.character(f$subject_id),
      cluster_id = f$first_careunit,
      outcome = num(f$hospital_mortality_after_landmark),
      age = num(f$anchor_age),
      male = ifelse(f$gender == "M", 1, 0),
      weight_kg = num(f$weight_kg),
      severity_primary = num(f$first_day_sofa),
      severity_secondary = num(f$sofa2_total),
      index_hour_from_icu = as.numeric(difftime(t_index, t_icu, units = "hours")),
      baseline_map_mean = num(f$baseline_map_mean),
      response_map_mean = num(f$response_map_mean),
      baseline_vaso_burden = num(f$baseline_neq_avg),
      response_vaso_burden = num(f$response_neq_avg),
      baseline_uo_ml_kg_h = num(f$baseline_uo_ml_kg_h),
      response_uo_ml_kg_h = num(f$response_uo_ml_kg_h),
      baseline_creatinine_mean = num(f$baseline_creatinine_mean),
      response_creatinine_mean = num(f$response_creatinine_mean),
      baseline_lactate_mean = num(f$baseline_lactate_mean),
      response_lactate_mean = num(f$response_lactate_mean),
      mcs_any_icu = num(f$mcs_any_icu),
      mcs_any_0_24h = num(f$mcs_any_0_24h),
      cardiogenic_shock_icd = num(f$cardiogenic_shock_icd),
      stringsAsFactors = FALSE
    )
    h <- data.frame(id = as.character(s$stay_id), hrc_core_z = num(s$hrc_core_z), low_hrc_q1 = num(s$low_hrc_q1))
  } else if (database == "eICU") {
    f <- read.csv("outputs/eicu_formal/eicu_formal_cohort.csv", stringsAsFactors = FALSE)
    s <- read.csv("outputs/eicu_formal/eicu_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    p <- read.csv(file.path(outdir, "eicu_post24_mechanism.csv"), stringsAsFactors = FALSE)
    d <- data.frame(
      database = database,
      id = as.character(f$patientunitstayid),
      patient_id = ifelse(f$uniquepid == "", as.character(f$patientunitstayid), f$uniquepid),
      cluster_id = as.character(f$hospitalid),
      outcome = num(f$hospital_mortality_after_landmark),
      age = num(f$age),
      male = ifelse(f$gender == "Male", 1, 0),
      weight_kg = num(f$weight_kg),
      severity_primary = num(f$acutephysiologyscore),
      severity_secondary = num(f$apachescore),
      index_hour_from_icu = num(f$index_offset) / 60,
      baseline_map_mean = num(f$baseline_map_mean),
      response_map_mean = num(f$response_map_mean),
      baseline_vaso_burden = num(f$baseline_vaso_burden),
      response_vaso_burden = num(f$response_vaso_burden),
      baseline_uo_ml_kg_h = num(f$baseline_uo_ml_kg_h),
      response_uo_ml_kg_h = num(f$response_uo_ml_kg_h),
      baseline_creatinine_mean = num(f$baseline_creatinine_mean),
      response_creatinine_mean = num(f$response_creatinine_mean),
      baseline_lactate_mean = num(f$baseline_lactate_mean),
      response_lactate_mean = num(f$response_lactate_mean),
      mcs_any_icu = NA_real_,
      mcs_any_0_24h = NA_real_,
      cardiogenic_shock_icd = NA_real_,
      stringsAsFactors = FALSE
    )
    h <- data.frame(id = as.character(s$patientunitstayid), hrc_core_z = num(s$hrc_core_z), low_hrc_q1 = num(s$low_hrc_q1))
  } else if (database == "SICdb") {
    f <- read.csv("outputs/sicdb_formal/sicdb_formal_cohort.csv", stringsAsFactors = FALSE)
    s <- read.csv("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    p <- read.csv(file.path(outdir, "sicdb_post24_mechanism.csv"), stringsAsFactors = FALSE)
    d <- data.frame(
      database = database,
      id = as.character(f$caseid),
      patient_id = as.character(f$patientid),
      cluster_id = as.character(f$hospitalunit),
      outcome = num(f$hospital_mortality_after_landmark),
      age = num(f$age),
      male = num(f$male),
      weight_kg = num(f$weight_kg),
      severity_primary = num(f$saps3),
      severity_secondary = NA_real_,
      index_hour_from_icu = (num(f$index_offset_sec) - num(f$icu_offset_sec)) / 3600,
      baseline_map_mean = num(f$baseline_map_mean),
      response_map_mean = num(f$response_map_mean),
      baseline_vaso_burden = num(f$baseline_vaso_burden),
      response_vaso_burden = num(f$response_vaso_burden),
      baseline_uo_ml_kg_h = num(f$baseline_uo_ml_kg_h),
      response_uo_ml_kg_h = num(f$response_uo_ml_kg_h),
      baseline_creatinine_mean = num(f$baseline_creatinine_mean),
      response_creatinine_mean = num(f$response_creatinine_mean),
      baseline_lactate_mean = num(f$baseline_lactate_mean),
      response_lactate_mean = num(f$response_lactate_mean),
      mcs_any_icu = NA_real_,
      mcs_any_0_24h = NA_real_,
      cardiogenic_shock_icd = NA_real_,
      stringsAsFactors = FALSE
    )
    h <- data.frame(id = as.character(s$caseid), hrc_core_z = num(s$hrc_core_z), low_hrc_q1 = num(s$low_hrc_q1))
  } else {
    stop("Unknown database")
  }
  p$id <- as.character(p$stay_id)
  for (v in setdiff(names(p), c("database", "patient_id", "encounter_id", "stay_id", "id"))) p[[v]] <- num(p[[v]])
  d$weight_kg[d$weight_kg <= 0 | d$weight_kg > 300] <- NA_real_
  d <- merge(d, h, by = "id", all = FALSE)
  d <- merge(d, p[, setdiff(names(p), c("database", "patient_id", "encounter_id", "stay_id")), drop = FALSE], by = "id", all.x = TRUE)
  d$baseline_vaso_log <- log1p(pmax(d$baseline_vaso_burden, 0))
  d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))
  d$post24_vaso_burden[is.na(d$post24_vaso_burden)] <- 0
  d$persistent_vaso_24_72 <- ifelse(d$post24_vaso_burden > 0.001, 1, 0)
  d$oliguria_24_72 <- ifelse(
    is.finite(d$post24_uo_ml_kg_h) & d$post24_followup_hours >= 6 & d$post24_uo_n > 0,
    ifelse(d$post24_uo_ml_kg_h < 0.5, 1, 0),
    NA_real_
  )
  d$creatinine_worsening_24_72 <- ifelse(
    is.finite(d$baseline_creatinine_mean) & is.finite(d$post24_creatinine_max),
    ifelse(d$post24_creatinine_max >= d$baseline_creatinine_mean + 0.3 |
             d$post24_creatinine_max >= 1.5 * d$baseline_creatinine_mean, 1, 0),
    NA_real_
  )
  d$lactate_nonclearance_24_72 <- ifelse(
    is.finite(d$baseline_lactate_mean) & d$baseline_lactate_mean >= 2 & is.finite(d$post24_lactate_min),
    ifelse(d$post24_lactate_min >= 0.9 * d$baseline_lactate_mean, 1, 0),
    NA_real_
  )
  d$persistent_hyperlactatemia_24_72 <- ifelse(
    is.finite(d$post24_lactate_min),
    ifelse(d$post24_lactate_min >= 2, 1, 0),
    NA_real_
  )
  organ_mat <- cbind(d$oliguria_24_72, d$creatinine_worsening_24_72, d$lactate_nonclearance_24_72)
  d$organ_nonrecovery_24_72 <- ifelse(rowSums(!is.na(organ_mat)) > 0, ifelse(rowSums(organ_mat == 1, na.rm = TRUE) > 0, 1, 0), NA_real_)
  support_mat <- cbind(d$persistent_vaso_24_72, d$oliguria_24_72, d$creatinine_worsening_24_72, d$lactate_nonclearance_24_72)
  d$support_or_organ_nonrecovery_24_72 <- ifelse(rowSums(!is.na(support_mat)) > 0, ifelse(rowSums(support_mat == 1, na.rm = TRUE) > 0, 1, 0), NA_real_)
  d$post24_rrt_any <- ifelse(is.finite(d$post24_rrt_any), d$post24_rrt_any, NA_real_)
  d$post24_invasive_vent_any <- ifelse(is.finite(d$post24_invasive_vent_any), d$post24_invasive_vent_any, NA_real_)

  d$creatinine_worsening_0_24 <- ifelse(
    is.finite(d$baseline_creatinine_mean) & is.finite(d$response_creatinine_mean),
    ifelse(d$response_creatinine_mean >= d$baseline_creatinine_mean + 0.3 |
             d$response_creatinine_mean >= 1.5 * d$baseline_creatinine_mean, 1, 0),
    NA_real_
  )
  d$lactate_nonclearance_0_24 <- ifelse(
    is.finite(d$baseline_lactate_mean) & d$baseline_lactate_mean >= 2 & is.finite(d$response_lactate_mean),
    ifelse(d$response_lactate_mean >= 0.9 * d$baseline_lactate_mean, 1, 0),
    NA_real_
  )
  early_org <- cbind(ifelse(is.finite(d$response_uo_ml_kg_h), d$response_uo_ml_kg_h < 0.5, NA), d$creatinine_worsening_0_24, d$lactate_nonclearance_0_24)
  d$organ_nonrecovery_0_24 <- ifelse(rowSums(!is.na(early_org)) > 0, ifelse(rowSums(early_org == 1, na.rm = TRUE) > 0, 1, 0), NA_real_)
  d$map_restored_0_24 <- ifelse(is.finite(d$response_map_mean), ifelse(d$response_map_mean >= 65, 1, 0), NA_real_)
  d$discordant_group <- NA_character_
  d$discordant_group[d$map_restored_0_24 == 0] <- "MAP_not_restored"
  d$discordant_group[d$map_restored_0_24 == 1 & d$organ_nonrecovery_0_24 == 0] <- "MAP_restored_organ_recovered"
  d$discordant_group[d$map_restored_0_24 == 1 & d$organ_nonrecovery_0_24 == 1] <- "MAP_restored_organ_not_recovered"
  d
}

adjust_terms <- function(dat, exposure = "hrc_core_z") {
  terms <- c(exposure, term_for_numeric(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, term_for_numeric(dat, "index_hour_from_icu", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  terms <- unique(terms[nzchar(terms)])
  terms
}

fit_binary <- function(dat, outcome, database, analysis, exposure = "hrc_core_z") {
  x <- dat[complete.cases(dat[, c(outcome, exposure)]), , drop = FALSE]
  x <- x[is.finite(x[[outcome]]) & is.finite(x[[exposure]]), , drop = FALSE]
  if (nrow(x) < 300 || sum(x[[outcome]] == 1, na.rm = TRUE) < 30 || !has_variation(x[[outcome]])) return(data.frame())
  for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log")) {
    if (v %in% names(x) && is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  terms <- adjust_terms(x, exposure)
  fit <- glm(as.formula(paste(outcome, "~", paste(terms, collapse = " + "))), data = x, family = binomial())
  tab <- summary(fit)$coefficients
  ci <- confint.default(fit)
  data.frame(
    database = database,
    endpoint = outcome,
    analysis = analysis,
    n = nrow(x),
    events = sum(x[[outcome]] == 1, na.rm = TRUE),
    beta = tab[exposure, "Estimate"],
    se = tab[exposure, "Std. Error"],
    odds_ratio = exp(tab[exposure, "Estimate"]),
    ci_low = exp(ci[exposure, 1]),
    ci_high = exp(ci[exposure, 2]),
    p_value = tab[exposure, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

net_benefit <- function(y, p, thresholds) {
  out <- lapply(thresholds, function(pt) {
    pred <- p >= pt
    tp <- sum(pred & y == 1, na.rm = TRUE)
    fp <- sum(pred & y == 0, na.rm = TRUE)
    n <- sum(!is.na(y) & !is.na(p))
    data.frame(threshold = pt, net_benefit = tp / n - fp / n * pt / (1 - pt))
  })
  do.call(rbind, out)
}

fit_dca <- function(dat, database, outcome, label, thresholds) {
  x <- dat[complete.cases(dat[, c(outcome, "hrc_core_z")]), , drop = FALSE]
  x <- x[is.finite(x[[outcome]]) & is.finite(x$hrc_core_z), , drop = FALSE]
  for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log")) {
    if (v %in% names(x) && is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  base_terms <- setdiff(adjust_terms(x, "hrc_core_z"), "hrc_core_z")
  base_fit <- glm(as.formula(paste(outcome, "~", paste(base_terms, collapse = " + "))), data = x, family = binomial())
  hrc_fit <- glm(as.formula(paste(outcome, "~", paste(c("hrc_core_z", base_terms), collapse = " + "))), data = x, family = binomial())
  y <- x[[outcome]]
  base_nb <- net_benefit(y, predict(base_fit, type = "response"), thresholds)
  hrc_nb <- net_benefit(y, predict(hrc_fit, type = "response"), thresholds)
  treat_all_nb <- sapply(thresholds, function(pt) mean(y == 1) - mean(y == 0) * pt / (1 - pt))
  rows <- data.frame(
    database = database,
    outcome = outcome,
    label = label,
    threshold = thresholds,
    baseline_net_benefit = base_nb$net_benefit,
    hrc_net_benefit = hrc_nb$net_benefit,
    delta_net_benefit = hrc_nb$net_benefit - base_nb$net_benefit,
    treat_all_net_benefit = treat_all_nb,
    treat_none_net_benefit = 0,
    stringsAsFactors = FALSE
  )
  data.frame(
    database = database,
    outcome = outcome,
    label = label,
    n = nrow(x),
    events = sum(y == 1, na.rm = TRUE),
    threshold_min = min(thresholds),
    threshold_max = max(thresholds),
    mean_delta_net_benefit = mean(rows$delta_net_benefit, na.rm = TRUE),
    positive_delta_thresholds = sum(rows$delta_net_benefit > 0, na.rm = TRUE),
    total_thresholds = length(thresholds),
    stringsAsFactors = FALSE
  ) -> summary
  list(curve = rows, summary = summary)
}

summarize_phenotype <- function(dat, database) {
  endpoints <- c(
    "outcome", "persistent_vaso_24_72", "oliguria_24_72",
    "creatinine_worsening_24_72", "lactate_nonclearance_24_72",
    "organ_nonrecovery_24_72", "support_or_organ_nonrecovery_24_72"
  )
  bind_rows(lapply(c(0, 1), function(low) {
    x <- dat[dat$low_hrc_q1 == low, , drop = FALSE]
    vals <- lapply(endpoints, function(ep) {
      y <- x[[ep]]
      data.frame(endpoint = ep, evaluable = sum(!is.na(y)), events = sum(y == 1, na.rm = TRUE), rate = mean(y == 1, na.rm = TRUE))
    })
    out <- do.call(rbind, vals)
    out$database <- database
    out$low_hrc_q1 <- low
    out$n <- nrow(x)
    out[, c("database", "low_hrc_q1", "n", "endpoint", "evaluable", "events", "rate")]
  }))
}

discordant_summary <- function(dat, database) {
  x <- dat[!is.na(dat$discordant_group), , drop = FALSE]
  groups <- c("MAP_restored_organ_recovered", "MAP_restored_organ_not_recovered", "MAP_not_restored")
  bind_rows(lapply(groups, function(g) {
    y <- x[x$discordant_group == g, , drop = FALSE]
    data.frame(
      database = database,
      discordant_group = g,
      n = nrow(y),
      deaths = sum(y$outcome == 1, na.rm = TRUE),
      mortality = mean(y$outcome == 1, na.rm = TRUE),
      post24_organ_nonrecovery_events = sum(y$organ_nonrecovery_24_72 == 1, na.rm = TRUE),
      post24_organ_nonrecovery_rate = mean(y$organ_nonrecovery_24_72 == 1, na.rm = TRUE),
      median_hrc = median(y$hrc_core_z, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

discordant_or <- function(dat, database) {
  x <- dat[!is.na(dat$discordant_group), , drop = FALSE]
  x$discordant_group <- factor(
    x$discordant_group,
    levels = c("MAP_restored_organ_recovered", "MAP_restored_organ_not_recovered", "MAP_not_restored")
  )
  x <- x[complete.cases(x[, c("outcome", "discordant_group")]), , drop = FALSE]
  if (nrow(x) < 300) return(data.frame())
  for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log")) {
    if (v %in% names(x) && is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  terms <- c("discordant_group", setdiff(adjust_terms(x, "hrc_core_z"), "hrc_core_z"))
  fit <- glm(as.formula(paste("outcome ~", paste(terms, collapse = " + "))), data = x, family = binomial())
  tab <- summary(fit)$coefficients
  ci <- confint.default(fit)
  terms_keep <- grep("^discordant_group", rownames(tab), value = TRUE)
  data.frame(
    database = database,
    term = terms_keep,
    n = nrow(x),
    deaths = sum(x$outcome == 1, na.rm = TRUE),
    odds_ratio = exp(tab[terms_keep, "Estimate"]),
    ci_low = exp(ci[terms_keep, 1]),
    ci_high = exp(ci[terms_keep, 2]),
    p_value = tab[terms_keep, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

cvp_mpp_analysis <- function(dat) {
  x <- dat[dat$database == "MIMIC-IV", , drop = FALSE]
  x$baseline_mpp <- x$baseline_map_mean - x$baseline_cvp_mean
  x$response_mpp <- x$response_map_mean - x$response_cvp_mean
  x$high_response_cvp <- ifelse(is.finite(x$response_cvp_mean), ifelse(x$response_cvp_mean >= 12, 1, 0), NA_real_)
  x$low_response_mpp <- ifelse(is.finite(x$response_mpp), ifelse(x$response_mpp < 60, 1, 0), NA_real_)
  x$low_hrc <- ifelse(x$hrc_core_z <= quantile(x$hrc_core_z, 0.25, na.rm = TRUE), 1, 0)
  x$congestive_low_mpp <- ifelse(x$high_response_cvp == 1 | x$low_response_mpp == 1, 1, 0)
  x$joint_cvp_hrc <- NA_character_
  x$joint_cvp_hrc[x$low_hrc == 0 & x$congestive_low_mpp == 0] <- "nonlow_HRC_no_congestion_lowMPP"
  x$joint_cvp_hrc[x$low_hrc == 1 & x$congestive_low_mpp == 0] <- "low_HRC_only"
  x$joint_cvp_hrc[x$low_hrc == 0 & x$congestive_low_mpp == 1] <- "congestion_lowMPP_only"
  x$joint_cvp_hrc[x$low_hrc == 1 & x$congestive_low_mpp == 1] <- "low_HRC_plus_congestion_lowMPP"

  sub <- x[is.finite(x$response_cvp_mean) & is.finite(x$response_mpp), , drop = FALSE]
  hrc_or <- fit_binary(sub, "outcome", "MIMIC-IV", "CVP/MPP subcohort mortality HRC association")
  org_or <- fit_binary(sub, "organ_nonrecovery_24_72", "MIMIC-IV", "CVP/MPP subcohort post24 organ nonrecovery HRC association")
  joint <- do.call(rbind, lapply(unique(na.omit(sub$joint_cvp_hrc)), function(g) {
    y <- sub[sub$joint_cvp_hrc == g, , drop = FALSE]
    data.frame(
      group = g,
      n = nrow(y),
      deaths = sum(y$outcome == 1, na.rm = TRUE),
      mortality = mean(y$outcome == 1, na.rm = TRUE),
      post24_organ_nonrecovery = mean(y$organ_nonrecovery_24_72 == 1, na.rm = TRUE),
      median_response_cvp = median(y$response_cvp_mean, na.rm = TRUE),
      median_response_mpp = median(y$response_mpp, na.rm = TRUE),
      median_hrc = median(y$hrc_core_z, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  list(subcohort = sub, or = rbind(hrc_or, org_or), joint = joint)
}

mcs_analysis <- function(dat) {
  x <- dat[dat$database == "MIMIC-IV" & dat$mcs_any_0_24h == 1, , drop = FALSE]
  rows <- rbind(
    fit_binary(x, "outcome", "MIMIC-IV", "MCS 0-24h mortality HRC association"),
    fit_binary(x, "organ_nonrecovery_24_72", "MIMIC-IV", "MCS 0-24h post24 organ nonrecovery HRC association"),
    fit_binary(x, "support_or_organ_nonrecovery_24_72", "MIMIC-IV", "MCS 0-24h support/organ nonrecovery HRC association")
  )
  disc <- discordant_summary(x, "MIMIC-IV MCS 0-24h")
  list(or = rows, discordant = disc)
}

databases <- c("MIMIC-IV", "eICU", "SICdb")
data_list <- setNames(lapply(databases, read_dataset), databases)
all_data <- bind_rows(data_list)
write.csv(all_data, file.path(outdir, "mechanism_analysis_dataset.csv"), row.names = FALSE)

endpoints <- c(
  "persistent_vaso_24_72",
  "oliguria_24_72",
  "creatinine_worsening_24_72",
  "lactate_nonclearance_24_72",
  "persistent_hyperlactatemia_24_72",
  "organ_nonrecovery_24_72",
  "support_or_organ_nonrecovery_24_72",
  "post24_rrt_any",
  "post24_invasive_vent_any"
)

endpoint_or <- bind_rows(lapply(databases, function(db) {
  d <- data_list[[db]]
  bind_rows(lapply(endpoints, function(ep) {
    if (!(ep %in% names(d))) return(data.frame())
    fit_binary(d, ep, db, paste0("Post-24h endpoint: ", ep))
  }))
}))
endpoint_meta <- bind_rows(lapply(split(endpoint_or, endpoint_or$endpoint), function(x) {
  meta_summary(x, unique(x$endpoint))
}))

phenotype_table <- bind_rows(lapply(databases, function(db) summarize_phenotype(data_list[[db]], db)))
discordant_rates <- bind_rows(lapply(databases, function(db) discordant_summary(data_list[[db]], db)))
discordant_models <- bind_rows(lapply(databases, function(db) discordant_or(data_list[[db]], db)))

dca_results <- lapply(databases, function(db) {
  d <- data_list[[db]]
  mortality <- fit_dca(d, db, "outcome", "Hospital mortality", seq(0.05, 0.40, by = 0.05))
  organ <- fit_dca(d, db, "organ_nonrecovery_24_72", "Post-24h organ nonrecovery", seq(0.10, 0.70, by = 0.05))
  list(curve = rbind(mortality$curve, organ$curve), summary = rbind(mortality$summary, organ$summary))
})
dca_curve <- bind_rows(lapply(dca_results, `[[`, "curve"))
dca_summary <- bind_rows(lapply(dca_results, `[[`, "summary"))

cvp <- cvp_mpp_analysis(all_data)
mcs <- mcs_analysis(all_data)

write.csv(endpoint_or, file.path(outdir, "post24_endpoint_hrc_or.csv"), row.names = FALSE)
write.csv(endpoint_meta, file.path(outdir, "post24_endpoint_hrc_meta.csv"), row.names = FALSE)
write.csv(phenotype_table, file.path(outdir, "low_hrc_post24_phenotype_table.csv"), row.names = FALSE)
write.csv(discordant_rates, file.path(outdir, "discordant_pressure_organ_recovery_rates.csv"), row.names = FALSE)
write.csv(discordant_models, file.path(outdir, "discordant_pressure_organ_recovery_or.csv"), row.names = FALSE)
write.csv(dca_curve, file.path(outdir, "decision_curve_net_benefit.csv"), row.names = FALSE)
write.csv(dca_summary, file.path(outdir, "decision_curve_summary.csv"), row.names = FALSE)
write.csv(cvp$or, file.path(outdir, "mimic_cvp_mpp_hrc_or.csv"), row.names = FALSE)
write.csv(cvp$joint, file.path(outdir, "mimic_cvp_mpp_joint_phenotype.csv"), row.names = FALSE)
write.csv(mcs$or, file.path(outdir, "mimic_mcs_mechanism_hrc_or.csv"), row.names = FALSE)
write.csv(mcs$discordant, file.path(outdir, "mimic_mcs_discordant_rates.csv"), row.names = FALSE)

endpoint_md <- endpoint_or
endpoint_md$OR_95CI <- fmt_or(endpoint_md$odds_ratio, endpoint_md$ci_low, endpoint_md$ci_high)
endpoint_md$p_value <- fmt_p(endpoint_md$p_value)
write_md_table(
  endpoint_md[, c("database", "endpoint", "n", "events", "OR_95CI", "p_value")],
  file.path(outdir, "post24_endpoint_hrc_or.md"),
  "Post-24h organ dysfunction endpoint associations"
)

meta_md <- endpoint_meta
meta_md$random_OR_95CI <- fmt_or(meta_md$random_or, meta_md$random_ci_low, meta_md$random_ci_high)
meta_md$i2_percent <- sprintf("%.1f", meta_md$i2_percent)
write_md_table(
  meta_md[, c("analysis", "k", "random_OR_95CI", "i2_percent")],
  file.path(outdir, "post24_endpoint_hrc_meta.md"),
  "Cross-database meta-analysis for post-24h endpoints"
)

disc_md <- discordant_rates
disc_md$mortality <- sprintf("%.3f", disc_md$mortality)
disc_md$post24_organ_nonrecovery_rate <- sprintf("%.3f", disc_md$post24_organ_nonrecovery_rate)
disc_md$median_hrc <- sprintf("%.2f", disc_md$median_hrc)
write_md_table(
  disc_md,
  file.path(outdir, "discordant_pressure_organ_recovery_rates.md"),
  "MAP restoration versus organ recovery discordance"
)

dca_md <- dca_summary
dca_md$mean_delta_net_benefit <- sprintf("%.4f", dca_md$mean_delta_net_benefit)
write_md_table(
  dca_md,
  file.path(outdir, "decision_curve_summary.md"),
  "Decision curve summary"
)

cvp_md <- cvp$joint
cvp_md$mortality <- sprintf("%.3f", cvp_md$mortality)
cvp_md$post24_organ_nonrecovery <- sprintf("%.3f", cvp_md$post24_organ_nonrecovery)
cvp_md$median_response_cvp <- sprintf("%.1f", cvp_md$median_response_cvp)
cvp_md$median_response_mpp <- sprintf("%.1f", cvp_md$median_response_mpp)
cvp_md$median_hrc <- sprintf("%.2f", cvp_md$median_hrc)
write_md_table(
  cvp_md,
  file.path(outdir, "mimic_cvp_mpp_joint_phenotype.md"),
  "MIMIC CVP/MPP joint phenotype"
)

mcs_md <- mcs$or
mcs_md$OR_95CI <- fmt_or(mcs_md$odds_ratio, mcs_md$ci_low, mcs_md$ci_high)
mcs_md$p_value <- fmt_p(mcs_md$p_value)
write_md_table(
  mcs_md[, c("database", "endpoint", "analysis", "n", "events", "OR_95CI", "p_value")],
  file.path(outdir, "mimic_mcs_mechanism_hrc_or.md"),
  "MIMIC MCS mechanism endpoint associations"
)

key_lines <- c(
  "# Mechanism and clinical utility key results",
  "",
  "## Post-24h organ dysfunction",
  "",
  paste(apply(endpoint_md[endpoint_md$endpoint %in% c("organ_nonrecovery_24_72", "support_or_organ_nonrecovery_24_72"), ], 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["endpoint"]], ": n=", row[["n"]], ", events=", row[["events"]], ", OR=", row[["OR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n"),
  "",
  "## Cross-database endpoint meta-analysis",
  "",
  paste(apply(meta_md[meta_md$analysis %in% c("organ_nonrecovery_24_72", "support_or_organ_nonrecovery_24_72", "creatinine_worsening_24_72", "oliguria_24_72"), ], 1, function(row) {
    paste0("- ", row[["analysis"]], ": random-effects OR=", row[["random_OR_95CI"]], ", I2=", row[["i2_percent"]], "%")
  }), collapse = "\n"),
  "",
  "## Discordant pressure-organ recovery",
  "",
  paste(apply(disc_md[disc_md$discordant_group == "MAP_restored_organ_not_recovered", ], 1, function(row) {
    paste0("- ", row[["database"]], ": MAP restored but organ not recovered, n=", row[["n"]], ", mortality=", row[["mortality"]], ", post24 organ nonrecovery=", row[["post24_organ_nonrecovery_rate"]])
  }), collapse = "\n"),
  "",
  "## Clinical utility",
  "",
  paste(apply(dca_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["label"]], ": mean delta net benefit=", row[["mean_delta_net_benefit"]], " across ", row[["positive_delta_thresholds"]], "/", row[["total_thresholds"]], " thresholds")
  }), collapse = "\n"),
  "",
  "## MIMIC CVP/MPP mechanism",
  "",
  paste(apply(cvp_md, 1, function(row) {
    paste0("- ", row[["group"]], ": n=", row[["n"]], ", mortality=", row[["mortality"]], ", post24 organ nonrecovery=", row[["post24_organ_nonrecovery"]], ", median CVP=", row[["median_response_cvp"]], ", median MPP=", row[["median_response_mpp"]])
  }), collapse = "\n"),
  "",
  "## MIMIC MCS mechanism",
  "",
  paste(apply(mcs_md, 1, function(row) {
    paste0("- ", row[["analysis"]], " / ", row[["endpoint"]], ": n=", row[["n"]], ", events=", row[["events"]], ", OR=", row[["OR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n")
)
writeLines(key_lines, file.path(outdir, "mechanism_clinical_utility_key_results.md"))

cat("Wrote mechanism and clinical utility outputs to", outdir, "\n")
