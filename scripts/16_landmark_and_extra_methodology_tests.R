#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

outdir <- "outputs/landmark_sensitivity"
extra_dir <- "outputs/extra_methodology_tests"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(extra_dir, recursive = TRUE, showWarnings = FALSE)
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

read_source <- function(db) {
  file <- switch(
    db,
    "MIMIC-IV" = "outputs/landmark_sensitivity/mimic_landmark_source.csv",
    "eICU" = "outputs/landmark_sensitivity/eicu_landmark_source.csv",
    "SICdb" = "outputs/landmark_sensitivity/sicdb_landmark_source.csv"
  )
  d <- read.csv(file, stringsAsFactors = FALSE)
  for (v in c(
    "age", "male", "weight_kg", "severity_primary", "severity_secondary",
    "severity_cardiovascular", "severity_renal", "index_hour_from_icu",
    "icu_los_hours", "landmark_eligible", "early_death_before_landmark",
    "early_icu_discharge_alive_before_landmark", "hospital_mortality_anytime",
    "patient_early_support_rank", "heart_surgery"
  )) {
    if (v %in% names(d)) d[[v]] <- num(d[[v]])
  }
  d$stay_key <- as.character(d$stay_key)
  d
}

read_scores <- function(db) {
  if (db == "MIMIC-IV") {
    d <- read.csv("outputs/mimic_formal/mimic_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = db,
      stay_key = as.character(d$stay_id),
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      hrc_map_residual_z = num(d$hrc_map_residual_z),
      hrc_vaso_residual_z = num(d$hrc_vaso_residual_z),
      hrc_uo_residual_z = num(d$hrc_uo_residual_z),
      age = num(d$anchor_age),
      male = ifelse(d$gender == "M", 1, 0),
      severity_primary = num(d$first_day_sofa),
      severity_secondary = NA_real_,
      baseline_vaso_log = log1p(pmax(num(d$baseline_neq_avg), 0)),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      baseline_vaso_raw = num(d$baseline_neq_avg),
      stringsAsFactors = FALSE
    )
  } else if (db == "eICU") {
    d <- read.csv("outputs/eicu_formal/eicu_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = db,
      stay_key = as.character(d$patientunitstayid),
      outcome = num(d$hospital_mortality_after_landmark),
      hrc_core_z = num(d$hrc_core_z),
      hrc_map_residual_z = num(d$hrc_map_residual_z),
      hrc_vaso_residual_z = num(d$hrc_vaso_residual_z),
      hrc_uo_residual_z = num(d$hrc_uo_residual_z),
      age = num(d$age),
      male = ifelse(d$gender == "Male", 1, 0),
      severity_primary = num(d$acutephysiologyscore),
      severity_secondary = num(d$apachescore),
      baseline_vaso_log = log1p(pmax(num(d$baseline_vaso_burden), 0)),
      baseline_map_mean = num(d$baseline_map_mean),
      baseline_uo_log = log1p(pmax(num(d$baseline_uo_ml_kg_h), 0)),
      baseline_vaso_raw = num(d$baseline_vaso_burden),
      stringsAsFactors = FALSE
    )
  } else if (db == "SICdb") {
    d <- read.csv("outputs/sicdb_formal/sicdb_hrc_formal_scores.csv", stringsAsFactors = FALSE)
    out <- data.frame(
      database = db,
      stay_key = as.character(d$caseid),
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
      baseline_vaso_raw = num(d$baseline_vaso_burden),
      stringsAsFactors = FALSE
    )
  }
  out
}

read_formal <- function(db) {
  file <- switch(
    db,
    "MIMIC-IV" = "outputs/mimic_formal/mimic_formal_cohort.csv",
    "eICU" = "outputs/eicu_formal/eicu_formal_cohort.csv",
    "SICdb" = "outputs/sicdb_formal/sicdb_formal_cohort.csv"
  )
  d <- read.csv(file, stringsAsFactors = FALSE)
  if (db == "MIMIC-IV") {
    d$stay_key <- as.character(d$stay_id)
    d$vaso_reduction <- num(d$log_neq_reduction)
  } else if (db == "eICU") {
    d$stay_key <- as.character(d$patientunitstayid)
    d$vaso_reduction <- num(d$log_vaso_burden_reduction)
  } else {
    d$stay_key <- as.character(d$caseid)
    d$vaso_reduction <- num(d$log_vaso_burden_reduction)
  }
  count_cols <- grep("(_n$|baseline_.*_n$|response_.*_n$)", names(d), value = TRUE)
  for (v in c(
    count_cols, "delta_lactate_reduction", "delta_creatinine_reduction",
    "log_uo_recovery", "vaso_reduction", "delta_map",
    "baseline_lactate_mean", "response_lactate_mean",
    "baseline_creatinine_mean", "response_creatinine_mean"
  )) {
    if (v %in% names(d)) d[[v]] <- num(d[[v]])
  }
  d
}

model_terms <- function(dat, exposure = "hrc_core_z", include_measurement = FALSE) {
  terms <- c(exposure, "ns(age, df = 3)")
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, term_for_numeric(dat, "severity_primary", 4))
  terms <- c(terms, term_for_numeric(dat, "severity_secondary", 4))
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) terms <- c(terms, "heart_surgery")
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  if (include_measurement) terms <- c(terms, term_for_numeric(dat, "measurement_count_log", 4))
  terms
}

fit_or <- function(dat, exposure = "hrc_core_z", analysis = "", db = "", weights = NULL, include_measurement = FALSE, simplified = FALSE) {
  candidate_cols <- c(
    "outcome", exposure, "age", "male", "severity_primary", "severity_secondary",
    "heart_surgery", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log",
    "measurement_count_log"
  )
  cols <- unique(candidate_cols)
  cols <- cols[cols %in% names(dat)]
  x <- dat[complete.cases(dat[, c("outcome", exposure)]), cols, drop = FALSE]
  x <- x[is.finite(x$outcome) & is.finite(x[[exposure]]), , drop = FALSE]
  for (v in setdiff(names(x), c("outcome", exposure))) {
    if (is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  if (simplified) {
    terms <- c(exposure, "ns(age, df = 3)")
    if ("male" %in% names(x) && has_variation(x$male)) terms <- c(terms, "male")
    terms <- c(terms, term_for_numeric(x, "severity_primary", 4))
    terms <- c(terms, term_for_numeric(x, "severity_secondary", 4))
  } else {
    terms <- model_terms(x, exposure, include_measurement)
  }
  form <- as.formula(paste("outcome ~", paste(terms, collapse = " + ")))
  if (!is.null(weights)) {
    x$w <- weights[match(rownames(x), names(weights))]
    x$w <- impute_median(x$w)
    fit <- suppressWarnings(glm(form, data = x, weights = w, family = quasibinomial()))
  } else {
    fit <- glm(form, data = x, family = binomial())
  }
  tab <- summary(fit)$coefficients
  ci <- confint.default(fit)
  p_col <- if ("Pr(>|t|)" %in% colnames(tab)) "Pr(>|t|)" else "Pr(>|z|)"
  data.frame(
    database = db,
    analysis = analysis,
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

fit_logit_event <- function(src, event, db) {
  x <- src
  x$event <- num(x[[event]])
  terms <- c("ns(age, df = 3)")
  if ("male" %in% names(x) && has_variation(x$male)) terms <- c(terms, "male")
  if ("severity_primary" %in% names(x) && has_variation(x$severity_primary)) terms <- c(terms, "ns(severity_primary, df = 4)")
  if ("severity_secondary" %in% names(x) && has_variation(x$severity_secondary)) terms <- c(terms, "ns(severity_secondary, df = 4)")
  if ("heart_surgery" %in% names(x) && has_variation(x$heart_surgery)) terms <- c(terms, "heart_surgery")
  if ("index_hour_from_icu" %in% names(x) && has_variation(x$index_hour_from_icu)) terms <- c(terms, "ns(index_hour_from_icu, df = 3)")
  cols <- unique(c("event", all.vars(as.formula(paste("~", paste(terms, collapse = "+"))))))
  x <- x[complete.cases(x[, c("event", "age")]), cols, drop = FALSE]
  for (v in setdiff(names(x), "event")) if (is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  fit <- glm(as.formula(paste("event ~", paste(terms, collapse = "+"))), data = x, family = binomial())
  p <- predict(fit, type = "response")
  data.frame(
    database = db,
    event = event,
    n = nrow(x),
    events = sum(x$event == 1),
    event_rate = mean(x$event == 1),
    model_auc = auc_rank(x$event, p),
    stringsAsFactors = FALSE
  )
}

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

meta_summary <- function(rows, analysis) {
  rows <- rows[is.finite(rows$beta) & is.finite(rows$se) & rows$se > 0, ]
  w <- 1 / rows$se^2
  beta_fe <- sum(w * rows$beta) / sum(w)
  q <- if (nrow(rows) > 1) sum(w * (rows$beta - beta_fe)^2) else NA_real_
  df <- nrow(rows) - 1
  tau2 <- if (nrow(rows) > 1) max(0, (q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (rows$se^2 + tau2)
  beta <- sum(wr * rows$beta) / sum(wr)
  se <- sqrt(1 / sum(wr))
  data.frame(
    analysis = analysis,
    k = nrow(rows),
    random_or = exp(beta),
    random_ci_low = exp(beta - 1.96 * se),
    random_ci_high = exp(beta + 1.96 * se),
    random_p = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    i2_percent = ifelse(is.finite(q) && q > 0, max(0, (q - df) / q) * 100, 0),
    stringsAsFactors = FALSE
  )
}

databases <- c("MIMIC-IV", "eICU", "SICdb")
sources <- setNames(lapply(databases, read_source), databases)
scores <- setNames(lapply(databases, read_scores), databases)
formal <- setNames(lapply(databases, read_formal), databases)

landmark_profile <- do.call(rbind, lapply(databases, function(db) {
  s <- sources[[db]]
  groups <- list(
    landmark_eligible = s[s$landmark_eligible == 1, ],
    early_death = s[s$early_death_before_landmark == 1, ],
    early_discharge_alive = s[s$early_icu_discharge_alive_before_landmark == 1, ]
  )
  do.call(rbind, lapply(names(groups), function(g) {
    x <- groups[[g]]
    data.frame(
      database = db,
      group = g,
      n = nrow(x),
      percent_of_early_support = 100 * nrow(x) / nrow(s),
      mortality_anytime = mean(x$hospital_mortality_anytime == 1, na.rm = TRUE),
      mean_age = mean(x$age, na.rm = TRUE),
      mean_severity = mean(x$severity_primary, na.rm = TRUE),
      median_index_hour = median(x$index_hour_from_icu, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(landmark_profile, file.path(outdir, "landmark_early_event_profile.csv"), row.names = FALSE)
write_md_table(landmark_profile, file.path(outdir, "landmark_early_event_profile.md"), "Landmark early death/discharge profile")

early_event_models <- do.call(rbind, lapply(databases, function(db) {
  rbind(
    fit_logit_event(sources[[db]], "landmark_eligible", db),
    fit_logit_event(sources[[db]], "early_death_before_landmark", db),
    fit_logit_event(sources[[db]], "early_icu_discharge_alive_before_landmark", db)
  )
}))
write.csv(early_event_models, file.path(outdir, "early_death_discharge_prediction_models.csv"), row.names = FALSE)
write_md_table(early_event_models, file.path(outdir, "early_death_discharge_prediction_models.md"), "Early death/discharge separate prediction models")

ipcw_rows <- do.call(rbind, lapply(databases, function(db) {
  src <- sources[[db]]
  sc <- scores[[db]]
  src$baseline_prob_row <- seq_len(nrow(src))
  terms <- c("ns(age, df = 3)")
  if (has_variation(src$male)) terms <- c(terms, "male")
  if (has_variation(src$severity_primary)) terms <- c(terms, "ns(severity_primary, df = 4)")
  if ("severity_secondary" %in% names(src) && has_variation(src$severity_secondary)) terms <- c(terms, "ns(severity_secondary, df = 4)")
  if ("heart_surgery" %in% names(src) && has_variation(src$heart_surgery)) terms <- c(terms, "heart_surgery")
  terms <- c(terms, "ns(index_hour_from_icu, df = 3)")
  cols <- unique(c("landmark_eligible", "stay_key", all.vars(as.formula(paste("~", paste(terms, collapse = "+"))))))
  x <- src[, cols, drop = FALSE]
  for (v in setdiff(names(x), c("landmark_eligible", "stay_key"))) if (is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  fit <- glm(as.formula(paste("landmark_eligible ~", paste(terms, collapse = "+"))), data = x, family = binomial())
  src$p_landmark <- pmin(pmax(predict(fit, newdata = x, type = "response"), 0.02), 0.98)
  w_lookup <- setNames(mean(src$landmark_eligible == 1) / src$p_landmark, src$stay_key)
  w <- w_lookup[sc$stay_key]
  lo <- quantile(w, 0.01, na.rm = TRUE)
  hi <- quantile(w, 0.99, na.rm = TRUE)
  w <- pmin(pmax(w, lo), hi)
  names(w) <- seq_len(nrow(sc))
  rownames(sc) <- names(w)
  fit_or(sc, "hrc_core_z", "Landmark IPCW HRC", db, weights = w)
}))
write.csv(ipcw_rows, file.path(outdir, "landmark_ipcw_hrc_or.csv"), row.names = FALSE)
write_md_table(ipcw_rows, file.path(outdir, "landmark_ipcw_hrc_or.md"), "Landmark IPCW HRC sensitivity")

composite_rows <- do.call(rbind, lapply(databases, function(db) {
  src <- sources[[db]]
  sc <- scores[[db]]
  q_worst <- quantile(sc$hrc_core_z, 0.01, na.rm = TRUE) - 0.1
  q_best <- quantile(sc$hrc_core_z, 0.99, na.rm = TRUE) + 0.1
  early_death <- src[src$early_death_before_landmark == 1, ]
  death_rows <- data.frame(
    outcome = 1,
    hrc_core_z = q_worst,
    age = early_death$age,
    male = early_death$male,
    severity_primary = early_death$severity_primary,
    severity_secondary = if ("severity_secondary" %in% names(early_death)) early_death$severity_secondary else NA_real_
  )
  early_discharge <- src[src$early_icu_discharge_alive_before_landmark == 1, ]
  discharge_rows <- data.frame(
    outcome = 0,
    hrc_core_z = q_best,
    age = early_discharge$age,
    male = early_discharge$male,
    severity_primary = early_discharge$severity_primary,
    severity_secondary = if ("severity_secondary" %in% names(early_discharge)) early_discharge$severity_secondary else NA_real_
  )
  core_rows <- sc[, c("outcome", "hrc_core_z", "age", "male", "severity_primary", "severity_secondary")]
  a <- rbind(core_rows, death_rows)
  b <- rbind(core_rows, death_rows, discharge_rows)
  rbind(
    fit_or(a, "hrc_core_z", "Early deaths assigned worst observed HRC", db, simplified = TRUE),
    fit_or(b, "hrc_core_z", "Early deaths worst + early discharges best HRC", db, simplified = TRUE)
  )
}))
write.csv(composite_rows, file.path(outdir, "early_event_composite_worst_response_or.csv"), row.names = FALSE)
write_md_table(composite_rows, file.path(outdir, "early_event_composite_worst_response_or.md"), "Early event composite worst-response sensitivity")

pca_rows <- do.call(rbind, lapply(databases, function(db) {
  sc <- scores[[db]]
  comp <- sc[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")]
  keep <- complete.cases(comp)
  pc <- prcomp(comp[keep, ], center = TRUE, scale. = TRUE)
  sc$hrc_pca_z <- NA_real_
  pc1 <- pc$x[, 1]
  if (cor(pc1, sc$hrc_core_z[keep], use = "complete.obs") < 0) pc1 <- -pc1
  sc$hrc_pca_z[keep] <- zscore(pc1)
  row <- fit_or(sc, "hrc_pca_z", "PCA-weighted HRC", db)
  loadings <- paste(names(pc$rotation[, 1]), sprintf("%.3f", pc$rotation[, 1]), collapse = "; ")
  row$pca_variance_explained <- pc$sdev[1]^2 / sum(pc$sdev^2)
  row$pca_loadings <- loadings
  row
}))
write.csv(pca_rows, file.path(extra_dir, "pca_weighted_hrc_or.csv"), row.names = FALSE)
write_md_table(pca_rows, file.path(extra_dir, "pca_weighted_hrc_or.md"), "PCA-weighted HRC sensitivity")

measurement_rows <- do.call(rbind, lapply(databases, function(db) {
  sc <- scores[[db]]
  fm <- formal[[db]]
  counts <- grep("_n$", names(fm), value = TRUE)
  fm$measurement_count_log <- log1p(rowSums(fm[, counts, drop = FALSE], na.rm = TRUE))
  merged <- merge(sc, fm[, c("stay_key", "measurement_count_log")], by = "stay_key", all.x = TRUE)
  rbind(
    fit_or(merged, "hrc_core_z", "No measurement-frequency adjustment", db),
    fit_or(merged, "hrc_core_z", "Adjusted for measurement frequency", db, include_measurement = TRUE)
  )
}))
write.csv(measurement_rows, file.path(extra_dir, "measurement_frequency_adjusted_hrc_or.csv"), row.names = FALSE)
write_md_table(measurement_rows, file.path(extra_dir, "measurement_frequency_adjusted_hrc_or.md"), "Measurement-frequency adjusted HRC sensitivity")

construct_rows <- do.call(rbind, lapply(databases, function(db) {
  sc <- scores[[db]]
  fm <- formal[[db]]
  merged <- merge(sc, fm, by = "stay_key", all.x = TRUE)
  targets <- c("delta_lactate_reduction", "delta_creatinine_reduction", "vaso_reduction", "log_uo_recovery", "delta_map")
  do.call(rbind, lapply(targets[targets %in% names(merged)], function(t) {
    keep <- complete.cases(merged[, c("hrc_core_z", t)])
    rho <- suppressWarnings(cor(merged$hrc_core_z[keep], merged[[t]][keep], method = "spearman"))
    data.frame(
      database = db,
      construct_target = t,
      n = sum(keep),
      spearman_rho = rho,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(construct_rows, file.path(extra_dir, "construct_validation_correlations.csv"), row.names = FALSE)
write_md_table(construct_rows, file.path(extra_dir, "construct_validation_correlations.md"), "Construct validation correlations")

vaso_strata_rows <- do.call(rbind, lapply(databases, function(db) {
  sc <- scores[[db]]
  sc$vaso_stratum <- NA_character_
  idx <- which(!is.na(sc$baseline_vaso_raw))
  idx <- idx[order(sc$baseline_vaso_raw[idx], seq_along(idx))]
  if (length(idx) >= 3) {
    grp <- cut(
      seq_along(idx),
      breaks = c(0, floor(length(idx) / 3), floor(2 * length(idx) / 3), length(idx)),
      labels = c("low", "mid", "high"),
      include.lowest = TRUE
    )
    sc$vaso_stratum[idx] <- as.character(grp)
  }
  sc$vaso_stratum <- factor(sc$vaso_stratum, levels = c("low", "mid", "high"))
  do.call(rbind, lapply(levels(sc$vaso_stratum), function(level) {
    x <- sc[sc$vaso_stratum == level, ]
    fit_or(x, "hrc_core_z", paste0("Baseline vaso stratum: ", level), db)
  }))
}))
write.csv(vaso_strata_rows, file.path(extra_dir, "baseline_vaso_stratified_hrc_or.csv"), row.names = FALSE)
write_md_table(vaso_strata_rows, file.path(extra_dir, "baseline_vaso_stratified_hrc_or.md"), "Baseline vasopressor-burden stratified HRC sensitivity")

key <- c(
  "# Landmark and extra methodology key results",
  "",
  "## 24h landmark problem",
  paste(apply(ipcw_rows, 1, function(r) paste0("- ", r[["database"]], " landmark-IPCW OR: ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["ci_low"]]), as.numeric(r[["ci_high"]])), ".")), collapse = "\n"),
  "",
  "## Early event composite sensitivity",
  paste(apply(composite_rows, 1, function(r) paste0("- ", r[["database"]], " / ", r[["analysis"]], ": OR ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["ci_low"]]), as.numeric(r[["ci_high"]])), ".")), collapse = "\n"),
  "",
  "## PCA weighting",
  paste(apply(pca_rows, 1, function(r) paste0("- ", r[["database"]], " PCA-HRC OR: ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["ci_low"]]), as.numeric(r[["ci_high"]])), "; PC1 variance=", sprintf("%.1f%%", 100 * as.numeric(r[["pca_variance_explained"]])), ".")), collapse = "\n"),
  "",
  "## Measurement-frequency adjustment",
  paste(apply(measurement_rows[measurement_rows$analysis == "Adjusted for measurement frequency", ], 1, function(r) paste0("- ", r[["database"]], ": OR ", fmt_or(as.numeric(r[["odds_ratio"]]), as.numeric(r[["ci_low"]]), as.numeric(r[["ci_high"]])), ".")), collapse = "\n"),
  "",
  "## Construct validation",
  paste(apply(construct_rows[construct_rows$construct_target %in% c("delta_lactate_reduction", "delta_creatinine_reduction"), ], 1, function(r) paste0("- ", r[["database"]], " ", r[["construct_target"]], ": Spearman rho=", sprintf("%.3f", as.numeric(r[["spearman_rho"]])), ", n=", r[["n"]], ".")), collapse = "\n")
)
writeLines(key, file.path(extra_dir, "landmark_and_extra_methodology_key_results.md"))

cat("Landmark outputs:", outdir, "\n")
cat("Extra methodology outputs:", extra_dir, "\n")
cat("Key results:", file.path(extra_dir, "landmark_and_extra_methodology_key_results.md"), "\n")
