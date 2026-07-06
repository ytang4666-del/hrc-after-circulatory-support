#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

outdir <- "outputs/time_window_sensitivity"
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

read_database <- function(database) {
  file <- switch(
    database,
    "MIMIC-IV" = file.path(outdir, "mimic_time_window_cohort.csv"),
    "eICU" = file.path(outdir, "eicu_time_window_cohort.csv"),
    "SICdb" = file.path(outdir, "sicdb_time_window_cohort.csv")
  )
  d <- read.csv(file, stringsAsFactors = FALSE)
  numeric_cols <- c(
    "baseline_start_h", "baseline_end_h", "response_start_h", "response_end_h",
    "baseline_hours", "response_hours", "age", "male", "weight_kg",
    "severity_primary", "severity_secondary", "severity_respiration",
    "severity_cardiovascular", "severity_renal", "heart_surgery",
    "hospital_mortality_after_landmark", "baseline_map_mean",
    "response_map_mean", "delta_map", "baseline_map_n", "response_map_n",
    "baseline_vaso_burden", "response_vaso_burden",
    "delta_vaso_burden_reduction", "log_vaso_burden_reduction",
    "baseline_vaso_observed_units", "response_vaso_observed_units",
    "baseline_uo_ml", "response_uo_ml", "baseline_uo_ml_kg_h",
    "response_uo_ml_kg_h", "delta_uo_ml_kg_h", "log_uo_recovery",
    "baseline_uo_n", "response_uo_n"
  )
  for (col in intersect(numeric_cols, names(d))) d[[col]] <- num(d[[col]])
  d$database <- database
  d$weight_kg[d$weight_kg <= 0 | d$weight_kg > 300] <- NA_real_
  d$baseline_vaso_log <- log1p(pmax(d$baseline_vaso_burden, 0))
  d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))
  d
}

expected_terms <- function(dat) {
  terms <- c(term_for_numeric(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, term_for_numeric(dat, "weight_kg", 3))
  terms <- c(terms, term_for_numeric(dat, "severity_primary", 4))
  terms <- c(terms, term_for_numeric(dat, "severity_secondary", 4))
  terms <- c(terms, term_for_numeric(dat, "severity_respiration", 3))
  terms <- c(terms, term_for_numeric(dat, "severity_cardiovascular", 3))
  terms <- c(terms, term_for_numeric(dat, "severity_renal", 3))
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) terms <- c(terms, "heart_surgery")
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  unique(terms[nzchar(terms)])
}

mortality_terms <- function(dat) {
  terms <- c(term_for_numeric(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, term_for_numeric(dat, "severity_primary", 4))
  terms <- c(terms, term_for_numeric(dat, "severity_secondary", 4))
  if ("heart_surgery" %in% names(dat) && has_variation(dat$heart_surgery)) terms <- c(terms, "heart_surgery")
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  unique(terms[nzchar(terms)])
}

crossfit_component <- function(dat, outcome, folds, terms) {
  expected <- rep(NA_real_, nrow(dat))
  fold_rmse <- rep(NA_real_, max(folds))
  form <- as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
  for (fold in sort(unique(folds))) {
    train <- dat[folds != fold, , drop = FALSE]
    test <- dat[folds == fold, , drop = FALSE]
    fit <- lm(form, data = train)
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

fit_window <- function(d, database, window_strategy) {
  x <- d[d$window_strategy == window_strategy, , drop = FALSE]
  required <- c(
    "delta_map", "log_vaso_burden_reduction", "log_uo_recovery",
    "hospital_mortality_after_landmark", "baseline_map_mean",
    "baseline_vaso_burden", "baseline_uo_ml_kg_h"
  )
  x <- x[complete.cases(x[, required]), , drop = FALSE]
  x <- x[
    is.finite(x$delta_map) &
      is.finite(x$log_vaso_burden_reduction) &
      is.finite(x$log_uo_recovery) &
      is.finite(x$hospital_mortality_after_landmark),
    ,
    drop = FALSE
  ]
  if (nrow(x) < 100) stop(paste("Too few complete rows:", database, window_strategy))

  covars <- c(
    "age", "male", "weight_kg", "severity_primary", "severity_secondary",
    "severity_respiration", "severity_cardiovascular", "severity_renal",
    "heart_surgery", "baseline_map_mean", "baseline_vaso_log", "baseline_uo_log",
    "baseline_vaso_burden", "baseline_uo_ml_kg_h"
  )
  for (col in intersect(covars, names(x))) x[[col]] <- impute_median(x[[col]])

  folds <- sample(rep(1:5, length.out = nrow(x)))
  e_terms <- expected_terms(x)
  map_cf <- crossfit_component(x, "delta_map", folds, e_terms)
  vaso_cf <- crossfit_component(x, "log_vaso_burden_reduction", folds, e_terms)
  uo_cf <- crossfit_component(x, "log_uo_recovery", folds, e_terms)

  x$expected_delta_map <- map_cf$expected
  x$expected_log_vaso_burden_reduction <- vaso_cf$expected
  x$expected_log_uo_recovery <- uo_cf$expected
  x$hrc_map_residual_z <- map_cf$residual_z
  x$hrc_vaso_residual_z <- vaso_cf$residual_z
  x$hrc_uo_residual_z <- uo_cf$residual_z
  x$hrc_core_raw <- rowMeans(
    x[, c("hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")],
    na.rm = FALSE
  )
  x$hrc_core_z <- zscore(x$hrc_core_raw)

  m_terms <- c("hrc_core_z", mortality_terms(x))
  m_terms <- unique(m_terms[nzchar(m_terms)])
  fit <- glm(
    as.formula(paste("hospital_mortality_after_landmark ~", paste(m_terms, collapse = " + "))),
    data = x,
    family = binomial()
  )
  tab <- summary(fit)$coefficients
  ci <- confint.default(fit)

  quartile_cut <- quantile(x$hrc_core_z, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  x$hrc_quartile <- cut(
    x$hrc_core_z,
    breaks = c(-Inf, quartile_cut, Inf),
    labels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
    include.lowest = TRUE
  )
  qsum <- do.call(rbind, lapply(levels(x$hrc_quartile), function(q) {
    y <- x[x$hrc_quartile == q, , drop = FALSE]
    data.frame(
      database = database,
      window_strategy = window_strategy,
      hrc_quartile = q,
      n = nrow(y),
      deaths = sum(y$hospital_mortality_after_landmark == 1, na.rm = TRUE),
      mortality_rate = mean(y$hospital_mortality_after_landmark == 1, na.rm = TRUE),
      median_hrc = median(y$hrc_core_z, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(
    scores = x,
    or = data.frame(
      database = database,
      window_strategy = window_strategy,
      n = nrow(x),
      deaths = sum(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
      beta = tab["hrc_core_z", "Estimate"],
      se = tab["hrc_core_z", "Std. Error"],
      odds_ratio = exp(tab["hrc_core_z", "Estimate"]),
      ci_low = exp(ci["hrc_core_z", 1]),
      ci_high = exp(ci["hrc_core_z", 2]),
      p_value = tab["hrc_core_z", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    ),
    metrics = data.frame(
      database = database,
      window_strategy = window_strategy,
      component = c("MAP recovery", "Vasopressor-burden reduction", "Urine output recovery"),
      observed_variable = c("delta_map", "log_vaso_burden_reduction", "log_uo_recovery"),
      n = nrow(x),
      oof_rmse = c(map_cf$rmse, vaso_cf$rmse, uo_cf$rmse),
      oof_r2 = c(map_cf$r2_oof, vaso_cf$r2_oof, uo_cf$r2_oof),
      stringsAsFactors = FALSE
    ),
    quartile = qsum
  )
}

databases <- c("MIMIC-IV", "eICU", "SICdb")
window_levels <- c("primary_0_6_6_24", "early_0_6_6_12", "late_0_12_12_24")

all_or <- list()
all_metrics <- list()
all_quartiles <- list()
all_scores <- list()

for (db in databases) {
  message("Reading ", db, " time-window cohort...")
  d <- read_database(db)
  for (w in window_levels) {
    message("Fitting ", db, " / ", w, "...")
    res <- fit_window(d, db, w)
    all_or[[paste(db, w, sep = "_")]] <- res$or
    all_metrics[[paste(db, w, sep = "_")]] <- res$metrics
    all_quartiles[[paste(db, w, sep = "_")]] <- res$quartile
    all_scores[[paste(db, w, sep = "_")]] <- res$scores
  }
}

or_rows <- do.call(rbind, all_or)
metrics_rows <- do.call(rbind, all_metrics)
quartile_rows <- do.call(rbind, all_quartiles)
score_rows <- do.call(rbind, all_scores)

write.csv(or_rows, file.path(outdir, "time_window_hrc_or.csv"), row.names = FALSE)
write.csv(metrics_rows, file.path(outdir, "time_window_component_model_metrics.csv"), row.names = FALSE)
write.csv(quartile_rows, file.path(outdir, "time_window_quartile_mortality.csv"), row.names = FALSE)
write.csv(score_rows, file.path(outdir, "time_window_hrc_scores.csv"), row.names = FALSE)

meta_rows <- do.call(rbind, lapply(window_levels, function(w) {
  meta_summary(or_rows[or_rows$window_strategy == w, , drop = FALSE], w)
}))
write.csv(meta_rows, file.path(outdir, "time_window_hrc_meta.csv"), row.names = FALSE)

or_md <- or_rows
or_md$OR_95CI <- fmt_or(or_md$odds_ratio, or_md$ci_low, or_md$ci_high)
or_md$p_value <- signif(or_md$p_value, 3)
or_md <- or_md[, c("database", "window_strategy", "n", "deaths", "OR_95CI", "p_value")]
write_md_table(or_md, file.path(outdir, "time_window_hrc_or.md"), "Time-window HRC mortality association")

meta_md <- meta_rows
meta_md$fixed_OR_95CI <- fmt_or(meta_md$fixed_or, meta_md$fixed_ci_low, meta_md$fixed_ci_high)
meta_md$random_OR_95CI <- fmt_or(meta_md$random_or, meta_md$random_ci_low, meta_md$random_ci_high)
meta_md$i2_percent <- sprintf("%.1f", meta_md$i2_percent)
meta_md <- meta_md[, c("analysis", "k", "fixed_OR_95CI", "random_OR_95CI", "i2_percent", "q_p")]
meta_md$q_p <- signif(meta_md$q_p, 3)
write_md_table(meta_md, file.path(outdir, "time_window_hrc_meta.md"), "Time-window cross-database meta-analysis")

quartile_md <- quartile_rows
quartile_md$mortality_rate <- sprintf("%.3f", quartile_md$mortality_rate)
quartile_md$median_hrc <- sprintf("%.2f", quartile_md$median_hrc)
write_md_table(
  quartile_md[, c("database", "window_strategy", "hrc_quartile", "n", "deaths", "mortality_rate", "median_hrc")],
  file.path(outdir, "time_window_quartile_mortality.md"),
  "Time-window HRC quartile mortality"
)

primary <- or_rows[or_rows$window_strategy == "primary_0_6_6_24", ]
early <- or_rows[or_rows$window_strategy == "early_0_6_6_12", ]
late <- or_rows[or_rows$window_strategy == "late_0_12_12_24", ]

key_lines <- c(
  "# Time-window sensitivity key results",
  "",
  "## Method",
  "",
  "- The formal 24-hour landmark cohorts were kept fixed.",
  "- HRC was re-estimated from raw MAP, vasoactive-burden, and urine-output measurements for each alternative window.",
  "- Expected recovery models were recalibrated within each database and each window using 5-fold cross-fitting.",
  "- Primary inference remains continuous HRC per 1 SD increase; lower OR means better recovery capacity is associated with lower post-landmark hospital mortality.",
  "",
  "## Continuous HRC OR by database",
  "",
  paste(
    apply(or_md, 1, function(row) {
      paste0(
        "- ", row[["database"]], " / ", row[["window_strategy"]],
        ": n=", row[["n"]], ", deaths=", row[["deaths"]],
        ", OR=", row[["OR_95CI"]], ", p=", row[["p_value"]]
      )
    }),
    collapse = "\n"
  ),
  "",
  "## Cross-database meta-analysis",
  "",
  paste(
    apply(meta_md, 1, function(row) {
      paste0(
        "- ", row[["analysis"]],
        ": fixed-effect OR=", row[["fixed_OR_95CI"]],
        "; random-effects OR=", row[["random_OR_95CI"]],
        "; I2=", row[["i2_percent"]], "%"
      )
    }),
    collapse = "\n"
  )
)
writeLines(key_lines, file.path(outdir, "time_window_key_results.md"))

cat("Wrote time-window sensitivity outputs to", outdir, "\n")
