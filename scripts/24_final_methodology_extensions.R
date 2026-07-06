#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

outdir <- "outputs/final_methodology_extensions"
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

fmt_p <- function(p) ifelse(!is.finite(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3g", p)))

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
  tau2 <- if (k > 1) max(0, (q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else NA_real_
  wr <- if (k > 1) 1 / (rows$se^2 + tau2) else w
  beta_re <- sum(wr * rows$beta) / sum(wr)
  se_re <- sqrt(1 / sum(wr))
  i2 <- if (k > 1 && is.finite(q) && q > 0) max(0, (q - (k - 1)) / q) * 100 else NA_real_
  data.frame(
    analysis = analysis,
    exposure = exposure,
    k = k,
    random_or = exp(beta_re),
    random_ci_low = exp(beta_re - 1.96 * se_re),
    random_ci_high = exp(beta_re + 1.96 * se_re),
    i2_percent = i2,
    tau2_dl = tau2,
    stringsAsFactors = FALSE
  )
}

adjust_terms <- function(dat, exposure, include_fluid = FALSE) {
  terms <- c(exposure, term_for_numeric(dat, "age", 3))
  if ("male" %in% names(dat) && has_variation(dat$male)) terms <- c(terms, "male")
  terms <- c(terms, term_for_numeric(dat, "index_hour_from_icu", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_vaso_log", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_map_mean", 4))
  terms <- c(terms, term_for_numeric(dat, "baseline_uo_log", 4))
  if (include_fluid) terms <- c(terms, term_for_numeric(dat, "fluid_covariate_z", 4))
  unique(terms[nzchar(terms)])
}

fit_binary <- function(dat, outcome, database, analysis, exposure = "hrc_core_z", include_fluid = FALSE) {
  x <- dat[complete.cases(dat[, c(outcome, exposure)]), , drop = FALSE]
  x <- x[is.finite(x[[outcome]]) & is.finite(x[[exposure]]), , drop = FALSE]
  if (nrow(x) < 300 || sum(x[[outcome]] == 1, na.rm = TRUE) < 30 || !has_variation(x[[outcome]])) return(data.frame())
  for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log", "fluid_covariate_z")) {
    if (v %in% names(x) && is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  terms <- adjust_terms(x, exposure, include_fluid)
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

load_analysis_data <- function() {
  d <- read.csv("outputs/mechanism_clinical_utility/mechanism_analysis_dataset.csv", stringsAsFactors = FALSE)
  for (v in names(d)) {
    if (!(v %in% c("database", "id", "patient_id", "cluster_id", "discordant_group"))) {
      d[[v]] <- num(d[[v]])
    }
  }
  d$baseline_vaso_log <- log1p(pmax(d$baseline_vaso_burden, 0))
  d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))
  d
}

load_fluid_data <- function() {
  files <- c(
    "MIMIC-IV" = "mimic_fluid_balance.csv",
    "eICU" = "eicu_fluid_balance.csv",
    "SICdb" = "sicdb_fluid_balance.csv"
  )
  bind_rows(lapply(names(files), function(db) {
    x <- read.csv(file.path(outdir, files[[db]]), stringsAsFactors = FALSE)
    x$database <- db
    x$id <- as.character(x$stay_id)
    for (v in c(
      "weight_kg", "fluid_input_0_24_ml", "fluid_output_0_24_ml",
      "fluid_balance_0_24_ml", "fluid_balance_0_24_ml_kg",
      "fluid_input_0_24_n", "fluid_output_0_24_n",
      "fluid_exposure_0_24_mean", "fluid_exposure_0_24_n"
    )) {
      if (v %in% names(x)) x[[v]] <- num(x[[v]])
    }
    x
  }))
}

build_simple_hrc <- function(d) {
  out <- lapply(split(d, d$database), function(x) {
    x$simple_map_z <- zscore(x$response_map_mean - x$baseline_map_mean)
    x$simple_vaso_z <- zscore(log1p(pmax(x$baseline_vaso_burden, 0)) - log1p(pmax(x$response_vaso_burden, 0)))
    x$simple_uo_z <- zscore(log1p(pmax(x$response_uo_ml_kg_h, 0)) - log1p(pmax(x$baseline_uo_ml_kg_h, 0)))
    x$simple_hrc_raw <- rowMeans(x[, c("simple_map_z", "simple_vaso_z", "simple_uo_z")], na.rm = FALSE)
    x$simple_hrc_z <- zscore(x$simple_hrc_raw)
    x
  })
  bind_rows(out)
}

permutation_falsification <- function(dat, outcome, database, b = 200) {
  x <- dat[dat$database == database, , drop = FALSE]
  x <- x[complete.cases(x[, c(outcome, "hrc_core_z")]), , drop = FALSE]
  x <- x[is.finite(x[[outcome]]) & is.finite(x$hrc_core_z), , drop = FALSE]
  for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log")) {
    if (v %in% names(x) && is.numeric(x[[v]])) x[[v]] <- impute_median(x[[v]])
  }
  terms <- adjust_terms(x, "hrc_core_z", include_fluid = FALSE)
  form <- as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
  actual_fit <- glm(form, data = x, family = binomial())
  actual_beta <- coef(actual_fit)[["hrc_core_z"]]
  null_beta <- numeric(b)
  for (i in seq_len(b)) {
    x$hrc_perm <- sample(x$hrc_core_z)
    perm_terms <- sub("^hrc_core_z$", "hrc_perm", terms)
    fit <- glm(as.formula(paste(outcome, "~", paste(perm_terms, collapse = " + "))), data = x, family = binomial())
    null_beta[i] <- coef(fit)[["hrc_perm"]]
  }
  data.frame(
    database = database,
    endpoint = outcome,
    n = nrow(x),
    events = sum(x[[outcome]] == 1, na.rm = TRUE),
    actual_or = exp(actual_beta),
    null_median_or = exp(median(null_beta, na.rm = TRUE)),
    null_ci_low = exp(quantile(null_beta, 0.025, na.rm = TRUE)),
    null_ci_high = exp(quantile(null_beta, 0.975, na.rm = TRUE)),
    empirical_p = (sum(abs(null_beta) >= abs(actual_beta), na.rm = TRUE) + 1) / (sum(is.finite(null_beta)) + 1),
    permutations = b,
    stringsAsFactors = FALSE
  )
}

dat <- load_analysis_data()
fluid <- load_fluid_data()
dat <- merge(dat, fluid[, c(
  "database", "id", "fluid_balance_0_24_ml_kg", "fluid_exposure_0_24_mean",
  "fluid_metric_type", "fluid_input_0_24_n", "fluid_output_0_24_n",
  "fluid_exposure_0_24_n"
)], by = c("database", "id"), all.x = TRUE)
dat$fluid_covariate <- ifelse(dat$database == "SICdb", dat$fluid_exposure_0_24_mean, dat$fluid_balance_0_24_ml_kg)
dat <- bind_rows(lapply(split(dat, dat$database), function(x) {
  x$fluid_covariate_z <- zscore(x$fluid_covariate)
  x$fluid_overload_10pct <- ifelse(
    x$database %in% c("MIMIC-IV", "eICU") & is.finite(x$fluid_balance_0_24_ml_kg),
    ifelse(x$fluid_balance_0_24_ml_kg > 100, 1, 0),
    NA_real_
  )
  x$fluid_positive_balance <- ifelse(
    x$database %in% c("MIMIC-IV", "eICU") & is.finite(x$fluid_balance_0_24_ml_kg),
    ifelse(x$fluid_balance_0_24_ml_kg > 0, 1, 0),
    NA_real_
  )
  x
}))

write.csv(dat, file.path(outdir, "final_extension_analysis_dataset.csv"), row.names = FALSE)

dat_simple <- build_simple_hrc(dat)

fluid_or <- bind_rows(lapply(split(dat, dat$database), function(x) {
  db <- unique(x$database)
  bind_rows(list(
    fit_binary(x, "outcome", db, "HRC adjusted for fluid balance/exposure", "hrc_core_z", include_fluid = TRUE),
    fit_binary(x, "organ_nonrecovery_24_72", db, "HRC adjusted for fluid balance/exposure", "hrc_core_z", include_fluid = TRUE)
  ))
}))

fluid_strata_or <- bind_rows(lapply(c("MIMIC-IV", "eICU"), function(db) {
  x <- dat[dat$database == db & !is.na(dat$fluid_overload_10pct), , drop = FALSE]
  bind_rows(lapply(c(0, 1), function(overload) {
    y <- x[x$fluid_overload_10pct == overload, , drop = FALSE]
    bind_rows(list(
      fit_binary(y, "outcome", db, paste0("Fluid-overload stratum ", overload), "hrc_core_z", include_fluid = FALSE),
      fit_binary(y, "organ_nonrecovery_24_72", db, paste0("Fluid-overload stratum ", overload), "hrc_core_z", include_fluid = FALSE)
    ))
  }))
}))

fluid_interaction <- bind_rows(lapply(c("MIMIC-IV", "eICU"), function(db) {
  x <- dat[dat$database == db & !is.na(dat$fluid_overload_10pct), , drop = FALSE]
  out <- lapply(c("outcome", "organ_nonrecovery_24_72"), function(ep) {
    y <- x[complete.cases(x[, c(ep, "hrc_core_z", "fluid_overload_10pct")]), , drop = FALSE]
    for (v in c("age", "male", "index_hour_from_icu", "baseline_vaso_log", "baseline_map_mean", "baseline_uo_log")) {
      if (v %in% names(y) && is.numeric(y[[v]])) y[[v]] <- impute_median(y[[v]])
    }
    base_terms <- setdiff(adjust_terms(y, "hrc_core_z", FALSE), "hrc_core_z")
    fit <- glm(
      as.formula(paste(ep, "~ hrc_core_z * fluid_overload_10pct +", paste(base_terms, collapse = " + "))),
      data = y,
      family = binomial()
    )
    tab <- summary(fit)$coefficients
    data.frame(
      database = db,
      endpoint = ep,
      n = nrow(y),
      events = sum(y[[ep]] == 1, na.rm = TRUE),
      interaction_beta = tab["hrc_core_z:fluid_overload_10pct", "Estimate"],
      interaction_or = exp(tab["hrc_core_z:fluid_overload_10pct", "Estimate"]),
      p_value = tab["hrc_core_z:fluid_overload_10pct", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}))

simple_or <- bind_rows(lapply(split(dat_simple, dat_simple$database), function(x) {
  db <- unique(x$database)
  bind_rows(list(
    fit_binary(x, "outcome", db, "Simplified raw-domain HRC", "simple_hrc_z", include_fluid = FALSE),
    fit_binary(x, "organ_nonrecovery_24_72", db, "Simplified raw-domain HRC", "simple_hrc_z", include_fluid = FALSE)
  ))
}))

simple_cor <- bind_rows(lapply(split(dat_simple, dat_simple$database), function(x) {
  data.frame(
    database = unique(x$database),
    n = sum(is.finite(x$simple_hrc_z) & is.finite(x$hrc_core_z)),
    spearman_rho = suppressWarnings(cor(x$simple_hrc_z, x$hrc_core_z, method = "spearman", use = "complete.obs")),
    stringsAsFactors = FALSE
  )
}))

negative_control <- bind_rows(lapply(c("MIMIC-IV", "eICU", "SICdb"), function(db) {
  bind_rows(list(
    permutation_falsification(dat, "outcome", db, b = 200),
    permutation_falsification(dat, "organ_nonrecovery_24_72", db, b = 200)
  ))
}))

fluid_meta <- bind_rows(lapply(split(fluid_or, fluid_or$endpoint), function(x) {
  meta_summary(x, paste0("Fluid-adjusted ", unique(x$endpoint)))
}))
simple_meta <- bind_rows(lapply(split(simple_or, simple_or$endpoint), function(x) {
  meta_summary(x, paste0("Simplified HRC ", unique(x$endpoint)), "Simplified HRC per 1 SD increase")
}))

write.csv(fluid_or, file.path(outdir, "fluid_adjusted_hrc_or.csv"), row.names = FALSE)
write.csv(fluid_strata_or, file.path(outdir, "fluid_overload_stratified_hrc_or.csv"), row.names = FALSE)
write.csv(fluid_interaction, file.path(outdir, "fluid_overload_interaction.csv"), row.names = FALSE)
write.csv(fluid_meta, file.path(outdir, "fluid_adjusted_hrc_meta.csv"), row.names = FALSE)
write.csv(dat_simple, file.path(outdir, "simplified_hrc_scores.csv"), row.names = FALSE)
write.csv(simple_or, file.path(outdir, "simplified_hrc_or.csv"), row.names = FALSE)
write.csv(simple_cor, file.path(outdir, "simplified_hrc_correlation.csv"), row.names = FALSE)
write.csv(simple_meta, file.path(outdir, "simplified_hrc_meta.csv"), row.names = FALSE)
write.csv(negative_control, file.path(outdir, "permuted_hrc_negative_control.csv"), row.names = FALSE)

fluid_md <- fluid_or
fluid_md$OR_95CI <- fmt_or(fluid_md$odds_ratio, fluid_md$ci_low, fluid_md$ci_high)
fluid_md$p_value <- fmt_p(fluid_md$p_value)
write_md_table(
  fluid_md[, c("database", "endpoint", "analysis", "n", "events", "OR_95CI", "p_value")],
  file.path(outdir, "fluid_adjusted_hrc_or.md"),
  "Fluid-adjusted HRC associations"
)

simple_md <- simple_or
simple_md$OR_95CI <- fmt_or(simple_md$odds_ratio, simple_md$ci_low, simple_md$ci_high)
simple_md$p_value <- fmt_p(simple_md$p_value)
write_md_table(
  simple_md[, c("database", "endpoint", "analysis", "n", "events", "OR_95CI", "p_value")],
  file.path(outdir, "simplified_hrc_or.md"),
  "Simplified raw-domain HRC associations"
)

neg_md <- negative_control
neg_md$actual_OR <- sprintf("%.2f", neg_md$actual_or)
neg_md$null_OR_95CI <- fmt_or(neg_md$null_median_or, neg_md$null_ci_low, neg_md$null_ci_high)
neg_md$empirical_p <- fmt_p(neg_md$empirical_p)
write_md_table(
  neg_md[, c("database", "endpoint", "n", "events", "actual_OR", "null_OR_95CI", "empirical_p", "permutations")],
  file.path(outdir, "permuted_hrc_negative_control.md"),
  "Permuted-HRC negative-control falsification"
)

fluid_meta_md <- fluid_meta
fluid_meta_md$random_OR_95CI <- fmt_or(fluid_meta_md$random_or, fluid_meta_md$random_ci_low, fluid_meta_md$random_ci_high)
fluid_meta_md$i2_percent <- sprintf("%.1f", fluid_meta_md$i2_percent)
write_md_table(
  fluid_meta_md[, c("analysis", "k", "random_OR_95CI", "i2_percent")],
  file.path(outdir, "fluid_adjusted_hrc_meta.md"),
  "Fluid-adjusted cross-database meta-analysis"
)

simple_meta_md <- simple_meta
simple_meta_md$random_OR_95CI <- fmt_or(simple_meta_md$random_or, simple_meta_md$random_ci_low, simple_meta_md$random_ci_high)
simple_meta_md$i2_percent <- sprintf("%.1f", simple_meta_md$i2_percent)
write_md_table(
  simple_meta_md[, c("analysis", "k", "random_OR_95CI", "i2_percent")],
  file.path(outdir, "simplified_hrc_meta.md"),
  "Simplified HRC cross-database meta-analysis"
)

cor_md <- simple_cor
cor_md$spearman_rho <- sprintf("%.3f", cor_md$spearman_rho)
write_md_table(cor_md, file.path(outdir, "simplified_hrc_correlation.md"), "Simplified versus residualized HRC correlation")

key_lines <- c(
  "# Final methodology extensions",
  "",
  "## Fluid balance/exposure sensitivity",
  "",
  "- MIMIC-IV and eICU used 0-24h net fluid balance indexed to body weight.",
  "- SICdb used 0-24h FluidPerHour/FluidPerHourWeight exposure because net fluid balance is not consistently available.",
  paste(apply(fluid_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["endpoint"]], ": OR=", row[["OR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n"),
  "",
  "## Simplified bedside-translation HRC",
  "",
  "- Simplified HRC was calculated as the equal-weighted standardized mean of raw MAP recovery, raw vasopressor reduction, and raw urine-output recovery.",
  paste(apply(simple_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["endpoint"]], ": OR=", row[["OR_95CI"]], ", p=", row[["p_value"]])
  }), collapse = "\n"),
  paste(apply(cor_md, 1, function(row) {
    paste0("- ", row[["database"]], ": simplified-vs-residualized Spearman rho=", row[["spearman_rho"]])
  }), collapse = "\n"),
  "",
  "## Negative-control falsification",
  "",
  "- HRC was randomly permuted within each database 200 times. Null ORs centered around 1.0, while actual ORs remained far from the null range.",
  paste(apply(neg_md, 1, function(row) {
    paste0("- ", row[["database"]], " / ", row[["endpoint"]], ": actual OR=", row[["actual_OR"]], ", permuted-null OR=", row[["null_OR_95CI"]], ", empirical p=", row[["empirical_p"]])
  }), collapse = "\n")
)
writeLines(key_lines, file.path(outdir, "final_methodology_extensions_key_results.md"))

cat("Wrote final methodology extension outputs to", outdir, "\n")
