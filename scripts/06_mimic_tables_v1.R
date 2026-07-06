#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(splines))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "outputs/mimic_formal"
table_dir <- file.path(outdir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cohort <- read.csv(file.path(outdir, "mimic_formal_cohort.csv"), stringsAsFactors = FALSE)
scores <- read.csv(file.path(outdir, "mimic_hrc_formal_scores.csv"), stringsAsFactors = FALSE)

num_cols <- c(
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
for (col in intersect(num_cols, names(cohort))) cohort[[col]] <- suppressWarnings(as.numeric(cohort[[col]]))
for (col in intersect(num_cols, names(scores))) scores[[col]] <- suppressWarnings(as.numeric(scores[[col]]))

score_keep <- c("stay_id", "hrc_core_z", "hrc_quartile", "low_hrc_q1",
                "hrc_map_residual_z", "hrc_vaso_residual_z", "hrc_uo_residual_z")
d <- merge(cohort, scores[, score_keep], by = "stay_id", all.x = FALSE)
d$gender_male <- ifelse(d$gender == "M", 1, 0)
d$baseline_neq_log <- log1p(pmax(d$baseline_neq_avg, 0))
d$baseline_uo_log <- log1p(pmax(d$baseline_uo_ml_kg_h, 0))
d$hrc_group <- ifelse(d$low_hrc_q1 == 1, "Low HRC", "Non-low HRC")
d$hrc_group <- factor(d$hrc_group, levels = c("Non-low HRC", "Low HRC"))

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "", format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE))
}

fmt_n_pct <- function(x, denom = length(x)) {
  n <- sum(x == 1, na.rm = TRUE)
  pct <- 100 * n / denom
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), pct)
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f-%.1f)", q[[2]], q[[1]], q[[3]])
}

smd_cont <- function(x, g) {
  x0 <- x[g == "Non-low HRC"]
  x1 <- x[g == "Low HRC"]
  den <- sqrt((var(x0, na.rm = TRUE) + var(x1, na.rm = TRUE)) / 2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / den
}

smd_bin <- function(x, g) {
  p0 <- mean(x[g == "Non-low HRC"] == 1, na.rm = TRUE)
  p1 <- mean(x[g == "Low HRC"] == 1, na.rm = TRUE)
  den <- sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)
  if (!is.finite(den) || den == 0) return(NA_real_)
  (p1 - p0) / den
}

write_md_table <- function(df, path, title = NULL) {
  lines <- character()
  if (!is.null(title)) {
    lines <- c(lines, paste0("# ", title), "")
  }
  cols <- names(df)
  lines <- c(lines, paste(cols, collapse = " | "))
  lines <- c(lines, paste(rep("---", length(cols)), collapse = " | "))
  for (i in seq_len(nrow(df))) {
    vals <- vapply(df[i, , drop = FALSE], as.character, character(1))
    vals <- gsub("\\|", "/", vals)
    lines <- c(lines, paste(vals, collapse = " | "))
  }
  writeLines(lines, con = path)
}

make_table1 <- function(data) {
  groups <- list(
    Overall = data,
    `Non-low HRC` = data[data$hrc_group == "Non-low HRC", ],
    `Low HRC` = data[data$hrc_group == "Low HRC", ]
  )
  variables <- list(
    list(label = "Age, years", var = "anchor_age", type = "continuous"),
    list(label = "Male sex", var = "gender_male", type = "binary"),
    list(label = "Weight, kg", var = "weight_kg", type = "continuous"),
    list(label = "First-day SOFA", var = "first_day_sofa", type = "continuous"),
    list(label = "SOFA cardiovascular", var = "sofa_cardiovascular", type = "continuous"),
    list(label = "SOFA renal", var = "sofa_renal", type = "continuous"),
    list(label = "Baseline MAP, mmHg", var = "baseline_map_mean", type = "continuous"),
    list(label = "Response MAP, mmHg", var = "response_map_mean", type = "continuous"),
    list(label = "Baseline NE-equivalent", var = "baseline_neq_avg", type = "continuous"),
    list(label = "Response NE-equivalent", var = "response_neq_avg", type = "continuous"),
    list(label = "Baseline urine output, ml/kg/h", var = "baseline_uo_ml_kg_h", type = "continuous"),
    list(label = "Response urine output, ml/kg/h", var = "response_uo_ml_kg_h", type = "continuous"),
    list(label = "Baseline creatinine, mg/dl", var = "baseline_creatinine_mean", type = "continuous"),
    list(label = "Response creatinine, mg/dl", var = "response_creatinine_mean", type = "continuous"),
    list(label = "Baseline lactate, mmol/l", var = "baseline_lactate_mean", type = "continuous"),
    list(label = "Response lactate, mmol/l", var = "response_lactate_mean", type = "continuous"),
    list(label = "Cardiogenic shock ICD", var = "cardiogenic_shock_icd", type = "binary"),
    list(label = "Any ICU MCS", var = "mcs_any_icu", type = "binary"),
    list(label = "IABP", var = "iabp_any_icu", type = "binary"),
    list(label = "Impella", var = "impella_any_icu", type = "binary"),
    list(label = "ECMO", var = "ecmo_any_icu", type = "binary"),
    list(label = "Post-landmark hospital death", var = "hospital_mortality_after_landmark", type = "binary")
  )
  rows <- lapply(variables, function(spec) {
    vals <- lapply(groups, function(gd) {
      if (spec$type == "continuous") fmt_median_iqr(gd[[spec$var]]) else fmt_n_pct(gd[[spec$var]], nrow(gd))
    })
    smd <- if (spec$type == "continuous") {
      smd_cont(data[[spec$var]], data$hrc_group)
    } else {
      smd_bin(data[[spec$var]], data$hrc_group)
    }
    data.frame(
      Characteristic = spec$label,
      Overall = vals$Overall,
      `Non-low HRC` = vals$`Non-low HRC`,
      `Low HRC` = vals$`Low HRC`,
      SMD = fmt_num(smd, 2),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

table1 <- make_table1(d)
write.csv(table1, file.path(table_dir, "table1_baseline_low_hrc.csv"), row.names = FALSE)
write_md_table(table1, file.path(table_dir, "table1_baseline_low_hrc.md"), "Table 1. Baseline characteristics by HRC phenotype")

extract_or <- function(fit, term) {
  s <- summary(fit)$coefficients
  ci <- confint.default(fit)
  if (!(term %in% rownames(s))) {
    return(data.frame(beta = NA, se = NA, odds_ratio = NA, ci_low = NA, ci_high = NA, p_value = NA))
  }
  data.frame(
    beta = s[term, "Estimate"],
    se = s[term, "Std. Error"],
    odds_ratio = exp(s[term, "Estimate"]),
    ci_low = exp(ci[term, 1]),
    ci_high = exp(ci[term, 2]),
    p_value = s[term, "Pr(>|z|)"],
    row.names = NULL
  )
}

fmt_or <- function(row) {
  sprintf("%.2f (%.2f-%.2f)", row$odds_ratio, row$ci_low, row$ci_high)
}

model_data <- d[complete.cases(d[, c(
  "hospital_mortality_after_landmark", "hrc_core_z", "low_hrc_q1",
  "anchor_age", "gender_male", "first_day_sofa",
  "baseline_neq_log", "baseline_map_mean", "baseline_uo_log",
  "cardiogenic_shock_icd", "mcs_any_icu"
)]), ]

model_specs <- list(
  Crude = "hospital_mortality_after_landmark ~ EXPOSURE",
  Demographic = "hospital_mortality_after_landmark ~ EXPOSURE + ns(anchor_age, df = 3) + gender_male",
  Severity = "hospital_mortality_after_landmark ~ EXPOSURE + ns(anchor_age, df = 3) + gender_male + ns(first_day_sofa, df = 4)",
  `Baseline physiology` = "hospital_mortality_after_landmark ~ EXPOSURE + ns(anchor_age, df = 3) + gender_male + ns(first_day_sofa, df = 4) + ns(baseline_neq_log, df = 4) + ns(baseline_map_mean, df = 4) + ns(baseline_uo_log, df = 4)",
  `Expanded CS/MCS sensitivity` = "hospital_mortality_after_landmark ~ EXPOSURE + ns(anchor_age, df = 3) + gender_male + ns(first_day_sofa, df = 4) + ns(baseline_neq_log, df = 4) + ns(baseline_map_mean, df = 4) + ns(baseline_uo_log, df = 4) + cardiogenic_shock_icd + mcs_any_icu"
)

table2_rows <- list()
for (nm in names(model_specs)) {
  f1 <- as.formula(gsub("EXPOSURE", "hrc_core_z", model_specs[[nm]], fixed = TRUE))
  fit1 <- glm(f1, data = model_data, family = binomial())
  or1 <- extract_or(fit1, "hrc_core_z")
  table2_rows[[length(table2_rows) + 1]] <- data.frame(
    Exposure = "HRC per 1 SD increase",
    Model = nm,
    N = nrow(model_data),
    Events = sum(model_data$hospital_mortality_after_landmark == 1),
    `OR (95% CI)` = fmt_or(or1),
    `P value` = signif(or1$p_value, 3),
    check.names = FALSE
  )

  f2 <- as.formula(gsub("EXPOSURE", "low_hrc_q1", model_specs[[nm]], fixed = TRUE))
  fit2 <- glm(f2, data = model_data, family = binomial())
  or2 <- extract_or(fit2, "low_hrc_q1")
  table2_rows[[length(table2_rows) + 1]] <- data.frame(
    Exposure = "Low-HRC phenotype",
    Model = nm,
    N = nrow(model_data),
    Events = sum(model_data$hospital_mortality_after_landmark == 1),
    `OR (95% CI)` = fmt_or(or2),
    `P value` = signif(or2$p_value, 3),
    check.names = FALSE
  )
}
table2 <- do.call(rbind, table2_rows)
write.csv(table2, file.path(table_dir, "table2_primary_models.csv"), row.names = FALSE)
write_md_table(table2, file.path(table_dir, "table2_primary_models.md"), "Table 2. Association of HRC with post-landmark hospital mortality")

or_2x2 <- function(a, b, c, d0) {
  if (any(c(a, b, c, d0) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d0 <- d0 + 0.5
  }
  or <- (a / b) / (c / d0)
  se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d0)
  c(or = or, low = exp(log(or) - 1.96 * se), high = exp(log(or) + 1.96 * se))
}

subgroups <- list(
  Overall = d,
  `No CS/MCS` = d[d$cs_or_mcs_subgroup == 0, ],
  `Cardiogenic shock ICD` = d[d$cardiogenic_shock_icd == 1, ],
  `Any MCS during ICU` = d[d$mcs_any_icu == 1, ],
  `MCS during 0-24h` = d[d$mcs_any_0_24h == 1, ],
  IABP = d[d$iabp_any_icu == 1, ],
  Impella = d[d$impella_any_icu == 1, ],
  ECMO = d[d$ecmo_any_icu == 1, ],
  `CS or MCS` = d[d$cs_or_mcs_subgroup == 1, ]
)

table3 <- do.call(rbind, lapply(names(subgroups), function(nm) {
  x <- subgroups[[nm]]
  low <- x[x$low_hrc_q1 == 1, ]
  non <- x[x$low_hrc_q1 == 0, ]
  a <- sum(low$hospital_mortality_after_landmark == 1, na.rm = TRUE)
  b <- sum(low$hospital_mortality_after_landmark == 0, na.rm = TRUE)
  c0 <- sum(non$hospital_mortality_after_landmark == 1, na.rm = TRUE)
  d0 <- sum(non$hospital_mortality_after_landmark == 0, na.rm = TRUE)
  orci <- or_2x2(a, b, c0, d0)
  data.frame(
    Subgroup = nm,
    N = nrow(x),
    Events = sum(x$hospital_mortality_after_landmark == 1, na.rm = TRUE),
    `Overall mortality` = fmt_n_pct(x$hospital_mortality_after_landmark, nrow(x)),
    `Low-HRC n` = nrow(low),
    `Low-HRC mortality` = fmt_n_pct(low$hospital_mortality_after_landmark, nrow(low)),
    `Non-low-HRC mortality` = fmt_n_pct(non$hospital_mortality_after_landmark, nrow(non)),
    `Low vs non-low OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", orci[["or"]], orci[["low"]], orci[["high"]]),
    check.names = FALSE
  )
}))
write.csv(table3, file.path(table_dir, "table3_cs_mcs_subgroups.csv"), row.names = FALSE)
write_md_table(table3, file.path(table_dir, "table3_cs_mcs_subgroups.md"), "Table 3. Low-HRC phenotype in cardiogenic shock and MCS subgroups")

complete_rows <- data.frame(
  Variable = c(
    "MAP baseline and response",
    "Vasopressor burden baseline and response",
    "Urine output + weight baseline and response",
    "Creatinine baseline and response",
    "Lactate baseline and response",
    "Core HRC complete"
  ),
  N = c(
    sum(!is.na(cohort$baseline_map_mean) & !is.na(cohort$response_map_mean)),
    sum(!is.na(cohort$baseline_neq_avg) & !is.na(cohort$response_neq_avg)),
    sum(!is.na(cohort$baseline_uo_ml_kg_h) & !is.na(cohort$response_uo_ml_kg_h)),
    sum(!is.na(cohort$baseline_creatinine_mean) & !is.na(cohort$response_creatinine_mean)),
    sum(!is.na(cohort$baseline_lactate_mean) & !is.na(cohort$response_lactate_mean)),
    nrow(d)
  ),
  Denominator = nrow(cohort),
  Percent = c(
    mean(!is.na(cohort$baseline_map_mean) & !is.na(cohort$response_map_mean)),
    mean(!is.na(cohort$baseline_neq_avg) & !is.na(cohort$response_neq_avg)),
    mean(!is.na(cohort$baseline_uo_ml_kg_h) & !is.na(cohort$response_uo_ml_kg_h)),
    mean(!is.na(cohort$baseline_creatinine_mean) & !is.na(cohort$response_creatinine_mean)),
    mean(!is.na(cohort$baseline_lactate_mean) & !is.na(cohort$response_lactate_mean)),
    nrow(d) / nrow(cohort)
  )
)
complete_rows$Percent <- sprintf("%.1f%%", 100 * complete_rows$Percent)
write.csv(complete_rows, file.path(table_dir, "supp_table_variable_completeness.csv"), row.names = FALSE)
write_md_table(complete_rows, file.path(table_dir, "supp_table_variable_completeness.md"), "Supplementary Table. Variable completeness")

mcs_audit <- read.csv(file.path(outdir, "mimic_mcs_procedureevents_audit.csv"), stringsAsFactors = FALSE)
cs_audit <- read.csv(file.path(outdir, "mimic_cardiogenic_shock_icd_audit.csv"), stringsAsFactors = FALSE)
write_md_table(mcs_audit, file.path(table_dir, "supp_table_mcs_itemid_audit.md"), "Supplementary Table. MCS procedureevents definitions")
write_md_table(cs_audit, file.path(table_dir, "supp_table_cardiogenic_shock_icd_audit.md"), "Supplementary Table. Cardiogenic shock ICD definitions")
write.csv(mcs_audit, file.path(table_dir, "supp_table_mcs_itemid_audit.csv"), row.names = FALSE)
write.csv(cs_audit, file.path(table_dir, "supp_table_cardiogenic_shock_icd_audit.csv"), row.names = FALSE)

index <- data.frame(
  Table = c(
    "Table 1",
    "Table 2",
    "Table 3",
    "Supplementary variable completeness",
    "Supplementary MCS itemid audit",
    "Supplementary cardiogenic shock ICD audit"
  ),
  CSV = c(
    "table1_baseline_low_hrc.csv",
    "table2_primary_models.csv",
    "table3_cs_mcs_subgroups.csv",
    "supp_table_variable_completeness.csv",
    "supp_table_mcs_itemid_audit.csv",
    "supp_table_cardiogenic_shock_icd_audit.csv"
  ),
  Markdown = c(
    "table1_baseline_low_hrc.md",
    "table2_primary_models.md",
    "table3_cs_mcs_subgroups.md",
    "supp_table_variable_completeness.md",
    "supp_table_mcs_itemid_audit.md",
    "supp_table_cardiogenic_shock_icd_audit.md"
  )
)
write.csv(index, file.path(table_dir, "table_index.csv"), row.names = FALSE)
write_md_table(index, file.path(table_dir, "table_index.md"), "MIMIC table outputs")

cat("Tables written to:", table_dir, "\n")
cat("Table 1 rows:", nrow(table1), "\n")
cat("Table 2 rows:", nrow(table2), "\n")
cat("Table 3 rows:", nrow(table3), "\n")
