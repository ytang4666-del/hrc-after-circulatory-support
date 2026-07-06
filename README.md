# Hemodynamic Recovery Capacity After Early Circulatory Support

This repository contains the analysis code and aggregate outputs for the manuscript:

**Hemodynamic Recovery Capacity After Early Circulatory Support Across Three Critical Care Databases**

Target journal: Annals of Intensive Care

## Study Overview

The study derives and validates Hemodynamic Recovery Capacity (HRC), a model-based residualized measure of early physiologic recovery after initiation of early vasoactive circulatory support. Analyses were performed across:

- MIMIC-IV v3.1
- eICU Collaborative Research Database v2.0
- SICdb v1.0.8

The primary outcome was post-landmark hospital mortality. The key secondary outcome was post-landmark organ nonrecovery.

## Important Data Restriction

This repository does **not** contain raw or patient-level derived MIMIC-IV, eICU, or SICdb data. Access to those databases is governed by the original database licenses and credentialing requirements:

- MIMIC-IV v3.1: https://physionet.org/content/mimiciv/3.1/
- eICU-CRD v2.0: https://physionet.org/content/eicu-crd/2.0/
- SICdb v1.0.8: https://physionet.org/content/sicdb/1.0.8/

Users must obtain independent database access and run the extraction scripts in their own database environment.

## Repository Structure

```text
scripts/
  Cohort extraction, HRC construction, modeling, sensitivity analyses, and figure/table scripts.

metadata/
  Variable harmonization, expected-recovery model specification, and STROBE checklist.

outputs_aggregate/
  Aggregate figure source data and supplementary table source data.

manuscript_files/
  Submission-ready manuscript draft, title page, cover letter, and compliance checklist.
```

## Database Connections

The public release uses environment-variable connection templates:

- `scripts/pgadmin_conn.py` for MIMIC-IV/eICU PostgreSQL connections.
- `scripts/sicdb_conn.py` for SICdb PostgreSQL connections.

Set environment variables before running extraction scripts. Do not store credentials in this repository.

Example:

```bash
export MIMIC_DB_HOST="localhost"
export MIMIC_DB_PORT="5432"
export MIMIC_DB_NAME="mimiciv"
export MIMIC_DB_USER="your_user"
export MIMIC_DB_PASSWORD="your_password"

export EICU_DB_HOST="localhost"
export EICU_DB_PORT="5432"
export EICU_DB_NAME="eicu"
export EICU_DB_USER="your_user"
export EICU_DB_PASSWORD="your_password"

export SICDB_HOST="localhost"
export SICDB_PORT="55432"
export SICDB_DBNAME="sicdb"
export SICDB_USER="your_user"
export SICDB_PASSWORD="your_password"
```

## Reproducibility Notes

1. Run database extraction scripts only after obtaining database access.
2. Do not commit raw database exports or patient-level derived analytic files.
3. Use aggregate output tables and figure source data for manuscript figure verification.
4. Regenerate software version files if analyses are rerun in a different environment.

Some supplementary generation scripts require a locally generated patient-level intermediate analysis file. In the public release this file is not included. If needed, generate it locally and point scripts to it with:

```bash
export HRC_PATIENT_LEVEL_ANALYSIS_DATA="outputs_patient_level/analysis_dataset_core_hrc.csv"
```

The `outputs_patient_level/` folder is ignored by `.gitignore` and should not be committed.

## Software

The manuscript workspace used:

- R 4.5.1
- Python 3.9.6

Key R packages included `ggplot2`, `patchwork`, `svglite`, `ragg`, `stringr`, `dplyr`, `readr`, and `tidyr`.

## Citation

If using this code, cite the associated manuscript and the archived Zenodo release DOI once available.
