from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from pgadmin_conn import connect_pgadmin_server
from sicdb_conn import connect_sicdb


DEFAULT_OUTDIR = Path("outputs/high_priority_methodology")
NUMERIC_RE_SQL = "^-?[0-9]+(\\.[0-9]+)?$"


def load_script_module(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_DIR / file_name)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mimic_formal = load_script_module("mimic_formal_extract_surv", "03_mimic_formal_cohort_extract.py")
eicu_formal = load_script_module("eicu_formal_extract_surv", "07_eicu_formal_cohort_extract.py")
sicdb_formal = load_script_module("sicdb_formal_extract_surv", "10_sicdb_formal_cohort_extract.py")


def log(message: str) -> None:
    print(message, flush=True)


def copy_query_to_csv(cur, query: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with cur.copy(f"COPY ({query}) TO STDOUT WITH CSV HEADER") as copy:
        with out_path.open("wb") as f:
            for data in copy:
                f.write(data)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fetch_scalar(cur, sql: str) -> Any:
    row = cur.execute(sql).fetchone()
    if isinstance(row, dict):
        return next(iter(row.values()))
    return row[0]


def extract_mimic(outdir: Path, server: str) -> None:
    log("MIMIC-IV: extracting post-landmark hospital death/discharge times...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '60min'")
            cur.execute("SET work_mem = '512MB'")
            mimic_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table mimic_survival_source as
                select
                  c.subject_id::text as patient_id,
                  c.hadm_id::text as encounter_id,
                  c.stay_id::text as stay_id,
                  c.index_time,
                  c.landmark_time,
                  c.outtime as icu_outtime,
                  a.dischtime as hospital_discharge_time,
                  a.deathtime as death_time,
                  a.hospital_expire_flag,
                  case
                    when a.hospital_expire_flag = 1 and a.deathtime >= c.landmark_time then 1
                    else 0
                  end as death_event,
                  case
                    when a.hospital_expire_flag = 0 and a.dischtime >= c.landmark_time then 1
                    else 0
                  end as alive_discharge_event,
                  extract(epoch from (
                    case
                      when a.hospital_expire_flag = 1 and a.deathtime >= c.landmark_time then a.deathtime
                      when a.dischtime >= c.landmark_time then a.dischtime
                      else null
                    end - c.landmark_time
                  )) / 86400.0 as time_to_hospital_exit_days,
                  extract(epoch from (c.outtime - c.landmark_time)) / 86400.0 as time_to_icu_exit_days
                from formal_index c
                join mimiciv_hosp.admissions a
                  on c.subject_id = a.subject_id
                 and c.hadm_id = a.hadm_id
                """
            )
            copy_query_to_csv(
                cur,
                "select * from mimic_survival_source order by patient_id, index_time, stay_id",
                outdir / "mimic_survival_source.csv",
            )
            log(f"MIMIC-IV survival rows: {fetch_scalar(cur, 'select count(*) from mimic_survival_source')}")


def extract_eicu(outdir: Path, server: str) -> None:
    log("eICU: extracting post-landmark hospital death/discharge offsets...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '60min'")
            cur.execute("SET work_mem = '512MB'")
            eicu_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table eicu_survival_source as
                select
                  coalesce(nullif(c.uniquepid, ''), c.patientunitstayid::text) as patient_id,
                  c.patienthealthsystemstayid::text as encounter_id,
                  c.patientunitstayid::text as stay_id,
                  c.index_offset,
                  c.landmark_offset,
                  p.unitdischargeoffset,
                  p.hospitaldischargeoffset,
                  p.hospitaldischargestatus,
                  case when p.hospitaldischargestatus = 'Expired' then 1 else 0 end as death_event,
                  case when p.hospitaldischargestatus = 'Alive' then 1 else 0 end as alive_discharge_event,
                  (p.hospitaldischargeoffset - c.landmark_offset) / 1440.0 as time_to_hospital_exit_days,
                  (p.unitdischargeoffset - c.landmark_offset) / 1440.0 as time_to_icu_exit_days
                from eicu_formal_index c
                join eicuii.patient p
                  on c.patientunitstayid = p.patientunitstayid
                """
            )
            copy_query_to_csv(
                cur,
                "select * from eicu_survival_source order by patient_id, index_offset, stay_id",
                outdir / "eicu_survival_source.csv",
            )
            log(f"eICU survival rows: {fetch_scalar(cur, 'select count(*) from eicu_survival_source')}")


def extract_sicdb(outdir: Path) -> None:
    log("SICdb: extracting post-landmark hospital death/discharge offsets...")
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            sicdb_formal.create_source_tables(cur)
            cur.execute(
                f"""
                create temp table sicdb_survival_source as
                with hospital_times as (
                  select
                    c.caseid,
                    case
                      when rawc.hospitalstaydays ~ '{NUMERIC_RE_SQL}'
                      then nullif(rawc.hospitalstaydays, '')::numeric * 86400.0
                      else null
                    end as hospital_discharge_offset_sec
                  from sicdb_formal_index c
                  join raw.cases rawc
                    on c.caseid = rawc.caseid
                )
                select
                  c.patientid::text as patient_id,
                  c.caseid::text as encounter_id,
                  c.caseid::text as stay_id,
                  c.index_offset_sec,
                  c.landmark_offset_sec,
                  c.icu_out_offset_sec,
                  ht.hospital_discharge_offset_sec,
                  c.offset_of_death_sec as death_offset_sec,
                  c.hospital_mortality,
                  case
                    when c.hospital_mortality = 1
                     and c.offset_of_death_sec >= c.landmark_offset_sec
                    then 1 else 0
                  end as death_event,
                  case
                    when c.hospital_mortality = 0
                     and ht.hospital_discharge_offset_sec >= c.landmark_offset_sec
                    then 1 else 0
                  end as alive_discharge_event,
                  case
                    when c.hospital_mortality = 1
                     and c.offset_of_death_sec >= c.landmark_offset_sec
                    then (c.offset_of_death_sec - c.landmark_offset_sec) / 86400.0
                    when ht.hospital_discharge_offset_sec >= c.landmark_offset_sec
                    then (ht.hospital_discharge_offset_sec - c.landmark_offset_sec) / 86400.0
                    else null
                  end as time_to_hospital_exit_days,
                  (c.icu_out_offset_sec - c.landmark_offset_sec) / 86400.0 as time_to_icu_exit_days
                from sicdb_formal_index c
                left join hospital_times ht
                  on c.caseid = ht.caseid
                """
            )
            copy_query_to_csv(
                cur,
                "select * from sicdb_survival_source order by patient_id, index_offset_sec, stay_id",
                outdir / "sicdb_survival_source.csv",
            )
            log(f"SICdb survival rows: {fetch_scalar(cur, 'select count(*) from sicdb_survival_source')}")


def write_summary(outdir: Path) -> None:
    rows: list[dict[str, Any]] = []
    for database, file_name in [
        ("MIMIC-IV", "mimic_survival_source.csv"),
        ("eICU", "eicu_survival_source.csv"),
        ("SICdb", "sicdb_survival_source.csv"),
    ]:
        path = outdir / file_name
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as f:
            data = list(csv.DictReader(f))
        analyzable = [
            r for r in data
            if r.get("time_to_hospital_exit_days") not in ("", "NA")
            and float(r["time_to_hospital_exit_days"]) > 0
        ]
        rows.append(
            {
                "database": database,
                "rows": len(data),
                "hospital_time_analyzable": len(analyzable),
                "death_events": sum(1 for r in analyzable if r.get("death_event") == "1"),
                "alive_discharge_events": sum(1 for r in analyzable if r.get("alive_discharge_event") == "1"),
            }
        )
    write_csv(
        outdir / "survival_source_summary.csv",
        rows,
        ["database", "rows", "hospital_time_analyzable", "death_events", "alive_discharge_events"],
    )
    md = [
        "# Survival/competing-risk source extraction",
        "",
        "database | rows | hospital_time_analyzable | death_events | alive_discharge_events",
        "--- | ---: | ---: | ---: | ---:",
    ]
    for row in rows:
        md.append(
            f"{row['database']} | {row['rows']} | {row['hospital_time_analyzable']} | "
            f"{row['death_events']} | {row['alive_discharge_events']}"
        )
    (outdir / "survival_source_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract post-landmark survival/competing-risk source fields.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--mimic-server", default="mimiciv")
    parser.add_argument("--eicu-server", default="eicu")
    parser.add_argument("--database", choices=["all", "mimic", "eicu", "sicdb"], default="all")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.database in ("all", "mimic"):
        extract_mimic(args.outdir, args.mimic_server)
    if args.database in ("all", "eicu"):
        extract_eicu(args.outdir, args.eicu_server)
    if args.database in ("all", "sicdb"):
        extract_sicdb(args.outdir)
    write_summary(args.outdir)
    log(f"Summary: {args.outdir / 'survival_source_summary.md'}")


if __name__ == "__main__":
    main()
