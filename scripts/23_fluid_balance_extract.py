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


DEFAULT_OUTDIR = Path("outputs/final_methodology_extensions")


def load_script_module(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_DIR / file_name)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mimic_formal = load_script_module("mimic_formal_extract_fluid", "03_mimic_formal_cohort_extract.py")
eicu_formal = load_script_module("eicu_formal_extract_fluid", "07_eicu_formal_cohort_extract.py")
sicdb_formal = load_script_module("sicdb_formal_extract_fluid", "10_sicdb_formal_cohort_extract.py")


def log(message: str) -> None:
    print(message, flush=True)


def fetch_scalar(cur, sql: str) -> Any:
    row = cur.execute(sql).fetchone()
    if isinstance(row, dict):
        return next(iter(row.values()))
    return row[0]


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


def extract_mimic(outdir: Path, server: str) -> None:
    log("MIMIC-IV: extracting 0-24h fluid balance...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            mimic_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table mimic_fluid_input_0_24 as
                with input_pieces as (
                  select
                    c.stay_id,
                    ie.amount::double precision as amount_ml,
                    ie.starttime,
                    ie.endtime,
                    c.index_time,
                    c.landmark_time,
                    extract(epoch from (ie.endtime - ie.starttime)) / 3600.0 as duration_h,
                    greatest(
                      0,
                      extract(epoch from (
                        least(coalesce(ie.endtime, ie.starttime), c.landmark_time)
                        - greatest(ie.starttime, c.index_time)
                      )) / 3600.0
                    ) as overlap_h
                  from formal_index c
                  join mimiciv_icu.inputevents ie
                    on c.stay_id = ie.stay_id
                   and coalesce(ie.endtime, ie.starttime) > c.index_time
                   and ie.starttime < c.landmark_time
                   and lower(ie.amountuom) = 'ml'
                   and ie.amount > 0
                   and ie.amount < 100000
                   and coalesce(ie.statusdescription, '') not in ('Rewritten')
                )
                select
                  stay_id,
                  sum(
                    case
                      when duration_h > 0 and overlap_h > 0 then amount_ml * overlap_h / duration_h
                      when duration_h <= 0 and starttime >= index_time and starttime < landmark_time then amount_ml
                      else 0
                    end
                  ) as fluid_input_0_24_ml,
                  count(*) as fluid_input_0_24_n
                from input_pieces
                group by stay_id
                """
            )
            cur.execute("create index on mimic_fluid_input_0_24(stay_id)")
            cur.execute(
                """
                create temp table mimic_fluid_output_0_24 as
                select
                  c.stay_id,
                  sum(oe.value::double precision) as fluid_output_0_24_ml,
                  count(*) as fluid_output_0_24_n
                from formal_index c
                left join mimiciv_icu.outputevents oe
                  on c.stay_id = oe.stay_id
                 and oe.charttime >= c.index_time
                 and oe.charttime < c.landmark_time
                 and lower(oe.valueuom) = 'ml'
                 and oe.value > 0
                 and oe.value < 100000
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_fluid_output_0_24(stay_id)")
            cur.execute(
                """
                create temp table mimic_fluid_balance as
                select
                  'MIMIC-IV'::text as database,
                  c.subject_id::text as patient_id,
                  c.hadm_id::text as encounter_id,
                  c.stay_id::text as stay_id,
                  w.weight::double precision as weight_kg,
                  coalesce(i.fluid_input_0_24_ml, 0) as fluid_input_0_24_ml,
                  coalesce(o.fluid_output_0_24_ml, 0) as fluid_output_0_24_ml,
                  coalesce(i.fluid_input_0_24_ml, 0) - coalesce(o.fluid_output_0_24_ml, 0) as fluid_balance_0_24_ml,
                  (coalesce(i.fluid_input_0_24_ml, 0) - coalesce(o.fluid_output_0_24_ml, 0)) / nullif(w.weight, 0) as fluid_balance_0_24_ml_kg,
                  coalesce(i.fluid_input_0_24_n, 0) as fluid_input_0_24_n,
                  coalesce(o.fluid_output_0_24_n, 0) as fluid_output_0_24_n,
                  null::double precision as fluid_exposure_0_24_mean,
                  null::double precision as fluid_exposure_0_24_n,
                  'net_balance'::text as fluid_metric_type
                from formal_index c
                left join mimiciv_derived.first_day_weight w on c.stay_id = w.stay_id
                left join mimic_fluid_input_0_24 i on c.stay_id = i.stay_id
                left join mimic_fluid_output_0_24 o on c.stay_id = o.stay_id
                """
            )
            copy_query_to_csv(
                cur,
                "select * from mimic_fluid_balance order by patient_id, stay_id",
                outdir / "mimic_fluid_balance.csv",
            )
            log(f"MIMIC-IV fluid rows: {fetch_scalar(cur, 'select count(*) from mimic_fluid_balance')}")


def extract_eicu(outdir: Path, server: str) -> None:
    log("eICU: extracting 0-24h fluid balance from intakeoutput...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            eicu_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table eicu_io_events_0_24 as
                select distinct
                  c.patientunitstayid,
                  io.intakeoutputoffset,
                  io.intaketotal::double precision as intaketotal,
                  io.outputtotal::double precision as outputtotal,
                  io.nettotal::double precision as nettotal
                from eicu_formal_index c
                join eicuii.intakeoutput io
                  on c.patientunitstayid = io.patientunitstayid
                 and io.intakeoutputoffset >= c.index_offset
                 and io.intakeoutputoffset < c.landmark_offset
                 and (
                      io.intaketotal is not null
                   or io.outputtotal is not null
                   or io.nettotal is not null
                 )
                where abs(coalesce(io.intaketotal, 0)) < 100000
                  and abs(coalesce(io.outputtotal, 0)) < 100000
                  and abs(coalesce(io.nettotal, 0)) < 100000
                """
            )
            cur.execute("create index on eicu_io_events_0_24(patientunitstayid)")
            cur.execute(
                """
                create temp table eicu_fluid_balance as
                select
                  'eICU'::text as database,
                  coalesce(nullif(c.uniquepid, ''), c.patientunitstayid::text) as patient_id,
                  c.patienthealthsystemstayid::text as encounter_id,
                  c.patientunitstayid::text as stay_id,
                  w.weight_kg::double precision as weight_kg,
                  coalesce(sum(greatest(e.intaketotal, 0)), 0) as fluid_input_0_24_ml,
                  coalesce(sum(greatest(e.outputtotal, 0)), 0) as fluid_output_0_24_ml,
                  coalesce(sum(e.nettotal), 0) as fluid_balance_0_24_ml,
                  coalesce(sum(e.nettotal), 0) / nullif(w.weight_kg, 0) as fluid_balance_0_24_ml_kg,
                  count(e.intaketotal) filter (where e.intaketotal is not null) as fluid_input_0_24_n,
                  count(e.outputtotal) filter (where e.outputtotal is not null) as fluid_output_0_24_n,
                  null::double precision as fluid_exposure_0_24_mean,
                  null::double precision as fluid_exposure_0_24_n,
                  'net_balance'::text as fluid_metric_type
                from eicu_formal_index c
                left join eicu_weight w on c.patientunitstayid = w.patientunitstayid
                left join eicu_io_events_0_24 e on c.patientunitstayid = e.patientunitstayid
                group by c.uniquepid, c.patientunitstayid, c.patienthealthsystemstayid, w.weight_kg
                """
            )
            copy_query_to_csv(
                cur,
                "select * from eicu_fluid_balance order by patient_id, stay_id",
                outdir / "eicu_fluid_balance.csv",
            )
            log(f"eICU fluid rows: {fetch_scalar(cur, 'select count(*) from eicu_fluid_balance')}")


def extract_sicdb(outdir: Path) -> None:
    log("SICdb: extracting 0-24h fluid exposure signal...")
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            sicdb_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table sicdb_fluid_exposure as
                select
                  'SICdb'::text as database,
                  c.patientid::text as patient_id,
                  c.caseid::text as encounter_id,
                  c.caseid::text as stay_id,
                  c.weight_kg::double precision as weight_kg,
                  null::double precision as fluid_input_0_24_ml,
                  null::double precision as fluid_output_0_24_ml,
                  null::double precision as fluid_balance_0_24_ml,
                  null::double precision as fluid_balance_0_24_ml_kg,
                  null::double precision as fluid_input_0_24_n,
                  null::double precision as fluid_output_0_24_n,
                  avg(nullif(s.val, '')::numeric) filter (
                    where s.dataid in ('2200', '2201')
                      and nullif(s.val, '')::numeric between 0 and 10000
                  ) as fluid_exposure_0_24_mean,
                  count(*) filter (
                    where s.dataid in ('2200', '2201')
                      and nullif(s.val, '')::numeric between 0 and 10000
                  ) as fluid_exposure_0_24_n,
                  'fluid_exposure_not_net_balance'::text as fluid_metric_type
                from sicdb_formal_index c
                left join lateral (
                  select s.caseid, s.dataid, s.offset, s.val
                  from raw.data_float_h s
                  where s.caseid = c.caseid
                    and s.dataid in ('2200', '2201')
                    and s.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and s.val ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and nullif(s.offset, '')::numeric >= c.index_offset_sec
                    and nullif(s.offset, '')::numeric < c.landmark_offset_sec
                ) s on true
                group by c.patientid, c.caseid, c.weight_kg
                """
            )
            copy_query_to_csv(
                cur,
                "select * from sicdb_fluid_exposure order by patient_id, stay_id",
                outdir / "sicdb_fluid_balance.csv",
            )
            log(f"SICdb fluid rows: {fetch_scalar(cur, 'select count(*) from sicdb_fluid_exposure')}")


def write_summary(outdir: Path) -> None:
    rows: list[dict[str, Any]] = []
    for database, file_name in [
        ("MIMIC-IV", "mimic_fluid_balance.csv"),
        ("eICU", "eicu_fluid_balance.csv"),
        ("SICdb", "sicdb_fluid_balance.csv"),
    ]:
        path = outdir / file_name
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as f:
            data = list(csv.DictReader(f))
        rows.append(
            {
                "database": database,
                "rows": len(data),
                "net_balance_available": sum(1 for r in data if r.get("fluid_balance_0_24_ml_kg") not in ("", "NA")),
                "fluid_exposure_available": sum(1 for r in data if r.get("fluid_exposure_0_24_mean") not in ("", "NA")),
                "metric_type": data[0].get("fluid_metric_type", "") if data else "",
            }
        )
    write_csv(
        outdir / "fluid_extraction_summary.csv",
        rows,
        ["database", "rows", "net_balance_available", "fluid_exposure_available", "metric_type"],
    )
    md = [
        "# Fluid extraction summary",
        "",
        "database | rows | net_balance_available | fluid_exposure_available | metric_type",
        "--- | ---: | ---: | ---: | ---",
    ]
    for row in rows:
        md.append(
            f"{row['database']} | {row['rows']} | {row['net_balance_available']} | "
            f"{row['fluid_exposure_available']} | {row['metric_type']}"
        )
    (outdir / "fluid_extraction_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract fluid balance/exposure sensitivity variables.")
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
    log(f"Summary: {args.outdir / 'fluid_extraction_summary.md'}")


if __name__ == "__main__":
    main()
