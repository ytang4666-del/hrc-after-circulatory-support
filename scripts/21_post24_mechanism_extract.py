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


DEFAULT_OUTDIR = Path("outputs/mechanism_clinical_utility")


def load_script_module(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_DIR / file_name)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mimic_formal = load_script_module("mimic_formal_extract_post24", "03_mimic_formal_cohort_extract.py")
eicu_formal = load_script_module("eicu_formal_extract_post24", "07_eicu_formal_cohort_extract.py")
sicdb_formal = load_script_module("sicdb_formal_extract_post24", "10_sicdb_formal_cohort_extract.py")


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
    log("MIMIC-IV: extracting 24-72h organ dysfunction and CVP/MPP mechanism variables...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '120min'")
            cur.execute("SET work_mem = '768MB'")
            mimic_formal.create_source_tables(cur)

            cur.execute(
                """
                create temp table mimic_post24_window as
                select
                  c.*,
                  c.landmark_time as post24_start,
                  least(c.outtime, c.index_time + interval '72 hours') as post24_end,
                  extract(epoch from (
                    least(c.outtime, c.index_time + interval '72 hours') - c.landmark_time
                  )) / 3600.0 as post24_followup_hours
                from formal_index c
                where least(c.outtime, c.index_time + interval '72 hours') > c.landmark_time
                """
            )
            cur.execute("create index on mimic_post24_window(stay_id)")

            log("MIMIC-IV: post-24h vasopressor burden...")
            cur.execute(
                """
                create temp table mimic_post24_vaso as
                with va as (
                  select
                    c.stay_id,
                    (
                      coalesce(v.norepinephrine, 0)
                      + coalesce(v.epinephrine, 0)
                      + coalesce(v.phenylephrine, 0) / 10.0
                      + coalesce(v.dopamine, 0) / 100.0
                      + coalesce(v.vasopressin, 0) * 2.5
                    )::double precision as ne_equiv,
                    greatest(v.starttime, c.post24_start) as overlap_start,
                    least(coalesce(v.endtime, v.starttime), c.post24_end) as overlap_end
                  from mimic_post24_window c
                  join mimiciv_derived.vasoactive_agent v
                    on c.stay_id = v.stay_id
                   and coalesce(v.endtime, v.starttime) > c.post24_start
                   and v.starttime < c.post24_end
                ),
                pieces as (
                  select
                    c.stay_id,
                    va.ne_equiv,
                    greatest(0, extract(epoch from (va.overlap_end - va.overlap_start)) / 3600.0) as overlap_hours
                  from mimic_post24_window c
                  left join va on c.stay_id = va.stay_id
                )
                select
                  c.stay_id,
                  coalesce(sum(p.ne_equiv * p.overlap_hours), 0) / nullif(max(c.post24_followup_hours), 0) as post24_vaso_burden,
                  coalesce(sum(p.overlap_hours), 0) as post24_vaso_hours
                from mimic_post24_window c
                left join pieces p on c.stay_id = p.stay_id
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_post24_vaso(stay_id)")

            log("MIMIC-IV: post-24h urine output...")
            cur.execute(
                """
                create temp table mimic_post24_urine as
                select
                  c.stay_id,
                  sum(u.urineoutput::double precision) as post24_uo_ml,
                  count(*) as post24_uo_n
                from mimic_post24_window c
                left join mimiciv_derived.urine_output u
                  on c.stay_id = u.stay_id
                 and u.charttime >= c.post24_start
                 and u.charttime < c.post24_end
                 and u.urineoutput between 0 and 5000
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_post24_urine(stay_id)")

            log("MIMIC-IV: post-24h creatinine/lactate...")
            cur.execute(
                """
                create temp table mimic_post24_labs as
                select
                  c.stay_id,
                  max(ch.creatinine::double precision) filter (where ch.creatinine between 0.1 and 30) as post24_creatinine_max,
                  avg(ch.creatinine::double precision) filter (where ch.creatinine between 0.1 and 30) as post24_creatinine_mean,
                  count(ch.creatinine) filter (where ch.creatinine between 0.1 and 30) as post24_creatinine_n,
                  min(bg.lactate::double precision) filter (where bg.lactate between 0.1 and 30) as post24_lactate_min,
                  avg(bg.lactate::double precision) filter (where bg.lactate between 0.1 and 30) as post24_lactate_mean,
                  count(bg.lactate) filter (where bg.lactate between 0.1 and 30) as post24_lactate_n
                from mimic_post24_window c
                left join mimiciv_derived.chemistry ch
                  on c.hadm_id = ch.hadm_id
                 and ch.charttime >= c.post24_start
                 and ch.charttime < c.post24_end
                left join mimiciv_derived.bg bg
                  on c.hadm_id = bg.hadm_id
                 and bg.charttime >= c.post24_start
                 and bg.charttime < c.post24_end
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_post24_labs(stay_id)")

            log("MIMIC-IV: post-24h RRT and mechanical ventilation...")
            cur.execute(
                """
                create temp table mimic_post24_support as
                select
                  c.stay_id,
                  max(case when r.dialysis_active = 1 or r.dialysis_present = 1 then 1 else 0 end) as post24_rrt_any,
                  coalesce(sum(
                    greatest(
                      0,
                      extract(epoch from (
                        least(v.endtime, c.post24_end) - greatest(v.starttime, c.post24_start)
                      )) / 3600.0
                    )
                  ) filter (where v.ventilation_status in ('InvasiveVent', 'Tracheostomy')), 0) as post24_invasive_vent_hours
                from mimic_post24_window c
                left join mimiciv_derived.rrt r
                  on c.stay_id = r.stay_id
                 and r.charttime >= c.post24_start
                 and r.charttime < c.post24_end
                left join mimiciv_derived.ventilation v
                  on c.stay_id = v.stay_id
                 and v.endtime > c.post24_start
                 and v.starttime < c.post24_end
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_post24_support(stay_id)")

            log("MIMIC-IV: CVP and MPP windows...")
            cur.execute(
                """
                create temp table mimic_cvp as
                select
                  c.stay_id,
                  avg(ce.valuenum::double precision) filter (
                    where ce.charttime >= c.index_time
                      and ce.charttime < c.index_time + interval '6 hours'
                      and ce.valuenum between 0 and 40
                  ) as baseline_cvp_mean,
                  count(*) filter (
                    where ce.charttime >= c.index_time
                      and ce.charttime < c.index_time + interval '6 hours'
                      and ce.valuenum between 0 and 40
                  ) as baseline_cvp_n,
                  avg(ce.valuenum::double precision) filter (
                    where ce.charttime >= c.index_time + interval '6 hours'
                      and ce.charttime < c.index_time + interval '24 hours'
                      and ce.valuenum between 0 and 40
                  ) as response_cvp_mean,
                  count(*) filter (
                    where ce.charttime >= c.index_time + interval '6 hours'
                      and ce.charttime < c.index_time + interval '24 hours'
                      and ce.valuenum between 0 and 40
                  ) as response_cvp_n
                from formal_index c
                left join mimiciv_icu.chartevents ce
                  on c.stay_id = ce.stay_id
                 and ce.itemid = 220074
                 and ce.charttime >= c.index_time
                 and ce.charttime < c.index_time + interval '24 hours'
                group by c.stay_id
                """
            )
            cur.execute("create index on mimic_cvp(stay_id)")

            cur.execute(
                """
                create temp table mimic_post24_mechanism as
                select
                  'MIMIC-IV'::text as database,
                  c.subject_id::text as patient_id,
                  c.hadm_id::text as encounter_id,
                  c.stay_id::text as stay_id,
                  c.post24_followup_hours,
                  pv.post24_vaso_burden,
                  pv.post24_vaso_hours,
                  pu.post24_uo_ml,
                  pu.post24_uo_n,
                  pu.post24_uo_ml / nullif(w.weight, 0) / nullif(c.post24_followup_hours, 0) as post24_uo_ml_kg_h,
                  pl.post24_creatinine_max,
                  pl.post24_creatinine_mean,
                  pl.post24_creatinine_n,
                  pl.post24_lactate_min,
                  pl.post24_lactate_mean,
                  pl.post24_lactate_n,
                  ps.post24_rrt_any,
                  ps.post24_invasive_vent_hours,
                  case when ps.post24_invasive_vent_hours > 0 then 1 else 0 end as post24_invasive_vent_any,
                  cvp.baseline_cvp_mean,
                  cvp.baseline_cvp_n,
                  cvp.response_cvp_mean,
                  cvp.response_cvp_n
                from mimic_post24_window c
                left join mimiciv_derived.first_day_weight w on c.stay_id = w.stay_id
                left join mimic_post24_vaso pv on c.stay_id = pv.stay_id
                left join mimic_post24_urine pu on c.stay_id = pu.stay_id
                left join mimic_post24_labs pl on c.stay_id = pl.stay_id
                left join mimic_post24_support ps on c.stay_id = ps.stay_id
                left join mimic_cvp cvp on c.stay_id = cvp.stay_id
                """
            )
            copy_query_to_csv(
                cur,
                "select * from mimic_post24_mechanism order by patient_id, stay_id",
                outdir / "mimic_post24_mechanism.csv",
            )
            log(f"MIMIC-IV mechanism rows: {fetch_scalar(cur, 'select count(*) from mimic_post24_mechanism')}")


def extract_eicu(outdir: Path, server: str) -> None:
    log("eICU: extracting 24-72h organ dysfunction variables...")
    active_expr = eicu_formal.active_agent_expression("pi")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            eicu_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table eicu_post24_window as
                select
                  c.*,
                  c.landmark_offset as post24_start_offset,
                  least(c.unitdischargeoffset, c.index_offset + 4320) as post24_end_offset,
                  (least(c.unitdischargeoffset, c.index_offset + 4320) - c.landmark_offset) / 60.0 as post24_followup_hours
                from eicu_formal_index c
                where least(c.unitdischargeoffset, c.index_offset + 4320) > c.landmark_offset
                """
            )
            cur.execute("create index on eicu_post24_window(patientunitstayid)")

            log("eICU: post-24h vasoactive burden...")
            cur.execute(
                f"""
                create temp table eicu_post24_vaso as
                with burden as (
                  select
                    c.patientunitstayid,
                    pi.chartoffset,
                    ({active_expr})::double precision as active_agent_count
                  from eicu_post24_window c
                  join eicuii.pivoted_infusion pi
                    on c.patientunitstayid = pi.patientunitstayid
                   and pi.chartoffset >= c.post24_start_offset
                   and pi.chartoffset < c.post24_end_offset
                )
                select
                  c.patientunitstayid,
                  coalesce(avg(b.active_agent_count), 0) as post24_vaso_burden,
                  coalesce(sum(case when b.active_agent_count > 0 then 1 else 0 end), 0) as post24_vaso_observed_units
                from eicu_post24_window c
                left join burden b on c.patientunitstayid = b.patientunitstayid
                group by c.patientunitstayid
                """
            )
            cur.execute("create index on eicu_post24_vaso(patientunitstayid)")

            log("eICU: post-24h urine and labs...")
            cur.execute(
                """
                create temp table eicu_post24_urine as
                select
                  c.patientunitstayid,
                  sum(u.urineoutput::double precision) as post24_uo_ml,
                  count(*) filter (where u.urineoutput is not null) as post24_uo_n
                from eicu_post24_window c
                left join eicuii.pivoted_uo u
                  on c.patientunitstayid = u.patientunitstayid
                 and u.chartoffset >= c.post24_start_offset
                 and u.chartoffset < c.post24_end_offset
                 and u.urineoutput between 0 and 5000
                group by c.patientunitstayid
                """
            )
            cur.execute("create index on eicu_post24_urine(patientunitstayid)")

            cur.execute(
                """
                create temp table eicu_post24_labs as
                select
                  c.patientunitstayid,
                  max(l.creatinine::double precision) filter (where l.creatinine between 0.1 and 30) as post24_creatinine_max,
                  avg(l.creatinine::double precision) filter (where l.creatinine between 0.1 and 30) as post24_creatinine_mean,
                  count(l.creatinine) filter (where l.creatinine between 0.1 and 30) as post24_creatinine_n,
                  min(l.lactate::double precision) filter (where l.lactate between 0.1 and 30) as post24_lactate_min,
                  avg(l.lactate::double precision) filter (where l.lactate between 0.1 and 30) as post24_lactate_mean,
                  count(l.lactate) filter (where l.lactate between 0.1 and 30) as post24_lactate_n
                from eicu_post24_window c
                left join eicuii.pivoted_lab l
                  on c.patientunitstayid = l.patientunitstayid
                 and l.chartoffset >= c.post24_start_offset
                 and l.chartoffset < c.post24_end_offset
                group by c.patientunitstayid
                """
            )
            cur.execute("create index on eicu_post24_labs(patientunitstayid)")

            cur.execute(
                """
                create temp table eicu_post24_mechanism as
                select
                  'eICU'::text as database,
                  coalesce(nullif(c.uniquepid, ''), c.patientunitstayid::text) as patient_id,
                  c.patienthealthsystemstayid::text as encounter_id,
                  c.patientunitstayid::text as stay_id,
                  c.post24_followup_hours,
                  pv.post24_vaso_burden,
                  pv.post24_vaso_observed_units,
                  pu.post24_uo_ml,
                  pu.post24_uo_n,
                  pu.post24_uo_ml / nullif(w.weight_kg, 0) / nullif(c.post24_followup_hours, 0) as post24_uo_ml_kg_h,
                  pl.post24_creatinine_max,
                  pl.post24_creatinine_mean,
                  pl.post24_creatinine_n,
                  pl.post24_lactate_min,
                  pl.post24_lactate_mean,
                  pl.post24_lactate_n,
                  null::double precision as post24_rrt_any,
                  null::double precision as post24_invasive_vent_hours,
                  null::double precision as post24_invasive_vent_any,
                  null::double precision as baseline_cvp_mean,
                  null::double precision as baseline_cvp_n,
                  null::double precision as response_cvp_mean,
                  null::double precision as response_cvp_n
                from eicu_post24_window c
                left join eicu_weight w on c.patientunitstayid = w.patientunitstayid
                left join eicu_post24_vaso pv on c.patientunitstayid = pv.patientunitstayid
                left join eicu_post24_urine pu on c.patientunitstayid = pu.patientunitstayid
                left join eicu_post24_labs pl on c.patientunitstayid = pl.patientunitstayid
                """
            )
            copy_query_to_csv(
                cur,
                "select * from eicu_post24_mechanism order by patient_id, stay_id",
                outdir / "eicu_post24_mechanism.csv",
            )
            log(f"eICU mechanism rows: {fetch_scalar(cur, 'select count(*) from eicu_post24_mechanism')}")


def extract_sicdb(outdir: Path) -> None:
    log("SICdb: extracting 24-72h organ dysfunction variables...")
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '120min'")
            cur.execute("SET work_mem = '768MB'")
            sicdb_formal.create_source_tables(cur)
            cur.execute(
                """
                create temp table sicdb_post24_window as
                select
                  c.*,
                  c.landmark_offset_sec as post24_start_sec,
                  least(c.icu_out_offset_sec, c.index_offset_sec + 259200) as post24_end_sec,
                  (least(c.icu_out_offset_sec, c.index_offset_sec + 259200) - c.landmark_offset_sec) / 3600.0 as post24_followup_hours
                from sicdb_formal_index c
                where least(c.icu_out_offset_sec, c.index_offset_sec + 259200) > c.landmark_offset_sec
                """
            )
            cur.execute("create index on sicdb_post24_window(caseid)")

            log("SICdb: post-24h vasopressor burden...")
            cur.execute(
                """
                create temp table sicdb_post24_vaso as
                with pieces as (
                  select
                    c.caseid,
                    greatest(
                      0,
                      least(v.end_offset_sec, c.post24_end_sec)
                      - greatest(v.start_offset_sec, c.post24_start_sec)
                    ) / 3600.0 as overlap_hours
                  from sicdb_post24_window c
                  left join sicdb_vaso_events v
                    on c.caseid = v.caseid
                   and v.end_offset_sec > c.post24_start_sec
                   and v.start_offset_sec < c.post24_end_sec
                )
                select
                  c.caseid,
                  coalesce(sum(p.overlap_hours), 0) / nullif(max(c.post24_followup_hours), 0) as post24_vaso_burden,
                  coalesce(sum(p.overlap_hours), 0) as post24_vaso_hours
                from sicdb_post24_window c
                left join pieces p on c.caseid = p.caseid
                group by c.caseid
                """
            )
            cur.execute("create index on sicdb_post24_vaso(caseid)")

            log("SICdb: scanning post-24h MAP/urine subset...")
            cur.execute(
                """
                create temp table sicdb_post24_signal_subset as
                select
                  s.caseid,
                  s.dataid,
                  nullif(s.offset, '')::numeric as offset_sec,
                  nullif(s.val, '')::numeric as val_num
                from sicdb_post24_window c
                join lateral (
                  select s.caseid, s.dataid, s.offset, s.val
                  from raw.data_float_h s
                  where s.caseid = c.caseid
                    and s.dataid = any(%s)
                    and s.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and s.val ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and nullif(s.offset, '')::numeric >= c.post24_start_sec
                    and nullif(s.offset, '')::numeric < c.post24_end_sec
                ) s on true
                """,
                (sicdb_formal.URINE_IDS,),
            )
            cur.execute("create index on sicdb_post24_signal_subset(caseid, dataid, offset_sec)")
            cur.execute("analyze sicdb_post24_signal_subset")
            log(f"SICdb: post-24h signal rows kept: {fetch_scalar(cur, 'select count(*) from sicdb_post24_signal_subset')}")

            cur.execute(
                """
                create temp table sicdb_post24_urine as
                select
                  c.caseid,
                  sum(s.val_num) filter (where s.dataid = any(%s) and s.val_num between 0 and 5000) as post24_uo_ml,
                  count(*) filter (where s.dataid = any(%s) and s.val_num between 0 and 5000) as post24_uo_n
                from sicdb_post24_window c
                left join sicdb_post24_signal_subset s on c.caseid = s.caseid
                group by c.caseid
                """,
                (sicdb_formal.URINE_IDS, sicdb_formal.URINE_IDS),
            )
            cur.execute("create index on sicdb_post24_urine(caseid)")

            log("SICdb: post-24h creatinine/lactate...")
            cur.execute(
                """
                create temp table sicdb_post24_labs as
                select
                  c.caseid,
                  max(nullif(l.laboratoryvalue, '')::numeric) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_creatinine_max,
                  avg(nullif(l.laboratoryvalue, '')::numeric) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_creatinine_mean,
                  count(*) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_creatinine_n,
                  min(nullif(l.laboratoryvalue, '')::numeric) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_lactate_min,
                  avg(nullif(l.laboratoryvalue, '')::numeric) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_lactate_mean,
                  count(*) filter (
                    where l.laboratoryid = any(%s)
                      and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
                  ) as post24_lactate_n
                from sicdb_post24_window c
                left join raw.laboratory l
                  on c.caseid = l.caseid
                 and l.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
                 and l.laboratoryvalue ~ '^-?[0-9]+(\\.[0-9]+)?$'
                 and nullif(l.offset, '')::numeric >= c.post24_start_sec
                 and nullif(l.offset, '')::numeric < c.post24_end_sec
                 and l.laboratoryid = any(%s)
                group by c.caseid
                """,
                (
                    sicdb_formal.CREATININE_IDS,
                    sicdb_formal.CREATININE_IDS,
                    sicdb_formal.CREATININE_IDS,
                    sicdb_formal.LACTATE_IDS,
                    sicdb_formal.LACTATE_IDS,
                    sicdb_formal.LACTATE_IDS,
                    sicdb_formal.CREATININE_IDS + sicdb_formal.LACTATE_IDS,
                ),
            )
            cur.execute("create index on sicdb_post24_labs(caseid)")

            cur.execute(
                """
                create temp table sicdb_post24_mechanism as
                select
                  'SICdb'::text as database,
                  c.patientid::text as patient_id,
                  c.caseid::text as encounter_id,
                  c.caseid::text as stay_id,
                  c.post24_followup_hours,
                  pv.post24_vaso_burden,
                  pv.post24_vaso_hours,
                  pu.post24_uo_ml,
                  pu.post24_uo_n,
                  pu.post24_uo_ml / nullif(c.weight_kg, 0) / nullif(c.post24_followup_hours, 0) as post24_uo_ml_kg_h,
                  pl.post24_creatinine_max,
                  pl.post24_creatinine_mean,
                  pl.post24_creatinine_n,
                  pl.post24_lactate_min,
                  pl.post24_lactate_mean,
                  pl.post24_lactate_n,
                  null::double precision as post24_rrt_any,
                  null::double precision as post24_invasive_vent_hours,
                  null::double precision as post24_invasive_vent_any,
                  null::double precision as baseline_cvp_mean,
                  null::double precision as baseline_cvp_n,
                  null::double precision as response_cvp_mean,
                  null::double precision as response_cvp_n
                from sicdb_post24_window c
                left join sicdb_post24_vaso pv on c.caseid = pv.caseid
                left join sicdb_post24_urine pu on c.caseid = pu.caseid
                left join sicdb_post24_labs pl on c.caseid = pl.caseid
                """
            )
            copy_query_to_csv(
                cur,
                "select * from sicdb_post24_mechanism order by patient_id, stay_id",
                outdir / "sicdb_post24_mechanism.csv",
            )
            log(f"SICdb mechanism rows: {fetch_scalar(cur, 'select count(*) from sicdb_post24_mechanism')}")


def write_summary(outdir: Path) -> None:
    rows: list[dict[str, Any]] = []
    for database, file_name in [
        ("MIMIC-IV", "mimic_post24_mechanism.csv"),
        ("eICU", "eicu_post24_mechanism.csv"),
        ("SICdb", "sicdb_post24_mechanism.csv"),
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
                "with_post24_uo": sum(1 for r in data if r.get("post24_uo_ml_kg_h") not in ("", "NA")),
                "with_post24_creatinine": sum(1 for r in data if r.get("post24_creatinine_max") not in ("", "NA")),
                "with_post24_lactate": sum(1 for r in data if r.get("post24_lactate_min") not in ("", "NA")),
                "with_cvp": sum(1 for r in data if r.get("baseline_cvp_mean") not in ("", "NA")),
            }
        )
    write_csv(
        outdir / "post24_mechanism_extraction_summary.csv",
        rows,
        ["database", "rows", "with_post24_uo", "with_post24_creatinine", "with_post24_lactate", "with_cvp"],
    )
    md = [
        "# Post-24h mechanism extraction summary",
        "",
        "database | rows | with_post24_uo | with_post24_creatinine | with_post24_lactate | with_cvp",
        "--- | ---: | ---: | ---: | ---: | ---:",
    ]
    for row in rows:
        md.append(
            f"{row['database']} | {row['rows']} | {row['with_post24_uo']} | "
            f"{row['with_post24_creatinine']} | {row['with_post24_lactate']} | {row['with_cvp']}"
        )
    (outdir / "post24_mechanism_extraction_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract post-24h organ dysfunction and mechanism variables.")
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
    log(f"Summary: {args.outdir / 'post24_mechanism_extraction_summary.md'}")


if __name__ == "__main__":
    main()
