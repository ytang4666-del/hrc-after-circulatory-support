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


DEFAULT_OUTDIR = Path("outputs/time_window_sensitivity")


def load_script_module(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_DIR / file_name)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {file_name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mimic_formal = load_script_module("mimic_formal_extract", "03_mimic_formal_cohort_extract.py")
eicu_formal = load_script_module("eicu_formal_extract", "07_eicu_formal_cohort_extract.py")
sicdb_formal = load_script_module("sicdb_formal_extract", "10_sicdb_formal_cohort_extract.py")


def log(message: str) -> None:
    print(message, flush=True)


def fetch_scalar(cur, sql: str) -> Any:
    row = cur.execute(sql).fetchone()
    if isinstance(row, dict):
        return next(iter(row.values()))
    return row[0]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def copy_query_to_csv(cur, query: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with cur.copy(f"COPY ({query}) TO STDOUT WITH CSV HEADER") as copy:
        with out_path.open("wb") as f:
            for data in copy:
                f.write(data)


def create_window_strategies(cur) -> None:
    cur.execute(
        """
        create temp table window_strategies (
          window_strategy text primary key,
          baseline_start_h int not null,
          baseline_end_h int not null,
          response_start_h int not null,
          response_end_h int not null,
          baseline_hours double precision not null,
          response_hours double precision not null
        )
        """
    )
    cur.execute(
        """
        insert into window_strategies values
          ('primary_0_6_6_24', 0, 6, 6, 24, 6.0, 18.0),
          ('early_0_6_6_12', 0, 6, 6, 12, 6.0, 6.0),
          ('late_0_12_12_24', 0, 12, 12, 24, 12.0, 12.0)
        """
    )


def extract_mimic(outdir: Path, server: str) -> None:
    out_path = outdir / "mimic_time_window_cohort.csv"
    log("MIMIC-IV: creating source cohort...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            mimic_formal.create_source_tables(cur)
            create_window_strategies(cur)

            log("MIMIC-IV: aggregating MAP by alternative windows...")
            cur.execute(
                """
                create temp table mimic_tw_map as
                select
                  c.stay_id,
                  w.window_strategy,
                  avg(coalesce(v.mbp, v.mbp_ni)::double precision) filter (
                    where v.charttime >= c.index_time + w.baseline_start_h * interval '1 hour'
                      and v.charttime < c.index_time + w.baseline_end_h * interval '1 hour'
                  ) as baseline_map_mean,
                  count(*) filter (
                    where v.charttime >= c.index_time + w.baseline_start_h * interval '1 hour'
                      and v.charttime < c.index_time + w.baseline_end_h * interval '1 hour'
                      and coalesce(v.mbp, v.mbp_ni) is not null
                  ) as baseline_map_n,
                  avg(coalesce(v.mbp, v.mbp_ni)::double precision) filter (
                    where v.charttime >= c.index_time + w.response_start_h * interval '1 hour'
                      and v.charttime < c.index_time + w.response_end_h * interval '1 hour'
                  ) as response_map_mean,
                  count(*) filter (
                    where v.charttime >= c.index_time + w.response_start_h * interval '1 hour'
                      and v.charttime < c.index_time + w.response_end_h * interval '1 hour'
                      and coalesce(v.mbp, v.mbp_ni) is not null
                  ) as response_map_n
                from formal_index c
                cross join window_strategies w
                left join mimiciv_derived.vitalsign v
                  on c.stay_id = v.stay_id
                 and v.charttime >= c.index_time
                 and v.charttime < c.index_time + interval '24 hours'
                 and coalesce(v.mbp, v.mbp_ni) between 20 and 200
                group by c.stay_id, w.window_strategy
                """
            )
            cur.execute("create index on mimic_tw_map(stay_id, window_strategy)")

            log("MIMIC-IV: aggregating vasopressor burden by alternative windows...")
            cur.execute(
                """
                create temp table mimic_tw_vaso as
                with va as (
                  select
                    c.stay_id,
                    greatest(v.starttime, c.index_time) as overlap_start,
                    least(coalesce(v.endtime, v.starttime), c.index_time + interval '24 hours') as overlap_end,
                    (
                      coalesce(v.norepinephrine, 0)
                      + coalesce(v.epinephrine, 0)
                      + coalesce(v.phenylephrine, 0) / 10.0
                      + coalesce(v.dopamine, 0) / 100.0
                      + coalesce(v.vasopressin, 0) * 2.5
                    )::double precision as ne_equiv
                  from formal_index c
                  join mimiciv_derived.vasoactive_agent v
                    on c.stay_id = v.stay_id
                   and coalesce(v.endtime, v.starttime) > c.index_time
                   and v.starttime < c.index_time + interval '24 hours'
                ),
                pieces as (
                  select
                    c.stay_id,
                    w.window_strategy,
                    w.baseline_hours as denominator_baseline_hours,
                    w.response_hours as denominator_response_hours,
                    va.ne_equiv,
                    greatest(
                      0,
                      extract(epoch from (
                        least(va.overlap_end, c.index_time + w.baseline_end_h * interval '1 hour')
                        - greatest(va.overlap_start, c.index_time + w.baseline_start_h * interval '1 hour')
                      )) / 3600.0
                    ) as baseline_hours,
                    greatest(
                      0,
                      extract(epoch from (
                        least(va.overlap_end, c.index_time + w.response_end_h * interval '1 hour')
                        - greatest(va.overlap_start, c.index_time + w.response_start_h * interval '1 hour')
                      )) / 3600.0
                    ) as response_hours
                  from formal_index c
                  cross join window_strategies w
                  left join va
                    on c.stay_id = va.stay_id
                )
                select
                  stay_id,
                  window_strategy,
                  coalesce(sum(ne_equiv * baseline_hours), 0) / max(denominator_baseline_hours) as baseline_vaso_burden,
                  coalesce(sum(ne_equiv * response_hours), 0) / max(denominator_response_hours) as response_vaso_burden,
                  coalesce(sum(baseline_hours), 0) as baseline_vaso_observed_units,
                  coalesce(sum(response_hours), 0) as response_vaso_observed_units
                from pieces
                group by stay_id, window_strategy
                """
            )
            cur.execute("create index on mimic_tw_vaso(stay_id, window_strategy)")

            log("MIMIC-IV: aggregating urine output by alternative windows...")
            cur.execute(
                """
                create temp table mimic_tw_urine as
                select
                  c.stay_id,
                  w.window_strategy,
                  sum(u.urineoutput::double precision) filter (
                    where u.charttime >= c.index_time + w.baseline_start_h * interval '1 hour'
                      and u.charttime < c.index_time + w.baseline_end_h * interval '1 hour'
                  ) as baseline_uo_ml,
                  count(*) filter (
                    where u.charttime >= c.index_time + w.baseline_start_h * interval '1 hour'
                      and u.charttime < c.index_time + w.baseline_end_h * interval '1 hour'
                  ) as baseline_uo_n,
                  sum(u.urineoutput::double precision) filter (
                    where u.charttime >= c.index_time + w.response_start_h * interval '1 hour'
                      and u.charttime < c.index_time + w.response_end_h * interval '1 hour'
                  ) as response_uo_ml,
                  count(*) filter (
                    where u.charttime >= c.index_time + w.response_start_h * interval '1 hour'
                      and u.charttime < c.index_time + w.response_end_h * interval '1 hour'
                  ) as response_uo_n
                from formal_index c
                cross join window_strategies w
                left join mimiciv_derived.urine_output u
                  on c.stay_id = u.stay_id
                 and u.charttime >= c.index_time
                 and u.charttime < c.index_time + interval '24 hours'
                 and u.urineoutput between 0 and 5000
                group by c.stay_id, w.window_strategy
                """
            )
            cur.execute("create index on mimic_tw_urine(stay_id, window_strategy)")

            log("MIMIC-IV: materializing time-window cohort...")
            cur.execute(
                """
                create temp table mimic_time_window_cohort as
                select
                  'MIMIC-IV'::text as database,
                  w.window_strategy,
                  c.subject_id::text as patient_id,
                  c.hadm_id::text as encounter_id,
                  c.stay_id::text as stay_id,
                  c.first_careunit::text as cluster_id,
                  c.index_time::text as index_time,
                  c.landmark_time::text as landmark_time,
                  w.baseline_start_h,
                  w.baseline_end_h,
                  w.response_start_h,
                  w.response_end_h,
                  w.baseline_hours,
                  w.response_hours,
                  c.anchor_age::double precision as age,
                  case when c.gender = 'M' then 1.0 else 0.0 end as male,
                  wt.weight::double precision as weight_kg,
                  s.sofa::double precision as severity_primary,
                  s2.sofa2_total::double precision as severity_secondary,
                  s.respiration::double precision as severity_respiration,
                  s.cardiovascular::double precision as severity_cardiovascular,
                  s.renal::double precision as severity_renal,
                  null::double precision as heart_surgery,
                  case
                    when c.hospital_expire_flag = 1
                     and (c.deathtime is null or c.deathtime >= c.landmark_time)
                    then 1 else 0
                  end as hospital_mortality_after_landmark,
                  mp.baseline_map_mean,
                  mp.response_map_mean,
                  mp.response_map_mean - mp.baseline_map_mean as delta_map,
                  mp.baseline_map_n,
                  mp.response_map_n,
                  vb.baseline_vaso_burden,
                  vb.response_vaso_burden,
                  vb.baseline_vaso_burden - vb.response_vaso_burden as delta_vaso_burden_reduction,
                  log(1 + vb.baseline_vaso_burden) - log(1 + vb.response_vaso_burden) as log_vaso_burden_reduction,
                  vb.baseline_vaso_observed_units,
                  vb.response_vaso_observed_units,
                  u.baseline_uo_ml,
                  u.response_uo_ml,
                  u.baseline_uo_ml / nullif(wt.weight, 0) / w.baseline_hours as baseline_uo_ml_kg_h,
                  u.response_uo_ml / nullif(wt.weight, 0) / w.response_hours as response_uo_ml_kg_h,
                  (u.response_uo_ml / nullif(wt.weight, 0) / w.response_hours)
                    - (u.baseline_uo_ml / nullif(wt.weight, 0) / w.baseline_hours) as delta_uo_ml_kg_h,
                  log(1 + u.response_uo_ml / nullif(wt.weight, 0) / w.response_hours)
                    - log(1 + u.baseline_uo_ml / nullif(wt.weight, 0) / w.baseline_hours) as log_uo_recovery,
                  u.baseline_uo_n,
                  u.response_uo_n
                from formal_index c
                cross join window_strategies w
                left join mimiciv_derived.first_day_weight wt
                  on c.stay_id = wt.stay_id
                left join mimiciv_derived.first_day_sofa s
                  on c.stay_id = s.stay_id
                left join mimiciv_derived.first_day_sofa2 s2
                  on c.stay_id = s2.stay_id
                left join mimic_tw_map mp
                  on c.stay_id = mp.stay_id and w.window_strategy = mp.window_strategy
                left join mimic_tw_vaso vb
                  on c.stay_id = vb.stay_id and w.window_strategy = vb.window_strategy
                left join mimic_tw_urine u
                  on c.stay_id = u.stay_id and w.window_strategy = u.window_strategy
                """
            )
            copy_query_to_csv(
                cur,
                "select * from mimic_time_window_cohort order by stay_id, window_strategy",
                out_path,
            )
            log(
                "MIMIC-IV: wrote "
                f"{fetch_scalar(cur, 'select count(*) from mimic_time_window_cohort')} rows to {out_path}"
            )


def extract_eicu(outdir: Path, server: str) -> None:
    out_path = outdir / "eicu_time_window_cohort.csv"
    active_expr = eicu_formal.active_agent_expression("pi")
    log("eICU: creating source cohort...")
    with connect_pgadmin_server(server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '90min'")
            cur.execute("SET work_mem = '512MB'")
            eicu_formal.create_source_tables(cur)
            create_window_strategies(cur)

            log("eICU: aggregating MAP by alternative windows...")
            cur.execute(
                """
                create temp table eicu_tw_map as
                select
                  c.patientunitstayid,
                  w.window_strategy,
                  avg(coalesce(v.ibp_mean, v.nibp_mean)::double precision) filter (
                    where v.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and v.chartoffset < c.index_offset + w.baseline_end_h * 60
                  ) as baseline_map_mean,
                  count(*) filter (
                    where v.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and v.chartoffset < c.index_offset + w.baseline_end_h * 60
                      and coalesce(v.ibp_mean, v.nibp_mean) is not null
                  ) as baseline_map_n,
                  avg(coalesce(v.ibp_mean, v.nibp_mean)::double precision) filter (
                    where v.chartoffset >= c.index_offset + w.response_start_h * 60
                      and v.chartoffset < c.index_offset + w.response_end_h * 60
                  ) as response_map_mean,
                  count(*) filter (
                    where v.chartoffset >= c.index_offset + w.response_start_h * 60
                      and v.chartoffset < c.index_offset + w.response_end_h * 60
                      and coalesce(v.ibp_mean, v.nibp_mean) is not null
                  ) as response_map_n
                from eicu_formal_index c
                cross join window_strategies w
                left join eicuii.pivoted_vital v
                  on c.patientunitstayid = v.patientunitstayid
                 and v.chartoffset >= c.index_offset
                 and v.chartoffset < c.index_offset + 1440
                 and coalesce(v.ibp_mean, v.nibp_mean) between 20 and 200
                group by c.patientunitstayid, w.window_strategy
                """
            )
            cur.execute("create index on eicu_tw_map(patientunitstayid, window_strategy)")

            log("eICU: aggregating vasoactive burden by alternative windows...")
            cur.execute(
                f"""
                create temp table eicu_tw_vaso as
                with burden as (
                  select
                    c.patientunitstayid,
                    pi.chartoffset,
                    ({active_expr})::double precision as active_agent_count
                  from eicu_formal_index c
                  join eicuii.pivoted_infusion pi
                    on c.patientunitstayid = pi.patientunitstayid
                   and pi.chartoffset >= c.index_offset
                   and pi.chartoffset < c.index_offset + 1440
                )
                select
                  c.patientunitstayid,
                  w.window_strategy,
                  avg(b.active_agent_count) filter (
                    where b.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and b.chartoffset < c.index_offset + w.baseline_end_h * 60
                  ) as baseline_vaso_burden,
                  count(*) filter (
                    where b.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and b.chartoffset < c.index_offset + w.baseline_end_h * 60
                  ) as baseline_vaso_observed_units,
                  avg(b.active_agent_count) filter (
                    where b.chartoffset >= c.index_offset + w.response_start_h * 60
                      and b.chartoffset < c.index_offset + w.response_end_h * 60
                  ) as response_vaso_burden,
                  count(*) filter (
                    where b.chartoffset >= c.index_offset + w.response_start_h * 60
                      and b.chartoffset < c.index_offset + w.response_end_h * 60
                  ) as response_vaso_observed_units
                from eicu_formal_index c
                cross join window_strategies w
                left join burden b
                  on c.patientunitstayid = b.patientunitstayid
                group by c.patientunitstayid, w.window_strategy
                """
            )
            cur.execute("create index on eicu_tw_vaso(patientunitstayid, window_strategy)")

            log("eICU: aggregating urine output by alternative windows...")
            cur.execute(
                """
                create temp table eicu_tw_urine as
                select
                  c.patientunitstayid,
                  w.window_strategy,
                  sum(u.urineoutput::double precision) filter (
                    where u.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and u.chartoffset < c.index_offset + w.baseline_end_h * 60
                  ) as baseline_uo_ml,
                  count(*) filter (
                    where u.chartoffset >= c.index_offset + w.baseline_start_h * 60
                      and u.chartoffset < c.index_offset + w.baseline_end_h * 60
                      and u.urineoutput is not null
                  ) as baseline_uo_n,
                  sum(u.urineoutput::double precision) filter (
                    where u.chartoffset >= c.index_offset + w.response_start_h * 60
                      and u.chartoffset < c.index_offset + w.response_end_h * 60
                  ) as response_uo_ml,
                  count(*) filter (
                    where u.chartoffset >= c.index_offset + w.response_start_h * 60
                      and u.chartoffset < c.index_offset + w.response_end_h * 60
                      and u.urineoutput is not null
                  ) as response_uo_n
                from eicu_formal_index c
                cross join window_strategies w
                left join eicuii.pivoted_uo u
                  on c.patientunitstayid = u.patientunitstayid
                 and u.chartoffset >= c.index_offset
                 and u.chartoffset < c.index_offset + 1440
                 and u.urineoutput between 0 and 5000
                group by c.patientunitstayid, w.window_strategy
                """
            )
            cur.execute("create index on eicu_tw_urine(patientunitstayid, window_strategy)")

            log("eICU: materializing time-window cohort...")
            cur.execute(
                """
                create temp table eicu_time_window_cohort as
                select
                  'eICU'::text as database,
                  w.window_strategy,
                  coalesce(nullif(c.uniquepid, ''), c.patientunitstayid::text) as patient_id,
                  c.patienthealthsystemstayid::text as encounter_id,
                  c.patientunitstayid::text as stay_id,
                  c.hospitalid::text as cluster_id,
                  c.index_offset::text as index_time,
                  c.landmark_offset::text as landmark_time,
                  w.baseline_start_h,
                  w.baseline_end_h,
                  w.response_start_h,
                  w.response_end_h,
                  w.baseline_hours,
                  w.response_hours,
                  c.age_int::double precision as age,
                  case when c.gender = 'Male' then 1.0 else 0.0 end as male,
                  wt.weight_kg::double precision as weight_kg,
                  ap.acutephysiologyscore::double precision as severity_primary,
                  ap.apachescore::double precision as severity_secondary,
                  null::double precision as severity_respiration,
                  null::double precision as severity_cardiovascular,
                  null::double precision as severity_renal,
                  null::double precision as heart_surgery,
                  c.hospital_mortality as hospital_mortality_after_landmark,
                  mp.baseline_map_mean,
                  mp.response_map_mean,
                  mp.response_map_mean - mp.baseline_map_mean as delta_map,
                  mp.baseline_map_n,
                  mp.response_map_n,
                  vb.baseline_vaso_burden,
                  vb.response_vaso_burden,
                  vb.baseline_vaso_burden - vb.response_vaso_burden as delta_vaso_burden_reduction,
                  log(1 + vb.baseline_vaso_burden) - log(1 + vb.response_vaso_burden) as log_vaso_burden_reduction,
                  vb.baseline_vaso_observed_units,
                  vb.response_vaso_observed_units,
                  u.baseline_uo_ml,
                  u.response_uo_ml,
                  u.baseline_uo_ml / nullif(wt.weight_kg, 0) / w.baseline_hours as baseline_uo_ml_kg_h,
                  u.response_uo_ml / nullif(wt.weight_kg, 0) / w.response_hours as response_uo_ml_kg_h,
                  (u.response_uo_ml / nullif(wt.weight_kg, 0) / w.response_hours)
                    - (u.baseline_uo_ml / nullif(wt.weight_kg, 0) / w.baseline_hours) as delta_uo_ml_kg_h,
                  log(1 + u.response_uo_ml / nullif(wt.weight_kg, 0) / w.response_hours)
                    - log(1 + u.baseline_uo_ml / nullif(wt.weight_kg, 0) / w.baseline_hours) as log_uo_recovery,
                  u.baseline_uo_n,
                  u.response_uo_n
                from eicu_formal_index c
                cross join window_strategies w
                left join eicu_weight wt
                  on c.patientunitstayid = wt.patientunitstayid
                left join eicu_apache ap
                  on c.patientunitstayid = ap.patientunitstayid
                left join eicu_tw_map mp
                  on c.patientunitstayid = mp.patientunitstayid and w.window_strategy = mp.window_strategy
                left join eicu_tw_vaso vb
                  on c.patientunitstayid = vb.patientunitstayid and w.window_strategy = vb.window_strategy
                left join eicu_tw_urine u
                  on c.patientunitstayid = u.patientunitstayid and w.window_strategy = u.window_strategy
                """
            )
            copy_query_to_csv(
                cur,
                "select * from eicu_time_window_cohort order by stay_id, window_strategy",
                out_path,
            )
            log(
                "eICU: wrote "
                f"{fetch_scalar(cur, 'select count(*) from eicu_time_window_cohort')} rows to {out_path}"
            )


def extract_sicdb(outdir: Path) -> None:
    out_path = outdir / "sicdb_time_window_cohort.csv"
    log("SICdb: creating source cohort...")
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '120min'")
            cur.execute("SET work_mem = '768MB'")
            sicdb_formal.create_source_tables(cur)
            create_window_strategies(cur)

            log("SICdb: aggregating vasoactive burden by alternative windows...")
            cur.execute(
                """
                create temp table sicdb_tw_vaso as
                with pieces as (
                  select
                    c.caseid,
                    w.window_strategy,
                    w.baseline_hours as denominator_baseline_hours,
                    w.response_hours as denominator_response_hours,
                    v.drugid,
                    greatest(
                      0,
                      least(v.end_offset_sec, c.index_offset_sec + w.baseline_end_h * 3600)
                      - greatest(v.start_offset_sec, c.index_offset_sec + w.baseline_start_h * 3600)
                    ) / 3600.0 as baseline_hours,
                    greatest(
                      0,
                      least(v.end_offset_sec, c.index_offset_sec + w.response_end_h * 3600)
                      - greatest(v.start_offset_sec, c.index_offset_sec + w.response_start_h * 3600)
                    ) / 3600.0 as response_hours
                  from sicdb_formal_index c
                  cross join window_strategies w
                  left join sicdb_vaso_events v
                    on c.caseid = v.caseid
                   and v.end_offset_sec > c.index_offset_sec
                   and v.start_offset_sec < c.index_offset_sec + 86400
                )
                select
                  caseid,
                  window_strategy,
                  coalesce(sum(baseline_hours), 0) / max(denominator_baseline_hours) as baseline_vaso_burden,
                  coalesce(sum(response_hours), 0) / max(denominator_response_hours) as response_vaso_burden,
                  count(distinct drugid) filter (where baseline_hours > 0) as baseline_vaso_observed_units,
                  count(distinct drugid) filter (where response_hours > 0) as response_vaso_observed_units
                from pieces
                group by caseid, window_strategy
                """
            )
            cur.execute("create index on sicdb_tw_vaso(caseid, window_strategy)")

            log("SICdb: scanning MAP/urine signal table once for 0-24h...")
            cur.execute(
                """
                create temp table sicdb_signal_subset as
                select
                  s.caseid,
                  s.dataid,
                  nullif(s.offset, '')::numeric as offset_sec,
                  nullif(s.val, '')::numeric as val_num
                from sicdb_formal_index c
                join lateral (
                  select s.caseid, s.dataid, s.offset, s.val
                  from raw.data_float_h s
                  where s.caseid = c.caseid
                    and s.dataid = any(%s)
                    and s.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and s.val ~ '^-?[0-9]+(\\.[0-9]+)?$'
                    and nullif(s.offset, '')::numeric >= c.index_offset_sec
                    and nullif(s.offset, '')::numeric < c.index_offset_sec + 86400
                ) s on true
                """,
                (sicdb_formal.MAP_IDS + sicdb_formal.URINE_IDS,),
            )
            cur.execute("create index on sicdb_signal_subset(caseid, dataid, offset_sec)")
            cur.execute("analyze sicdb_signal_subset")
            log(f"SICdb: signal rows kept: {fetch_scalar(cur, 'select count(*) from sicdb_signal_subset')}")

            log("SICdb: aggregating MAP by alternative windows...")
            cur.execute(
                """
                create temp table sicdb_tw_map as
                select
                  c.caseid,
                  w.window_strategy,
                  avg(s.val_num) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.baseline_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.baseline_end_h * 3600
                      and s.val_num between 20 and 200
                  ) as baseline_map_mean,
                  count(*) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.baseline_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.baseline_end_h * 3600
                      and s.val_num between 20 and 200
                  ) as baseline_map_n,
                  avg(s.val_num) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.response_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.response_end_h * 3600
                      and s.val_num between 20 and 200
                  ) as response_map_mean,
                  count(*) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.response_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.response_end_h * 3600
                      and s.val_num between 20 and 200
                  ) as response_map_n
                from sicdb_formal_index c
                cross join window_strategies w
                left join sicdb_signal_subset s
                  on c.caseid = s.caseid
                group by c.caseid, w.window_strategy
                """,
                (sicdb_formal.MAP_IDS, sicdb_formal.MAP_IDS, sicdb_formal.MAP_IDS, sicdb_formal.MAP_IDS),
            )
            cur.execute("create index on sicdb_tw_map(caseid, window_strategy)")

            log("SICdb: aggregating urine output by alternative windows...")
            cur.execute(
                """
                create temp table sicdb_tw_urine as
                select
                  c.caseid,
                  w.window_strategy,
                  sum(s.val_num) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.baseline_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.baseline_end_h * 3600
                      and s.val_num between 0 and 5000
                  ) as baseline_uo_ml,
                  count(*) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.baseline_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.baseline_end_h * 3600
                      and s.val_num between 0 and 5000
                  ) as baseline_uo_n,
                  sum(s.val_num) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.response_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.response_end_h * 3600
                      and s.val_num between 0 and 5000
                  ) as response_uo_ml,
                  count(*) filter (
                    where s.dataid = any(%s)
                      and s.offset_sec >= c.index_offset_sec + w.response_start_h * 3600
                      and s.offset_sec < c.index_offset_sec + w.response_end_h * 3600
                      and s.val_num between 0 and 5000
                  ) as response_uo_n
                from sicdb_formal_index c
                cross join window_strategies w
                left join sicdb_signal_subset s
                  on c.caseid = s.caseid
                group by c.caseid, w.window_strategy
                """,
                (sicdb_formal.URINE_IDS, sicdb_formal.URINE_IDS, sicdb_formal.URINE_IDS, sicdb_formal.URINE_IDS),
            )
            cur.execute("create index on sicdb_tw_urine(caseid, window_strategy)")

            log("SICdb: materializing time-window cohort...")
            cur.execute(
                """
                create temp table sicdb_time_window_cohort as
                select
                  'SICdb'::text as database,
                  w.window_strategy,
                  c.patientid::text as patient_id,
                  c.caseid::text as encounter_id,
                  c.caseid::text as stay_id,
                  c.hospitalunit::text as cluster_id,
                  c.index_offset_sec::text as index_time,
                  c.landmark_offset_sec::text as landmark_time,
                  w.baseline_start_h,
                  w.baseline_end_h,
                  w.response_start_h,
                  w.response_end_h,
                  w.baseline_hours,
                  w.response_hours,
                  c.age::double precision as age,
                  c.male::double precision as male,
                  c.weight_kg::double precision as weight_kg,
                  c.saps3::double precision as severity_primary,
                  null::double precision as severity_secondary,
                  null::double precision as severity_respiration,
                  null::double precision as severity_cardiovascular,
                  null::double precision as severity_renal,
                  case when c.heartsurgeryadditionaldata = '740' then 1.0 else 0.0 end as heart_surgery,
                  case
                    when c.hospital_mortality = 1
                     and (c.offset_of_death_sec is null or c.offset_of_death_sec >= c.landmark_offset_sec)
                    then 1 else 0
                  end as hospital_mortality_after_landmark,
                  mp.baseline_map_mean,
                  mp.response_map_mean,
                  mp.response_map_mean - mp.baseline_map_mean as delta_map,
                  mp.baseline_map_n,
                  mp.response_map_n,
                  vb.baseline_vaso_burden,
                  vb.response_vaso_burden,
                  vb.baseline_vaso_burden - vb.response_vaso_burden as delta_vaso_burden_reduction,
                  log(1 + vb.baseline_vaso_burden) - log(1 + vb.response_vaso_burden) as log_vaso_burden_reduction,
                  vb.baseline_vaso_observed_units,
                  vb.response_vaso_observed_units,
                  u.baseline_uo_ml,
                  u.response_uo_ml,
                  u.baseline_uo_ml / nullif(c.weight_kg, 0) / w.baseline_hours as baseline_uo_ml_kg_h,
                  u.response_uo_ml / nullif(c.weight_kg, 0) / w.response_hours as response_uo_ml_kg_h,
                  (u.response_uo_ml / nullif(c.weight_kg, 0) / w.response_hours)
                    - (u.baseline_uo_ml / nullif(c.weight_kg, 0) / w.baseline_hours) as delta_uo_ml_kg_h,
                  log(1 + u.response_uo_ml / nullif(c.weight_kg, 0) / w.response_hours)
                    - log(1 + u.baseline_uo_ml / nullif(c.weight_kg, 0) / w.baseline_hours) as log_uo_recovery,
                  u.baseline_uo_n,
                  u.response_uo_n
                from sicdb_formal_index c
                cross join window_strategies w
                left join sicdb_tw_map mp
                  on c.caseid = mp.caseid and w.window_strategy = mp.window_strategy
                left join sicdb_tw_vaso vb
                  on c.caseid = vb.caseid and w.window_strategy = vb.window_strategy
                left join sicdb_tw_urine u
                  on c.caseid = u.caseid and w.window_strategy = u.window_strategy
                """
            )
            copy_query_to_csv(
                cur,
                "select * from sicdb_time_window_cohort order by stay_id, window_strategy",
                out_path,
            )
            log(
                "SICdb: wrote "
                f"{fetch_scalar(cur, 'select count(*) from sicdb_time_window_cohort')} rows to {out_path}"
            )


def write_summary(outdir: Path) -> None:
    summary_rows: list[dict[str, Any]] = []
    for database, file_name in [
        ("MIMIC-IV", "mimic_time_window_cohort.csv"),
        ("eICU", "eicu_time_window_cohort.csv"),
        ("SICdb", "sicdb_time_window_cohort.csv"),
    ]:
        path = outdir / file_name
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        for window in sorted({r["window_strategy"] for r in rows}):
            x = [r for r in rows if r["window_strategy"] == window]
            complete = [
                r for r in x
                if r["delta_map"] not in ("", "NA")
                and r["log_vaso_burden_reduction"] not in ("", "NA")
                and r["log_uo_recovery"] not in ("", "NA")
            ]
            deaths = sum(1 for r in x if r["hospital_mortality_after_landmark"] == "1")
            summary_rows.append(
                {
                    "database": database,
                    "window_strategy": window,
                    "rows": len(x),
                    "core_complete_rows": len(complete),
                    "hospital_deaths": deaths,
                    "death_rate": deaths / len(x) if x else None,
                }
            )
    write_csv(
        outdir / "time_window_extraction_summary.csv",
        summary_rows,
        ["database", "window_strategy", "rows", "core_complete_rows", "hospital_deaths", "death_rate"],
    )
    md = [
        "# Time-window sensitivity extraction",
        "",
        "Same formal landmark-eligible cohorts were reused. Only baseline/response windows were changed.",
        "",
        "database | window_strategy | rows | core_complete_rows | hospital_deaths | death_rate",
        "--- | --- | ---: | ---: | ---: | ---:",
    ]
    for row in summary_rows:
        md.append(
            f"{row['database']} | {row['window_strategy']} | {row['rows']} | "
            f"{row['core_complete_rows']} | {row['hospital_deaths']} | {row['death_rate']:.3f}"
        )
    (outdir / "time_window_extraction_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract alternative-window HRC components across MIMIC/eICU/SICdb.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--mimic-server", default="mimiciv")
    parser.add_argument("--eicu-server", default="eicu")
    parser.add_argument(
        "--database",
        choices=["all", "mimic", "eicu", "sicdb"],
        default="all",
        help="Run one database or all databases.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.database in ("all", "mimic"):
        extract_mimic(args.outdir, args.mimic_server)
    if args.database in ("all", "eicu"):
        extract_eicu(args.outdir, args.eicu_server)
    if args.database in ("all", "sicdb"):
        extract_sicdb(args.outdir)
    write_summary(args.outdir)
    log(f"Extraction summary: {args.outdir / 'time_window_extraction_summary.md'}")


if __name__ == "__main__":
    main()
