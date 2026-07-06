from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

from pgadmin_conn import connect_pgadmin_server


DEFAULT_OUTDIR = Path("outputs/mimic_feasibility")
ANALYSIS_CSV = "mimic_hrc_feasibility_dataset.csv"
SUMMARY_CSV = "mimic_hrc_feasibility_summary.csv"


def fetch_scalar(cur, sql: str, params: tuple[Any, ...] | None = None) -> Any:
    row = cur.execute(sql, params).fetchone()
    if isinstance(row, dict):
        return next(iter(row.values()))
    return row[0]


def copy_query_to_csv(cur, query: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with cur.copy(f"COPY ({query}) TO STDOUT WITH CSV HEADER") as copy:
        with out_path.open("wb") as f:
            for data in copy:
                f.write(data)


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(rows)


def create_cohort_tables(cur) -> None:
    cur.execute("SET statement_timeout = '30min'")
    cur.execute("SET work_mem = '256MB'")

    cur.execute(
        """
        create temp table hrc_adult_icu as
        select
          i.subject_id,
          i.hadm_id,
          i.stay_id,
          i.first_careunit,
          i.last_careunit,
          i.intime,
          i.outtime,
          i.los,
          p.gender,
          p.anchor_age,
          a.deathtime,
          a.hospital_expire_flag
        from mimiciv_icu.icustays i
        join mimiciv_hosp.patients p
          on i.subject_id = p.subject_id
        join mimiciv_hosp.admissions a
          on i.subject_id = a.subject_id
         and i.hadm_id = a.hadm_id
        where p.anchor_age >= 18
          and i.outtime > i.intime
        """
    )
    cur.execute("create index on hrc_adult_icu(stay_id)")
    cur.execute("create index on hrc_adult_icu(subject_id)")

    cur.execute(
        """
        create temp table hrc_qualifying_support as
        select
          i.*,
          min(v.starttime) as index_time
        from hrc_adult_icu i
        join mimiciv_derived.vasoactive_agent v
          on i.stay_id = v.stay_id
         and v.starttime >= i.intime
         and v.starttime < least(i.outtime, i.intime + interval '24 hours')
        group by
          i.subject_id, i.hadm_id, i.stay_id, i.first_careunit, i.last_careunit,
          i.intime, i.outtime, i.los, i.gender, i.anchor_age,
          i.deathtime, i.hospital_expire_flag
        """
    )
    cur.execute("create index on hrc_qualifying_support(stay_id)")
    cur.execute("create index on hrc_qualifying_support(subject_id)")

    cur.execute(
        """
        create temp table hrc_index as
        select *
        from (
          select
            q.*,
            q.index_time + interval '24 hours' as landmark_time,
            row_number() over (partition by q.subject_id order by q.index_time, q.stay_id) as rn
          from hrc_qualifying_support q
          where q.outtime >= q.index_time + interval '24 hours'
            and (q.deathtime is null or q.deathtime >= q.index_time + interval '24 hours')
        ) x
        where rn = 1
        """
    )
    cur.execute("create index on hrc_index(stay_id)")
    cur.execute("create index on hrc_index(hadm_id)")
    cur.execute("create index on hrc_index(index_time)")
    cur.execute("analyze hrc_index")


def create_component_tables(cur) -> None:
    cur.execute(
        """
        create temp table hrc_map as
        select
          c.stay_id,
          avg(coalesce(v.mbp, v.mbp_ni)::double precision)
            filter (
              where v.charttime >= c.index_time
                and v.charttime < c.index_time + interval '6 hours'
            ) as baseline_map_mean,
          count(*) filter (
              where v.charttime >= c.index_time
                and v.charttime < c.index_time + interval '6 hours'
                and coalesce(v.mbp, v.mbp_ni) is not null
            ) as baseline_map_n,
          avg(coalesce(v.mbp, v.mbp_ni)::double precision)
            filter (
              where v.charttime >= c.index_time + interval '6 hours'
                and v.charttime < c.index_time + interval '24 hours'
            ) as response_map_mean,
          count(*) filter (
              where v.charttime >= c.index_time + interval '6 hours'
                and v.charttime < c.index_time + interval '24 hours'
                and coalesce(v.mbp, v.mbp_ni) is not null
            ) as response_map_n
        from hrc_index c
        left join mimiciv_derived.vitalsign v
          on c.stay_id = v.stay_id
         and v.charttime >= c.index_time
         and v.charttime < c.index_time + interval '24 hours'
         and coalesce(v.mbp, v.mbp_ni) between 20 and 200
        group by c.stay_id
        """
    )
    cur.execute("create index on hrc_map(stay_id)")

    cur.execute(
        """
        create temp table hrc_vaso as
        with va as (
          select
            c.stay_id,
            greatest(v.starttime, c.index_time) as overlap_start_0_24,
            least(coalesce(v.endtime, v.starttime), c.index_time + interval '24 hours') as overlap_end_0_24,
            (
              coalesce(v.norepinephrine, 0)
              + coalesce(v.epinephrine, 0)
              + coalesce(v.phenylephrine, 0) / 10.0
              + coalesce(v.dopamine, 0) / 100.0
              + coalesce(v.vasopressin, 0) * 2.5
            )::double precision as ne_equiv
          from hrc_index c
          join mimiciv_derived.vasoactive_agent v
            on c.stay_id = v.stay_id
           and coalesce(v.endtime, v.starttime) > c.index_time
           and v.starttime < c.index_time + interval '24 hours'
        ),
        pieces as (
          select
            c.stay_id,
            va.ne_equiv,
            greatest(
              0,
              extract(epoch from (
                least(va.overlap_end_0_24, c.index_time + interval '6 hours')
                - greatest(va.overlap_start_0_24, c.index_time)
              )) / 3600.0
            ) as baseline_hours,
            greatest(
              0,
              extract(epoch from (
                least(va.overlap_end_0_24, c.index_time + interval '24 hours')
                - greatest(va.overlap_start_0_24, c.index_time + interval '6 hours')
              )) / 3600.0
            ) as response_hours
          from hrc_index c
          left join va
            on c.stay_id = va.stay_id
        )
        select
          stay_id,
          coalesce(sum(ne_equiv * baseline_hours), 0) / 6.0 as baseline_neq_avg,
          coalesce(sum(ne_equiv * response_hours), 0) / 18.0 as response_neq_avg,
          coalesce(sum(baseline_hours), 0) as baseline_vaso_observed_hours,
          coalesce(sum(response_hours), 0) as response_vaso_observed_hours
        from pieces
        group by stay_id
        """
    )
    cur.execute("create index on hrc_vaso(stay_id)")

    cur.execute(
        """
        create temp table hrc_urine as
        select
          c.stay_id,
          sum(u.urineoutput::double precision) filter (
            where u.charttime >= c.index_time
              and u.charttime < c.index_time + interval '6 hours'
          ) as baseline_uo_ml,
          count(*) filter (
            where u.charttime >= c.index_time
              and u.charttime < c.index_time + interval '6 hours'
          ) as baseline_uo_n,
          sum(u.urineoutput::double precision) filter (
            where u.charttime >= c.index_time + interval '6 hours'
              and u.charttime < c.index_time + interval '24 hours'
          ) as response_uo_ml,
          count(*) filter (
            where u.charttime >= c.index_time + interval '6 hours'
              and u.charttime < c.index_time + interval '24 hours'
          ) as response_uo_n
        from hrc_index c
        left join mimiciv_derived.urine_output u
          on c.stay_id = u.stay_id
         and u.charttime >= c.index_time
         and u.charttime < c.index_time + interval '24 hours'
         and u.urineoutput between 0 and 5000
        group by c.stay_id
        """
    )
    cur.execute("create index on hrc_urine(stay_id)")

    cur.execute(
        """
        create temp table hrc_creatinine as
        select
          c.stay_id,
          avg(ch.creatinine::double precision) filter (
            where ch.charttime >= c.index_time
              and ch.charttime < c.index_time + interval '6 hours'
          ) as baseline_creatinine_mean,
          count(*) filter (
            where ch.charttime >= c.index_time
              and ch.charttime < c.index_time + interval '6 hours'
              and ch.creatinine is not null
          ) as baseline_creatinine_n,
          avg(ch.creatinine::double precision) filter (
            where ch.charttime >= c.index_time + interval '6 hours'
              and ch.charttime < c.index_time + interval '24 hours'
          ) as response_creatinine_mean,
          count(*) filter (
            where ch.charttime >= c.index_time + interval '6 hours'
              and ch.charttime < c.index_time + interval '24 hours'
              and ch.creatinine is not null
          ) as response_creatinine_n
        from hrc_index c
        left join mimiciv_derived.chemistry ch
          on c.hadm_id = ch.hadm_id
         and ch.charttime >= c.index_time
         and ch.charttime < c.index_time + interval '24 hours'
         and ch.creatinine between 0.1 and 30
        group by c.stay_id
        """
    )
    cur.execute("create index on hrc_creatinine(stay_id)")

    cur.execute(
        """
        create temp table hrc_lactate as
        select
          c.stay_id,
          avg(bg.lactate::double precision) filter (
            where bg.charttime >= c.index_time
              and bg.charttime < c.index_time + interval '6 hours'
          ) as baseline_lactate_mean,
          count(*) filter (
            where bg.charttime >= c.index_time
              and bg.charttime < c.index_time + interval '6 hours'
              and bg.lactate is not null
          ) as baseline_lactate_n,
          avg(bg.lactate::double precision) filter (
            where bg.charttime >= c.index_time + interval '6 hours'
              and bg.charttime < c.index_time + interval '24 hours'
          ) as response_lactate_mean,
          count(*) filter (
            where bg.charttime >= c.index_time + interval '6 hours'
              and bg.charttime < c.index_time + interval '24 hours'
              and bg.lactate is not null
          ) as response_lactate_n
        from hrc_index c
        left join mimiciv_derived.bg bg
          on c.hadm_id = bg.hadm_id
         and bg.charttime >= c.index_time
         and bg.charttime < c.index_time + interval '24 hours'
         and bg.lactate between 0.1 and 30
        group by c.stay_id
        """
    )
    cur.execute("create index on hrc_lactate(stay_id)")


FINAL_QUERY = """
select
  c.subject_id,
  c.hadm_id,
  c.stay_id,
  c.first_careunit,
  c.last_careunit,
  c.intime as icu_intime,
  c.outtime as icu_outtime,
  c.index_time,
  c.landmark_time,
  c.los as icu_los_days,
  c.anchor_age,
  c.gender,
  w.weight as weight_kg,
  s.sofa as first_day_sofa,
  s.respiration as sofa_respiration,
  s.cardiovascular as sofa_cardiovascular,
  s.renal as sofa_renal,
  s2.sofa2_total,
  c.hospital_expire_flag,
  case
    when c.hospital_expire_flag = 1
     and (c.deathtime is null or c.deathtime >= c.landmark_time)
    then 1 else 0
  end as hospital_mortality_after_landmark,
  m.baseline_map_mean,
  m.response_map_mean,
  m.response_map_mean - m.baseline_map_mean as delta_map,
  m.baseline_map_n,
  m.response_map_n,
  v.baseline_neq_avg,
  v.response_neq_avg,
  v.baseline_neq_avg - v.response_neq_avg as delta_neq_reduction,
  v.baseline_vaso_observed_hours,
  v.response_vaso_observed_hours,
  u.baseline_uo_ml,
  u.response_uo_ml,
  u.baseline_uo_ml / nullif(w.weight, 0) / 6.0 as baseline_uo_ml_kg_h,
  u.response_uo_ml / nullif(w.weight, 0) / 18.0 as response_uo_ml_kg_h,
  (u.response_uo_ml / nullif(w.weight, 0) / 18.0)
    - (u.baseline_uo_ml / nullif(w.weight, 0) / 6.0) as delta_uo_ml_kg_h,
  u.baseline_uo_n,
  u.response_uo_n,
  cr.baseline_creatinine_mean,
  cr.response_creatinine_mean,
  cr.baseline_creatinine_mean - cr.response_creatinine_mean as delta_creatinine_reduction,
  cr.baseline_creatinine_n,
  cr.response_creatinine_n,
  la.baseline_lactate_mean,
  la.response_lactate_mean,
  la.baseline_lactate_mean - la.response_lactate_mean as delta_lactate_reduction,
  la.baseline_lactate_n,
  la.response_lactate_n
from hrc_index c
left join mimiciv_derived.first_day_weight w
  on c.stay_id = w.stay_id
left join mimiciv_derived.first_day_sofa s
  on c.stay_id = s.stay_id
left join mimiciv_derived.first_day_sofa2 s2
  on c.stay_id = s2.stay_id
left join hrc_map m
  on c.stay_id = m.stay_id
left join hrc_vaso v
  on c.stay_id = v.stay_id
left join hrc_urine u
  on c.stay_id = u.stay_id
left join hrc_creatinine cr
  on c.stay_id = cr.stay_id
left join hrc_lactate la
  on c.stay_id = la.stay_id
order by c.subject_id, c.index_time, c.stay_id
"""


def collect_summary(cur) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = {
        "adult_icu_stays": "select count(*) from hrc_adult_icu",
        "adult_icu_stays_with_vasoactive_start_in_first_24h": "select count(*) from hrc_qualifying_support",
        "first_patient_level_landmark_eligible_stays": "select count(*) from hrc_index",
        "analysis_rows": "select count(*) from hrc_final",
        "complete_map_baseline_response": (
            "select count(*) from hrc_final "
            "where baseline_map_mean is not null and response_map_mean is not null"
        ),
        "complete_urine_weight_baseline_response": (
            "select count(*) from hrc_final "
            "where weight_kg is not null "
            "and baseline_uo_ml_kg_h is not null "
            "and response_uo_ml_kg_h is not null"
        ),
        "complete_creatinine_baseline_response": (
            "select count(*) from hrc_final "
            "where baseline_creatinine_mean is not null and response_creatinine_mean is not null"
        ),
        "complete_lactate_baseline_response": (
            "select count(*) from hrc_final "
            "where baseline_lactate_mean is not null and response_lactate_mean is not null"
        ),
        "hospital_mortality_after_landmark_events": (
            "select count(*) from hrc_final "
            "where hospital_mortality_after_landmark = 1"
        ),
    }
    for metric, sql in metrics.items():
        rows.append({"metric": metric, "value": fetch_scalar(cur, sql)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract a MIMIC-IV feasibility dataset for HRC component modeling."
    )
    parser.add_argument("--server", default="mimiciv", help="pgAdmin saved server name")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with connect_pgadmin_server(args.server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating MIMIC HRC cohort tables...")
            create_cohort_tables(cur)
            print("Aggregating HRC component windows...")
            create_component_tables(cur)
            print("Materializing final analysis table...")
            cur.execute(f"create temp table hrc_final as {FINAL_QUERY}")
            cur.execute("create index on hrc_final(stay_id)")
            print("Writing analysis dataset...")
            analysis_path = args.outdir / ANALYSIS_CSV
            copy_query_to_csv(
                cur,
                """
                select *
                from hrc_final
                order by subject_id, index_time, stay_id
                """,
                analysis_path,
            )
            print("Writing summary...")
            summary_rows = collect_summary(cur)
            summary_path = args.outdir / SUMMARY_CSV
            write_summary(summary_path, summary_rows)

    print(f"Analysis dataset: {analysis_path}")
    print(f"Summary: {summary_path}")
    for row in summary_rows:
        print(f"{row['metric']}: {row['value']}")


if __name__ == "__main__":
    main()
