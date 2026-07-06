from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

from pgadmin_conn import connect_pgadmin_server


DEFAULT_OUTDIR = Path("outputs/eicu_formal")

VASO_COLUMNS = [
    "dopamine",
    "dobutamine",
    "norepinephrine",
    "phenylephrine",
    "epinephrine",
    "vasopressin",
    "milrinone",
]


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


def active_agent_expression(alias: str = "pi") -> str:
    parts = [
        f"case when coalesce({alias}.{col}, 0) > 0 then 1 else 0 end"
        for col in VASO_COLUMNS
    ]
    return " + ".join(parts)


def create_source_tables(cur) -> None:
    cur.execute("SET statement_timeout = '45min'")
    cur.execute("SET work_mem = '256MB'")

    cur.execute(
        """
        create temp table eicu_adult_valid as
        select
          p.patientunitstayid,
          p.patienthealthsystemstayid,
          p.uniquepid,
          p.hospitalid,
          p.wardid,
          p.unittype,
          p.unitvisitnumber,
          p.unitstaytype,
          p.age_int,
          p.gender,
          p.ethnicity,
          p.apacheadmissiondx,
          p.admissionweight,
          p.unitdischargeoffset,
          p.unitdischargestatus,
          p.hospitaldischargestatus,
          case when p.hospitaldischargestatus = 'Expired' then 1 else 0 end as hospital_mortality,
          case when p.unitdischargestatus = 'Expired' then 1 else 0 end as icu_mortality
        from eicuii.patient p
        where p.age_int >= 18
          and p.unitdischargeoffset > 0
          and p.hospitaldischargestatus in ('Alive', 'Expired')
        """
    )
    cur.execute("create index on eicu_adult_valid(patientunitstayid)")
    cur.execute("create index on eicu_adult_valid(uniquepid)")

    active_expr = active_agent_expression("pi")
    cur.execute(
        f"""
        create temp table eicu_support_start as
        select
          p.patientunitstayid,
          min(pi.chartoffset) as index_offset
        from eicu_adult_valid p
        join eicuii.pivoted_infusion pi
          on p.patientunitstayid = pi.patientunitstayid
         and pi.chartoffset >= 0
         and pi.chartoffset < 1440
        where ({active_expr}) > 0
        group by p.patientunitstayid
        """
    )
    cur.execute("create index on eicu_support_start(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_early_support as
        select
          p.*,
          s.index_offset,
          s.index_offset + 1440 as landmark_offset
        from eicu_adult_valid p
        join eicu_support_start s
          on p.patientunitstayid = s.patientunitstayid
        """
    )
    cur.execute("create index on eicu_early_support(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_landmark_eligible as
        select *
        from eicu_early_support
        where unitdischargeoffset >= landmark_offset
        """
    )
    cur.execute("create index on eicu_landmark_eligible(patientunitstayid)")
    cur.execute("create index on eicu_landmark_eligible(uniquepid)")

    cur.execute(
        """
        create temp table eicu_formal_index as
        select *
        from (
          select
            l.*,
            row_number() over (
              partition by coalesce(nullif(l.uniquepid, ''), l.patientunitstayid::text)
              order by l.index_offset, l.patientunitstayid
            ) as patient_support_rank
          from eicu_landmark_eligible l
        ) x
        where patient_support_rank = 1
        """
    )
    cur.execute("create index on eicu_formal_index(patientunitstayid)")
    cur.execute("create index on eicu_formal_index(index_offset)")

    cur.execute(
        """
        create temp table eicu_apache as
        select
          patientunitstayid,
          max(acutephysiologyscore) as acutephysiologyscore,
          max(apachescore) as apachescore,
          max(actualventdays) as actualventdays
        from eicuii.apachepatientresult
        group by patientunitstayid
        """
    )
    cur.execute("create index on eicu_apache(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_apsvar as
        select
          patientunitstayid,
          max(meanbp) as apache_meanbp,
          max(creatinine) as apache_creatinine,
          max(urine) as apache_urine,
          max(vent) as apache_vent,
          max(dialysis) as apache_dialysis
        from eicuii.apacheapsvar
        group by patientunitstayid
        """
    )
    cur.execute("create index on eicu_apsvar(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_weight as
        select
          p.patientunitstayid,
          coalesce(
            nullif(p.admissionweight, 0),
            percentile_cont(0.5) within group (order by w.weight)
              filter (where w.weight between 20 and 300)
          )::double precision as weight_kg
        from eicu_formal_index p
        left join eicuii.pivoted_weight w
          on p.patientunitstayid = w.patientunitstayid
         and w.chartoffset >= -360
         and w.chartoffset < 1440
        group by p.patientunitstayid, p.admissionweight
        """
    )
    cur.execute("create index on eicu_weight(patientunitstayid)")


def create_component_tables(cur) -> None:
    cur.execute(
        """
        create temp table eicu_map as
        select
          c.patientunitstayid,
          avg(coalesce(v.ibp_mean, v.nibp_mean)::double precision)
            filter (
              where v.chartoffset >= c.index_offset
                and v.chartoffset < c.index_offset + 360
            ) as baseline_map_mean,
          count(*) filter (
              where v.chartoffset >= c.index_offset
                and v.chartoffset < c.index_offset + 360
                and coalesce(v.ibp_mean, v.nibp_mean) is not null
            ) as baseline_map_n,
          avg(coalesce(v.ibp_mean, v.nibp_mean)::double precision)
            filter (
              where v.chartoffset >= c.index_offset + 360
                and v.chartoffset < c.index_offset + 1440
            ) as response_map_mean,
          count(*) filter (
              where v.chartoffset >= c.index_offset + 360
                and v.chartoffset < c.index_offset + 1440
                and coalesce(v.ibp_mean, v.nibp_mean) is not null
            ) as response_map_n
        from eicu_formal_index c
        left join eicuii.pivoted_vital v
          on c.patientunitstayid = v.patientunitstayid
         and v.chartoffset >= c.index_offset
         and v.chartoffset < c.index_offset + 1440
         and coalesce(v.ibp_mean, v.nibp_mean) between 20 and 200
        group by c.patientunitstayid
        """
    )
    cur.execute("create index on eicu_map(patientunitstayid)")

    active_expr = active_agent_expression("pi")
    cur.execute(
        f"""
        create temp table eicu_vaso_burden as
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
          avg(b.active_agent_count) filter (
            where b.chartoffset >= c.index_offset
              and b.chartoffset < c.index_offset + 360
          ) as baseline_vaso_burden,
          count(*) filter (
            where b.chartoffset >= c.index_offset
              and b.chartoffset < c.index_offset + 360
          ) as baseline_vaso_n,
          avg(b.active_agent_count) filter (
            where b.chartoffset >= c.index_offset + 360
              and b.chartoffset < c.index_offset + 1440
          ) as response_vaso_burden,
          count(*) filter (
            where b.chartoffset >= c.index_offset + 360
              and b.chartoffset < c.index_offset + 1440
          ) as response_vaso_n
        from eicu_formal_index c
        left join burden b
          on c.patientunitstayid = b.patientunitstayid
        group by c.patientunitstayid
        """
    )
    cur.execute("create index on eicu_vaso_burden(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_urine as
        select
          c.patientunitstayid,
          sum(u.urineoutput::double precision) filter (
            where u.chartoffset >= c.index_offset
              and u.chartoffset < c.index_offset + 360
          ) as baseline_uo_ml,
          count(*) filter (
            where u.chartoffset >= c.index_offset
              and u.chartoffset < c.index_offset + 360
              and u.urineoutput is not null
          ) as baseline_uo_n,
          sum(u.urineoutput::double precision) filter (
            where u.chartoffset >= c.index_offset + 360
              and u.chartoffset < c.index_offset + 1440
          ) as response_uo_ml,
          count(*) filter (
            where u.chartoffset >= c.index_offset + 360
              and u.chartoffset < c.index_offset + 1440
              and u.urineoutput is not null
          ) as response_uo_n
        from eicu_formal_index c
        left join eicuii.pivoted_uo u
          on c.patientunitstayid = u.patientunitstayid
         and u.chartoffset >= c.index_offset
         and u.chartoffset < c.index_offset + 1440
         and u.urineoutput between 0 and 5000
        group by c.patientunitstayid
        """
    )
    cur.execute("create index on eicu_urine(patientunitstayid)")

    cur.execute(
        """
        create temp table eicu_labs as
        select
          c.patientunitstayid,
          avg(l.creatinine::double precision) filter (
            where l.chartoffset >= c.index_offset
              and l.chartoffset < c.index_offset + 360
          ) as baseline_creatinine_mean,
          count(*) filter (
            where l.chartoffset >= c.index_offset
              and l.chartoffset < c.index_offset + 360
              and l.creatinine is not null
          ) as baseline_creatinine_n,
          avg(l.creatinine::double precision) filter (
            where l.chartoffset >= c.index_offset + 360
              and l.chartoffset < c.index_offset + 1440
          ) as response_creatinine_mean,
          count(*) filter (
            where l.chartoffset >= c.index_offset + 360
              and l.chartoffset < c.index_offset + 1440
              and l.creatinine is not null
          ) as response_creatinine_n,
          avg(l.lactate::double precision) filter (
            where l.chartoffset >= c.index_offset
              and l.chartoffset < c.index_offset + 360
          ) as baseline_lactate_mean,
          count(*) filter (
            where l.chartoffset >= c.index_offset
              and l.chartoffset < c.index_offset + 360
              and l.lactate is not null
          ) as baseline_lactate_n,
          avg(l.lactate::double precision) filter (
            where l.chartoffset >= c.index_offset + 360
              and l.chartoffset < c.index_offset + 1440
          ) as response_lactate_mean,
          count(*) filter (
            where l.chartoffset >= c.index_offset + 360
              and l.chartoffset < c.index_offset + 1440
              and l.lactate is not null
          ) as response_lactate_n
        from eicu_formal_index c
        left join eicuii.pivoted_lab l
          on c.patientunitstayid = l.patientunitstayid
         and l.chartoffset >= c.index_offset
         and l.chartoffset < c.index_offset + 1440
         and (
              (l.creatinine between 0.1 and 30)
           or (l.lactate between 0.1 and 30)
         )
        group by c.patientunitstayid
        """
    )
    cur.execute("create index on eicu_labs(patientunitstayid)")


FINAL_QUERY = """
select
  c.patientunitstayid,
  c.patienthealthsystemstayid,
  c.uniquepid,
  c.hospitalid,
  c.wardid,
  c.unittype,
  c.unitvisitnumber,
  c.unitstaytype,
  c.age_int as age,
  c.gender,
  c.ethnicity,
  c.apacheadmissiondx,
  c.index_offset,
  c.landmark_offset,
  c.unitdischargeoffset,
  c.unitdischargeoffset / 60.0 as icu_los_hours,
  c.hospitaldischargestatus,
  c.unitdischargestatus,
  c.hospital_mortality as hospital_mortality_after_landmark,
  c.icu_mortality,
  w.weight_kg,
  ap.acutephysiologyscore,
  ap.apachescore,
  ap.actualventdays,
  av.apache_meanbp,
  av.apache_creatinine,
  av.apache_urine,
  av.apache_vent,
  av.apache_dialysis,
  m.baseline_map_mean,
  m.response_map_mean,
  m.response_map_mean - m.baseline_map_mean as delta_map,
  m.baseline_map_n,
  m.response_map_n,
  vb.baseline_vaso_burden,
  vb.response_vaso_burden,
  vb.baseline_vaso_burden - vb.response_vaso_burden as delta_vaso_burden_reduction,
  log(1 + vb.baseline_vaso_burden) - log(1 + vb.response_vaso_burden) as log_vaso_burden_reduction,
  vb.baseline_vaso_n,
  vb.response_vaso_n,
  u.baseline_uo_ml,
  u.response_uo_ml,
  u.baseline_uo_ml / nullif(w.weight_kg, 0) / 6.0 as baseline_uo_ml_kg_h,
  u.response_uo_ml / nullif(w.weight_kg, 0) / 18.0 as response_uo_ml_kg_h,
  (u.response_uo_ml / nullif(w.weight_kg, 0) / 18.0)
    - (u.baseline_uo_ml / nullif(w.weight_kg, 0) / 6.0) as delta_uo_ml_kg_h,
  log(1 + u.response_uo_ml / nullif(w.weight_kg, 0) / 18.0)
    - log(1 + u.baseline_uo_ml / nullif(w.weight_kg, 0) / 6.0) as log_uo_recovery,
  u.baseline_uo_n,
  u.response_uo_n,
  l.baseline_creatinine_mean,
  l.response_creatinine_mean,
  l.baseline_creatinine_mean - l.response_creatinine_mean as delta_creatinine_reduction,
  l.baseline_creatinine_n,
  l.response_creatinine_n,
  l.baseline_lactate_mean,
  l.response_lactate_mean,
  l.baseline_lactate_mean - l.response_lactate_mean as delta_lactate_reduction,
  l.baseline_lactate_n,
  l.response_lactate_n
from eicu_formal_index c
left join eicu_weight w
  on c.patientunitstayid = w.patientunitstayid
left join eicu_apache ap
  on c.patientunitstayid = ap.patientunitstayid
left join eicu_apsvar av
  on c.patientunitstayid = av.patientunitstayid
left join eicu_map m
  on c.patientunitstayid = m.patientunitstayid
left join eicu_vaso_burden vb
  on c.patientunitstayid = vb.patientunitstayid
left join eicu_urine u
  on c.patientunitstayid = u.patientunitstayid
left join eicu_labs l
  on c.patientunitstayid = l.patientunitstayid
"""


def write_summary(cur, outdir: Path) -> None:
    metrics = {
        "adult_valid_icu_stays": "select count(*) from eicu_adult_valid",
        "early_vasoactive_support_stays": "select count(*) from eicu_early_support",
        "landmark_eligible_stays": "select count(*) from eicu_landmark_eligible",
        "formal_first_patient_support_stays": "select count(*) from eicu_formal_index",
        "formal_dataset_rows": "select count(*) from eicu_formal_cohort",
        "core_hrc_complete_rows": (
            "select count(*) from eicu_formal_cohort "
            "where delta_map is not null and log_vaso_burden_reduction is not null and log_uo_recovery is not null"
        ),
        "map_baseline_response_complete": (
            "select count(*) from eicu_formal_cohort "
            "where baseline_map_mean is not null and response_map_mean is not null"
        ),
        "vaso_burden_baseline_response_complete": (
            "select count(*) from eicu_formal_cohort "
            "where baseline_vaso_burden is not null and response_vaso_burden is not null"
        ),
        "urine_weight_baseline_response_complete": (
            "select count(*) from eicu_formal_cohort "
            "where weight_kg is not null and baseline_uo_ml_kg_h is not null and response_uo_ml_kg_h is not null"
        ),
        "creatinine_baseline_response_complete": (
            "select count(*) from eicu_formal_cohort "
            "where baseline_creatinine_mean is not null and response_creatinine_mean is not null"
        ),
        "lactate_baseline_response_complete": (
            "select count(*) from eicu_formal_cohort "
            "where baseline_lactate_mean is not null and response_lactate_mean is not null"
        ),
        "hospital_mortality_after_landmark_events": (
            "select count(*) from eicu_formal_cohort where hospital_mortality_after_landmark = 1"
        ),
    }
    rows = [{"metric": metric, "value": fetch_scalar(cur, sql)} for metric, sql in metrics.items()]
    write_csv(outdir / "eicu_formal_cohort_summary.csv", rows, ["metric", "value"])

    md = [
        "# eICU formal cohort extraction",
        "",
        "## Cohort definition",
        "",
        "- Adult valid ICU stay: age >=18 years, positive ICU discharge offset, known hospital discharge status.",
        "- Main support cohort: first patient-level ICU stay with pivoted vasoactive infusion record within 24 hours of ICU admission.",
        "- Index time: first vasoactive infusion chart offset in the first ICU day.",
        "- Baseline window: 0-6 hours after index.",
        "- Response window: 6-24 hours after index.",
        "- Landmark eligibility: ICU discharge offset at least 24 hours after index.",
        "- Vasopressor domain uses active vasoactive-agent count/log burden because eICU infusion harmonization is not a stable NE-equivalent dose system.",
        "",
        "## Counts",
        "",
    ]
    md.extend(f"- {row['metric']}: {row['value']}" for row in rows)
    (outdir / "eicu_formal_cohort_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract eICU HRC external validation cohort.")
    parser.add_argument("--server", default="eicu", help="pgAdmin saved server name")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with connect_pgadmin_server(args.server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating eICU source cohort...")
            create_source_tables(cur)
            print("Aggregating eICU HRC windows...")
            create_component_tables(cur)
            print("Materializing eICU cohort...")
            cur.execute(f"create temp table eicu_formal_cohort as {FINAL_QUERY}")
            cur.execute("create index on eicu_formal_cohort(patientunitstayid)")
            print("Writing outputs...")
            copy_query_to_csv(
                cur,
                "select * from eicu_formal_cohort order by uniquepid, index_offset, patientunitstayid",
                args.outdir / "eicu_formal_cohort.csv",
            )
            write_summary(cur, args.outdir)

    print(f"Formal cohort: {args.outdir / 'eicu_formal_cohort.csv'}")
    print(f"Summary: {args.outdir / 'eicu_formal_cohort_summary.md'}")


if __name__ == "__main__":
    main()
