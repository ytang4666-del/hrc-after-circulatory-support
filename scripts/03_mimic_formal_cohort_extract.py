from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

from pgadmin_conn import connect_pgadmin_server


DEFAULT_OUTDIR = Path("outputs/mimic_formal")

CS_ICD_CODES = {
    9: ["78551", "99801"],
    10: ["R570", "T8111", "T8111XA", "T8111XD", "T8111XS"],
}

MCS_ITEMIDS = {
    "iabp": [224272],
    "impella": [228169],
    "ecmo": [229529, 229530],
}


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


def create_source_tables(cur) -> None:
    cur.execute("SET statement_timeout = '45min'")
    cur.execute("SET work_mem = '256MB'")

    cur.execute(
        """
        create temp table adult_valid_icu as
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
    cur.execute("create index on adult_valid_icu(stay_id)")
    cur.execute("create index on adult_valid_icu(hadm_id)")
    cur.execute("create index on adult_valid_icu(subject_id)")

    cur.execute(
        """
        create temp table early_vaso_support as
        select
          i.*,
          min(v.starttime) as index_time
        from adult_valid_icu i
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
    cur.execute("create index on early_vaso_support(stay_id)")
    cur.execute("create index on early_vaso_support(hadm_id)")
    cur.execute("create index on early_vaso_support(subject_id)")

    cur.execute(
        """
        create temp table landmark_eligible_stays as
        select
          q.*,
          q.index_time + interval '24 hours' as landmark_time
        from early_vaso_support q
        where q.outtime >= q.index_time + interval '24 hours'
          and (q.deathtime is null or q.deathtime >= q.index_time + interval '24 hours')
        """
    )
    cur.execute("create index on landmark_eligible_stays(stay_id)")
    cur.execute("create index on landmark_eligible_stays(hadm_id)")
    cur.execute("create index on landmark_eligible_stays(subject_id)")

    cur.execute(
        """
        create temp table formal_index as
        select *
        from (
          select
            l.*,
            row_number() over (partition by l.subject_id order by l.index_time, l.stay_id) as patient_support_rank
          from landmark_eligible_stays l
        ) x
        where patient_support_rank = 1
        """
    )
    cur.execute("create index on formal_index(stay_id)")
    cur.execute("create index on formal_index(hadm_id)")
    cur.execute("create index on formal_index(subject_id)")
    cur.execute("create index on formal_index(index_time)")

    all_mcs_itemids = [itemid for values in MCS_ITEMIDS.values() for itemid in values]
    cur.execute(
        """
        create temp table mcs_events as
        select
          pe.subject_id,
          pe.hadm_id,
          pe.stay_id,
          pe.starttime,
          pe.endtime,
          pe.itemid,
          d.label,
          case
            when pe.itemid = any(%s) then 'IABP'
            when pe.itemid = any(%s) then 'Impella'
            when pe.itemid = any(%s) then 'ECMO'
            else 'Other'
          end as mcs_type
        from mimiciv_icu.procedureevents pe
        left join mimiciv_icu.d_items d
          on pe.itemid = d.itemid
        where pe.itemid = any(%s)
        """,
        (
            MCS_ITEMIDS["iabp"],
            MCS_ITEMIDS["impella"],
            MCS_ITEMIDS["ecmo"],
            all_mcs_itemids,
        ),
    )
    cur.execute("create index on mcs_events(stay_id)")
    cur.execute("create index on mcs_events(hadm_id)")

    cs_codes_flat = [(version, code) for version, codes in CS_ICD_CODES.items() for code in codes]
    values_sql = ",".join(
        f"({version}, '{code}')" for version, code in cs_codes_flat
    )
    cur.execute(
        f"""
        create temp table cs_diagnoses as
        with cs_codes(icd_version, icd_code_clean) as (
          values {values_sql}
        )
        select
          d.subject_id,
          d.hadm_id,
          max(case when c.icd_code_clean is not null then 1 else 0 end) as cardiogenic_shock_icd,
          string_agg(distinct d.icd_version::text || ':' || regexp_replace(d.icd_code, '\\s+', '', 'g'), ';') as cs_icd_codes
        from mimiciv_hosp.diagnoses_icd d
        join cs_codes c
          on d.icd_version = c.icd_version
         and regexp_replace(d.icd_code, '\\s+', '', 'g') = c.icd_code_clean
        group by d.subject_id, d.hadm_id
        """
    )
    cur.execute("create index on cs_diagnoses(hadm_id)")


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
        from formal_index c
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
            va.ne_equiv,
            greatest(
              0,
              extract(epoch from (
                least(va.overlap_end, c.index_time + interval '6 hours')
                - greatest(va.overlap_start, c.index_time)
              )) / 3600.0
            ) as baseline_hours,
            greatest(
              0,
              extract(epoch from (
                least(va.overlap_end, c.index_time + interval '24 hours')
                - greatest(va.overlap_start, c.index_time + interval '6 hours')
              )) / 3600.0
            ) as response_hours
          from formal_index c
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
        from formal_index c
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
        from formal_index c
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
        from formal_index c
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
  s.coagulation as sofa_coagulation,
  s.liver as sofa_liver,
  s.cardiovascular as sofa_cardiovascular,
  s.cns as sofa_cns,
  s.renal as sofa_renal,
  s2.sofa2_total,
  c.hospital_expire_flag,
  case
    when c.hospital_expire_flag = 1
     and (c.deathtime is null or c.deathtime >= c.landmark_time)
    then 1 else 0
  end as hospital_mortality_after_landmark,
  coalesce(cs.cardiogenic_shock_icd, 0) as cardiogenic_shock_icd,
  cs.cs_icd_codes,
  coalesce(mc.mcs_any_icu, 0) as mcs_any_icu,
  coalesce(mc.mcs_any_0_24h, 0) as mcs_any_0_24h,
  coalesce(mc.mcs_pre_index, 0) as mcs_pre_index,
  coalesce(mc.iabp_any_icu, 0) as iabp_any_icu,
  coalesce(mc.impella_any_icu, 0) as impella_any_icu,
  coalesce(mc.ecmo_any_icu, 0) as ecmo_any_icu,
  coalesce(mc.iabp_0_24h, 0) as iabp_0_24h,
  coalesce(mc.impella_0_24h, 0) as impella_0_24h,
  coalesce(mc.ecmo_0_24h, 0) as ecmo_0_24h,
  mc.first_mcs_time,
  mc.first_mcs_type,
  mc.mcs_types_icu,
  case
    when coalesce(cs.cardiogenic_shock_icd, 0) = 1 or coalesce(mc.mcs_any_icu, 0) = 1
    then 1 else 0
  end as cs_or_mcs_subgroup,
  mp.baseline_map_mean,
  mp.response_map_mean,
  mp.response_map_mean - mp.baseline_map_mean as delta_map,
  mp.baseline_map_n,
  mp.response_map_n,
  v.baseline_neq_avg,
  v.response_neq_avg,
  v.baseline_neq_avg - v.response_neq_avg as delta_neq_reduction,
  log(1 + v.baseline_neq_avg) - log(1 + v.response_neq_avg) as log_neq_reduction,
  v.baseline_vaso_observed_hours,
  v.response_vaso_observed_hours,
  u.baseline_uo_ml,
  u.response_uo_ml,
  u.baseline_uo_ml / nullif(w.weight, 0) / 6.0 as baseline_uo_ml_kg_h,
  u.response_uo_ml / nullif(w.weight, 0) / 18.0 as response_uo_ml_kg_h,
  (u.response_uo_ml / nullif(w.weight, 0) / 18.0)
    - (u.baseline_uo_ml / nullif(w.weight, 0) / 6.0) as delta_uo_ml_kg_h,
  log(1 + u.response_uo_ml / nullif(w.weight, 0) / 18.0)
    - log(1 + u.baseline_uo_ml / nullif(w.weight, 0) / 6.0) as log_uo_recovery,
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
from formal_index c
left join mimiciv_derived.first_day_weight w
  on c.stay_id = w.stay_id
left join mimiciv_derived.first_day_sofa s
  on c.stay_id = s.stay_id
left join mimiciv_derived.first_day_sofa2 s2
  on c.stay_id = s2.stay_id
left join cs_diagnoses cs
  on c.hadm_id = cs.hadm_id
left join (
  select
    c.stay_id,
    max(case when e.stay_id is not null then 1 else 0 end) as mcs_any_icu,
    max(case when e.starttime >= c.index_time and e.starttime < c.landmark_time then 1 else 0 end) as mcs_any_0_24h,
    max(case when e.starttime < c.index_time then 1 else 0 end) as mcs_pre_index,
    max(case when e.mcs_type = 'IABP' then 1 else 0 end) as iabp_any_icu,
    max(case when e.mcs_type = 'Impella' then 1 else 0 end) as impella_any_icu,
    max(case when e.mcs_type = 'ECMO' then 1 else 0 end) as ecmo_any_icu,
    max(case when e.mcs_type = 'IABP' and e.starttime >= c.index_time and e.starttime < c.landmark_time then 1 else 0 end) as iabp_0_24h,
    max(case when e.mcs_type = 'Impella' and e.starttime >= c.index_time and e.starttime < c.landmark_time then 1 else 0 end) as impella_0_24h,
    max(case when e.mcs_type = 'ECMO' and e.starttime >= c.index_time and e.starttime < c.landmark_time then 1 else 0 end) as ecmo_0_24h,
    min(e.starttime) as first_mcs_time,
    (array_agg(e.mcs_type order by e.starttime, e.itemid))[1] as first_mcs_type,
    string_agg(distinct e.mcs_type, ';' order by e.mcs_type) as mcs_types_icu
  from formal_index c
  left join mcs_events e
    on c.stay_id = e.stay_id
   and e.starttime >= c.intime
   and e.starttime < c.outtime
  group by c.stay_id
) mc
  on c.stay_id = mc.stay_id
left join hrc_map mp
  on c.stay_id = mp.stay_id
left join hrc_vaso v
  on c.stay_id = v.stay_id
left join hrc_urine u
  on c.stay_id = u.stay_id
left join hrc_creatinine cr
  on c.stay_id = cr.stay_id
left join hrc_lactate la
  on c.stay_id = la.stay_id
"""


def audit_mcs_procedureevents(cur) -> list[dict[str, Any]]:
    return cur.execute(
        """
        select
          e.mcs_type,
          e.itemid,
          e.label,
          count(*) as n_records,
          count(distinct e.stay_id) as n_stays,
          min(e.starttime) as first_record_time,
          max(e.starttime) as last_record_time
        from mcs_events e
        group by e.mcs_type, e.itemid, e.label
        order by e.mcs_type, e.itemid
        """
    ).fetchall()


def audit_cs_icd(cur) -> list[dict[str, Any]]:
    return cur.execute(
        """
        select
          d.icd_version,
          regexp_replace(d.icd_code, '\\s+', '', 'g') as icd_code,
          dx.long_title,
          count(*) as n_diagnoses,
          count(distinct d.hadm_id) as n_hadm
        from mimiciv_hosp.diagnoses_icd d
        join mimiciv_hosp.d_icd_diagnoses dx
          on d.icd_version = dx.icd_version
         and d.icd_code = dx.icd_code
        where (d.icd_version = 9 and regexp_replace(d.icd_code, '\\s+', '', 'g') in ('78551', '99801'))
           or (d.icd_version = 10 and regexp_replace(d.icd_code, '\\s+', '', 'g') in ('R570', 'T8111', 'T8111XA', 'T8111XD', 'T8111XS'))
        group by d.icd_version, regexp_replace(d.icd_code, '\\s+', '', 'g'), dx.long_title
        order by d.icd_version, icd_code
        """
    ).fetchall()


def write_summary(cur, outdir: Path) -> None:
    metrics = {
        "adult_valid_icu_stays": "select count(*) from adult_valid_icu",
        "early_vaso_support_stays": "select count(*) from early_vaso_support",
        "landmark_eligible_stays": "select count(*) from landmark_eligible_stays",
        "formal_first_patient_support_stays": "select count(*) from formal_index",
        "formal_dataset_rows": "select count(*) from formal_cohort",
        "core_hrc_complete_rows": (
            "select count(*) from formal_cohort "
            "where delta_map is not null and log_neq_reduction is not null and log_uo_recovery is not null"
        ),
        "cardiogenic_shock_icd_rows": "select count(*) from formal_cohort where cardiogenic_shock_icd = 1",
        "mcs_any_icu_rows": "select count(*) from formal_cohort where mcs_any_icu = 1",
        "mcs_any_0_24h_rows": "select count(*) from formal_cohort where mcs_any_0_24h = 1",
        "cs_or_mcs_subgroup_rows": "select count(*) from formal_cohort where cs_or_mcs_subgroup = 1",
        "iabp_any_icu_rows": "select count(*) from formal_cohort where iabp_any_icu = 1",
        "impella_any_icu_rows": "select count(*) from formal_cohort where impella_any_icu = 1",
        "ecmo_any_icu_rows": "select count(*) from formal_cohort where ecmo_any_icu = 1",
        "hospital_mortality_after_landmark_events": (
            "select count(*) from formal_cohort where hospital_mortality_after_landmark = 1"
        ),
    }
    rows = [{"metric": metric, "value": fetch_scalar(cur, sql)} for metric, sql in metrics.items()]
    write_csv(outdir / "mimic_formal_cohort_summary.csv", rows, ["metric", "value"])

    md = [
        "# MIMIC formal cohort extraction",
        "",
        "## Cohort definition",
        "",
        "- Adult valid ICU stay: age >=18 years and ICU outtime after intime.",
        "- Main support cohort: first patient-level ICU stay with vasoactive-agent initiation within 24 hours of ICU admission.",
        "- Index time: first vasoactive-agent start time within the first ICU day.",
        "- Baseline window: 0-6 hours after index.",
        "- Response window: 6-24 hours after index.",
        "- Landmark eligibility: alive and still in ICU at 24 hours after index.",
        "- Mechanistic subgroup: cardiogenic shock ICD diagnosis or ICU MCS procedureevent evidence.",
        "",
        "## Counts",
        "",
    ]
    md.extend(f"- {row['metric']}: {row['value']}" for row in rows)
    (outdir / "mimic_formal_cohort_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract formal MIMIC HRC cohort with CS/MCS flags.")
    parser.add_argument("--server", default="mimiciv", help="pgAdmin saved server name")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with connect_pgadmin_server(args.server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating formal source cohort and CS/MCS flags...")
            create_source_tables(cur)
            print("Aggregating HRC component windows...")
            create_component_tables(cur)
            print("Materializing formal cohort...")
            cur.execute(f"create temp table formal_cohort as {FINAL_QUERY}")
            cur.execute("create index on formal_cohort(stay_id)")
            cur.execute("create index on formal_cohort(hadm_id)")
            cur.execute("create index on formal_cohort(cardiogenic_shock_icd)")
            cur.execute("create index on formal_cohort(mcs_any_icu)")
            print("Writing outputs...")
            copy_query_to_csv(
                cur,
                "select * from formal_cohort order by subject_id, index_time, stay_id",
                args.outdir / "mimic_formal_cohort.csv",
            )
            write_summary(cur, args.outdir)
            mcs_audit = audit_mcs_procedureevents(cur)
            cs_audit = audit_cs_icd(cur)

    write_csv(
        args.outdir / "mimic_mcs_procedureevents_audit.csv",
        mcs_audit,
        ["mcs_type", "itemid", "label", "n_records", "n_stays", "first_record_time", "last_record_time"],
    )
    write_csv(
        args.outdir / "mimic_cardiogenic_shock_icd_audit.csv",
        cs_audit,
        ["icd_version", "icd_code", "long_title", "n_diagnoses", "n_hadm"],
    )

    print(f"Formal cohort: {args.outdir / 'mimic_formal_cohort.csv'}")
    print(f"Summary: {args.outdir / 'mimic_formal_cohort_summary.md'}")


if __name__ == "__main__":
    main()
