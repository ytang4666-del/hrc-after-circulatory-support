from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

from sicdb_conn import connect_sicdb


DEFAULT_OUTDIR = Path("outputs/sicdb_formal")
NUMERIC_RE_SQL = "^-?[0-9]+(\\.[0-9]+)?$"

VASO_DRUGS = {
    "norepinephrine": ["1562"],
    "epinephrine": ["1502"],
    "vasopressin": ["1550"],
    "dobutamine": ["1559"],
    "milrinone": ["1560"],
    "dopamine": ["1618"],
    "phenylephrine": ["1593"],
}
VASO_DRUG_IDS = [drug_id for ids in VASO_DRUGS.values() for drug_id in ids]

MAP_IDS = ["703", "706"]
URINE_IDS = ["725"]
LACTATE_IDS = ["454", "657", "465"]
CREATININE_IDS = ["367", "368", "369", "339"]


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


def create_source_tables(cur) -> None:
    cur.execute("SET statement_timeout = '90min'")
    cur.execute("SET work_mem = '512MB'")

    cur.execute(
        """
        create temp table sicdb_adult_valid as
        select
          c.caseid,
          c.patientid,
          nullif(c.icuoffset, '')::numeric as icu_offset_sec,
          nullif(c.timeofstay, '')::numeric as time_of_stay_sec,
          nullif(c.offsetafterfirstadmission, '')::numeric as offset_after_first_admission_sec,
          nullif(c.offsetofdeath, '')::numeric as offset_of_death_sec,
          nullif(c.ageonadmission, '')::numeric as age,
          case
            when lower(sex_ref.referencevalue) in ('male', 'mann', 'm') then 1
            when lower(sex_ref.referencevalue) in ('female', 'frau', 'f') then 0
            else null
          end as male,
          case
            when nullif(c.weightonadmission, '')::numeric > 300
            then nullif(c.weightonadmission, '')::numeric / 1000.0
            else nullif(c.weightonadmission, '')::numeric
          end as weight_kg,
          nullif(c.saps3, '')::numeric as saps3,
          c.hospitaldischargetype,
          hosp_ref.referencevalue as hospitaldischargetype_label,
          c.dischargestate,
          discharge_ref.referencevalue as dischargestate_label,
          c.hospitalunit,
          unit_ref.referencevalue as hospitalunit_label,
          c.surgicalsite,
          surg_ref.referencevalue as surgicalsite_label,
          c.heartsurgeryadditionaldata,
          nullif(c.heartsurgerycpbtime, '')::numeric as heart_surgery_cpb_min,
          nullif(c.heartsurgeryendoffset, '')::numeric as heart_surgery_end_sec,
          case when c.hospitaldischargetype = '3130' then 1 else 0 end as hospital_mortality
        from raw.cases c
        left join raw.d_references sex_ref on c.sex = sex_ref.referenceglobalid
        left join raw.d_references hosp_ref on c.hospitaldischargetype = hosp_ref.referenceglobalid
        left join raw.d_references discharge_ref on c.dischargestate = discharge_ref.referenceglobalid
        left join raw.d_references unit_ref on c.hospitalunit = unit_ref.referenceglobalid
        left join raw.d_references surg_ref on c.surgicalsite = surg_ref.referenceglobalid
        where c.ageonadmission ~ '^-?[0-9]+(\\.[0-9]+)?$'
          and nullif(c.ageonadmission, '')::numeric >= 18
          and c.icuoffset ~ '^-?[0-9]+(\\.[0-9]+)?$'
          and c.timeofstay ~ '^-?[0-9]+(\\.[0-9]+)?$'
          and nullif(c.timeofstay, '')::numeric > 0
        """
    )
    cur.execute("create index on sicdb_adult_valid(caseid)")
    cur.execute("create index on sicdb_adult_valid(patientid)")

    cur.execute(
        """
        create temp table sicdb_vaso_events as
        select
          m.caseid,
          m.patientid,
          m.drugid,
          ref.referencevalue as drug_label,
          nullif(m.offset, '')::numeric as start_offset_sec,
          case
            when m.offsetdrugend ~ '^-?[0-9]+(\\.[0-9]+)?$'
             and nullif(m.offsetdrugend, '')::numeric > nullif(m.offset, '')::numeric
            then nullif(m.offsetdrugend, '')::numeric
            else nullif(m.offset, '')::numeric + 3600
          end as end_offset_sec
        from raw.medication m
        left join raw.d_references ref on m.drugid = ref.referenceglobalid
        where m.drugid = any(%s)
          and m.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
        """,
        (VASO_DRUG_IDS,),
    )
    cur.execute("create index on sicdb_vaso_events(caseid)")
    cur.execute("create index on sicdb_vaso_events(caseid, start_offset_sec)")

    cur.execute(
        """
        create temp table sicdb_support_start as
        select
          c.caseid,
          min(v.start_offset_sec) as index_offset_sec
        from sicdb_adult_valid c
        join sicdb_vaso_events v
          on c.caseid = v.caseid
         and v.start_offset_sec >= c.icu_offset_sec
         and v.start_offset_sec < c.icu_offset_sec + 86400
        group by c.caseid
        """
    )
    cur.execute("create index on sicdb_support_start(caseid)")

    cur.execute(
        """
        create temp table sicdb_early_support as
        select
          c.*,
          s.index_offset_sec,
          s.index_offset_sec + 86400 as landmark_offset_sec,
          c.icu_offset_sec + c.time_of_stay_sec as icu_out_offset_sec
        from sicdb_adult_valid c
        join sicdb_support_start s using (caseid)
        """
    )
    cur.execute("create index on sicdb_early_support(caseid)")
    cur.execute("create index on sicdb_early_support(patientid)")

    cur.execute(
        """
        create temp table sicdb_landmark_eligible as
        select *
        from sicdb_early_support
        where icu_out_offset_sec >= landmark_offset_sec
          and (offset_of_death_sec is null or offset_of_death_sec >= landmark_offset_sec)
        """
    )
    cur.execute("create index on sicdb_landmark_eligible(caseid)")
    cur.execute("create index on sicdb_landmark_eligible(patientid)")

    cur.execute(
        """
        create temp table sicdb_formal_index as
        select *
        from (
          select
            l.*,
            row_number() over (
              partition by l.patientid
              order by l.index_offset_sec, l.caseid
            ) as patient_support_rank
          from sicdb_landmark_eligible l
        ) x
        where patient_support_rank = 1
        """
    )
    cur.execute("create index on sicdb_formal_index(caseid)")
    cur.execute("create index on sicdb_formal_index(index_offset_sec)")
    log(
        "  source counts: "
        f"adult={fetch_scalar(cur, 'select count(*) from sicdb_adult_valid')}, "
        f"early_support={fetch_scalar(cur, 'select count(*) from sicdb_early_support')}, "
        f"landmark={fetch_scalar(cur, 'select count(*) from sicdb_landmark_eligible')}, "
        f"formal={fetch_scalar(cur, 'select count(*) from sicdb_formal_index')}"
    )


def ensure_performance_indexes() -> None:
    """Create narrow persistent indexes needed for SICdb HRC extraction.

    SICdb was imported without indexes. The signal table is large, so these are
    partial indexes limited to HRC variables rather than broad raw-table indexes.
    """
    statements = [
        (
            "cases caseid",
            "create index if not exists idx_cases_hrc_caseid on raw.cases (caseid)",
        ),
        (
            "cases patientid",
            "create index if not exists idx_cases_hrc_patientid on raw.cases (patientid)",
        ),
        (
            "vaso medication case-offset",
            f"""
            create index if not exists idx_medication_hrc_vaso_case_offset
            on raw.medication (caseid, ((nullif("offset", '')::numeric)), drugid)
            where drugid in ({','.join("'" + x + "'" for x in VASO_DRUG_IDS)})
              and "offset" ~ '{NUMERIC_RE_SQL}'
            """,
        ),
        (
            "vaso medication patient",
            f"""
            create index if not exists idx_medication_hrc_vaso_patient
            on raw.medication (patientid)
            where drugid in ({','.join("'" + x + "'" for x in VASO_DRUG_IDS)})
              and "offset" ~ '{NUMERIC_RE_SQL}'
            """,
        ),
        (
            "MAP/urine signal case-offset",
            f"""
            create index if not exists idx_data_float_h_hrc_signal_case_offset
            on raw.data_float_h (caseid, ((nullif("offset", '')::numeric)), dataid)
            where dataid in ({','.join("'" + x + "'" for x in MAP_IDS + URINE_IDS)})
              and "offset" ~ '{NUMERIC_RE_SQL}'
              and val ~ '{NUMERIC_RE_SQL}'
            """,
        ),
        (
            "laboratory case-offset",
            f"""
            create index if not exists idx_laboratory_hrc_case_offset
            on raw.laboratory (caseid, ((nullif("offset", '')::numeric)), laboratoryid)
            where laboratoryid in ({','.join("'" + x + "'" for x in LACTATE_IDS + CREATININE_IDS)})
              and "offset" ~ '{NUMERIC_RE_SQL}'
              and laboratoryvalue ~ '{NUMERIC_RE_SQL}'
            """,
        ),
    ]

    with connect_sicdb() as conn:
        conn.autocommit = True
        with conn.cursor() as cur:
            cur.execute("SET statement_timeout = '6h'")
            cur.execute("SET maintenance_work_mem = '1GB'")
            for label, sql in statements:
                log(f"Creating/checking SICdb HRC index: {label}...")
                cur.execute(sql)
            for table in ("raw.cases", "raw.medication", "raw.data_float_h", "raw.laboratory"):
                log(f"Analyzing {table}...")
                cur.execute(f"analyze {table}")


def create_component_tables(cur) -> None:
    cur.execute(
        """
        create temp table sicdb_vaso_burden as
        with pieces as (
          select
            c.caseid,
            v.drugid,
            greatest(
              0,
              extract(epoch from (
                to_timestamp(least(v.end_offset_sec, c.index_offset_sec + 21600))
                - to_timestamp(greatest(v.start_offset_sec, c.index_offset_sec))
              )) / 3600.0
            ) as baseline_hours,
            greatest(
              0,
              extract(epoch from (
                to_timestamp(least(v.end_offset_sec, c.index_offset_sec + 86400))
                - to_timestamp(greatest(v.start_offset_sec, c.index_offset_sec + 21600))
              )) / 3600.0
            ) as response_hours
          from sicdb_formal_index c
          left join sicdb_vaso_events v
            on c.caseid = v.caseid
           and v.end_offset_sec > c.index_offset_sec
           and v.start_offset_sec < c.index_offset_sec + 86400
        )
        select
          caseid,
          coalesce(sum(baseline_hours), 0) / 6.0 as baseline_vaso_burden,
          coalesce(sum(response_hours), 0) / 18.0 as response_vaso_burden,
          count(distinct drugid) filter (where baseline_hours > 0) as baseline_vaso_agents_n,
          count(distinct drugid) filter (where response_hours > 0) as response_vaso_agents_n
        from pieces
        group by caseid
        """
    )
    cur.execute("create index on sicdb_vaso_burden(caseid)")

    log("  aggregating vasoactive burden...")
    log("  scanning SICdb signal table for MAP/urine subset...")
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
        (MAP_IDS + URINE_IDS,),
    )
    cur.execute("create index on sicdb_signal_subset(caseid, dataid, offset_sec)")
    cur.execute("analyze sicdb_signal_subset")
    log(f"  signal rows kept: {fetch_scalar(cur, 'select count(*) from sicdb_signal_subset')}")

    log("  aggregating MAP windows...")
    cur.execute(
        """
        create temp table sicdb_map as
        select
          c.caseid,
          avg(s.val_num) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec
              and s.offset_sec < c.index_offset_sec + 21600
              and s.val_num between 20 and 200
          ) as baseline_map_mean,
          count(*) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec
              and s.offset_sec < c.index_offset_sec + 21600
              and s.val_num between 20 and 200
          ) as baseline_map_n,
          avg(s.val_num) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec + 21600
              and s.offset_sec < c.index_offset_sec + 86400
              and s.val_num between 20 and 200
          ) as response_map_mean,
          count(*) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec + 21600
              and s.offset_sec < c.index_offset_sec + 86400
              and s.val_num between 20 and 200
          ) as response_map_n
        from sicdb_formal_index c
        left join sicdb_signal_subset s on c.caseid = s.caseid
        group by c.caseid
        """,
        (MAP_IDS, MAP_IDS, MAP_IDS, MAP_IDS),
    )
    cur.execute("create index on sicdb_map(caseid)")

    log("  aggregating urine windows...")
    cur.execute(
        """
        create temp table sicdb_urine as
        select
          c.caseid,
          sum(s.val_num) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec
              and s.offset_sec < c.index_offset_sec + 21600
              and s.val_num between 0 and 5000
          ) as baseline_uo_ml,
          count(*) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec
              and s.offset_sec < c.index_offset_sec + 21600
              and s.val_num between 0 and 5000
          ) as baseline_uo_n,
          sum(s.val_num) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec + 21600
              and s.offset_sec < c.index_offset_sec + 86400
              and s.val_num between 0 and 5000
          ) as response_uo_ml,
          count(*) filter (
            where s.dataid = any(%s)
              and s.offset_sec >= c.index_offset_sec + 21600
              and s.offset_sec < c.index_offset_sec + 86400
              and s.val_num between 0 and 5000
          ) as response_uo_n
        from sicdb_formal_index c
        left join sicdb_signal_subset s on c.caseid = s.caseid
        group by c.caseid
        """,
        (URINE_IDS, URINE_IDS, URINE_IDS, URINE_IDS),
    )
    cur.execute("create index on sicdb_urine(caseid)")

    log("  aggregating creatinine/lactate windows...")
    cur.execute(
        """
        create temp table sicdb_labs as
        select
          c.caseid,
          avg(nullif(l.laboratoryvalue, '')::numeric) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 21600
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as baseline_lactate_mean,
          count(*) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 21600
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as baseline_lactate_n,
          avg(nullif(l.laboratoryvalue, '')::numeric) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec + 21600
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 86400
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as response_lactate_mean,
          count(*) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec + 21600
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 86400
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as response_lactate_n,
          avg(nullif(l.laboratoryvalue, '')::numeric) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 21600
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as baseline_creatinine_mean,
          count(*) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 21600
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as baseline_creatinine_n,
          avg(nullif(l.laboratoryvalue, '')::numeric) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec + 21600
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 86400
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as response_creatinine_mean,
          count(*) filter (
            where l.laboratoryid = any(%s)
              and nullif(l.offset, '')::numeric >= c.index_offset_sec + 21600
              and nullif(l.offset, '')::numeric < c.index_offset_sec + 86400
              and nullif(l.laboratoryvalue, '')::numeric between 0.1 and 30
          ) as response_creatinine_n
        from sicdb_formal_index c
        left join raw.laboratory l
          on c.caseid = l.caseid
         and l.offset ~ '^-?[0-9]+(\\.[0-9]+)?$'
         and l.laboratoryvalue ~ '^-?[0-9]+(\\.[0-9]+)?$'
         and nullif(l.offset, '')::numeric >= c.index_offset_sec
         and nullif(l.offset, '')::numeric < c.index_offset_sec + 86400
         and l.laboratoryid = any(%s)
        group by c.caseid
        """,
        (
            LACTATE_IDS, LACTATE_IDS, LACTATE_IDS, LACTATE_IDS,
            CREATININE_IDS, CREATININE_IDS, CREATININE_IDS, CREATININE_IDS,
            LACTATE_IDS + CREATININE_IDS,
        ),
    )
    cur.execute("create index on sicdb_labs(caseid)")


FINAL_QUERY = """
select
  c.caseid,
  c.patientid,
  c.icu_offset_sec,
  c.time_of_stay_sec,
  c.icu_out_offset_sec,
  c.index_offset_sec,
  c.landmark_offset_sec,
  c.age,
  c.male,
  c.weight_kg,
  c.saps3,
  c.hospitaldischargetype,
  c.hospitaldischargetype_label,
  c.dischargestate,
  c.dischargestate_label,
  c.hospitalunit,
  c.hospitalunit_label,
  c.surgicalsite,
  c.surgicalsite_label,
  c.heartsurgeryadditionaldata,
  c.heart_surgery_cpb_min,
  c.heart_surgery_end_sec,
  case
    when c.hospital_mortality = 1
     and (c.offset_of_death_sec is null or c.offset_of_death_sec >= c.landmark_offset_sec)
    then 1 else 0
  end as hospital_mortality_after_landmark,
  m.baseline_map_mean,
  m.response_map_mean,
  m.response_map_mean - m.baseline_map_mean as delta_map,
  m.baseline_map_n,
  m.response_map_n,
  vb.baseline_vaso_burden,
  vb.response_vaso_burden,
  vb.baseline_vaso_burden - vb.response_vaso_burden as delta_vaso_burden_reduction,
  log(1 + vb.baseline_vaso_burden) - log(1 + vb.response_vaso_burden) as log_vaso_burden_reduction,
  vb.baseline_vaso_agents_n,
  vb.response_vaso_agents_n,
  u.baseline_uo_ml,
  u.response_uo_ml,
  u.baseline_uo_ml / nullif(c.weight_kg, 0) / 6.0 as baseline_uo_ml_kg_h,
  u.response_uo_ml / nullif(c.weight_kg, 0) / 18.0 as response_uo_ml_kg_h,
  (u.response_uo_ml / nullif(c.weight_kg, 0) / 18.0)
    - (u.baseline_uo_ml / nullif(c.weight_kg, 0) / 6.0) as delta_uo_ml_kg_h,
  log(1 + u.response_uo_ml / nullif(c.weight_kg, 0) / 18.0)
    - log(1 + u.baseline_uo_ml / nullif(c.weight_kg, 0) / 6.0) as log_uo_recovery,
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
from sicdb_formal_index c
left join sicdb_map m using (caseid)
left join sicdb_vaso_burden vb using (caseid)
left join sicdb_urine u using (caseid)
left join sicdb_labs l using (caseid)
"""


def write_audits(cur, outdir: Path) -> None:
    metrics = {
        "adult_valid_icu_cases": "select count(*) from sicdb_adult_valid",
        "early_vasoactive_support_cases": "select count(*) from sicdb_early_support",
        "landmark_eligible_cases": "select count(*) from sicdb_landmark_eligible",
        "formal_first_patient_support_cases": "select count(*) from sicdb_formal_index",
        "formal_dataset_rows": "select count(*) from sicdb_formal_cohort",
        "core_hrc_complete_rows": (
            "select count(*) from sicdb_formal_cohort "
            "where delta_map is not null and log_vaso_burden_reduction is not null and log_uo_recovery is not null"
        ),
        "map_baseline_response_complete": (
            "select count(*) from sicdb_formal_cohort "
            "where baseline_map_mean is not null and response_map_mean is not null"
        ),
        "vaso_burden_baseline_response_complete": (
            "select count(*) from sicdb_formal_cohort "
            "where baseline_vaso_burden is not null and response_vaso_burden is not null"
        ),
        "urine_weight_baseline_response_complete": (
            "select count(*) from sicdb_formal_cohort "
            "where weight_kg is not null and baseline_uo_ml_kg_h is not null and response_uo_ml_kg_h is not null"
        ),
        "creatinine_baseline_response_complete": (
            "select count(*) from sicdb_formal_cohort "
            "where baseline_creatinine_mean is not null and response_creatinine_mean is not null"
        ),
        "lactate_baseline_response_complete": (
            "select count(*) from sicdb_formal_cohort "
            "where baseline_lactate_mean is not null and response_lactate_mean is not null"
        ),
        "hospital_mortality_after_landmark_events": (
            "select count(*) from sicdb_formal_cohort where hospital_mortality_after_landmark = 1"
        ),
        "heart_surgery_rows": (
            "select count(*) from sicdb_formal_cohort where heartsurgeryadditionaldata = '740'"
        ),
    }
    rows = [{"metric": metric, "value": fetch_scalar(cur, sql)} for metric, sql in metrics.items()]
    write_csv(outdir / "sicdb_formal_cohort_summary.csv", rows, ["metric", "value"])

    vaso_audit = cur.execute(
        """
        select
          v.drugid,
          v.drug_label,
          count(*) as n_records,
          count(distinct v.caseid) as n_cases
        from sicdb_vaso_events v
        group by v.drugid, v.drug_label
        order by n_cases desc
        """
    ).fetchall()
    write_csv(outdir / "sicdb_vaso_drug_audit.csv", vaso_audit, ["drugid", "drug_label", "n_records", "n_cases"])

    md = [
        "# SICdb formal cohort extraction",
        "",
        "## Cohort definition",
        "",
        "- Adult valid ICU case: age >=18 years, valid ICU offset and positive ICU length of stay.",
        "- Main support cohort: first patient-level ICU case with vasoactive medication initiation within 24 hours of ICU admission.",
        "- Index time: first vasoactive medication offset in the first ICU day.",
        "- Baseline window: 0-6 hours after index.",
        "- Response window: 6-24 hours after index.",
        "- Landmark eligibility: alive and still in ICU at 24 hours after index.",
        "- MAP and urine output were derived from `raw.data_float_h`; vasoactive burden was derived from medication intervals.",
        "",
        "## Counts",
        "",
    ]
    md.extend(f"- {row['metric']}: {row['value']}" for row in rows)
    (outdir / "sicdb_formal_cohort_summary.md").write_text("\n".join(md), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract SICdb HRC external validation cohort.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--create-indexes",
        action="store_true",
        help="Create narrow persistent indexes for SICdb HRC variables before extraction.",
    )
    args = parser.parse_args()

    if args.create_indexes:
        ensure_performance_indexes()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            log("Creating SICdb source cohort...")
            create_source_tables(cur)
            log("Aggregating SICdb HRC windows...")
            create_component_tables(cur)
            log("Materializing SICdb cohort...")
            cur.execute(f"create temp table sicdb_formal_cohort as {FINAL_QUERY}")
            cur.execute("create index on sicdb_formal_cohort(caseid)")
            log("Writing outputs...")
            copy_query_to_csv(
                cur,
                "select * from sicdb_formal_cohort order by patientid, index_offset_sec, caseid",
                args.outdir / "sicdb_formal_cohort.csv",
            )
            write_audits(cur, args.outdir)

    log(f"Formal cohort: {args.outdir / 'sicdb_formal_cohort.csv'}")
    log(f"Summary: {args.outdir / 'sicdb_formal_cohort_summary.md'}")


if __name__ == "__main__":
    main()
