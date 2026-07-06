from __future__ import annotations

from pathlib import Path
import sys

from psycopg.rows import dict_row

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pgadmin_conn import connect_pgadmin_server
from sicdb_conn import connect_sicdb
import importlib


OUTDIR = Path("outputs/landmark_sensitivity")


def copy_query_to_csv(cur, query: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with cur.copy(f"COPY ({query}) TO STDOUT WITH CSV HEADER") as copy:
        with out_path.open("wb") as f:
            for data in copy:
                f.write(data)


def extract_mimic() -> None:
    mimic = importlib.import_module("03_mimic_formal_cohort_extract")
    with connect_pgadmin_server("mimiciv") as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating MIMIC source tables...", flush=True)
            mimic.create_source_tables(cur)
            query = """
            select
              'MIMIC-IV' as database,
              e.subject_id::text as patient_key,
              e.hadm_id::text as admission_key,
              e.stay_id::text as stay_key,
              e.anchor_age as age,
              case when e.gender = 'M' then 1 else 0 end as male,
              w.weight as weight_kg,
              s.sofa as severity_primary,
              s.cardiovascular as severity_cardiovascular,
              s.renal as severity_renal,
              extract(epoch from (e.index_time - e.intime)) / 3600.0 as index_hour_from_icu,
              e.los * 24.0 as icu_los_hours,
              case
                when e.outtime >= e.index_time + interval '24 hours'
                 and (e.deathtime is null or e.deathtime >= e.index_time + interval '24 hours')
                then 1 else 0
              end as landmark_eligible,
              case
                when e.deathtime is not null
                 and e.deathtime < e.index_time + interval '24 hours'
                then 1 else 0
              end as early_death_before_landmark,
              case
                when e.outtime < e.index_time + interval '24 hours'
                 and (e.deathtime is null or e.deathtime >= e.index_time + interval '24 hours')
                then 1 else 0
              end as early_icu_discharge_alive_before_landmark,
              e.hospital_expire_flag as hospital_mortality_anytime,
              row_number() over (partition by e.subject_id order by e.index_time, e.stay_id) as patient_early_support_rank
            from early_vaso_support e
            left join mimiciv_derived.first_day_weight w
              on e.stay_id = w.stay_id
            left join mimiciv_derived.first_day_sofa s
              on e.stay_id = s.stay_id
            """
            copy_query_to_csv(cur, query, OUTDIR / "mimic_landmark_source.csv")


def extract_eicu() -> None:
    eicu = importlib.import_module("07_eicu_formal_cohort_extract")
    with connect_pgadmin_server("eicu") as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating eICU source tables...", flush=True)
            eicu.create_source_tables(cur)
            query = """
            select
              'eICU' as database,
              coalesce(nullif(e.uniquepid, ''), e.patientunitstayid::text) as patient_key,
              e.patienthealthsystemstayid::text as admission_key,
              e.patientunitstayid::text as stay_key,
              e.hospitalid::text as center_key,
              e.age_int as age,
              case when e.gender = 'Male' then 1 when e.gender = 'Female' then 0 else null end as male,
              nullif(e.admissionweight, 0) as weight_kg,
              ap.acutephysiologyscore as severity_primary,
              ap.apachescore as severity_secondary,
              e.index_offset / 60.0 as index_hour_from_icu,
              e.unitdischargeoffset / 60.0 as icu_los_hours,
              case when e.unitdischargeoffset >= e.landmark_offset then 1 else 0 end as landmark_eligible,
              case
                when e.unitdischargeoffset < e.landmark_offset
                 and e.unitdischargestatus = 'Expired'
                then 1 else 0
              end as early_death_before_landmark,
              case
                when e.unitdischargeoffset < e.landmark_offset
                 and e.unitdischargestatus <> 'Expired'
                then 1 else 0
              end as early_icu_discharge_alive_before_landmark,
              e.hospital_mortality as hospital_mortality_anytime,
              row_number() over (
                partition by coalesce(nullif(e.uniquepid, ''), e.patientunitstayid::text)
                order by e.index_offset, e.patientunitstayid
              ) as patient_early_support_rank
            from eicu_early_support e
            left join eicu_apache ap
              on e.patientunitstayid = ap.patientunitstayid
            """
            copy_query_to_csv(cur, query, OUTDIR / "eicu_landmark_source.csv")


def extract_sicdb() -> None:
    sicdb = importlib.import_module("10_sicdb_formal_cohort_extract")
    with connect_sicdb() as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            print("Creating SICdb source tables...", flush=True)
            sicdb.create_source_tables(cur)
            query = """
            select
              'SICdb' as database,
              e.patientid::text as patient_key,
              e.caseid::text as stay_key,
              e.hospitalunit::text as center_key,
              e.age,
              e.male,
              e.weight_kg,
              e.saps3 as severity_primary,
              case when e.heartsurgeryadditionaldata = '740' then 1 else 0 end as heart_surgery,
              (e.index_offset_sec - e.icu_offset_sec) / 3600.0 as index_hour_from_icu,
              e.time_of_stay_sec / 3600.0 as icu_los_hours,
              case
                when e.icu_out_offset_sec >= e.landmark_offset_sec
                 and (e.offset_of_death_sec is null or e.offset_of_death_sec >= e.landmark_offset_sec)
                then 1 else 0
              end as landmark_eligible,
              case
                when e.offset_of_death_sec is not null
                 and e.offset_of_death_sec < e.landmark_offset_sec
                then 1 else 0
              end as early_death_before_landmark,
              case
                when e.icu_out_offset_sec < e.landmark_offset_sec
                 and (e.offset_of_death_sec is null or e.offset_of_death_sec >= e.landmark_offset_sec)
                then 1 else 0
              end as early_icu_discharge_alive_before_landmark,
              e.hospital_mortality as hospital_mortality_anytime,
              row_number() over (
                partition by e.patientid
                order by e.index_offset_sec, e.caseid
              ) as patient_early_support_rank
            from sicdb_early_support e
            """
            copy_query_to_csv(cur, query, OUTDIR / "sicdb_landmark_source.csv")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    extract_mimic()
    extract_eicu()
    extract_sicdb()
    print(f"Landmark source outputs: {OUTDIR}", flush=True)


if __name__ == "__main__":
    main()
