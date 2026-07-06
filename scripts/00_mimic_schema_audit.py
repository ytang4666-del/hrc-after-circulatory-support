from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from psycopg.rows import dict_row

from pgadmin_conn import connect_pgadmin_server


DEFAULT_OUTDIR = Path("outputs/mimic_feasibility")

REQUIRED_TABLES = [
    ("mimiciv_icu", "icustays", "cohort backbone"),
    ("mimiciv_hosp", "admissions", "mortality/admission dates"),
    ("mimiciv_hosp", "patients", "age/sex"),
    ("mimiciv_derived", "vitalsign", "MAP and vital signs"),
    ("mimiciv_derived", "vasoactive_agent", "vasopressor support"),
    ("mimiciv_derived", "urine_output", "renal response"),
    ("mimiciv_derived", "first_day_sofa", "baseline severity"),
    ("mimiciv_derived", "first_day_sofa2", "baseline severity sensitivity"),
    ("mimiciv_derived", "bg", "lactate/metabolic optional"),
    ("mimiciv_derived", "chemistry", "creatinine optional"),
    ("mimiciv_derived", "ventilation", "mechanical ventilation covariate/outcome"),
    ("mimiciv_derived", "rrt", "RRT covariate/outcome"),
    ("mimiciv_icu", "inputevents", "fluid/MCS medication source if needed"),
    ("mimiciv_hosp", "labevents", "raw labs fallback"),
    ("mimiciv_hosp", "d_labitems", "raw lab labels fallback"),
]

REQUIRED_COLUMNS = {
    ("mimiciv_icu", "icustays"): [
        "subject_id",
        "hadm_id",
        "stay_id",
        "intime",
        "outtime",
        "los",
    ],
    ("mimiciv_hosp", "admissions"): [
        "subject_id",
        "hadm_id",
        "admittime",
        "dischtime",
        "deathtime",
        "hospital_expire_flag",
    ],
    ("mimiciv_hosp", "patients"): [
        "subject_id",
        "gender",
        "anchor_age",
        "anchor_year",
    ],
    ("mimiciv_derived", "vitalsign"): [
        "subject_id",
        "stay_id",
        "charttime",
        "mbp",
        "mbp_ni",
        "heart_rate",
    ],
    ("mimiciv_derived", "vasoactive_agent"): [
        "stay_id",
        "starttime",
        "endtime",
        "norepinephrine",
        "epinephrine",
        "dopamine",
        "dobutamine",
        "vasopressin",
        "phenylephrine",
    ],
    ("mimiciv_derived", "urine_output"): [
        "stay_id",
        "charttime",
        "urineoutput",
    ],
    ("mimiciv_derived", "first_day_sofa"): [
        "stay_id",
        "sofa",
    ],
    ("mimiciv_derived", "first_day_sofa2"): [
        "stay_id",
        "sofa2_total",
    ],
    ("mimiciv_derived", "bg"): [
        "subject_id",
        "hadm_id",
        "charttime",
        "lactate",
    ],
    ("mimiciv_derived", "chemistry"): [
        "subject_id",
        "hadm_id",
        "charttime",
        "creatinine",
    ],
}


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fetch_columns(cur, schema: str, table: str) -> list[str]:
    rows = cur.execute(
        """
        select column_name
        from information_schema.columns
        where table_schema = %s
          and table_name = %s
        order by ordinal_position
        """,
        (schema, table),
    ).fetchall()
    return [row["column_name"] for row in rows]


def table_audit(cur) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for schema, table, role in REQUIRED_TABLES:
        info = cur.execute(
            """
            select
              to_regclass(%s) is not null as exists,
              coalesce(c.reltuples::bigint, null) as approx_rows
            from (select 1) x
            left join pg_class c on c.oid = to_regclass(%s)
            """,
            (f"{schema}.{table}", f"{schema}.{table}"),
        ).fetchone()
        rows.append(
            {
                "schema": schema,
                "table": table,
                "role": role,
                "exists": info["exists"],
                "approx_rows": info["approx_rows"] if info["exists"] else "",
            }
        )
    return rows


def column_audit(cur) -> tuple[list[dict[str, Any]], dict[tuple[str, str], list[str]]]:
    rows: list[dict[str, Any]] = []
    all_columns: dict[tuple[str, str], list[str]] = {}
    for schema, table in sorted(REQUIRED_COLUMNS):
        columns = fetch_columns(cur, schema, table)
        all_columns[(schema, table)] = columns
        column_set = set(columns)
        for required in REQUIRED_COLUMNS[(schema, table)]:
            rows.append(
                {
                    "schema": schema,
                    "table": table,
                    "required_column": required,
                    "present": required in column_set,
                }
            )
    return rows, all_columns


def write_summary(
    path: Path,
    connection_row: dict[str, Any],
    tables: list[dict[str, Any]],
    columns: list[dict[str, Any]],
    all_columns: dict[tuple[str, str], list[str]],
) -> None:
    missing_tables = [row for row in tables if not row["exists"]]
    missing_columns = [row for row in columns if not row["present"]]
    lines = [
        "# MIMIC-IV HRC feasibility audit",
        "",
        "## Connection",
        "",
        f"- database: `{connection_row['database']}`",
        f"- user: `{connection_row['user']}`",
        f"- server_version: `{connection_row['server_version']}`",
        "",
        "## Required tables",
        "",
    ]
    for row in tables:
        status = "OK" if row["exists"] else "MISSING"
        approx = f", approx_rows={row['approx_rows']}" if row["approx_rows"] != "" else ""
        lines.append(f"- {status}: `{row['schema']}.{row['table']}` ({row['role']}{approx})")

    lines.extend(["", "## Required columns", ""])
    if missing_columns:
        for row in missing_columns:
            lines.append(
                f"- MISSING: `{row['schema']}.{row['table']}.{row['required_column']}`"
            )
    else:
        lines.append("- All required columns are present.")

    lines.extend(["", "## Available columns in key HRC tables", ""])
    for schema, table in sorted(all_columns):
        lines.append(f"### {schema}.{table}")
        lines.append("")
        lines.append(", ".join(f"`{column}`" for column in all_columns[(schema, table)]) or "_No columns found._")
        lines.append("")

    lines.extend(["", "## Verdict", ""])
    if missing_tables or missing_columns:
        lines.append(
            "MIMIC connection works, but the extraction specification needs adjustment before HRC extraction."
        )
    else:
        lines.append(
            "MIMIC connection and core HRC variables are available. Proceed to cohort and HRC component extraction."
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit MIMIC-IV tables needed for HRC v1.0.")
    parser.add_argument("--server", default="mimiciv", help="pgAdmin saved server name")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    with connect_pgadmin_server(args.server) as conn:
        conn.row_factory = dict_row
        with conn.cursor() as cur:
            connection_row = cur.execute(
                """
                select
                  current_database() as database,
                  current_user as "user",
                  split_part(version(), ' ', 2) as server_version
                """
            ).fetchone()
            tables = table_audit(cur)
            columns, all_columns = column_audit(cur)

    write_csv(
        args.outdir / "mimic_required_table_audit.csv",
        tables,
        ["schema", "table", "role", "exists", "approx_rows"],
    )
    write_csv(
        args.outdir / "mimic_required_column_audit.csv",
        columns,
        ["schema", "table", "required_column", "present"],
    )
    write_csv(
        args.outdir / "mimic_connection_audit.csv",
        [connection_row],
        ["database", "user", "server_version"],
    )
    write_summary(
        args.outdir / "mimic_schema_audit_summary.md",
        connection_row,
        tables,
        columns,
        all_columns,
    )

    missing_tables = sum(1 for row in tables if not row["exists"])
    missing_columns = sum(1 for row in columns if not row["present"])
    print(f"Connected to {connection_row['database']} as {connection_row['user']}.")
    print(f"Required tables checked: {len(tables)}; missing: {missing_tables}.")
    print(f"Required columns checked: {len(columns)}; missing: {missing_columns}.")
    print(f"Summary: {args.outdir / 'mimic_schema_audit_summary.md'}")


if __name__ == "__main__":
    main()
