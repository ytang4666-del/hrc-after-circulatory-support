from __future__ import annotations

import os
from typing import Any

import psycopg


SERVER_PREFIX = {
    "mimiciv": "MIMIC_DB",
    "eicu": "EICU_DB",
}


class DatabaseConnectionError(RuntimeError):
    """Raised when required database connection environment variables are missing."""


def get_pgadmin_conninfo(server_name: str = "mimiciv", connect_timeout: int = 20) -> dict[str, Any]:
    """Return PostgreSQL connection info from environment variables.

    This public-release helper intentionally does not read pgAdmin, keychains,
    local config files, or stored passwords.
    """
    key = server_name.lower()
    if key not in SERVER_PREFIX:
        raise DatabaseConnectionError(
            f"Unknown server '{server_name}'. Expected one of: {', '.join(SERVER_PREFIX)}"
        )

    prefix = SERVER_PREFIX[key]
    conninfo = {
        "host": os.environ.get(f"{prefix}_HOST", "localhost"),
        "port": int(os.environ.get(f"{prefix}_PORT", "5432")),
        "dbname": os.environ.get(f"{prefix}_NAME", key),
        "user": os.environ.get(f"{prefix}_USER"),
        "password": os.environ.get(f"{prefix}_PASSWORD"),
        "connect_timeout": connect_timeout,
    }

    missing = [name for name in ("user", "password") if not conninfo.get(name)]
    if missing:
        env_names = ", ".join(f"{prefix}_{name.upper()}" for name in missing)
        raise DatabaseConnectionError(f"Missing required environment variable(s): {env_names}")

    return conninfo


def connect_pgadmin_server(server_name: str = "mimiciv") -> psycopg.Connection:
    return psycopg.connect(**get_pgadmin_conninfo(server_name))
