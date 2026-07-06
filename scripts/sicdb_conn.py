from __future__ import annotations

import os

import psycopg


class DatabaseConnectionError(RuntimeError):
    """Raised when required SICdb connection environment variables are missing."""


def get_sicdb_conninfo() -> dict:
    conninfo = {
        "host": os.environ.get("SICDB_HOST", "localhost"),
        "port": int(os.environ.get("SICDB_PORT", "55432")),
        "dbname": os.environ.get("SICDB_DBNAME", "sicdb"),
        "user": os.environ.get("SICDB_USER"),
        "password": os.environ.get("SICDB_PASSWORD"),
        "connect_timeout": int(os.environ.get("SICDB_CONNECT_TIMEOUT", "20")),
    }
    missing = [name for name in ("user", "password") if not conninfo.get(name)]
    if missing:
        env_names = ", ".join(f"SICDB_{name.upper()}" for name in missing)
        raise DatabaseConnectionError(f"Missing required environment variable(s): {env_names}")
    return conninfo


def connect_sicdb() -> psycopg.Connection:
    return psycopg.connect(**get_sicdb_conninfo())
