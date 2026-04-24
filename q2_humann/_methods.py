# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
import shutil
import pandas as pd
import subprocess
import tempfile

from q2_humann._types_and_formats import HumannDatabaseDirFmt

EXTERNAL_CMD_WARNING = (
    "Running external command line application. This may print additional "
    "output below."
)

TRANSLATED_SEARCH_DATABASE_BUILDS = (
    "uniref90_diamond",
    "uniref90_ec_filtered_diamond",
    "uniref50_diamond",
    "uniref50_ec_filtered_diamond",
)


def duplicate_table(table: pd.DataFrame) -> pd.DataFrame:
    return table


def run_command(
    cmd: list[str], env=None, verbose: bool = True, pipe: bool = False,
    **kwargs
):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")

    if pipe:
        result = subprocess.run(
            cmd,
            env=env,
            check=True,
            capture_output=True,
            text=True,
            **kwargs,
        )
        return result

    if env:
        subprocess.run(cmd, env=env, check=True, **kwargs)
    else:
        subprocess.run(cmd, check=True, **kwargs)


def _run_humann_command(cmd: list[str]) -> None:
    try:
        run_command(cmd)
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip() or exc.stdout.strip()
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {detail}"
        ) from exc


def _copy_directory_contents(source_dir: Path, target_dir: Path) -> None:
    for path in source_dir.iterdir():
        destination = target_dir / path.name
        if path.is_dir():
            shutil.copytree(path, destination)
        else:
            shutil.copy2(path, destination)


def _write_database_metadata(
    artifact: HumannDatabaseDirFmt, database_kind: str, build: str
) -> None:
    metadata = {
        "database_kind": database_kind,
        "build": build,
    }
    with (artifact.path / "metadata.json").open("w") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)


def _download_humann_database(
    database: str, build: str, database_kind: str
) -> HumannDatabaseDirFmt:
    with tempfile.TemporaryDirectory(prefix="q2-humann-db-") as tmpdir:
        install_dir = Path(tmpdir) / "downloaded-database"
        install_dir.mkdir()

        cmd = [
            "humann_databases",
            "--download",
            database,
            build,
            str(install_dir),
            "--update-config",
            "no",
        ]

        try:
            _run_humann_command(cmd)
        except RuntimeError as exc:
            raise RuntimeError(
                f"humann_databases failed: {exc}"
            ) from exc

        if not any(install_dir.iterdir()):
            raise RuntimeError(
                "humann_databases completed successfully but did not produce "
                "any database files."
            )

        artifact = HumannDatabaseDirFmt()
        payload_dir = artifact.path / "data"
        payload_dir.mkdir()
        _copy_directory_contents(install_dir, payload_dir)
        _write_database_metadata(artifact, database_kind, build)
        return artifact


def download_chocophlan_database() -> HumannDatabaseDirFmt:
    return _download_humann_database(
        database="chocophlan",
        build="full",
        database_kind="chocophlan",
    )


def download_translated_search_database(
    build: str = "uniref90_diamond",
) -> HumannDatabaseDirFmt:
    return _download_humann_database(
        database="uniref",
        build=build,
        database_kind="translated-search",
    )
