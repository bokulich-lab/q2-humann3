# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
import tempfile

from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    MetaphlanDatabaseDirFmt,
)
from q2_humann.utils import copy_directory_contents, run_humann_command


TRANSLATED_SEARCH_DATABASE_BUILDS = (
    "uniref90_diamond",
    "uniref90_ec_filtered_diamond",
    "uniref50_diamond",
    "uniref50_ec_filtered_diamond",
)


def _infer_metaphlan_index(install_dir: Path) -> str:
    for path in sorted(install_dir.glob("*.pkl")):
        return path.stem
    raise RuntimeError(
        "Unable to infer the MetaPhlAn database index from the downloaded "
        "files."
    )


def _validate_metaphlan_database(
    install_dir: Path, index: str
) -> None:
    pkl_path = install_dir / f"{index}.pkl"
    if not pkl_path.exists():
        raise RuntimeError(
            "MetaPhlAn download did not produce the expected database pickle "
            f"for index {index!r}."
        )

    bt2_ext = None
    if (install_dir / f"{index}.1.bt2").exists():
        bt2_ext = "bt2"
    elif (install_dir / f"{index}.1.bt2l").exists():
        bt2_ext = "bt2l"

    if bt2_ext is None:
        raise RuntimeError(
            "MetaPhlAn download did not produce Bowtie2 index files for "
            f"index {index!r}."
        )

    required_suffixes = (
        f"1.{bt2_ext}",
        f"2.{bt2_ext}",
        f"3.{bt2_ext}",
        f"4.{bt2_ext}",
        f"rev.1.{bt2_ext}",
        f"rev.2.{bt2_ext}",
    )
    missing = [
        suffix for suffix in required_suffixes
        if not (install_dir / f"{index}.{suffix}").exists()
    ]
    if missing:
        missing_str = ", ".join(missing)
        raise RuntimeError(
            "MetaPhlAn download produced an incomplete Bowtie2 database for "
            f"index {index!r}. Missing files: {missing_str}."
        )


def _write_database_metadata(
    artifact: HumannDatabaseDirFmt, database_kind: str, build: str
) -> None:
    metadata = {
        "database_kind": database_kind,
        "build": build,
    }
    with (artifact.path / "metadata.json").open("w") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)


def _write_metaphlan_database_metadata(
    artifact: MetaphlanDatabaseDirFmt, index: str
) -> None:
    metadata = {
        "database_kind": "metaphlan",
        "index": index,
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
            run_humann_command(cmd)
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
        copy_directory_contents(install_dir, payload_dir)
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


def download_metaphlan_database(
    index: str = "latest",
    cpus: int = 1,
) -> MetaphlanDatabaseDirFmt:
    with tempfile.TemporaryDirectory(prefix="q2-metaphlan-db-") as tmpdir:
        install_dir = Path(tmpdir) / "downloaded-database"
        install_dir.mkdir()

        cmd = [
            "metaphlan",
            "--install",
            "--bowtie2db",
            str(install_dir),
            "-x",
            index,
            "--nproc",
            str(cpus),
        ]
        run_humann_command(cmd)

        if not any(install_dir.iterdir()):
            raise RuntimeError(
                "MetaPhlAn completed successfully but did not produce any "
                "database files."
            )

        resolved_index = _infer_metaphlan_index(install_dir)
        _validate_metaphlan_database(install_dir, resolved_index)
        artifact = MetaphlanDatabaseDirFmt()
        payload_dir = artifact.path / "data"
        payload_dir.mkdir()
        copy_directory_contents(install_dir, payload_dir)
        _write_metaphlan_database_metadata(artifact, resolved_index)
        return artifact
