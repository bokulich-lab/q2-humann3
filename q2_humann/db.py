# ----------------------------------------------------------------------------
# Copyright (c) 2026, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
from typing import Union

from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    MetaphlanDatabaseDirFmt,
)
from q2_humann.utils import run_humann_command


TRANSLATED_SEARCH_DATABASE_BUILDS = (
    "uniref90_diamond",
    "uniref90_ec_filtered_diamond",
    "uniref50_diamond",
    "uniref50_ec_filtered_diamond",
)


def _infer_metaphlan_index(install_dir: Path) -> str:
    """Infer the installed MetaPhlAn index name from its pickle file."""
    for path in sorted(install_dir.glob("*.pkl")):
        return path.stem
    raise RuntimeError(
        "Unable to infer the MetaPhlAn database index from the downloaded "
        "files."
    )


def _validate_metaphlan_database(
    install_dir: Path, index: str
) -> None:
    """Validate that a MetaPhlAn install contains pickle and Bowtie2 files."""
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
    artifact: Union[MetaphlanDatabaseDirFmt, HumannDatabaseDirFmt],
    **metadata,
) -> None:
    """Write database metadata fields into an artifact."""
    with (artifact.path / "metadata.json").open("w") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)


def _download_humann_database(
    database: str, build: str, database_kind: str
) -> HumannDatabaseDirFmt:
    """Download a HUMANN database and stage it as a QIIME 2 artifact."""
    artifact = HumannDatabaseDirFmt()
    install_dir = artifact.path

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

    _write_database_metadata(
        artifact,
        database_kind=database_kind,
        build=build,
    )
    return artifact


def download_chocophlan_database() -> HumannDatabaseDirFmt:
    """Download and stage the full HUMANN ChocoPhlAn database."""
    return _download_humann_database(
        database="chocophlan",
        build="full",
        database_kind="chocophlan",
    )


def download_translated_search_database(
    build: str = "uniref90_diamond",
) -> HumannDatabaseDirFmt:
    """Download and stage a HUMANN translated-search protein database."""
    return _download_humann_database(
        database="uniref",
        build=build,
        database_kind="translated-search",
    )


def download_metaphlan_database(
    index: str = "latest",
    cpus: int = 1,
) -> MetaphlanDatabaseDirFmt:
    """Download, validate, and stage a MetaPhlAn database."""
    artifact = MetaphlanDatabaseDirFmt()
    install_dir = artifact.path

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
    _write_database_metadata(
        artifact,
        database_kind="metaphlan",
        index=resolved_index,
    )
    return artifact
