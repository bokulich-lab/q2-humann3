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

from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    HumannGeneFamilyDirectoryFormat,
    HumannPathAbundanceDirectoryFormat,
    HumannReactionDirectoryFormat,
    MetaphlanDatabaseDirFmt,
    MetaphlanMergedAbundanceDirectoryFormat,
)

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
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        detail = stderr or stdout or str(exc)
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {detail}"
        ) from exc


def _stage_sample_input(
    sample_manifest: pd.DataFrame, sample_id: str, work_dir: Path
) -> Path:
    staged_input = work_dir / f"{sample_id}.fastq.gz"
    forward_reads = sample_manifest.loc[
        sample_manifest["direction"] == "forward", "absolute_path"
    ].tolist()
    reverse_reads = sample_manifest.loc[
        sample_manifest["direction"] == "reverse", "absolute_path"
    ].tolist()

    input_paths = forward_reads + reverse_reads
    if not input_paths:
        raise RuntimeError(f"No reads found for sample {sample_id!r}.")

    if len(input_paths) == 1:
        shutil.copy2(input_paths[0], staged_input)
        return staged_input

    with open(staged_input, "wb") as out_fh:
        for input_path in input_paths:
            with open(input_path, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

    return staged_input


def _copy_to_single_file_dirfmt(source: Path, destination) -> None:
    shutil.copy2(source, destination.path / "table.tsv")


def _join_humann_tables(
    run_output_dir: Path, file_name: str, output_path: Path
) -> None:
    cmd = [
        "humann_join_tables",
        "--input",
        str(run_output_dir),
        "--output",
        str(output_path),
        "--file_name",
        file_name,
    ]
    _run_humann_command(cmd)


def _read_database_build(database: HumannDatabaseDirFmt) -> str:
    with (database.path / "metadata.json").open() as fh:
        metadata = json.load(fh)
    return metadata["build"]


def _regroup_gene_families_to_reactions(
    gene_families_path: Path,
    translated_search_database: HumannDatabaseDirFmt,
    output_path: Path,
) -> None:
    build = _read_database_build(translated_search_database)
    if build.startswith("uniref50"):
        groups = "uniref50_rxn"
    elif build.startswith("uniref90"):
        groups = "uniref90_rxn"
    else:
        raise RuntimeError(
            f"Unsupported translated-search database build {build!r} for "
            "reaction regrouping."
        )

    cmd = [
        "humann_regroup_table",
        "--input",
        str(gene_families_path),
        "--groups",
        groups,
        "--output",
        str(output_path),
    ]
    _run_humann_command(cmd)


def _copy_directory_contents(source_dir: Path, target_dir: Path) -> None:
    for path in source_dir.iterdir():
        destination = target_dir / path.name
        if path.is_dir():
            shutil.copytree(path, destination)
        else:
            shutil.copy2(path, destination)


def _infer_metaphlan_index(install_dir: Path) -> str:
    for path in sorted(install_dir.glob("*.pkl")):
        return path.stem
    raise RuntimeError(
        "Unable to infer the MetaPhlAn database index from the downloaded "
        "files."
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


def download_metaphlan_database(
    index: str = "latest",
) -> MetaphlanDatabaseDirFmt:
    with tempfile.TemporaryDirectory(prefix="q2-metaphlan-db-") as tmpdir:
        install_dir = Path(tmpdir) / "downloaded-database"
        install_dir.mkdir()

        cmd = [
            "metaphlan",
            "--install",
            "--db_dir",
            str(install_dir),
            "-x",
            index,
        ]
        _run_humann_command(cmd)

        if not any(install_dir.iterdir()):
            raise RuntimeError(
                "MetaPhlAn completed successfully but did not produce any "
                "database files."
            )

        resolved_index = _infer_metaphlan_index(install_dir)
        artifact = MetaphlanDatabaseDirFmt()
        payload_dir = artifact.path / "data"
        payload_dir.mkdir()
        _copy_directory_contents(install_dir, payload_dir)
        _write_metaphlan_database_metadata(artifact, resolved_index)
        return artifact


def run_humann(
    reads: pd.DataFrame,
    nucleotide_database: HumannDatabaseDirFmt,
    translated_search_database: HumannDatabaseDirFmt,
    threads: int = 1,
) -> tuple[
    HumannGeneFamilyDirectoryFormat,
    HumannPathAbundanceDirectoryFormat,
    MetaphlanMergedAbundanceDirectoryFormat,
    HumannReactionDirectoryFormat,
]:
    with tempfile.TemporaryDirectory(prefix="q2-humann-run-") as tmpdir:
        tmpdir = Path(tmpdir)
        run_output_dir = tmpdir / "humann-output"
        run_output_dir.mkdir()

        for sample_id, sample_manifest in reads.groupby("sample-id"):
            sample_work_dir = tmpdir / sample_id
            sample_work_dir.mkdir()
            staged_input = _stage_sample_input(
                sample_manifest, sample_id, sample_work_dir
            )

            cmd = [
                "humann",
                "--input",
                str(staged_input),
                "--output",
                str(run_output_dir),
                "--output-basename",
                sample_id,
                "--threads",
                str(threads),
                "--nucleotide-database",
                str(nucleotide_database.path / "data"),
                "--protein-database",
                str(translated_search_database.path / "data"),
                "--output-format",
                "tsv",
            ]
            _run_humann_command(cmd)

        joined_gene_families = tmpdir / "joined_genefamilies.tsv"
        joined_path_abundance = tmpdir / "joined_pathabundance.tsv"
        joined_metaphlan_profile = tmpdir / "joined_metaphlan_profile.tsv"
        joined_reactions = tmpdir / "joined_reactions.tsv"
        _join_humann_tables(
            run_output_dir, "genefamilies", joined_gene_families
        )
        _join_humann_tables(
            run_output_dir, "pathabundance", joined_path_abundance
        )
        _join_humann_tables(
            run_output_dir, "metaphlan_profile", joined_metaphlan_profile
        )
        _regroup_gene_families_to_reactions(
            joined_gene_families,
            translated_search_database,
            joined_reactions,
        )

        gene_families = HumannGeneFamilyDirectoryFormat()
        path_abundance = HumannPathAbundanceDirectoryFormat()
        metaphlan_profile = MetaphlanMergedAbundanceDirectoryFormat()
        reactions = HumannReactionDirectoryFormat()
        _copy_to_single_file_dirfmt(joined_gene_families, gene_families)
        _copy_to_single_file_dirfmt(joined_path_abundance, path_abundance)
        _copy_to_single_file_dirfmt(
            joined_metaphlan_profile, metaphlan_profile
        )
        _copy_to_single_file_dirfmt(joined_reactions, reactions)
        return gene_families, path_abundance, metaphlan_profile, reactions
