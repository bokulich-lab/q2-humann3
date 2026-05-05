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
import tempfile

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
)

from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    HumannReactionDirectoryFormat,
    MetaphlanDatabaseDirFmt,
)
from q2_humann.utils import run_humann_command
from q2_sapienns.plugin_setup import (
    HumannGeneFamilyDirectoryFormat,
    HumannPathAbundanceDirectoryFormat,
    MetaphlanMergedAbundanceDirectoryFormat,
)


def _stage_sample_input(
    sample_manifest: pd.Series, sample_id: str, work_dir: Path
) -> Path:
    """Return a sample FASTQ path, concatenating paired reads if needed."""
    staged_input = work_dir / f"{sample_id}.fastq.gz"
    input_paths = [
        sample_manifest.get(direction)
        for direction in ("forward", "reverse")
        if pd.notna(sample_manifest.get(direction))
    ]
    if not input_paths:
        raise RuntimeError(f"No reads found for sample {sample_id!r}.")

    if len(input_paths) == 1:
        return Path(input_paths[0])

    # the docs recommend concatenating paired reads
    work_dir.mkdir()
    with open(staged_input, "wb") as out_fh:
        for input_path in input_paths:
            with open(input_path, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

    return staged_input


def _join_humann_tables(
    run_output_dir: Path, file_name: str, output_path: Path
) -> None:
    """Join per-sample HUMANN output tables for a given file stem."""
    cmd = [
        "humann_join_tables",
        "--input",
        str(run_output_dir),
        "--output",
        str(output_path),
        "--file_name",
        file_name,
    ]
    run_humann_command(cmd)


def _merge_metaphlan_profiles(
    run_output_dir: Path, output_path: Path
) -> None:
    """Merge MetaPhlAn bugs-list profiles into one abundance table."""
    profile_paths = sorted(
        run_output_dir.glob("*_humann_temp/*_metaphlan_bugs_list.tsv")
    )
    if not profile_paths:
        raise RuntimeError(
            "No MetaPhlAn bugs-list files were available to merge."
        )

    cmd = [
        "merge_metaphlan_tables.py",
        *[str(path) for path in profile_paths],
        "-o",
        str(output_path),
    ]
    run_humann_command(cmd)
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise RuntimeError(
            "MetaPhlAn table merge did not produce a merged profile table."
        )


def _read_database_metadata_value(
    database: HumannDatabaseDirFmt | MetaphlanDatabaseDirFmt, key: str
) -> str:
    """Read a value from staged database artifact metadata."""
    with (database.path / "metadata.json").open() as fh:
        metadata = json.load(fh)
    return metadata[key]


def _regroup_gene_families_to_reactions(
    gene_families_path: Path,
    translated_search_database: HumannDatabaseDirFmt,
    output_path: Path,
) -> None:
    """Regroup HUMANN gene families into reactions for the database build."""
    build = _read_database_metadata_value(
        translated_search_database, "build"
    )
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
    run_humann_command(cmd)


def run_humann(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    nucleotide_database: HumannDatabaseDirFmt,
    translated_search_database: HumannDatabaseDirFmt,
    metaphlan_database: MetaphlanDatabaseDirFmt,
    threads: int = 1,
    memory_use: str = "minimum",
    prescreen_threshold: float = 0.01,
    nucleotide_identity_threshold: float = 0.0,
    nucleotide_query_coverage_threshold: float = 90.0,
    nucleotide_subject_coverage_threshold: float = 50.0,
    translated_identity_threshold: float = None,
    translated_query_coverage_threshold: float = 90.0,
    translated_subject_coverage_threshold: float = 50.0,
    evalue: float = 1.0,
    gap_fill: bool = True,
    minpath: bool = True,
    pathways: str = "metacyc",
    output_max_decimals: int = 10,
    log_level: str = "DEBUG",
) -> (HumannGeneFamilyDirectoryFormat,
      HumannPathAbundanceDirectoryFormat,
      MetaphlanMergedAbundanceDirectoryFormat,
      HumannReactionDirectoryFormat):
    """Run HUMANN per sample and return merged functional profile tables."""
    with tempfile.TemporaryDirectory(prefix="q2-humann-run-") as tmpdir:
        tmpdir = Path(tmpdir)
        run_output_dir = tmpdir / "humann-output"
        run_output_dir.mkdir()
        metaphlan_index = _read_database_metadata_value(
            metaphlan_database, "index"
        )
        metaphlan_options = (
            f"--bowtie2db {metaphlan_database.path} -x {metaphlan_index}"
        )

        for sample_id, sample_manifest in reads.manifest.iterrows():
            sample_work_dir = tmpdir / sample_id
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
                str(nucleotide_database.path / "chocophlan"),
                "--protein-database",
                str(translated_search_database.path / "uniref"),
                "--metaphlan-options",
                metaphlan_options,
                "--output-format",
                "tsv",
                "--memory-use",
                memory_use,
                "--prescreen-threshold",
                str(prescreen_threshold),
                "--nucleotide-identity-threshold",
                str(nucleotide_identity_threshold),
                "--nucleotide-query-coverage-threshold",
                str(nucleotide_query_coverage_threshold),
                "--nucleotide-subject-coverage-threshold",
                str(nucleotide_subject_coverage_threshold),
                "--evalue",
                str(evalue),
                "--gap-fill",
                "on" if gap_fill else "off",
                "--minpath",
                "on" if minpath else "off",
                "--pathways",
                pathways,
                "--output-max-decimals",
                str(output_max_decimals),
                "--log-level",
                log_level,
            ]
            if translated_identity_threshold is not None:
                cmd.extend(
                    [
                        "--translated-identity-threshold",
                        str(translated_identity_threshold),
                    ]
                )
            cmd.extend(
                [
                    "--translated-query-coverage-threshold",
                    str(translated_query_coverage_threshold),
                    "--translated-subject-coverage-threshold",
                    str(translated_subject_coverage_threshold),
                ]
            )
            run_humann_command(cmd)

        gene_families = HumannGeneFamilyDirectoryFormat()
        path_abundance = HumannPathAbundanceDirectoryFormat()
        metaphlan_profile = MetaphlanMergedAbundanceDirectoryFormat()
        reactions = HumannReactionDirectoryFormat()

        gene_families_table = gene_families.path / "table.tsv"
        path_abundance_table = path_abundance.path / "table.tsv"
        metaphlan_profile_table = metaphlan_profile.path / "table.tsv"
        reactions_table = reactions.path / "table.tsv"

        _join_humann_tables(
            run_output_dir, "genefamilies", gene_families_table
        )
        _join_humann_tables(
            run_output_dir, "pathabundance", path_abundance_table
        )
        _merge_metaphlan_profiles(
            run_output_dir, metaphlan_profile_table
        )
        _regroup_gene_families_to_reactions(
            gene_families_table,
            translated_search_database,
            reactions_table,
        )

        return gene_families, path_abundance, metaphlan_profile, reactions
