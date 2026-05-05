# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from io import StringIO
from pathlib import Path
import shutil
import tempfile

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    JoinedSequencesWithQuality,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData

from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    HumannReactionTable,
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
    unit_labels = {
        "genefamilies": "Abundance-RPKs",
        "pathabundance": "Abundance",
    }
    if file_name in unit_labels:
        _normalize_table_sample_headers(
            output_path, file_name, unit_labels[file_name]
        )


def _normalize_table_sample_headers(
    table_path: Path, file_name: str, unit_label: str
) -> None:
    """Normalize HUMANN sample headers to include expected unit labels."""
    lines = table_path.read_text().splitlines()
    if not lines:
        raise RuntimeError(f"HUMANN table {table_path} is empty.")

    header = lines[0].split("\t")
    for idx, field in enumerate(header[1:], start=1):
        if field.endswith(unit_label):
            continue
        if field.endswith(f"_{file_name}"):
            sample_id = field.removesuffix(f"_{file_name}")
        else:
            sample_id = field
        header[idx] = f"{sample_id}_{unit_label}"

    lines[0] = "\t".join(header)
    table_path.write_text("\n".join(lines) + "\n")
    _assert_table_sample_headers_end_with(
        table_path, unit_label.rsplit("-", 1)[-1]
    )


def _assert_table_sample_headers_end_with(
    table_path: Path, unit_label: str
) -> None:
    """Confirm table sample headers satisfy downstream format validation."""
    lines = table_path.read_text().splitlines()
    if not lines:
        raise RuntimeError(f"HUMANN table {table_path} is empty.")

    header = lines[0].split("\t")
    invalid_headers = [
        field for field in header[1:] if not field.endswith(unit_label)
    ]
    if invalid_headers:
        raise RuntimeError(
            "HUMANN table sample headers were not normalized correctly. "
            f"Expected sample columns to end with {unit_label!r}, found: "
            + ", ".join(invalid_headers)
        )


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


def _read_table(path: Path) -> pd.DataFrame:
    """Read a HUMANN-style TSV table."""
    lines = path.read_text().splitlines()
    while lines and lines[0].startswith("#mpa_"):
        lines.pop(0)
    if not lines:
        raise RuntimeError(f"Table {path} is empty.")

    return pd.read_csv(StringIO("\n".join(lines)), sep="\t")


def _merge_table_dirfmts(
    tables: list,
    output_dirfmt,
    key_columns: list[str] = None,
    unit_label: str = None,
    file_name: str = None,
):
    """Merge table directory formats into one output directory format."""
    if not tables:
        raise RuntimeError("No table artifacts were available to collate.")

    merged = None
    for table_dirfmt in tables:
        table_path = table_dirfmt.path / "table.tsv"
        table = _read_table(table_path)
        if key_columns is None:
            key_columns = [table.columns[0]]

        if merged is None:
            merged = table
        else:
            merged = merged.merge(table, on=key_columns, how="outer")

    output = output_dirfmt()
    merged.fillna(0).to_csv(output.path / "table.tsv", sep="\t", index=False)
    if unit_label is not None:
        _normalize_table_sample_headers(
            output.path / "table.tsv", file_name, unit_label
        )
    return output


def collate_gene_families(
    tables: HumannGeneFamilyDirectoryFormat,
) -> HumannGeneFamilyDirectoryFormat:
    """Collate HUMANN gene-family tables across partitions."""
    return _merge_table_dirfmts(
        tables,
        HumannGeneFamilyDirectoryFormat,
        unit_label="Abundance-RPKs",
        file_name="genefamilies",
    )


def collate_path_abundance(
    tables: HumannPathAbundanceDirectoryFormat,
) -> HumannPathAbundanceDirectoryFormat:
    """Collate HUMANN pathway-abundance tables across partitions."""
    return _merge_table_dirfmts(
        tables,
        HumannPathAbundanceDirectoryFormat,
        unit_label="Abundance",
        file_name="pathabundance",
    )


def collate_metaphlan_profiles(
    tables: MetaphlanMergedAbundanceDirectoryFormat,
) -> MetaphlanMergedAbundanceDirectoryFormat:
    """Collate merged MetaPhlAn profile tables across partitions."""
    return _merge_table_dirfmts(
        tables,
        MetaphlanMergedAbundanceDirectoryFormat,
        key_columns=["clade_name", "NCBI_tax_id"],
    )


def collate_reactions(
    tables: HumannReactionDirectoryFormat,
) -> HumannReactionDirectoryFormat:
    """Collate HUMANN reaction tables across partitions."""
    return _merge_table_dirfmts(tables, HumannReactionDirectoryFormat)


def _run_humann(
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
        _assert_table_sample_headers_end_with(gene_families_table, "RPKs")
        _assert_table_sample_headers_end_with(
            path_abundance_table, "Abundance"
        )

        return gene_families, path_abundance, metaphlan_profile, reactions


def run_humann(
    ctx,
    reads,
    nucleotide_database,
    translated_search_database,
    metaphlan_database,
    threads=1,
    memory_use="minimum",
    prescreen_threshold=0.01,
    nucleotide_identity_threshold=0.0,
    nucleotide_query_coverage_threshold=90.0,
    nucleotide_subject_coverage_threshold=50.0,
    translated_identity_threshold=None,
    translated_query_coverage_threshold=90.0,
    translated_subject_coverage_threshold=50.0,
    evalue=1.0,
    gap_fill=True,
    minpath=True,
    pathways="metacyc",
    output_max_decimals=10,
    log_level="DEBUG",
    num_partitions=1,
):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["ctx", "reads", "num_partitions"]
    }

    if reads.type <= SampleData[
        SequencesWithQuality | JoinedSequencesWithQuality
    ]:
        _partition_reads = ctx.get_action("types", "partition_samples_single")
    elif reads.type <= SampleData[PairedEndSequencesWithQuality]:
        _partition_reads = ctx.get_action("types", "partition_samples_paired")
    else:
        raise NotImplementedError()

    _humann = ctx.get_action("humann", "_run_humann")
    _collate_gene_families = ctx.get_action(
        "humann", "collate_gene_families"
    )
    _collate_path_abundance = ctx.get_action(
        "humann", "collate_path_abundance"
    )
    _collate_metaphlan_profiles = ctx.get_action(
        "humann", "collate_metaphlan_profiles"
    )

    (partitioned_reads,) = _partition_reads(reads, num_partitions)

    genes_all, path_all, metaphlan_all = [], [], []
    for _reads in partitioned_reads.values():
        (gene_families, path_abundance, metaphlan_profile, reactions) = _humann(
            reads=_reads, **kwargs
        )
        genes_all.append(gene_families)
        path_all.append(path_abundance)
        metaphlan_all.append(metaphlan_profile)

    (gene_families,) = _collate_gene_families(genes_all)
    (path_abundance,) = _collate_path_abundance(path_all)
    (metaphlan_profile,) = _collate_metaphlan_profiles(metaphlan_all)

    reactions = HumannReactionDirectoryFormat()
    _regroup_gene_families_to_reactions(
        gene_families.view(
            HumannGeneFamilyDirectoryFormat
        ).path / "table.tsv",
        translated_search_database.view(HumannDatabaseDirFmt),
        reactions.path / "table.tsv",
    )

    return (
        gene_families,
        path_abundance,
        metaphlan_profile,
        ctx.make_artifact(HumannReactionTable, reactions),
    )
