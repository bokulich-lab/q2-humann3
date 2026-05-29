# ----------------------------------------------------------------------------
# Copyright (c) 2026, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Bool, Citations, Choices, Float, Int, List, Plugin, Range, Str
)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    JoinedSequencesWithQuality,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)

from q2_humann3 import __version__
from q2_humann3.db import (
    TRANSLATED_SEARCH_DATABASE_BUILDS,
    download_chocophlan_database,
    download_metaphlan_database,
    download_translated_search_database,
)
from q2_humann3.run import (
    collate_gene_families,
    collate_metaphlan_profiles,
    collate_path_abundance,
    collate_reactions,
    run_humann,
    _run_humann,
)
from q2_humann3._types_and_formats import (
    ChocoPhlAn,
    HumannDatabase,
    HumannDatabaseDirFmt,
    HumannDatabaseFileFormat,
    HumannReactionDirectoryFormat,
    HumannReactionFormat,
    HumannReactionTable,
    MetaphlanDatabase,
    MetaphlanDatabaseDirFmt,
    MetaphlanDatabaseFileFormat,
    MetaphlanDatabaseMetadataFormat,
    HumannDatabaseMetadataFormat,
    TranslatedSearch,
)
from q2_sapienns.plugin_setup import (
    HumannGeneFamilyTable,
    HumannPathAbundanceTable,
    MetaphlanMergedAbundanceTable,
)

citations = Citations.load("citations.bib", package="q2_humann3")

plugin = Plugin(
    name="humann3",
    version=__version__,
    website="https://github.com/biobakery/humann",
    package="q2_humann3",
    description=(
        "A QIIME 2 plugin for wrapping HUMANN 3 functionality for functional "
        "profiling of microbial communities."
    ),
    short_description="QIIME 2 wrapper for HUMANN 3.",
    citations=[citations['Caporaso-Bolyen-2024']]
)

plugin.register_semantic_types(
    HumannDatabase,
    ChocoPhlAn,
    TranslatedSearch,
    HumannReactionTable,
    MetaphlanDatabase,
)
plugin.register_formats(
    HumannDatabaseMetadataFormat,
    HumannDatabaseFileFormat,
    HumannDatabaseDirFmt,
    HumannReactionFormat,
    HumannReactionDirectoryFormat,
    MetaphlanDatabaseMetadataFormat,
    MetaphlanDatabaseFileFormat,
    MetaphlanDatabaseDirFmt,
)
plugin.register_artifact_class(
    HumannDatabase[ChocoPhlAn],
    directory_format=HumannDatabaseDirFmt,
    description=(
        "A HUMANN ChocoPhlAn nucleotide database."
    ),
)
plugin.register_artifact_class(
    HumannDatabase[TranslatedSearch],
    directory_format=HumannDatabaseDirFmt,
    description=(
        "A HUMANN translated-search database."
    ),
)
plugin.register_artifact_class(
    MetaphlanDatabase,
    directory_format=MetaphlanDatabaseDirFmt,
    description=(
        "A MetaPhlAn database."
    ),
)
plugin.register_semantic_type_to_format(
    HumannReactionTable,
    HumannReactionDirectoryFormat,
)

plugin.methods.register_function(
    function=download_chocophlan_database,
    inputs={},
    parameters={},
    outputs=[("database", HumannDatabase[ChocoPhlAn])],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        "database": (
            "The downloaded HUMANN 3 ChocoPhlAn nucleotide database."
        )
    },
    name="Download ChocoPhlAn database",
    description=(
        "Download the full HUMANN 3 ChocoPhlAn nucleotide database without "
        "updating the user's HUMANN 3 configuration."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=download_metaphlan_database,
    inputs={},
    parameters={
        "index": Str,
        "cpus": Int % Range(1, None),
    },
    outputs=[("database", MetaphlanDatabase)],
    input_descriptions={},
    parameter_descriptions={
        "index": (
            "The MetaPhlAn 3 database index to download, such as 'latest' or "
            "a specific mpa_v... identifier."
        ),
        "cpus": (
            "Number of CPUs to pass to MetaPhlAn's Bowtie2 index build "
            "step via --nproc."
        ),
    },
    output_descriptions={
        "database": "The downloaded MetaPhlAn database."
    },
    name="Download MetaPhlAn 3 database",
    description=(
        "Download a MetaPhlAn 3 database."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=download_translated_search_database,
    inputs={},
    parameters={
        "build": Str % Choices(TRANSLATED_SEARCH_DATABASE_BUILDS),
    },
    outputs=[("database", HumannDatabase[TranslatedSearch])],
    input_descriptions={},
    parameter_descriptions={
        "build": (
            "The HUMANN 3translated-search database build to download."
        )
    },
    output_descriptions={
        "database": (
            "The downloaded HUMANN 3 translated-search database."
        )
    },
    name="Download translated-search database",
    description=(
        "Download a HUMANN 3 translated-search database without updating the "
        "user's HUMANN 3 configuration."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=_run_humann,
    inputs={
        "reads": SampleData[
            SequencesWithQuality | PairedEndSequencesWithQuality
        ],
        "nucleotide_database": HumannDatabase[ChocoPhlAn],
        "translated_search_database": HumannDatabase[TranslatedSearch],
        "metaphlan_database": MetaphlanDatabase,
    },
    parameters={
        "threads": Int % Range(1, None),
        "memory_use": Str % Choices(["minimum", "maximum"]),
        "prescreen_threshold": Float % Range(0, None),
        "nucleotide_identity_threshold": Float % Range(0, 100),
        "nucleotide_query_coverage_threshold": Float % Range(0, 100),
        "nucleotide_subject_coverage_threshold": Float % Range(0, 100),
        "translated_identity_threshold": Float % Range(0, 100),
        "translated_query_coverage_threshold": Float % Range(0, 100),
        "translated_subject_coverage_threshold": Float % Range(0, 100),
        "evalue": Float % Range(0, None),
        "gap_fill": Bool,
        "minpath": Bool,
        "pathways": Str % Choices(["metacyc", "unipathway"]),
        "output_max_decimals": Int % Range(0, None),
        "log_level": Str % Choices(
            ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        ),
    },
    outputs=[
        ("gene_families", HumannGeneFamilyTable),
        ("path_abundance", HumannPathAbundanceTable),
        ("metaphlan_profile", MetaphlanMergedAbundanceTable),
        ("reactions", HumannReactionTable),
    ],
    input_descriptions={
        "reads": (
            "Single-end or paired-end demultiplexed reads to profile with "
            "HUMANN 3."
        ),
        "nucleotide_database": (
            "The staged ChocoPhlAn database artifact to use for nucleotide "
            "search."
        ),
        "translated_search_database": (
            "The staged translated-search database artifact to use for "
            "protein search."
        ),
        "metaphlan_database": (
            "The staged MetaPhlAn 3 database artifact to use during taxonomic "
            "prescreening."
        ),
    },
    parameter_descriptions={
        "threads": "Number of worker threads to pass to HUMANN.",
        "memory_use": (
            "Amount of memory HUMANN 3 should use for intermediate "
            "processing."
        ),
        "prescreen_threshold": (
            "Minimum percentage of reads that must match a species during "
            "MetaPhlAn 3 prescreening."
        ),
        "nucleotide_identity_threshold": (
            "Minimum identity threshold for nucleotide alignments."
        ),
        "nucleotide_query_coverage_threshold": (
            "Minimum query coverage threshold for nucleotide alignments."
        ),
        "nucleotide_subject_coverage_threshold": (
            "Minimum subject coverage threshold for nucleotide alignments."
        ),
        "translated_identity_threshold": (
            "Minimum identity threshold for translated alignments. Leave "
            "unset to let HUMANN 3 choose based on the translated-search "
            "database."
        ),
        "translated_query_coverage_threshold": (
            "Minimum query coverage threshold for translated alignments."
        ),
        "translated_subject_coverage_threshold": (
            "Minimum subject coverage threshold for translated alignments."
        ),
        "evalue": (
            "Maximum e-value threshold for the translated search."
        ),
        "gap_fill": "Enable HUMANN 3 pathway gap filling.",
        "minpath": "Enable HUMANN 3 MinPath pathway computation.",
        "pathways": "Pathway database to use for pathway computations.",
        "output_max_decimals": (
            "Maximum number of decimal places in HUMANN 3 output tables."
        ),
        "log_level": "HUMANN 3 log level.",
    },
    output_descriptions={
        "gene_families": "Merged HUMANN 3 gene family abundance table.",
        "path_abundance": "Merged HUMANN 3 pathway abundance table.",
        "metaphlan_profile": "Merged MetaPhlAn 3 taxonomic profile table.",
        "reactions": (
            "Reaction table derived from merged HUMANN 3 gene families."
        ),
    },
    name="Run HUMANN3 ",
    description=(
        "Run HUMANN 3 on each sample in a single-end or paired-end read "
        "artifact, then merge the resulting gene family, pathway abundance, "
        "and MetaPhlAn profile tables across samples and derive a reactions "
        "table from the merged gene families, using staged HUMANN 3 and "
        "MetaPhlAn 3 database artifacts."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=collate_gene_families,
    inputs={"tables": List[HumannGeneFamilyTable]},
    parameters={},
    outputs={"collated_table": HumannGeneFamilyTable},
    input_descriptions={
        "tables": "HUMANN 3 gene-family tables to collate."
    },
    parameter_descriptions={},
    output_descriptions={
        "collated_table": "The collated HUMANN 3 gene-family table."
    },
    name="Collate HUMANN 3 gene-family tables",
    description=(
        "Collate multiple HUMANN 3 gene-family tables into one table."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=collate_path_abundance,
    inputs={"tables": List[HumannPathAbundanceTable]},
    parameters={},
    outputs={"collated_table": HumannPathAbundanceTable},
    input_descriptions={
        "tables": "HUMANN 3 pathway-abundance tables to collate."
    },
    parameter_descriptions={},
    output_descriptions={
        "collated_table": "The collated HUMANN 3 pathway-abundance table."
    },
    name="Collate HUMANN 3 pathway-abundance tables",
    description=(
        "Collate multiple HUMANN 3 pathway-abundance tables into one table."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=collate_metaphlan_profiles,
    inputs={"tables": List[MetaphlanMergedAbundanceTable]},
    parameters={},
    outputs={"collated_table": MetaphlanMergedAbundanceTable},
    input_descriptions={
        "tables": "Merged MetaPhlAn 3 profile tables to collate."
    },
    parameter_descriptions={},
    output_descriptions={
        "collated_table": "The collated MetaPhlAn 3 profile table."
    },
    name="Collate MetaPhlAn 3 profile tables",
    description=(
        "Collate multiple merged MetaPhlAn 3 profile tables into one table."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=collate_reactions,
    inputs={"tables": List[HumannReactionTable]},
    parameters={},
    outputs={"collated_table": HumannReactionTable},
    input_descriptions={
        "tables": "HUMANN 3 reaction tables to collate."
    },
    parameter_descriptions={},
    output_descriptions={
        "collated_table": "The collated HUMANN 3 reaction table."
    },
    name="Collate HUMANN 3 reaction tables",
    description="Collate multiple HUMANN 3 reaction tables into one table.",
    citations=[citations["Beghini-etal-2021"]],
)

plugin.pipelines.register_function(
    function=run_humann,
    inputs={
        "reads": SampleData[
            SequencesWithQuality
            | JoinedSequencesWithQuality
            | PairedEndSequencesWithQuality
        ],
        "nucleotide_database": HumannDatabase[ChocoPhlAn],
        "translated_search_database": HumannDatabase[TranslatedSearch],
        "metaphlan_database": MetaphlanDatabase,
    },
    parameters={
        "threads": Int % Range(1, None),
        "memory_use": Str % Choices(["minimum", "maximum"]),
        "prescreen_threshold": Float % Range(0, None),
        "nucleotide_identity_threshold": Float % Range(0, 100),
        "nucleotide_query_coverage_threshold": Float % Range(0, 100),
        "nucleotide_subject_coverage_threshold": Float % Range(0, 100),
        "translated_identity_threshold": Float % Range(0, 100),
        "translated_query_coverage_threshold": Float % Range(0, 100),
        "translated_subject_coverage_threshold": Float % Range(0, 100),
        "evalue": Float % Range(0, None),
        "gap_fill": Bool,
        "minpath": Bool,
        "pathways": Str % Choices(["metacyc", "unipathway"]),
        "output_max_decimals": Int % Range(0, None),
        "log_level": Str % Choices(
            ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        ),
        "num_partitions": Int % Range(1, None),
    },
    outputs=[
        ("gene_families", HumannGeneFamilyTable),
        ("path_abundance", HumannPathAbundanceTable),
        ("metaphlan_profile", MetaphlanMergedAbundanceTable),
        ("reactions", HumannReactionTable),
    ],
    input_descriptions={
        "reads": (
            "Single-end, joined, or paired-end demultiplexed reads to "
            "profile with HUMANN 3."
        ),
        "nucleotide_database": (
            "The staged ChocoPhlAn database artifact to use for nucleotide "
            "search."
        ),
        "translated_search_database": (
            "The staged translated-search database artifact to use for "
            "protein search."
        ),
        "metaphlan_database": (
            "The staged MetaPhlAn database artifact to use during taxonomic "
            "prescreening."
        ),
    },
    parameter_descriptions={
        "threads": "Number of worker threads to pass to HUMANN 3.",
        "memory_use": (
            "Amount of memory HUMANN 3 should use for intermediate "
            "processing."
        ),
        "prescreen_threshold": (
            "Minimum percentage of reads that must match a species during "
            "MetaPhlAn 3 prescreening."
        ),
        "nucleotide_identity_threshold": (
            "Minimum identity threshold for nucleotide alignments."
        ),
        "nucleotide_query_coverage_threshold": (
            "Minimum query coverage threshold for nucleotide alignments."
        ),
        "nucleotide_subject_coverage_threshold": (
            "Minimum subject coverage threshold for nucleotide alignments."
        ),
        "translated_identity_threshold": (
            "Minimum identity threshold for translated alignments. Leave "
            "unset to let HUMANN 3 choose based on the translated-search "
            "database."
        ),
        "translated_query_coverage_threshold": (
            "Minimum query coverage threshold for translated alignments."
        ),
        "translated_subject_coverage_threshold": (
            "Minimum subject coverage threshold for translated alignments."
        ),
        "evalue": "Maximum e-value threshold for the translated search.",
        "gap_fill": "Enable HUMANN 3 pathway gap filling.",
        "minpath": "Enable HUMANN 3 MinPath pathway computation.",
        "pathways": "Pathway database to use for pathway computations.",
        "output_max_decimals": (
            "Maximum number of decimal places in HUMANN 3 output tables."
        ),
        "log_level": "HUMANN 3 log level.",
        "num_partitions": (
            "Number of read partitions to process before merging outputs."
        ),
    },
    output_descriptions={
        "gene_families": "Merged HUMANN 3 gene family abundance table.",
        "path_abundance": "Merged HUMANN 3 pathway abundance table.",
        "metaphlan_profile": "Merged MetaPhlAn 3 taxonomic profile table.",
        "reactions": (
            "Reaction table derived from merged HUMANN 3 gene families."
        ),
    },
    name="Run HUMANN 3",
    description=(
        "Partition a single-end, joined, or paired-end read artifact, run "
        "HUMANN 3 on each partition, then merge the resulting gene family, "
        "pathway abundance, and MetaPhlAn 3 profile tables across partitions "
        "and derive a reactions table from the merged gene families."
    ),
    citations=[citations["Beghini-etal-2021"]],
)
