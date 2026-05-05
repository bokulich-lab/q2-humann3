# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Bool, Citations, Choices, Float, Int, Plugin, Range, Str
)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)

from q2_humann import __version__
from q2_humann.db import (
    TRANSLATED_SEARCH_DATABASE_BUILDS,
    download_chocophlan_database,
    download_metaphlan_database,
    download_translated_search_database,
)
from q2_humann.run import run_humann
from q2_humann._types_and_formats import (
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

citations = Citations.load("citations.bib", package="q2_humann")

plugin = Plugin(
    name="humann",
    version=__version__,
    website="https://github.com/biobakery/humann",
    package="q2_humann",
    description=(
        "A QIIME 2 plugin for wrapping HUMANN functionality for functional "
        "profiling of microbial communities."
    ),
    short_description="QIIME 2 wrapper for HUMANN.",
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
            "The downloaded HUMANN ChocoPhlAn nucleotide database."
        )
    },
    name="Download ChocoPhlAn database",
    description=(
        "Download the full HUMANN ChocoPhlAn nucleotide database without "
        "updating the user's HUMANN configuration."
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
            "The MetaPhlAn database index to download, such as 'latest' or "
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
    name="Download MetaPhlAn database",
    description=(
        "Download a MetaPhlAn database."
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
            "The HUMANN translated-search database build to download."
        )
    },
    output_descriptions={
        "database": (
            "The downloaded HUMANN translated-search database."
        )
    },
    name="Download translated-search database",
    description=(
        "Download a HUMANN translated-search database without updating the "
        "user's HUMANN configuration."
    ),
    citations=[citations["Beghini-etal-2021"]],
)

plugin.methods.register_function(
    function=run_humann,
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
            "HUMANN."
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
        "threads": "Number of worker threads to pass to HUMANN.",
        "memory_use": (
            "Amount of memory HUMANN should use for intermediate "
            "processing."
        ),
        "prescreen_threshold": (
            "Minimum percentage of reads that must match a species during "
            "MetaPhlAn prescreening."
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
            "unset to let HUMANN choose based on the translated-search "
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
        "gap_fill": "Enable HUMANN pathway gap filling.",
        "minpath": "Enable HUMANN MinPath pathway computation.",
        "pathways": "Pathway database to use for pathway computations.",
        "output_max_decimals": (
            "Maximum number of decimal places in HUMANN output tables."
        ),
        "log_level": "HUMANN log level.",
    },
    output_descriptions={
        "gene_families": "Merged HUMANN gene family abundance table.",
        "path_abundance": "Merged HUMANN pathway abundance table.",
        "metaphlan_profile": "Merged MetaPhlAn taxonomic profile table.",
        "reactions": (
            "Reaction table derived from merged HUMANN gene families."
        ),
    },
    name="Run HUMANN",
    description=(
        "Run HUMANN on each sample in a single-end or paired-end read "
        "artifact, then merge the resulting gene family, pathway abundance, "
        "and MetaPhlAn profile tables across samples and derive a reactions "
        "table from the merged gene families, using staged HUMANN and "
        "MetaPhlAn database artifacts."
    ),
    citations=[citations["Beghini-etal-2021"]],
)
