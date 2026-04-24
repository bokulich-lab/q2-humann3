# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Choices, Plugin, Str
from q2_types.feature_table import FeatureTable, Frequency

from q2_humann import __version__
from q2_humann._methods import (
    TRANSLATED_SEARCH_DATABASE_BUILDS,
    download_chocophlan_database,
    download_translated_search_database,
    duplicate_table,
)
from q2_humann._transformers import humann_database_dirfmt_to_path
from q2_humann._types_and_formats import (
    ChocoPhlAn,
    HumannDatabase,
    HumannDatabaseDirFmt,
    HumannDatabaseFileFormat,
    HumannDatabaseMetadataFormat,
    TranslatedSearch,
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

plugin.register_semantic_types(HumannDatabase, ChocoPhlAn, TranslatedSearch)
plugin.register_formats(
    HumannDatabaseMetadataFormat,
    HumannDatabaseFileFormat,
    HumannDatabaseDirFmt,
)
plugin.register_artifact_class(
    HumannDatabase[ChocoPhlAn],
    directory_format=HumannDatabaseDirFmt,
    description=(
        "A staged HUMANN ChocoPhlAn nucleotide database directory."
    ),
)
plugin.register_artifact_class(
    HumannDatabase[TranslatedSearch],
    directory_format=HumannDatabaseDirFmt,
    description=(
        "A staged HUMANN translated-search database directory."
    ),
)

plugin.register_transformer(humann_database_dirfmt_to_path)

plugin.methods.register_function(
    function=duplicate_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('new_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be duplicated.'},
    parameter_descriptions={},
    output_descriptions={'new_table': 'The duplicated feature table.'},
    name='Duplicate table',
    description=("Create a copy of a feature table with a new uuid. "
                 "This is for demonstration purposes only. "),
    citations=[]
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
            "The downloaded HUMANN ChocoPhlAn nucleotide database directory."
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
            "The downloaded HUMANN translated-search database directory."
        )
    },
    name="Download translated-search database",
    description=(
        "Download a HUMANN translated-search database without updating the "
        "user's HUMANN configuration."
    ),
    citations=[citations["Beghini-etal-2021"]],
)
