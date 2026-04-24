# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin
from q2_types.feature_table import FeatureTable, Frequency

from q2_humann import __version__
from q2_humann._methods import duplicate_table

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
