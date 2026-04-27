# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_humann._types_and_formats import HumannDatabaseDirFmt


def humann_database_dirfmt_to_path(dirfmt: HumannDatabaseDirFmt) -> Path:
    return dirfmt.path / "data"


def _manifest_dirfmt_to_df(dirfmt) -> pd.DataFrame:
    manifest = pd.read_csv(dirfmt.path / "MANIFEST")
    manifest["absolute_path"] = manifest["filename"].apply(
        lambda value: str(dirfmt.path / value)
    )
    return manifest


def single_end_reads_to_manifest_df(
    dirfmt: SingleLanePerSampleSingleEndFastqDirFmt,
) -> pd.DataFrame:
    return _manifest_dirfmt_to_df(dirfmt)


def paired_end_reads_to_manifest_df(
    dirfmt: SingleLanePerSamplePairedEndFastqDirFmt,
) -> pd.DataFrame:
    return _manifest_dirfmt_to_df(dirfmt)
