# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path

from q2_humann._types_and_formats import HumannDatabaseDirFmt


def humann_database_dirfmt_to_path(dirfmt: HumannDatabaseDirFmt) -> Path:
    return dirfmt.path / "data"
