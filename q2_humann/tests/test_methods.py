# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
import pandas as pd
import pandas.testing as pdt
import subprocess
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.feature_table import BIOMV100Format

from q2_humann._methods import (
    download_chocophlan_database,
    download_translated_search_database,
    duplicate_table,
)
from q2_humann._types_and_formats import HumannDatabaseDirFmt


class DuplicateTableTests(TestPluginBase):
    package = 'q2_humann.tests'

    def test_simple1(self):
        in_table = pd.DataFrame(
            [[1, 2, 3, 4, 5], [9, 10, 11, 12, 13]],
            columns=['abc', 'def', 'jkl', 'mno', 'pqr'],
            index=['sample-1', 'sample-2'])
        observed = duplicate_table(in_table)

        expected = in_table

        pdt.assert_frame_equal(observed, expected)

    def test_simple2(self):
        # test table duplication with table loaded from file this time
        # (for demonstration purposes)
        in_table = transform(
            self.get_data_path('table-1.biom'),
            from_type=BIOMV100Format,
            to_type=pd.DataFrame)
        observed = duplicate_table(in_table)

        expected = in_table

        pdt.assert_frame_equal(observed, expected)


class DownloadDatabaseTests(TestPluginBase):
    package = "q2_humann.tests"

    def _mock_download(self, expected_database: str, expected_build: str):
        def _runner(cmd):
            self.assertEqual(
                cmd[:5],
                [
                    "humann_databases",
                    "--download",
                    expected_database,
                    expected_build,
                    cmd[4],
                ],
            )
            self.assertEqual(cmd[5:], ["--update-config", "no"])
            install_dir = Path(cmd[4])
            (install_dir / f"{expected_build}.1").write_bytes(b"db-data")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        return _runner

    def test_download_chocophlan_database(self):
        with patch(
            "q2_humann._methods.run_command",
            side_effect=self._mock_download("chocophlan", "full"),
        ):
            observed = download_chocophlan_database()

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        self.assertTrue((observed.path / "data" / "full.1").exists())
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {"build": "full", "database_kind": "chocophlan"},
        )

    def test_download_translated_search_database(self):
        build = "uniref50_ec_filtered_diamond"
        with patch(
            "q2_humann._methods.run_command",
            side_effect=self._mock_download("uniref", build),
        ):
            observed = download_translated_search_database(build=build)

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        self.assertTrue((observed.path / "data" / f"{build}.1").exists())
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {"build": build, "database_kind": "translated-search"},
        )

    def test_download_database_failure(self):
        with patch(
            "q2_humann._methods.run_command",
            side_effect=RuntimeError("Command failed with exit code 23: download failed"),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "download failed"
            ):
                download_translated_search_database()
