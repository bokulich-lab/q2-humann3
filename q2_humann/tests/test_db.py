# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
from pathlib import Path
import subprocess
import tempfile
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_humann.db import (
    download_chocophlan_database,
    download_metaphlan_database,
    download_translated_search_database,
    _infer_metaphlan_index,
    _validate_metaphlan_database,
)
from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    MetaphlanDatabaseDirFmt,
)


class DownloadDatabaseTests(TestPluginBase):
    package = "q2_humann.tests"

    def _mock_humann_download(
        self, expected_database: str, expected_build: str
    ):
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
            database_dir = install_dir / expected_database
            database_dir.mkdir()
            (database_dir / f"{expected_build}.1").write_bytes(b"db-data")
            return subprocess.CompletedProcess(cmd, 0, "", "")

        return _runner

    def _mock_metaphlan_download(
        self, expected_index: str, expected_cpus: int, resolved_index: str
    ):
        def _runner(cmd):
            self.assertEqual(
                cmd,
                [
                    "metaphlan",
                    "--install",
                    "--bowtie2db",
                    cmd[3],
                    "-x",
                    expected_index,
                    "--nproc",
                    str(expected_cpus),
                ],
            )
            install_dir = Path(cmd[3])
            self._write_complete_metaphlan_database(
                install_dir, resolved_index
            )
            return subprocess.CompletedProcess(cmd, 0, "", "")

        return _runner

    def _write_complete_metaphlan_database(
        self, install_dir: Path, index: str, suffix: str = "bt2"
    ) -> None:
        (install_dir / f"{index}.pkl").write_bytes(b"taxonomy")
        for required in (
            f"1.{suffix}",
            f"2.{suffix}",
            f"3.{suffix}",
            f"4.{suffix}",
            f"rev.1.{suffix}",
            f"rev.2.{suffix}",
        ):
            (install_dir / f"{index}.{required}").write_bytes(b"bowtie")

    def test_download_chocophlan_database(self):
        with patch(
            "q2_humann.db.run_humann_command",
            side_effect=self._mock_humann_download("chocophlan", "full"),
        ) as run_command:
            observed = download_chocophlan_database()

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        observed_cmd = run_command.call_args.args[0]
        self.assertEqual(
            Path(observed_cmd[observed_cmd.index("full") + 1]),
            observed.path,
        )
        self.assertTrue(
            (observed.path / "chocophlan" / "full.1").exists()
        )
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {"build": "full", "database_kind": "chocophlan"},
        )

    def test_download_translated_search_database(self):
        build = "uniref50_ec_filtered_diamond"
        with patch(
            "q2_humann.db.run_humann_command",
            side_effect=self._mock_humann_download("uniref", build),
        ) as run_command:
            observed = download_translated_search_database(build=build)

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        observed_cmd = run_command.call_args.args[0]
        self.assertEqual(
            Path(observed_cmd[observed_cmd.index(build) + 1]),
            observed.path,
        )
        self.assertTrue(
            (observed.path / "uniref" / f"{build}.1").exists()
        )
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {"build": build, "database_kind": "translated-search"},
        )

    def test_download_database_failure(self):
        with patch(
            "q2_humann.db.run_humann_command",
            side_effect=RuntimeError(
                "Command failed with exit code 23: download failed"
            ),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "download failed"
            ):
                download_translated_search_database()

    def test_download_database_empty_result(self):
        with patch(
            "q2_humann.db.run_humann_command",
            return_value=subprocess.CompletedProcess([], 0, "", ""),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "did not produce any database files"
            ):
                download_chocophlan_database()

    def test_download_metaphlan_database(self):
        with patch(
            "q2_humann.db.run_humann_command",
            side_effect=self._mock_metaphlan_download(
                expected_index="latest",
                expected_cpus=4,
                resolved_index="mpa_vJan21_CHOCOPhlAnSGB_202103",
            ),
        ) as run_command:
            observed = download_metaphlan_database(cpus=4)

        self.assertIsInstance(observed, MetaphlanDatabaseDirFmt)
        observed_cmd = run_command.call_args.args[0]
        self.assertEqual(
            Path(observed_cmd[observed_cmd.index("--bowtie2db") + 1]),
            observed.path,
        )
        self.assertTrue(
            (
                observed.path
                / "mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"
            ).exists()
        )
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {
                "database_kind": "metaphlan",
                "index": "mpa_vJan21_CHOCOPhlAnSGB_202103",
            },
        )

    def test_download_metaphlan_database_empty_result(self):
        with patch(
            "q2_humann.db.run_humann_command",
            return_value=subprocess.CompletedProcess([], 0, "", ""),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "did not produce any database files"
            ):
                download_metaphlan_database()

    def test_infer_metaphlan_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            install_dir = Path(tmpdir)
            (install_dir / "mpa_v2.pkl").write_bytes(b"db")
            (install_dir / "mpa_v1.pkl").write_bytes(b"db")

            observed = _infer_metaphlan_index(install_dir)

        self.assertEqual(observed, "mpa_v1")

    def test_infer_metaphlan_index_missing_pickle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(
                RuntimeError, "Unable to infer the MetaPhlAn database index"
            ):
                _infer_metaphlan_index(Path(tmpdir))

    def test_validate_metaphlan_database_accepts_bt2l(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            install_dir = Path(tmpdir)
            self._write_complete_metaphlan_database(
                install_dir, "mpa_vTest", suffix="bt2l"
            )

            _validate_metaphlan_database(install_dir, "mpa_vTest")

    def test_validate_metaphlan_database_missing_pickle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            install_dir = Path(tmpdir)
            (install_dir / "mpa_vTest.1.bt2").write_bytes(b"bowtie")

            with self.assertRaisesRegex(
                RuntimeError, "expected database pickle"
            ):
                _validate_metaphlan_database(install_dir, "mpa_vTest")

    def test_validate_metaphlan_database_missing_bowtie_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            install_dir = Path(tmpdir)
            (install_dir / "mpa_vTest.pkl").write_bytes(b"taxonomy")

            with self.assertRaisesRegex(
                RuntimeError, "did not produce Bowtie2 index files"
            ):
                _validate_metaphlan_database(install_dir, "mpa_vTest")

    def test_validate_metaphlan_database_incomplete_bowtie_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            install_dir = Path(tmpdir)
            (install_dir / "mpa_vTest.pkl").write_bytes(b"taxonomy")
            (install_dir / "mpa_vTest.1.bt2").write_bytes(b"bowtie")

            with self.assertRaisesRegex(
                RuntimeError, "Missing files: 2.bt2"
            ):
                _validate_metaphlan_database(install_dir, "mpa_vTest")
