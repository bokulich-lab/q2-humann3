# ----------------------------------------------------------------------------
# Copyright (c) 2026, Bokulich Lab.
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
        def _side_effect(cmd):
            (Path(cmd[4]) / "chocophlan").mkdir()
            (Path(cmd[4]) / "chocophlan" / "full.ffn.gz").touch()

        with patch(
            "q2_humann.db.run_humann_command", side_effect=_side_effect
        ) as run_command:
            observed = download_chocophlan_database()

        run_command.assert_called_once_with([
            "humann_databases", "--download", "chocophlan", "full",
            str(observed.path), "--update-config", "no"
        ])

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        self.assertTrue(
            (observed.path / "chocophlan" / "full.ffn.gz").exists()
        )
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {"build": "full", "database_kind": "chocophlan"},
        )

    def test_download_translated_search_database(self):
        build = "uniref50_ec_filtered_diamond"

        def _side_effect(cmd):
            (Path(cmd[4]) / "uniref").mkdir()
            (Path(cmd[4]) / "uniref" / f"{build}.dmnd").touch()

        with patch(
            "q2_humann.db.run_humann_command", side_effect=_side_effect
        ) as run_command:
            observed = download_translated_search_database(build=build)

        run_command.assert_called_once_with([
            "humann_databases", "--download", "uniref", build,
            str(observed.path), "--update-config", "no"
        ])

        self.assertIsInstance(observed, HumannDatabaseDirFmt)
        self.assertTrue(
            (observed.path / "uniref" / f"{build}.dmnd").exists()
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
        index = "mpa_vJan21_CHOCOPhlAnSGB_202103"

        def _side_effect(cmd):
            self._write_complete_metaphlan_database(Path(cmd[3]), index)

        with patch(
            "q2_humann.db.run_humann_command", side_effect=_side_effect
        ) as run_command:
            observed = download_metaphlan_database(cpus=4)

        run_command.assert_called_once_with([
            "metaphlan", "--install", "--bowtie2db", str(observed.path),
            "-x", "latest", "--nproc", "4"
        ])

        self.assertIsInstance(observed, MetaphlanDatabaseDirFmt)
        self.assertTrue((observed.path / f"{index}.pkl").exists())
        with open(observed.path / "metadata.json") as fh:
            metadata = json.load(fh)
        self.assertEqual(
            metadata,
            {
                "database_kind": "metaphlan",
                "index": index,
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

    def test_download_chocophlan_database_prunes_unexpected(self):
        def _side_effect(cmd):
            install_dir = Path(cmd[4])
            (install_dir / "chocophlan").mkdir()
            (install_dir / "chocophlan" / "full.ffn.gz").touch()
            (install_dir / "chocophlan" / "junk.txt").touch()
            (install_dir / "junk_dir").mkdir()
            (install_dir / "root_junk.txt").touch()

        with patch(
            "q2_humann.db.run_humann_command", side_effect=_side_effect
        ):
            observed = download_chocophlan_database()

        self.assertTrue(
            (observed.path / "chocophlan" / "full.ffn.gz").exists()
        )
        self.assertFalse(
            (observed.path / "chocophlan" / "junk.txt").exists()
        )
        self.assertFalse((observed.path / "junk_dir").exists())
        self.assertFalse((observed.path / "root_junk.txt").exists())

    def test_download_metaphlan_database_prunes_unexpected(self):
        index = "mpa_vJan21_CHOCOPhlAnSGB_202103"

        def _side_effect(cmd):
            install_dir = Path(cmd[3])
            self._write_complete_metaphlan_database(install_dir, index)
            (install_dir / "other_index.pkl").touch()
            (install_dir / "metaphlan_log.txt").touch()

        with patch(
            "q2_humann.db.run_humann_command", side_effect=_side_effect
        ):
            observed = download_metaphlan_database()

        self.assertTrue((observed.path / f"{index}.pkl").exists())
        self.assertFalse((observed.path / "other_index.pkl").exists())
        self.assertFalse((observed.path / "metaphlan_log.txt").exists())

    def test_humann_database_validation(self):
        # Valid ChocoPhlAn
        db = HumannDatabaseDirFmt()
        (db.path / "chocophlan").mkdir()
        (db.path / "chocophlan" / "data.ffn.gz").write_bytes(b"data")
        with (db.path / "metadata.json").open("w") as fh:
            json.dump({"database_kind": "chocophlan", "build": "full"}, fh)
        db.validate()

        # Valid Translated Search
        db = HumannDatabaseDirFmt()
        (db.path / "uniref").mkdir()
        (db.path / "uniref" / "data.dmnd").write_bytes(b"data")
        with (db.path / "metadata.json").open("w") as fh:
            json.dump(
                {"database_kind": "translated-search", "build": "uniref90"}, fh
            )
        db.validate()

    def test_metaphlan_database_validation(self):
        db = MetaphlanDatabaseDirFmt()
        index = "mpa_vTest"
        (db.path / f"{index}.pkl").write_bytes(b"taxonomy")
        for suffix in (
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ):
            (db.path / f"{index}.{suffix}").write_bytes(b"bowtie")

        with (db.path / "metadata.json").open("w") as fh:
            json.dump({"database_kind": "metaphlan", "index": index}, fh)
        db.validate()


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
