# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import json
import os
from pathlib import Path
import shutil
import tempfile
from unittest.mock import patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_humann.run import _run_humann
from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    MetaphlanDatabaseDirFmt,
)


class _ReadsDirFmt:
    def __init__(self, manifest: pd.DataFrame):
        self.manifest = manifest


class RunHumannIntegrationTests(TestPluginBase):
    package = "q2_humann.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.temp_path = Path(self.temp_dir.name)

    def _gzip_fixture(self, fixture_name: str) -> str:
        source_path = Path(self.get_data_path(fixture_name))
        destination_path = self.temp_path / f"{fixture_name}.gz"
        with source_path.open("rb") as source_fh:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=destination_path.open("wb"),
                mtime=0,
            ) as destination_fh:
                destination_fh.write(source_fh.read())
        return str(destination_path)

    def _make_reads_manifest_for_fixture(
        self, sample_id: str, fixture_name: str
    ) -> _ReadsDirFmt:
        return _ReadsDirFmt(
            pd.DataFrame(
                {
                    "forward": {
                        sample_id: self._gzip_fixture(fixture_name),
                    },
                    "reverse": {
                        sample_id: None,
                    },
                }
            )
        )

    def _make_humann_demo_database(
        self,
        package_dir: str,
        payload_dir: str,
        database_kind: str,
        build: str,
    ) -> HumannDatabaseDirFmt:
        try:
            import humann
        except ImportError:
            self.skipTest("HUMANN package is not importable.")

        source_dir = (
            Path(humann.__file__).resolve().parent / "data" / package_dir
        )
        if not source_dir.exists():
            self.skipTest(f"HUMANN demo database is missing: {source_dir}")

        database = HumannDatabaseDirFmt()
        shutil.copytree(source_dir, database.path / payload_dir)
        with (database.path / "metadata.json").open("w") as fh:
            json.dump(
                {"database_kind": database_kind, "build": build},
                fh,
            )
        return database

    def _make_toy_metaphlan_database(self) -> MetaphlanDatabaseDirFmt:
        database = MetaphlanDatabaseDirFmt()
        (database.path / "mpa_vTest.pkl").write_bytes(b"taxonomy")
        with (database.path / "metadata.json").open("w") as fh:
            json.dump(
                {"database_kind": "metaphlan", "index": "mpa_vTest"},
                fh,
            )
        return database

    def _write_fake_metaphlan_executable(self, bin_dir: Path) -> None:
        executable = bin_dir / "metaphlan"
        executable.write_text(
            """#!/usr/bin/env python
from pathlib import Path
import sys

if "--version" in sys.argv:
    print("MetaPhlAn version 3.0.0")
    raise SystemExit(0)

profile = \"\"\"#mpa_v30_CHOCOPhlAn_201901
#mpa_v30_CHOCOPhlAn_201901
#SampleID\tMetaphlan_Analysis
#clade_name\tNCBI_tax_id\trelative_abundance
k__Bacteria\t2\t100.0
g__Bacteroides|s__Bacteroides_dorei\t357276\t100.0
\"\"\"

Path(sys.argv[sys.argv.index("-o") + 1]).write_text(profile)
Path(sys.argv[sys.argv.index("--bowtie2out") + 1]).write_text("")
"""
        )
        executable.chmod(0o755)

    def test_run_humann_with_humann_demo_data(self):
        missing_executables = [
            executable
            for executable in (
                "humann",
                "humann_join_tables",
                "humann_regroup_table",
                "bowtie2",
                "diamond",
            )
            if shutil.which(executable) is None
        ]
        if missing_executables:
            self.skipTest(
                "Missing integration-test executables: "
                + ", ".join(missing_executables)
            )

        reads = self._make_reads_manifest_for_fixture(
            "demo", "humann-demo.fastq"
        )
        nucleotide_database = self._make_humann_demo_database(
            "chocophlan_DEMO",
            "chocophlan",
            "chocophlan",
            "full",
        )
        translated_search_database = self._make_humann_demo_database(
            "uniref_DEMO",
            "uniref",
            "translated-search",
            "uniref90_diamond",
        )
        metaphlan_database = self._make_toy_metaphlan_database()
        fake_bin_dir = self.temp_path / "bin"
        fake_bin_dir.mkdir()
        self._write_fake_metaphlan_executable(fake_bin_dir)

        with patch.dict(
            os.environ,
            {"PATH": f"{fake_bin_dir}{os.pathsep}{os.environ['PATH']}"},
        ):
            gene_families, path_abundance, metaphlan_profile, reactions = (
                _run_humann(
                    reads,
                    nucleotide_database,
                    translated_search_database,
                    metaphlan_database,
                )
            )

        gene_families_text = (gene_families.path / "table.tsv").read_text()
        path_abundance_text = (path_abundance.path / "table.tsv").read_text()
        metaphlan_profile_text = (
            metaphlan_profile.path / "table.tsv"
        ).read_text()
        reactions_text = (reactions.path / "table.tsv").read_text()
        self.assertIn("demo", gene_families_text)
        self.assertIn("demo", path_abundance_text)
        self.assertIn("Bacteroides_dorei", metaphlan_profile_text)
        self.assertIn("demo", reactions_text)
