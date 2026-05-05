# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import gzip
from pathlib import Path
import tempfile
from unittest.mock import patch

import pandas as pd
from q2_types.per_sample_sequences import SequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.plugin.testing import TestPluginBase

from q2_humann.run import (
    _merge_metaphlan_profiles,
    _regroup_gene_families_to_reactions,
    _run_humann,
    _stage_sample_input,
    collate_gene_families,
    collate_metaphlan_profiles,
    collate_path_abundance,
    run_humann,
)
from q2_humann._types_and_formats import (
    HumannDatabaseDirFmt,
    HumannReactionDirectoryFormat,
    MetaphlanDatabaseDirFmt,
)
from q2_sapienns.plugin_setup import (
    HumannGeneFamilyDirectoryFormat,
    HumannPathAbundanceDirectoryFormat,
    MetaphlanMergedAbundanceDirectoryFormat,
)


class _ReadsDirFmt:
    def __init__(self, manifest: pd.DataFrame):
        self.manifest = manifest


class _Artifact:
    def __init__(self, directory_format):
        self.directory_format = directory_format

    def view(self, directory_format):
        assert isinstance(self.directory_format, directory_format)
        return self.directory_format


class _ReadsArtifact:
    type = SampleData[SequencesWithQuality]


class _Context:
    def __init__(self, partition_outputs):
        self.partition_outputs = partition_outputs
        self.created_artifacts = []

    def get_action(self, plugin_name, action_name):
        if (plugin_name, action_name) == (
            "types",
            "partition_samples_single",
        ):
            return self._partition_samples
        if (plugin_name, action_name) == ("humann", "_run_humann"):
            return self._run_partition
        if (plugin_name, action_name) == (
            "humann",
            "collate_gene_families",
        ):
            return self._collate_gene_families
        if (plugin_name, action_name) == (
            "humann",
            "collate_path_abundance",
        ):
            return self._collate_path_abundance
        if (plugin_name, action_name) == (
            "humann",
            "collate_metaphlan_profiles",
        ):
            return self._collate_metaphlan_profiles
        raise AssertionError(f"Unexpected action: {plugin_name}.{action_name}")

    def make_artifact(self, semantic_type, directory_format):
        artifact = _Artifact(directory_format)
        self.created_artifacts.append((semantic_type, directory_format))
        return artifact

    def _partition_samples(self, reads, num_partitions):
        return (
            {
                f"partition-{idx}": object()
                for idx in range(num_partitions)
            },
        )

    def _run_partition(self, reads, **kwargs):
        return self.partition_outputs.pop(0)

    def _collate_gene_families(self, tables):
        return (
            _Artifact(collate_gene_families(
                [table.view(HumannGeneFamilyDirectoryFormat)
                 for table in tables]
            )),
        )

    def _collate_path_abundance(self, tables):
        return (
            _Artifact(collate_path_abundance(
                [table.view(HumannPathAbundanceDirectoryFormat)
                 for table in tables]
            )),
        )

    def _collate_metaphlan_profiles(self, tables):
        return (
            _Artifact(collate_metaphlan_profiles(
                [table.view(MetaphlanMergedAbundanceDirectoryFormat)
                 for table in tables]
            )),
        )


class RunHumannTests(TestPluginBase):
    package = "q2_humann.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.compressed_reads_dir = Path(self.temp_dir.name)

    def _gzip_fixture(self, fixture_name: str) -> str:
        source_path = Path(self.get_data_path(fixture_name))
        destination_path = self.compressed_reads_dir / f"{fixture_name}.gz"
        with source_path.open("rb") as source_fh:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=destination_path.open("wb"),
                mtime=0,
            ) as destination_fh:
                destination_fh.write(source_fh.read())
        return str(destination_path)

    def _make_reads_manifest(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "forward": {
                    "sample-a": self._gzip_fixture(
                        "sample-a-forward.fastq"
                    ),
                    "sample-b": self._gzip_fixture(
                        "sample-b-forward.fastq"
                    ),
                },
                "reverse": {
                    "sample-a": self._gzip_fixture(
                        "sample-a-reverse.fastq"
                    ),
                    "sample-b": None,
                },
            }
        )

    def _make_humann_database(
        self, database_kind: str, build: str
    ) -> HumannDatabaseDirFmt:
        database = HumannDatabaseDirFmt()
        (database.path / database_kind).mkdir()
        with (database.path / "metadata.json").open("w") as fh:
            json.dump(
                {"database_kind": database_kind, "build": build},
                fh,
            )
        return database

    def _make_metaphlan_database(self) -> MetaphlanDatabaseDirFmt:
        database = MetaphlanDatabaseDirFmt()
        (database.path / "mpa_vTest.pkl").write_bytes(b"taxonomy")
        with (database.path / "metadata.json").open("w") as fh:
            json.dump(
                {"database_kind": "metaphlan", "index": "mpa_vTest"},
                fh,
            )
        return database

    def _write_metaphlan_profile(
        self, path: Path, clade_name: str, tax_id: str, abundance: str
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "\n".join(
                [
                    "#mpa_vTest",
                    f"{clade_name}\t{tax_id}\t{abundance}",
                    "",
                ]
            )
        )

    def _make_table_artifact(self, directory_format, text: str) -> _Artifact:
        artifact_dirfmt = directory_format()
        (artifact_dirfmt.path / "table.tsv").write_text(text)
        return _Artifact(artifact_dirfmt)

    def test_stage_sample_input_reuses_single_end_reads(self):
        manifest = self._make_reads_manifest()
        sample_manifest = manifest.loc["sample-b"]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sample_work_dir = tmpdir / "sample-b"
            observed = _stage_sample_input(
                sample_manifest, "sample-b", sample_work_dir
            )

            self.assertEqual(observed, Path(sample_manifest["forward"]))
            self.assertFalse(sample_work_dir.exists())

    def test_stage_sample_input_concatenates_paired_end_reads(self):
        manifest = self._make_reads_manifest()
        sample_manifest = manifest.loc["sample-a"]

        with tempfile.TemporaryDirectory() as tmpdir:
            sample_work_dir = Path(tmpdir) / "sample-a"
            observed = _stage_sample_input(
                sample_manifest, "sample-a", sample_work_dir
            )

            expected = (
                Path(self.get_data_path("sample-a-forward.fastq")).read_text()
                + Path(
                    self.get_data_path("sample-a-reverse.fastq")
                ).read_text()
            )
            self.assertEqual(observed, sample_work_dir / "sample-a.fastq.gz")
            self.assertTrue(sample_work_dir.is_dir())
            self.assertEqual(
                gzip.decompress(observed.read_bytes()).decode(), expected
            )

    def test_stage_sample_input_fails_when_sample_has_no_reads(self):
        manifest = pd.Series({"forward": None, "reverse": None})

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(RuntimeError, "No reads found"):
                _stage_sample_input(manifest, "sample-c", Path(tmpdir))

    def test_collate_gene_families(self):
        first = HumannGeneFamilyDirectoryFormat()
        second = HumannGeneFamilyDirectoryFormat()
        (first.path / "table.tsv").write_text(
            "# Gene Family\tsample-a\nUniRef90_A\t1.0\n"
        )
        (second.path / "table.tsv").write_text(
            "# Gene Family\tsample-b\nUniRef90_B\t2.0\n"
        )

        observed = collate_gene_families([first, second])

        self.assertEqual(
            (observed.path / "table.tsv").read_text(),
            "\n".join(
                [
                    "# Gene Family\tsample-a\tsample-b",
                    "UniRef90_A\t1.0\t0.0",
                    "UniRef90_B\t0.0\t2.0",
                    "",
                ]
            ),
        )

    def test_collate_metaphlan_profiles_skips_version_header(self):
        first = MetaphlanMergedAbundanceDirectoryFormat()
        second = MetaphlanMergedAbundanceDirectoryFormat()
        (first.path / "table.tsv").write_text(
            "#mpa_vTest\n"
            "clade_name\tNCBI_tax_id\tsample-a\n"
            "k__Bacteria\t2\t42.0\n"
        )
        (second.path / "table.tsv").write_text(
            "#mpa_vTest\n"
            "clade_name\tNCBI_tax_id\tsample-b\n"
            "k__Archaea\t2157\t7.5\n"
        )

        observed = collate_metaphlan_profiles([first, second])

        self.assertEqual(
            (observed.path / "table.tsv").read_text(),
            "\n".join(
                [
                    "clade_name\tNCBI_tax_id\tsample-a\tsample-b",
                    "k__Archaea\t2157\t0.0\t7.5",
                    "k__Bacteria\t2\t42.0\t0.0",
                    "",
                ]
            ),
        )

    def test_run_humann_pipeline_merges_partition_artifacts(self):
        translated_search_database = self._make_humann_database(
            "translated-search", "uniref90_diamond"
        )
        partition_outputs = [
            (
                self._make_table_artifact(
                    HumannGeneFamilyDirectoryFormat,
                    "# Gene Family\tsample-a\nUniRef90_A\t1.0\n",
                ),
                self._make_table_artifact(
                    HumannPathAbundanceDirectoryFormat,
                    "# Pathway\tsample-a\nPWY-A\t3.0\n",
                ),
                self._make_table_artifact(
                    MetaphlanMergedAbundanceDirectoryFormat,
                    "clade_name\tNCBI_tax_id\tsample-a\n"
                    "k__Bacteria\t2\t42.0\n",
                ),
                object(),
            ),
            (
                self._make_table_artifact(
                    HumannGeneFamilyDirectoryFormat,
                    "# Gene Family\tsample-b\nUniRef90_B\t2.0\n",
                ),
                self._make_table_artifact(
                    HumannPathAbundanceDirectoryFormat,
                    "# Pathway\tsample-b\nPWY-B\t4.0\n",
                ),
                self._make_table_artifact(
                    MetaphlanMergedAbundanceDirectoryFormat,
                    "clade_name\tNCBI_tax_id\tsample-b\n"
                    "k__Archaea\t2157\t7.5\n",
                ),
                object(),
            ),
        ]
        ctx = _Context(partition_outputs)

        def _mock_regroup(input_path, database, output_path):
            self.assertTrue(input_path.exists())
            self.assertEqual(
                database.path, translated_search_database.path
            )
            output_path.write_text(
                "reaction\tsample-a\tsample-b\nR1\t1.0\t2.0\n"
            )

        with patch(
            "q2_humann.run._regroup_gene_families_to_reactions",
            side_effect=_mock_regroup,
        ):
            gene_families, path_abundance, metaphlan_profile, reactions = (
                run_humann(
                    ctx,
                    _ReadsArtifact(),
                    object(),
                    _Artifact(translated_search_database),
                    object(),
                    num_partitions=2,
                )
            )

        self.assertIn(
            "UniRef90_A\t1.0\t0.0",
            (gene_families.view(
                HumannGeneFamilyDirectoryFormat
            ).path / "table.tsv").read_text(),
        )
        self.assertIn(
            "PWY-B\t0.0\t4.0",
            (path_abundance.view(
                HumannPathAbundanceDirectoryFormat
            ).path / "table.tsv").read_text(),
        )
        self.assertIn(
            "k__Archaea\t2157\t0.0\t7.5",
            (metaphlan_profile.view(
                MetaphlanMergedAbundanceDirectoryFormat
            ).path / "table.tsv").read_text(),
        )
        self.assertIn(
            "R1\t1.0\t2.0",
            (reactions.view(
                HumannReactionDirectoryFormat
            ).path / "table.tsv").read_text(),
        )

    def test_merge_metaphlan_profiles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            run_output_dir = tmpdir / "humann-output"
            self._write_metaphlan_profile(
                run_output_dir
                / "sample-a_humann_temp"
                / "sample-a_metaphlan_bugs_list.tsv",
                "k__Bacteria",
                "2",
                "42.0",
            )
            self._write_metaphlan_profile(
                run_output_dir
                / "sample-b_humann_temp"
                / "sample-b_metaphlan_bugs_list.tsv",
                "k__Archaea",
                "2157",
                "7.5",
            )
            output_path = tmpdir / "merged.tsv"

            def _mock_run_command(cmd):
                Path(cmd[cmd.index("-o") + 1]).write_text(
                    "#mpa_vTest_CHOCOPhlAn_202401\n"
                    "clade_name\tNCBI_tax_id\tsample-a\tsample-b\n"
                )

            with patch(
                "q2_humann.run.run_humann_command",
                side_effect=_mock_run_command,
            ) as run_command:
                _merge_metaphlan_profiles(run_output_dir, output_path)

            run_command.assert_called_once_with(
                [
                    "merge_metaphlan_tables.py",
                    str(
                        run_output_dir
                        / "sample-a_humann_temp"
                        / "sample-a_metaphlan_bugs_list.tsv"
                    ),
                    str(
                        run_output_dir
                        / "sample-b_humann_temp"
                        / "sample-b_metaphlan_bugs_list.tsv"
                    ),
                    "-o",
                    str(output_path),
                ]
            )

    def test_merge_metaphlan_profiles_fails_when_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            run_output_dir = tmpdir / "humann-output"
            run_output_dir.mkdir()

            with self.assertRaisesRegex(
                RuntimeError, "No MetaPhlAn bugs-list files"
            ):
                _merge_metaphlan_profiles(
                    run_output_dir, tmpdir / "merged.tsv"
                )

    def test_regroup_gene_families_to_reactions_uses_uniref50_groups(self):
        database = self._make_humann_database(
            "translated-search", "uniref50_diamond"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "genefamilies.tsv"
            output_path = tmpdir / "reactions.tsv"
            input_path.write_text("feature\tsample-a\nUniRef50_A\t1.0\n")

            with patch("q2_humann.run.run_humann_command") as run_command:
                _regroup_gene_families_to_reactions(
                    input_path, database, output_path
                )

        run_command.assert_called_once_with(
            [
                "humann_regroup_table",
                "--input",
                str(input_path),
                "--groups",
                "uniref50_rxn",
                "--output",
                str(output_path),
            ]
        )

    def test_regroup_gene_families_to_reactions_fails_on_unknown_build(self):
        database = self._make_humann_database(
            "translated-search", "custom_diamond"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(
                RuntimeError, "Unsupported translated-search"
            ):
                _regroup_gene_families_to_reactions(
                    Path(tmpdir) / "genefamilies.tsv",
                    database,
                    Path(tmpdir) / "reactions.tsv",
                )

    def test_run_humann(self):
        reads = _ReadsDirFmt(self._make_reads_manifest())
        nucleotide_database = self._make_humann_database(
            "chocophlan", "full"
        )
        translated_search_database = self._make_humann_database(
            "translated-search", "uniref90_diamond"
        )
        metaphlan_database = self._make_metaphlan_database()

        def _mock_run_command(cmd):
            if cmd[0] == "humann":
                run_output_dir = Path(cmd[cmd.index("--output") + 1])
                sample_id = cmd[cmd.index("--output-basename") + 1]
                self._write_metaphlan_profile(
                    run_output_dir
                    / f"{sample_id}_humann_temp"
                    / f"{sample_id}_metaphlan_bugs_list.tsv",
                    "k__Bacteria",
                    "2",
                    "42.0",
                )
            elif cmd[0] == "merge_metaphlan_tables.py":
                output_path = Path(cmd[cmd.index("-o") + 1])
                output_path.write_text(
                    "clade_name\tNCBI_tax_id\tsample-a\tsample-b\n"
                    "k__Bacteria\t2\t42.0\t42.0\n"
                )
            elif cmd[0] == "humann_join_tables":
                output_path = Path(cmd[cmd.index("--output") + 1])
                table_name = cmd[cmd.index("--file_name") + 1]
                output_path.write_text(
                    f"{table_name}\tsample-a\tsample-b\nfeature\t1\t2\n"
                )
            elif cmd[0] == "humann_regroup_table":
                output_path = Path(cmd[cmd.index("--output") + 1])
                output_path.write_text(
                    "reaction\tsample-a\tsample-b\nR1\t1\t2\n"
                )
            else:
                raise AssertionError(f"Unexpected command: {cmd!r}")

        with patch(
            "q2_humann.run.run_humann_command",
            side_effect=_mock_run_command,
        ) as run_command:
            observed = _run_humann(
                reads,
                nucleotide_database,
                translated_search_database,
                metaphlan_database,
                threads=3,
                memory_use="maximum",
                prescreen_threshold=0.02,
                nucleotide_identity_threshold=95.0,
                nucleotide_query_coverage_threshold=91.0,
                nucleotide_subject_coverage_threshold=51.0,
                translated_identity_threshold=80.0,
                translated_query_coverage_threshold=92.0,
                translated_subject_coverage_threshold=52.0,
                evalue=0.001,
                gap_fill=False,
                minpath=False,
                pathways="unipathway",
                output_max_decimals=5,
                log_level="INFO",
            )

        gene_families, path_abundance, metaphlan_profile, reactions = observed
        self.assertIsInstance(
            gene_families, HumannGeneFamilyDirectoryFormat
        )
        self.assertIsInstance(
            path_abundance, HumannPathAbundanceDirectoryFormat
        )
        self.assertIsInstance(
            metaphlan_profile, MetaphlanMergedAbundanceDirectoryFormat
        )
        self.assertIsInstance(reactions, HumannReactionDirectoryFormat)
        self.assertTrue((gene_families.path / "table.tsv").exists())
        self.assertTrue((path_abundance.path / "table.tsv").exists())
        self.assertTrue((metaphlan_profile.path / "table.tsv").exists())
        self.assertTrue((reactions.path / "table.tsv").exists())
        self.assertEqual(run_command.call_count, 6)
        first_humann_cmd = run_command.call_args_list[0].args[0]
        self.assertEqual(first_humann_cmd[0], "humann")
        self.assertEqual(
            first_humann_cmd[
                first_humann_cmd.index("--output-basename") + 1
            ],
            "sample-a",
        )
        self.assertIn(
            f"--bowtie2db {metaphlan_database.path} -x mpa_vTest",
            first_humann_cmd,
        )
        expected_options = {
            "--memory-use": "maximum",
            "--prescreen-threshold": "0.02",
            "--nucleotide-identity-threshold": "95.0",
            "--nucleotide-query-coverage-threshold": "91.0",
            "--nucleotide-subject-coverage-threshold": "51.0",
            "--translated-identity-threshold": "80.0",
            "--translated-query-coverage-threshold": "92.0",
            "--translated-subject-coverage-threshold": "52.0",
            "--evalue": "0.001",
            "--gap-fill": "off",
            "--minpath": "off",
            "--pathways": "unipathway",
            "--output-max-decimals": "5",
            "--log-level": "INFO",
        }
        for option, value in expected_options.items():
            self.assertEqual(
                first_humann_cmd[first_humann_cmd.index(option) + 1],
                value,
            )
        merge_cmd = run_command.call_args_list[4].args[0]
        self.assertEqual(merge_cmd[0], "merge_metaphlan_tables.py")
        self.assertEqual(merge_cmd[-2], "-o")
        self.assertEqual(Path(merge_cmd[-1]).name, "table.tsv")
        self.assertTrue(
            merge_cmd[1].endswith("sample-a_metaphlan_bugs_list.tsv")
        )
        self.assertTrue(
            merge_cmd[2].endswith("sample-b_metaphlan_bugs_list.tsv")
        )
        regroup_cmd = run_command.call_args_list[5].args[0]
        self.assertEqual(regroup_cmd[0], "humann_regroup_table")
        self.assertEqual(
            regroup_cmd[regroup_cmd.index("--groups") + 1],
            "uniref90_rxn",
        )
