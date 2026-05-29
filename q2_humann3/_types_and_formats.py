# ----------------------------------------------------------------------------
# Copyright (c) 2026, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

from qiime2.plugin import SemanticType, ValidationError, model
from q2_sapienns.plugin_setup import HumannTableFormat

HumannDatabase = SemanticType("HumannDatabase", field_names="kind")
ChocoPhlAn = SemanticType("ChocoPhlAn", variant_of=HumannDatabase.field["kind"])
TranslatedSearch = SemanticType(
    "TranslatedSearch", variant_of=HumannDatabase.field["kind"]
)

HumannReactionTable = SemanticType("HumannReactionTable")
MetaphlanDatabase = SemanticType("MetaphlanDatabase")


class HumannDatabaseMetadataFormat(model.TextFileFormat):
    _allowed_kinds = {"chocophlan", "translated-search"}

    def _validate_(self, level):
        with self.open() as fh:
            metadata = json.load(fh)

        missing = {"database_kind", "build"} - set(metadata)
        if missing:
            missing_fields = ", ".join(sorted(missing))
            raise ValidationError(
                f"Missing required HUMANN database metadata fields: "
                f"{missing_fields}."
            )

        if metadata["database_kind"] not in self._allowed_kinds:
            raise ValidationError(
                "Unsupported HUMANN database kind " f"{metadata['database_kind']!r}."
            )


class HumannDatabaseFileFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        if self.path.stat().st_size == 0:
            raise ValidationError(
                f"HUMANN database payload file {self.path.name!r} is empty."
            )


class HumannDatabaseDirFmt(model.DirectoryFormat):
    metadata = model.File("metadata.json", format=HumannDatabaseMetadataFormat)
    payload = model.FileCollection(
        r"(chocophlan|uniref)/.+\.(ffn\.gz|faa\.gz|fna\.gz|dmnd|fna)$",
        format=HumannDatabaseFileFormat,
    )

    @payload.set_path_maker
    def payload_path_maker(self, relative_path: str):
        return relative_path

    def _validate_(self, level):
        with (self.path / "metadata.json").open() as fh:
            metadata = json.load(fh)

        kind = metadata["database_kind"]
        if kind == "chocophlan":
            expected_dir = self.path / "chocophlan"
            if not expected_dir.is_dir():
                raise ValidationError(
                    "HUMANN chocophlan database must contain a 'chocophlan' "
                    "subdirectory."
                )
            if not any(expected_dir.glob("*.f*n.gz")):
                raise ValidationError(
                    "HUMANN chocophlan database must contain '.ffn.gz', "
                    "'.faa.gz', or '.fna.gz' files in the 'chocophlan' "
                    "subdirectory."
                )
        elif kind == "translated-search":
            expected_dir = self.path / "uniref"
            if not expected_dir.is_dir():
                raise ValidationError(
                    "HUMANN translated-search database must contain a "
                    "'uniref' subdirectory."
                )
            if not (
                list(expected_dir.glob("*.dmnd")) or list(expected_dir.glob("*.fna"))
            ):
                raise ValidationError(
                    "HUMANN translated-search database must contain '.dmnd' "
                    "or '.fna' files in the 'uniref' subdirectory."
                )


class MetaphlanDatabaseMetadataFormat(model.TextFileFormat):

    def _validate_(self, level):
        with self.open() as fh:
            metadata = json.load(fh)

        missing = {"database_kind", "index"} - set(metadata)
        if missing:
            missing_fields = ", ".join(sorted(missing))
            raise ValidationError(
                f"Missing required MetaPhlAn database metadata fields: "
                f"{missing_fields}."
            )

        if metadata["database_kind"] != "metaphlan":
            raise ValidationError(
                "Unsupported MetaPhlAn database kind " f"{metadata['database_kind']!r}."
            )


class MetaphlanDatabaseFileFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        if self.path.stat().st_size == 0:
            raise ValidationError(
                f"MetaPhlAn database payload file {self.path.name!r} is empty."
            )


class MetaphlanDatabaseDirFmt(model.DirectoryFormat):
    metadata = model.File("metadata.json", format=MetaphlanDatabaseMetadataFormat)
    payload = model.FileCollection(
        r".+\.(pkl|bt2|bt2l)$", format=MetaphlanDatabaseFileFormat
    )

    @payload.set_path_maker
    def payload_path_maker(self, relative_path: str):
        return relative_path

    def _validate_(self, level):
        with (self.path / "metadata.json").open() as fh:
            metadata = json.load(fh)

        index = metadata["index"]
        pkl_path = self.path / f"{index}.pkl"
        if not pkl_path.exists():
            raise ValidationError(
                f"MetaPhlAn database is missing the expected pickle file: "
                f"{index}.pkl"
            )

        bt2_ext = None
        if (self.path / f"{index}.1.bt2").exists():
            bt2_ext = "bt2"
        elif (self.path / f"{index}.1.bt2l").exists():
            bt2_ext = "bt2l"

        if bt2_ext is None:
            raise ValidationError(
                f"MetaPhlAn database is missing Bowtie2 index files for "
                f"index {index!r}."
            )

        required_suffixes = (
            f"1.{bt2_ext}",
            f"2.{bt2_ext}",
            f"3.{bt2_ext}",
            f"4.{bt2_ext}",
            f"rev.1.{bt2_ext}",
            f"rev.2.{bt2_ext}",
        )
        missing = [
            suffix
            for suffix in required_suffixes
            if not (self.path / f"{index}.{suffix}").exists()
        ]
        if missing:
            missing_str = ", ".join(missing)
            raise ValidationError(
                f"MetaPhlAn database is missing required Bowtie2 files for "
                f"index {index!r}: {missing_str}"
            )


class HumannReactionFormat(HumannTableFormat):
    # Regrouped from joined genefamilies; column names keep the RPKs suffix.
    _unit_label = "RPKs"


class MetaphlanMergedAbundanceFormat(model.TextFileFormat):

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline()
            while header_line.startswith("#"):
                header_line = fh.readline()
            n_header_fields = len(header_line.split("\t"))
            if n_header_fields < 3:
                raise ValidationError("No sample columns appear to be present.")
            for idx, line in enumerate(fh, 2):
                if n_lines is not None and idx > n_lines + 1:
                    break
                fields = line.strip().split("\t")
                n_fields = len(fields)
                if n_fields != n_header_fields:
                    raise ValidationError(
                        "Number of columns on line %d is inconsistent with "
                        "the header line." % line
                    )
                for value in fields[2:]:
                    try:
                        value = float(value)
                    except ValueError:
                        raise ValidationError(
                            "Values in table must be float-able. Found: %s" % value
                        )
                    if value > 100.0 or value < 0.0:
                        raise ValidationError(
                            "Values must be in range [0, 100]. Found: %f" % value
                        )

    def _validate_(self, level):
        level_to_n_lines = {"min": 5, "max": None}
        self._equal_number_of_columns(level_to_n_lines[level])


HumannReactionDirectoryFormat = model.SingleFileDirectoryFormat(
    "HumannReactionDirectoryFormat", "table.tsv", HumannReactionFormat
)

MetaphlanMergedAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    "MetaphlanMergedAbundanceDirectoryFormat",
    "table.tsv",
    MetaphlanMergedAbundanceFormat,
)
