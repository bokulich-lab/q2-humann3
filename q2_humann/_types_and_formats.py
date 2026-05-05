# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

from qiime2.plugin import SemanticType, ValidationError, model
from q2_sapienns.plugin_setup import HumannTableFormat

HumannDatabase = SemanticType("HumannDatabase", field_names="kind")
ChocoPhlAn = SemanticType(
    "ChocoPhlAn", variant_of=HumannDatabase.field["kind"]
)
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
                "Unsupported HUMANN database kind "
                f"{metadata['database_kind']!r}."
            )


class HumannDatabaseFileFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        if self.path.stat().st_size == 0:
            raise ValidationError(
                f"HUMANN database payload file {self.path.name!r} is empty."
            )


class HumannDatabaseDirFmt(model.DirectoryFormat):
    metadata = model.File(
        "metadata.json", format=HumannDatabaseMetadataFormat
    )
    payload = model.FileCollection(
        r"(?!metadata\.json$).+", format=HumannDatabaseFileFormat
    )

    @payload.set_path_maker
    def payload_path_maker(self, relative_path: str):
        return relative_path


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
                "Unsupported MetaPhlAn database kind "
                f"{metadata['database_kind']!r}."
            )


class MetaphlanDatabaseFileFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        if self.path.stat().st_size == 0:
            raise ValidationError(
                f"MetaPhlAn database payload file {self.path.name!r} is empty."
            )


class MetaphlanDatabaseDirFmt(model.DirectoryFormat):
    metadata = model.File(
        "metadata.json", format=MetaphlanDatabaseMetadataFormat
    )
    payload = model.FileCollection(
        r"(?!metadata\.json$).+", format=MetaphlanDatabaseFileFormat
    )

    @payload.set_path_maker
    def payload_path_maker(self, relative_path: str):
        return relative_path


class HumannReactionFormat(HumannTableFormat):
    # Regrouped from joined genefamilies; column names keep the RPKs suffix.
    _unit_label = 'RPKs'


class MetaphlanMergedAbundanceFormat(model.TextFileFormat):

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline()
            while header_line.startswith('#'):
                header_line = fh.readline()
            n_header_fields = len(header_line.split('\t'))
            if n_header_fields < 3:
                raise ValidationError('No sample columns appear to be present.')
            for idx, line in enumerate(fh, 2):
                if n_lines is not None and idx > n_lines + 1:
                    break
                fields = line.strip().split('\t')
                n_fields = len(fields)
                if n_fields != n_header_fields:
                    raise ValidationError(
                        'Number of columns on line %d is inconsistent with '
                        'the header line.' % line)
                for value in fields[2:]:
                    try:
                        value = float(value)
                    except ValueError:
                        raise ValidationError(
                            'Values in table must be float-able. Found: %s' %
                            value
                        )
                    if value > 100.0 or value < 0.0:
                        raise ValidationError(
                            'Values must be in range [0, 100]. Found: %f' %
                            value
                        )

    def _validate_(self, level):
        level_to_n_lines = {'min': 5, 'max': None}
        self._equal_number_of_columns(level_to_n_lines[level])


HumannReactionDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannReactionDirectoryFormat', 'table.tsv',
    HumannReactionFormat
)

MetaphlanMergedAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    'MetaphlanMergedAbundanceDirectoryFormat', 'table.tsv',
    MetaphlanMergedAbundanceFormat
)
