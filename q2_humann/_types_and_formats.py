# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path

import json

from qiime2.plugin import SemanticType, ValidationError, model

HumannDatabase = SemanticType("HumannDatabase", field_names="kind")
ChocoPhlAn = SemanticType(
    "ChocoPhlAn", variant_of=HumannDatabase.field["kind"]
)
TranslatedSearch = SemanticType(
    "TranslatedSearch", variant_of=HumannDatabase.field["kind"]
)
HumannPathAbundanceTable = SemanticType("HumannPathAbundanceTable")
HumannGeneFamilyTable = SemanticType("HumannGeneFamilyTable")
HumannReactionTable = SemanticType("HumannReactionTable")
MetaphlanMergedAbundanceTable = SemanticType("MetaphlanMergedAbundanceTable")


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
        r"data/.+", format=HumannDatabaseFileFormat
    )

    @payload.set_path_maker
    def payload_path_maker(self, relative_path: str):
        return str(Path("data") / relative_path)


class HumannTableFormat(model.TextFileFormat):
    _unit_label = None

    def _equal_number_of_columns(self, n_lines):
        with self.open() as fh:
            header_line = fh.readline().strip()
            header_fields = header_line.split('\t')
            n_header_fields = len(header_fields)
            if n_header_fields < 2:
                raise ValidationError('No sample columns appear to be present.')
            for sample_id in header_fields[1:]:
                if not sample_id.endswith(self._unit_label):
                    raise ValidationError(
                        'Expected sample ids (e.g., %s) to end with unit '
                        'descriptor %s' % (sample_id, self._unit_label)
                    )
            for idx, line in enumerate(fh, 2):
                if n_lines is not None and idx > n_lines + 1:
                    break
                n_fields = len(line.split('\t'))
                if n_fields != n_header_fields:
                    raise ValidationError(
                        'Number of columns on line %d is inconsistent with '
                        'the header line.' % line)

    def _validate_(self, level):
        level_to_n_lines = {'min': 5, 'max': None}
        self._equal_number_of_columns(level_to_n_lines[level])


class HumannPathAbundanceFormat(HumannTableFormat):
    _unit_label = 'Abundance'


class HumannGeneFamilyFormat(HumannTableFormat):
    _unit_label = 'RPKs'


class HumannReactionFormat(HumannTableFormat):
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


HumannPathAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannPathAbundanceDirectoryFormat', 'table.tsv',
    HumannPathAbundanceFormat
)

HumannGeneFamilyDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannGeneFamilyDirectoryFormat', 'table.tsv',
    HumannGeneFamilyFormat
)

HumannReactionDirectoryFormat = model.SingleFileDirectoryFormat(
    'HumannReactionDirectoryFormat', 'table.tsv',
    HumannReactionFormat
)

MetaphlanMergedAbundanceDirectoryFormat = model.SingleFileDirectoryFormat(
    'MetaphlanMergedAbundanceDirectoryFormat', 'table.tsv',
    MetaphlanMergedAbundanceFormat
)
