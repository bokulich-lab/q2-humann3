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
