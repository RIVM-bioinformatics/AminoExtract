"""This module contains dataclasses and enums to represent GFF files and their contents."""

import gzip
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, TextIO

from AminoExtract.file_utils import FileUtils


@dataclass(frozen=True)
class GFFColumn:
    """Dataclass to represent a GFF column"""

    index: int
    name: str
    dtype: Any


@dataclass
class GFFHeader:
    """Dataclass to represent the header of a GFF file"""

    raw_text: str

    @classmethod
    def from_file(cls, file_path: Path) -> "GFFHeader":
        """Load the header of a GFF file"""
        header_lines = []
        with cls._open_file(file_path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                header_lines.append(line)
        return cls("".join(header_lines))

    @staticmethod
    def _open_file(path: Path) -> TextIO:
        """
        Open a GFF file for reading

        Notes:
        -------
        This method should be used in a context manager to ensure the file is closed properly.

        """
        return gzip.open(path, "rt", encoding="utf-8") if FileUtils.is_gzipped(path) else open(path, "rt", encoding="utf-8")

    def export_header_text(self) -> str:
        """Export the header text of the GFF file, ending with a newline"""
        if not self.raw_text:
            return ""
        if not self.raw_text.endswith("\n"):
            self.raw_text += "\n"
        return self.raw_text


class GFFColumns(Enum):
    """Enum to represent the columns of a GFF file"""

    SEQID = GFFColumn(0, "seqid", str)
    SOURCE = GFFColumn(1, "source", str)
    TYPE = GFFColumn(2, "type", str)
    START = GFFColumn(3, "start", int)
    END = GFFColumn(4, "end", int)
    SCORE = GFFColumn(5, "score", object)
    STRAND = GFFColumn(6, "strand", str)
    PHASE = GFFColumn(7, "phase", str)
    ATTRIBUTES = GFFColumn(8, "attributes", str)

    @classmethod
    def get_names(cls) -> list[str]:
        """Get the names of the columns"""
        return [col.value.name for col in cls]

    @classmethod
    def get_dtypes(cls) -> dict[str, Any]:
        """Get the data types of the columns"""
        return {col.value.name: col.value.dtype for col in cls}

    @classmethod
    def get_indices(cls) -> dict[str, int]:
        """Get the indices of the columns"""
        return {col.value.name: col.value.index for col in cls}


@dataclass
class SplicingInfo:
    """Represents the splicing information of a gene"""

    cds_locations: list[tuple[int, int]]
    parent_ids: set[str]
    gene_id: str
