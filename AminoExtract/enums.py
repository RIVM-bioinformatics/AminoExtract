from dataclasses import dataclass
from enum import Enum
from typing import Any


@dataclass(frozen=True)
class GFFColumn:
    index: int
    name: str
    dtype: Any


class GFFColumns(Enum):
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
        return [col.value.name for col in cls]

    @classmethod
    def get_dtypes(cls) -> dict[str, Any]:
        return {col.value.name: col.value.dtype for col in cls}

    @classmethod
    def get_indices(cls) -> dict[str, int]:
        return {col.value.name: col.value.index for col in cls}
