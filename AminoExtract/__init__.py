"""
AminoExtract: Extract amino acid sequences from GFF annotations.
This package provides functionality to parse GFF3 and FASTA files and extract amino
acid sequences based on CDS features from genomic annotations.

Dependencies
-----------
- pandas
- biopython
- rich
- python-magic


See Also
--------
Documentation: https://github.com/RIVM-bioinformatics/AminoExtract

Notes
-----
License: MIT
"""

# pylint: disable=invalid-name
import contextlib
from pathlib import Path

__prog__ = "AminoExtract"
__version__ = "0.4.0"

ROOT_DIR = Path(__file__).resolve().parent.parent

# use contextlib to suppress the ImportError
# This may occur when this file is imported in downstream modules which may technically cause a circular import
with contextlib.suppress(ImportError):
    from AminoExtract.__main__ import (
        AminoAcidExtractor,
        get_feature_name_attribute,
        main,
    )
    from AminoExtract.filter import GFFRecordFilter
    from AminoExtract.reader import GFFDataFrame, SequenceReader
    from AminoExtract.sequences import SequenceExtractor
