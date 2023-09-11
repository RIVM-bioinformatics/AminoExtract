import contextlib

__prog__ = "AminoExtract"
__version__ = "0.3.0"

# use contextlib to suppress the ImportError
# This may occur when this file is imported in setup.py as dependencies are not yet installed
with contextlib.suppress(ImportError):
    from AminoExtract.__main__ import get_feature_name_attribute, main
    from AminoExtract.filter import filter_gff
    from AminoExtract.reader import GffDataFrame, read_gff, read_fasta
    from AminoExtract.sequences import extract_aminoacids
