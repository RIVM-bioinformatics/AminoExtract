import contextlib

__prog__ = "AminoExtract"
__version__ = "0.1.0"

# use contextlib to suppress the ImportError that may occur when this file is imported in setup.py as dependencies are not yet installed
with contextlib.suppress(ImportError):
    from AminoExtract.__main__ import get_feature_name_attribute, main
