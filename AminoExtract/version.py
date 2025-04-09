"""
Version information for AminoExtract.
This used to be in __init__.py, but that caused a circular import error.
"""

from pathlib import Path

__prog__ = "AminoExtract"
__version__ = "0.3.1"
ROOT_DIR = Path(__file__).resolve().parent.parent
