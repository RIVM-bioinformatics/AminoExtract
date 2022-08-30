__version__ = "0.1.0"
__prog__ = "AminoExtract"

import sys

from rich import print

from AminoExtract.args import validate_args
from AminoExtract.functions import log


def main():
    args = validate_args(sys.argv[1:])
    print(args)
