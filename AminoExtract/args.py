"""
Argument parsing for AminoExtract
"""

import argparse
import os
import re
import sys
from pathlib import Path, PurePath

from AminoExtract import __prog__, __version__
from AminoExtract.logging import QuickArgFormatter, RichParser, log

# When you log a message with f-strings, you always execute the f-string,
# even if the log level is not set to display the message.
# You can circumvent this by using ("lala %s", variable) instead.
# However, this is such a minor optimization that it's not worth it, imo.
# pylint: disable=logging-fstring-interpolation


def set_output_type(args: argparse.Namespace) -> argparse.Namespace:
    """Set the output type (dir or file) based on the given output path."""

    path = Path(args.output)

    if path.is_dir():
        log.info(
            "The given output seems to be a directory.\n"
            "All amino acid sequences will be written to individual files in this directory."
            f"\n([cyan]{args.output}[/cyan])"
        )
        args.outtype = 1
        return args
    log.info("The given output seems to be a file.\n" "All amino acid sequences will be written to this file.\n" f"([cyan]{args.output}[/cyan])")
    args.outtype = 0
    return args


def check_features(args: argparse.Namespace) -> argparse.Namespace:
    """Checks if the given feature type is valid."""

    if args.feature_type not in ["CDS", "gene", "all"]:
        log.error(
            f"[green]'{args.feature_type}'[/green] is not a valid feature type.\n"
            "Please use any of the following: [bold]'CDS'[/bold],[bold]'gene'[/bold],"
            " or [bold]'all'[/bold]\n"
            "These keywords are case-sensitive."
        )
        sys.exit(1)
    return args


def check_valid_output_filename(args: argparse.Namespace) -> argparse.Namespace:
    """Check if the output filename is valid."""
    output_ext = PurePath(args.output).name
    if not re.match(r"^[\w\-. ]+$", output_ext) or "/." in str(args.output):
        log.error(
            f"'[red]{output_ext}[/red]' does not seem to be a valid filename.\n" "Please use only alphanumeric characters, underscores, and dashes."
        )
        sys.exit(1)
    return args


def check_file_ext(fname: str | None = None, choices: list[str] | None = None, ftype: str | None = None) -> Path | None:
    """Check if the file exists and has a valid extension

    Parameters
    ----------
    fname
        The name of the file to check.
    choices
        A list of valid file extensions for the file type.
    ftype
        This is the type of file you're checking. It's used in the error message.

    Returns
    -------
        The absolute path of the file.

    """
    if fname is not None and os.path.isfile(fname):
        ext = "".join(Path(fname).suffixes)
        assert choices is not None, "No choices given"
        if not any(ext.endswith(c) for c in choices):
            log.error(
                f"{fname} does not have a valid {ftype} file extension.\n"
                "Please use any of the following extensions: [bold]{' '.join(choices)}[/bold]"
            )
            sys.exit(1)
        return Path(fname).resolve()
    log.error(f"[green]'{fname}'[/green] is not a file.\n Exiting...")
    sys.exit(1)


def get_args(givenargs: list[str] | None = None) -> argparse.Namespace:
    """This function takes in a list of arguments, parses them, and returns a Namespace object

    Parameters
    ----------
    givenargs
        The arguments given to the script.

    Returns
    -------
        The arguments that are being passed to the program.

    """

    parser = RichParser(
        prog=f"[bold]{__prog__}[/bold]",
        usage=rf"{__prog__} \[required options] \[optional options]",
        description=(f"[bold underline]{__prog__}[/bold underline]:" " A quick tool to extract amino acid sequences from a FASTA file."),
        formatter_class=QuickArgFormatter,
        add_help=False,
    )

    req_args = parser.add_argument_group(title="[bold underline]Required Arguments[/bold underline]")
    opt_args = parser.add_argument_group("[bold underline]Optional Arguments[/bold underline]")

    req_args.add_argument(
        "--input",
        "-i",
        type=lambda s: check_file_ext(s, [".fasta", ".fas", ".fna", ".fa"], "FASTA"),
        metavar="File",
        help="Input FASTA file with nucleotide sequences.",
        required=True,
    )

    req_args.add_argument(
        "--features",
        "-gff",
        metavar="File",
        type=lambda s: check_file_ext(s, [".gff", ".gff3"], "GFF"),
        help="GFF file containing the information of the amino acid sequences to extract.",
        required=True,
    )

    req_args.add_argument(
        "--output",
        "-o",
        type=lambda s: Path(s).absolute(),
        metavar="Path",
        help=(
            "Output path, either a [underline cyan]file[/underline cyan] or "
            "[underline magenta]directory[/underline magenta].\n"
            " * If a file path is given, then all amino acid sequences will be written to this "
            "file.\n"
            " * If a directory path is given, then each amino acid sequence will be written to a "
            "separate file in this directory.\n"
            " [underline]Please see the docs for more info[/underline]"
        ),
        required=True,
    )

    req_args.add_argument(
        "--name",
        "-n",
        type=str,
        metavar="Text",
        help=(
            "Name of the sample that is being processed.\n"
            " * This will be used to create the fasta headers when all amino acid sequences are "
            "written to a single file.\n"
            " * If the output is going to be written to individual files in an output-directory, "
            "then this name will be used as a prefix to create the output files.\n"
            " [underline]Please see the docs for more info[/underline]"
        ),
        required=True,
    )

    opt_args.add_argument(
        "--feature-type",
        "-ft",
        type=str,
        metavar="Text",
        default="CDS",
        help=(
            "Defines which feature types in the input gff will be processed to AA sequences. "
            "Defaults to 'CDS'.\n"
            "Options are 'CDS', 'gene', and 'all'"
        ),
        required=False,
    )

    opt_args.add_argument(
        "--keep-gaps",
        "-kg",
        action="store_true",
        default=False,
        help=(
            "If this flag is set then the AA translation will be done including gaps in the "
            "nucleotide sequence.\n"
            'This results in an "X" on gapped positions in the AA sequence as gap characters '
            '("-") will be replaced by "N" in the nucleotide sequence.\n'
            "[underline]By default, gaps are removed before translation.[/underline]"
        ),
        required=False,
    )

    opt_args.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help=f"Print the {__prog__} version and exit.",
    )

    opt_args.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )

    opt_args.add_argument(
        "--verbose",
        "-vb",
        action="store_true",
        default=False,
        help="Print out more information during the process.",
        required=False,
    )

    # TODO: add extra arg for extra gff filters such as source, strand, etc.

    return parser.parse_args(givenargs)


def validate_args(givenargs: list[str]) -> argparse.Namespace:
    """
    Validate the given arguments by setting the output type, checking the feature type,
    and checking the output filename.
    """
    parsed_args = set_output_type(get_args(givenargs))
    parsed_args = check_features(parsed_args)
    return check_valid_output_filename(parsed_args)
