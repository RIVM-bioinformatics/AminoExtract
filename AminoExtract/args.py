import argparse
import os
import pathlib
import re
import sys

from AminoExtract import __prog__, __version__
from AminoExtract.functions import QuickArgFormatter, RichParser, log


def set_output_type(args: argparse.Namespace) -> argparse.Namespace:
    args.outtype = 0

    if not pathlib.Path(args.output).suffixes:
        log.info(
            f"The given output seems to be a directory.\nAll amino acid sequences will be written to individual files in this directory.\n([cyan]{args.output}[/cyan])"
        )
        args.outtype = 1
        return args
    log.info(
        f"The given output seems to be a file.\nAll amino acid sequences will be written to this file.\n([cyan]{args.output}[/cyan])"
    )
    return args


def check_features(args: argparse.Namespace) -> argparse.Namespace:
    if args.feature_type not in ["CDS", "gene", "all"]:
        log.error(
            f"[green]'{args.feature_type}'[/green] is not a valid feature type.\nPlease use any of the following: [bold]'CDS'[/bold],[bold]'gene'[/bold], or [bold]'all'[/bold]\nThese keywords are case-sensitive."
        )
        sys.exit(1)
    return args


def check_valid_output_filename(args: argparse.Namespace) -> argparse.Namespace:
    output_ext = pathlib.PurePath(args.output).name
    if not re.match("^[\w\-. ]+$", output_ext) or "/." in str(args.output):
        log.error(
            f"'[red]{output_ext}[/red]' does not seem to be a valid filename.\nPlease use only alphanumeric characters, underscores, and dashes."
        )
        sys.exit(1)
    return args


def check_file_ext(
    fname: str | None = None, choices: list[str] | None = None, ftype: str | None = None
) -> pathlib.Path | None:
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
        ext = "".join(pathlib.Path(fname).suffixes)
        if not any(ext.endswith(c) for c in choices):
            log.error(
                f"{fname} does not have a valid {ftype} file extension.\nPlease use any of the following extensions: [bold]{' '.join(choices)}[/bold]"
            )
            sys.exit(1)
        return pathlib.Path(fname).resolve()
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
        usage=f"{__prog__} \[required options] \[optional options]",
        description=f"[bold underline]{__prog__}[/bold underline]: A quick tool to extract amino acid sequences from a FASTA file.",
        formatter_class=QuickArgFormatter,
        add_help=False,
    )

    req_args = parser.add_argument_group(
        title="[bold underline]Required Arguments[/bold underline]"
    )
    opt_args = parser.add_argument_group(
        "[bold underline]Optional Arguments[/bold underline]"
    )

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
        type=lambda s: pathlib.Path(s).absolute(),
        metavar="Path",
        help="Output path, either a [underline cyan]file[/underline cyan] or [underline magenta]directory[/underline magenta].\n * If a file path is given, then all amino acid sequences will be written to this file.\n * If a directory path is given, then each amino acid sequence will be written to a separate file in this directory.\n [underline]Please see the docs for more info[/underline]",
        required=True,
    )

    req_args.add_argument(
        "--name",
        "-n",
        type=str,
        metavar="Text",
        help="Name of the sample that is being processed.\n * This will be used to create the fasta headers when all amino acid sequences are written to a single file.\n * If the output is going to be written to individual files in an output-directory then this name will be used as a prefix to create the output files.\n [underline]Please see the docs for more info[/underline]",
        required=True,
    )

    opt_args.add_argument(
        "--feature-type",
        "-ft",
        type=str,
        metavar="Text",
        default="CDS",
        help="Defines which feature types in the input gff will be processed to amino acid sequences. Defaults to 'CDS'.\nOptions are 'CDS', 'gene', and 'all'",
        required=False,
    )

    opt_args.add_argument(
        "--keep-gaps",
        "-kg",
        action="store_true",
        default=False,
        help='If this flag is set then the amino acid translation will be done including gaps in the nucleotide sequence.\n This results in an "X" on gapped positions in the aminoacid sequence as gap characters ("-") will be replaced by "N" in the nucleotide sequence.\n  [underline]By default, gaps are removed before translation.[/underline]',
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

    # TODO: add extra arg for extra gff filters such as source, strand, etc.

    return parser.parse_args(givenargs)


def validate_args(givenargs: list[str]) -> argparse.Namespace:
    parsed_args = set_output_type(get_args(givenargs))
    parsed_args = check_features(parsed_args)
    return check_valid_output_filename(parsed_args)
