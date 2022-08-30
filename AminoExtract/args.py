import argparse
import os
import pathlib
import sys

from rich import print

from AminoExtract import __prog__, __version__
from AminoExtract.functions import QuickArgFormatter, RichParser, log


def validate_output_type(args, parser):
    outtype = None
    if args.output_type.upper() in ("COMBINED", "COMB"):
        if not pathlib.Path(args.output).suffixes:
            log.error(f"The given output-type is set to [green]'combined'[/green] as you provided [cyan]'{args.output_type}'[/cyan] as the output-type.\nHowever, the given output ([red]'{args.output}'[/red]) does not have a valid file extension such as [green]'.fasta'[/green] or [green]'.fa'[/green].\nPlease try again")
            sys.exit(1)
        outtype = 0
    if args.output_type.upper() in ("SEPARATE", "SEP"):
        if pathlib.Path(args.output).suffixes:
            log.error(f"The given output-type is set to [green]'separate'[/green] as you provided [cyan]'{args.output_type}'[/cyan] as the output-type.\nHowever, the given output([red]'{args.output}'[/red]) has a file extension while we expect a this to be a folder for this output-type.\nPlease try again")
            if args.yes is True:
                outtype = 1
            else:
                sys.exit(1)
        outtype = 1
    args.output_type = outtype

    if not isinstance(args.output_type, int) and args.output_type.upper() not in ("COMBINED", "SEPARATE", "COMB", "SEP"):

        log.error(f"'[red]{args.output_type}[/red]' is not a valid output type.\n Please use either 'Combined' or 'Separate', or their shortened versions 'comb' or 'sep'.\n See '[bold]{__prog__} --help[/bold]' for more info or check out the docs.")
        sys.exit(1)


    return args


def check_file_ext(fname, choices, ftype):
    '''> Check if the file exists and has a valid extension
    
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
    
    '''
    if os.path.isfile(fname):
        file_exts = "".join(pathlib.Path(fname).suffixes)
        if file_exts not in choices:
            log.error(
                f"{fname} does not have a valid {ftype} file extension.\nPlease use any of the following extensions: [bold]{' '.join(choices)}[/bold]"
            )
            sys.exit(1)
        return pathlib.Path(fname).absolute()
    log.error(f"[green]'{fname}'[/green] is not a file.\n Exiting...")
    sys.exit(1)


def get_args(givenargs):
    """> This function takes in a list of arguments, parses them, and returns a Namespace object

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
        type=lambda s: check_file_ext(s, (".fasta", ".fas", ".fna", ".fa"), "FASTA"),
        metavar="File",
        help="Input FASTA file with nucleotide sequences.",
        required=True,
    )

    req_args.add_argument(
        "--features",
        "-gff",
        metavar="File",
        type=lambda s: check_file_ext(s, (".gff", ".gff3"), "GFF"),
        help="GFF file containing the information of the amino acid sequences to extract.",
        required=True,
    )

    req_args.add_argument(
        "--output",
        "-o",
        type=lambda s: pathlib.Path(s).absolute(),
        metavar="Path",
        help="Output path, either a file or directory.\nSpecify what type of output you want with '--output-type'",
        required=True
    )

    req_args.add_argument(
        "--name",
        "-n",
        type=str,
        metavar="Text",
        help="Name of the sample that is being processed.\n * This will be used to create the fasta headers, if 'Separate' is specified as the output type then this name will also be used to create the output files.\n [underline]Please see the docs for more info[/underline]",
        required=True,
    )

    req_args.add_argument(
        "--output-type",
        "-ot",
        type=str,
        metavar="Combined/Separate",
        help=f"Output type; Should be either 'Combined' or 'Separate'.\n * If 'Combined' is specified then {__prog__} will output all found amino acid sequences to a single fasta file.\n * The output file should end with a valid fasta file extension.\n * If 'Separate' is specified then {__prog__} will output all found amino acid sequences to separate fasta files for each sequence.\n [underline]Please see the docs for more info[/underline]",
        required=True
    )
    
    opt_args.add_argument(
        '--feature-type',
        '-ft',
        type=str,
        metavar="Text",
        default="CDS",
        help="Defines which feature types in the input gff will be processed to amino acid sequences. Defaults to 'CDS'.\nOptions are 'CDS', 'gene', and 'all",
        required=False
    )
    
    opt_args.add_argument(
        "--yes",
        "-y",
        action="store_true",
        default=False,
        help="Do not ask for confirmation"
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

    return parser.parse_args(givenargs), parser


def validate_args(givenargs):
    parsed_args, parser = get_args(givenargs)
    parsed_args = validate_output_type(parsed_args, parser)
    # print(parsed_args)
    #
    # print(pathlib.Path.is_dir(parsed_args.output))
    # print(pathlib.Path(parsed_args.output).suffixes)
    #
    return parsed_args
