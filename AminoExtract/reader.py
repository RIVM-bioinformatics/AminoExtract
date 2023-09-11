import gzip
import sys

import magic
import pandas as pd
from Bio import SeqIO

from AminoExtract.functions import log


# It reads in a GFF file and stores its contents as a GFFdataframe object
class GffDataFrame(object):
    def __init__(
        self, logger=log, inputfile: str | None = None, verbose: bool = False, split_attributes: bool = False
    ) -> None:
        None if inputfile else sys.exit("Inputfile is not provided")
        if readable_file_type(inputfile):
            self.inputfile = inputfile
            self.log = logger
            self.verbose = verbose
            self.log.info(
                f"Parsing GFF input file: '[green]{inputfile}[/green]'"
            ) if verbose else None
            self._read()
            self._read_header()
            self.df = _split_attributes_column(self.df) if split_attributes else self.df
        else:
            self.log = log
            self.verbose = verbose
            self.log.error(
                f"Input file is not readable: {inputfile}"
            ) if self.verbose else None
            sys.exit(1)

    # def split_attributes_column(self) -> pd.DataFrame:
    #     """Takes a dataframe with a column called "attributes" that contains a string of attributes, and it
    #     returns a dataframe with the attributes split into separate columns

    #     Parameters
    #     ----------
    #     df : pd.DataFrame

    #     Returns
    #     -------
    #         A dataframe with the attributes column split into individual columns.

    #     """
    #     self.df["attributes"] = self.df["attributes"].apply(_attr_string_to_dict)
    #     df = self.df.join(pd.DataFrame(self.df["attributes"].to_dict()).T).drop(
    #         "attributes", axis=1
    #     )
    #     return df

    def _read(self) -> pd.DataFrame:
        if _is_gzipped(self.inputfile):
            return self._read_gff_gzipped()
        else:
            return self._read_gff_uncompressed()

    def _read_gff_gzipped(self) -> pd.DataFrame:
        self.df = pd.read_csv(
            self.inputfile,
            sep="\t",
            comment="#",
            names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            compression="gzip",
            keep_default_na=False,
        )
        return self.df

    def _read_gff_uncompressed(self) -> pd.DataFrame:
        self.df = pd.read_csv(
            self.inputfile,
            sep="\t",
            comment="#",
            names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
            keep_default_na=False,
        )
        return self.df

    def _read_header(self):
        self.header = ""
        if _is_gzipped(self.inputfile):
            with gzip.open(self.inputfile, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        self.header += line
                    else:
                        break
        else:
            with open(self.inputfile, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        self.header += line
                    else:
                        break
        return self.header

def _split_attributes_column(df: pd.DataFrame) -> pd.DataFrame:
        """
        Takes a dataframe with a column called "attributes" that contains a string of attributes, and it
        returns a dataframe with the attributes split into separate columns.

        Parameters
        ----------
        None

        Returns
        -------
        pd.DataFrame
            A dataframe with the attributes column split into individual columns.
        """
        df["attributes"] = df["attributes"].apply(_attr_string_to_dict)
        columns = df.columns.tolist()
        # remove key-value pair from the dictionary in the attributes column if the key is already a column
        df["attributes"] = df["attributes"].apply(
            lambda attr: {k: v for k, v in attr.items() if k not in columns}
        )
        df = df.join(pd.DataFrame(df["attributes"].to_dict()).T).drop(
            "attributes", axis=1
        )
        return df

def _attr_string_to_dict(string: str) -> dict:
    """
    Takes a string like "a=1;b=2;c=3" and returns a dictionary like {"a": "1", "b": "2", "c": "3"}

    Parameters
    ----------
    string : str
        The string to parse.

    Returns
    -------
    dict
        A dictionary with the key being the attribute name and the value being the attribute value.
    """
    return dict([x.split("=") for x in string.split(";") if "=" in x])


def _is_gzipped(infile: str) -> bool:
    """
    Returns `True` if the file is gzipped, and `False` otherwise.

    Parameters
    ----------
    infile : str
        The path to the file to be checked.

    Returns
    -------
    bool
        `True` if the file is gzipped, and `False` otherwise.

    """
    return magic.from_file(infile, mime=True) == "application/gzip"


def readable_file_type(infile: str) -> bool:
    """
    Returns `True` if the file is a plain text file or a gzip file, and `False` otherwise.

    Parameters
    ----------
    infile : str
        The file to check.

    Returns
    -------
    bool
        A boolean value indicating whether the file is a plain text file or a gzip file.
    """
    return bool(
        magic.from_file(infile, mime=True) == "text/plain" or "application/gzip"
    )


def read_gff(file: str, verbose: bool = False, split_attributes: bool = False) -> GffDataFrame:
    """
    Reads a GFF file and returns a GffDataFrame object.

    Parameters
    ----------
    file : str
        The path to the GFF file to be read.
    verbose : bool, optional
        If True, print out the number of lines read in.
    split_attributes : bool, optional
        If True, split the attributes column into separate columns.

    Returns
    -------
    GffDataFrame
        A GffDataFrame object containing the data from the GFF file.
    """
    return GffDataFrame(inputfile=file, verbose=verbose, split_attributes=split_attributes)


def read_fasta(file: str, verbose: bool = False) -> list:
    """
    Reads a FASTA file and returns a list of SeqRecord objects

    Parameters
    ----------
    file : str
        The path to the FASTA file to be read.
    verbose : bool, optional
        If True, log information about the file being parsed.

    Returns
    -------
    list
        A list of SeqRecord objects representing the sequences in the input file.
    """
    log.info(f"Parsing FASTA input file: '[green]{file}[/green]'") if verbose else None
    return list(SeqIO.parse(file, "fasta"))
