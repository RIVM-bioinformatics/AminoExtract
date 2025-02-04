import gzip
import re
import sys

import magic
import pandas as pd
from Bio import SeqIO

from AminoExtract.enums import GFFColumns
from AminoExtract.functions import log


# It reads in a GFF file and stores its contents as a GFFdataframe object
class GffDataFrame(object):
    def __init__(
        self,
        logger=log,
        inputfile: str | None = None,
        verbose: bool = False,
    ) -> None:
        if not inputfile:
            sys.exit("Inputfile is not provided")
        self.df: pd.DataFrame
        self.splicing_table: pd.DataFrame | None = None

        if readable_file_type(inputfile):
            self.inputfile = inputfile
            self.log = logger
            self.verbose = verbose
            if verbose:
                self.log.info(f"Parsing GFF input file: '[green]{inputfile}[/green]'")
            self._read()
            self._normalize_attributes()
            self._read_header()
            self.df = _split_attributes_column(self.df)
        else:
            self.log = log
            self.verbose = verbose
            if self.verbose:
                self.log.error(f"Input file is not readable: {inputfile}")
            sys.exit(1)

    def _read(self) -> pd.DataFrame:
        if _is_gzipped(self.inputfile):
            return self._read_gff(gzipped=True)
        return self._read_gff(gzipped=False)

    def _normalize_attributes(self) -> None:
        """
        Normalize the attributes in the GFF data frame.

        This function processes the 'attributes' column of the data frame,
        ensuring that any attribute containing the word 'name' is normalized
        to 'Name'.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        def _normalize(attr: str) -> str:
            r"""
            \b is a word boundary, so seperates non-word characters from word characters
            \w* matches any amount of word characters
            """
            return re.sub(r"\b\w*name\w*\b", "Name", attr, flags=re.IGNORECASE)

        self.df["attributes"] = self.df["attributes"].apply(_normalize)

    def _read_gff(self, gzipped: bool) -> pd.DataFrame:
        self.df = pd.read_csv(
            self.inputfile,
            sep="\t",
            comment="#",
            names=GFFColumns.get_names(),
            dtype=GFFColumns.get_dtypes(),
            compression="gzip" if gzipped else None,
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
    df = df.join(pd.DataFrame(df["attributes"].to_dict()).T).drop("attributes", axis=1)
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

    def _splitter(string: str, key: str) -> list[str]:
        return string.replace('"', "").replace("'", "").split(key)

    res_list = []
    for x in string.split(";"):
        if "=" in x:
            res_list.append(_splitter(x, "="))
        elif " " in x:
            res_list.append(_splitter(x, " "))
        else:
            raise ValueError("Attributes are not separated by '=' or ' '")
    return dict(res_list)


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


def read_gff(file: str, verbose: bool = False) -> GffDataFrame:
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
    return GffDataFrame(inputfile=file, verbose=verbose)


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
    if verbose:
        log.info(f"Parsing FASTA input file: '[green]{file}[/green]'")
    return list(SeqIO.parse(file, "fasta"))
