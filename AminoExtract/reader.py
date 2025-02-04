import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from AminoExtract.file_utils import FileUtils
from AminoExtract.functions import log
from AminoExtract.gff_data import GFFColumns, GFFHeader


# Methods:
# read_fasta, read_gff
class AttributeParser:
    """Handles GFF attributes"""

    @staticmethod
    def parse_attributes(attr_string: str) -> dict[str, str]:
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

        attr_string_without_quotes = (
            attr_string.replace('"', "").replace("'", "").strip()
        )

        try:
            pairs = [
                pair.split("=") if "=" in pair else pair.split(" ")
                for pair in attr_string_without_quotes.split(";")
            ]
        except ValueError:
            raise ValueError(f"{attr_string} is not separated by '=' or ' '")

        return {key: value for key, value in pairs}

    @staticmethod
    def normalize_name_attribute(attr: str) -> str:
        r"""
        \b is a word boundary, so seperates non-word characters from word characters
        \w* matches any amount of word characters
        """
        return re.sub(r"\b\w*name\w*\b", "Name", attr, flags=re.IGNORECASE)


class GffDataFrame(object):
    def __init__(
        self,
        inputfile: str | Path,
        logger=log,
        verbose: bool = False,
    ) -> None:
        self.file_path = Path(inputfile)
        self.logger = logger
        self.verbose = verbose
        self.header: GFFHeader | None = None
        self.df: pd.DataFrame | None = None
        self.splicing_table: pd.DataFrame | None = None

        if not self._validate_input():
            sys.exit(f"Input file is not readable: {self.file_path}")

        self._load_data()

    def _validate_input(self) -> bool:
        return self.file_path.exists() and FileUtils.is_readable(self.file_path)

    def _load_data(self) -> None:
        self.header = GFFHeader.from_file(self.file_path)
        self.df = self._read_gff_data()
        self._process_attributes()

    def _read_gff_data(self) -> pd.DataFrame:
        compression = "gzip" if FileUtils.is_gzipped(self.file_path) else None
        return pd.read_csv(
            self.file_path,
            sep="\t",
            comment="#",
            names=GFFColumns.get_names(),
            dtype=GFFColumns.get_dtypes(),
            compression=compression,
            keep_default_na=False,
        )

    def _process_attributes(self) -> None:
        if self.df is None:
            return
        self.df["attributes"] = (
            self.df["attributes"]
            .apply(AttributeParser.normalize_name_attribute)
            .apply(AttributeParser.parse_attributes)
        )
        self._expand_attributes()

    def _expand_attributes(self) -> None:
        """Expand attributes into separate columns, must be called after parsing attributes"""
        if self.df is None:
            return

        existing_cols = set(self.df.columns)
        attr_df = pd.json_normalize(self.df["attributes"])

        # Only add new columns
        new_cols = [col for col in attr_df.columns if col not in existing_cols]
        if new_cols:
            self.df = pd.concat([self.df, attr_df[new_cols]], axis=1)

        self.df.drop("attributes", axis=1, inplace=True)

    def validate_dataframe(self, feature_type: str | None = None) -> bool:
        """Returns True if the dataframe is not empty and the feature type is not None"""
        if self.df is None or feature_type is None:
            log.warning(
                f"The GFF file is empty after filtering.\nThis might mean that there are no records within the GFF that match the sequence ID(s) in the given Fasta file.\nThis could also mean that there are no records within the GFF that match the feature type '[cyan]{feature_type}[/cyan]'.\nPlease check your inputs and try again."
            )
            return False
        return True


class SequenceReader:
    """Handles reading sequence and gff files"""

    def __init__(self, verbose: bool = False) -> None:
        self.verbose = verbose

    def read_gff(self, file: str) -> GffDataFrame:
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
        return GffDataFrame(inputfile=file, verbose=self.verbose)

    def read_fasta(self, file: str) -> list:
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
        if self.verbose:
            log.info(f"Parsing FASTA input file: '[green]{file}[/green]'")
        return list(SeqIO.parse(file, "fasta"))
