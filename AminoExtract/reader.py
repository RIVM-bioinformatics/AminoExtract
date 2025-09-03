"""Module for reading GFF and FASTA files."""

import re
import sys
from logging import Logger
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SeqRecord

from AminoExtract.file_utils import FileUtils
from AminoExtract.gff_data import GFFColumns, GFFHeader, SplicingInfo
from AminoExtract.logging import log


class AttributeParser:
    """Handles GFF attributes"""

    @staticmethod
    def parse_attributes(attr_string: str) -> dict[str, str]:
        """
        Takes a string like "a=1;b=2;c=3",
        and returns a dictionary like {"a": "1", "b": "2", "c": "3"}

        Parameters
        ----------
        string : str
            The string to parse.

        Returns
        -------
        dict
            A dictionary with the key being the attribute name,
            and the value being the attribute value.
        """

        attr_string_without_quotes = attr_string.replace('"', "").replace("'", "").strip()

        try:
            pairs = [pair.split("=") if "=" in pair else pair.split(" ") for pair in attr_string_without_quotes.split(";")]
            for x in pairs:
                if len(x) == 1:
                    pairs.remove(x)
        except ValueError as e:
            raise ValueError(f"{attr_string} is not separated by '=' or ' '") from e

        return dict(pairs)

    @staticmethod
    def normalize_name_attribute(attr: str) -> str:
        """
        Checks if there is a 'Name' like attribute in the GFF attributes.
        If there is, it normalizes it to 'Name'.
        If there is not, copy the 'gene_name' like attribute to 'Name'.
        If there is also not a gene_name attribute, copy the 'id' like attribute to 'Name'.
        If there is still nothing, copy the 'seqid' like attribute to 'Name'.
        Else leave empty, it will turn in pd.NA later
        """
        if re.search(r"\bname\b", attr, flags=re.IGNORECASE):
            attr = re.sub(r"\bname\b", "Name", attr, flags=re.IGNORECASE)
        elif match := re.search(r"\bgene_name=([^;]+)", attr, flags=re.IGNORECASE):
            attr += f";Name={match.group(1)}"
        elif match := re.search(r"\bid=([^;]+)", attr, flags=re.IGNORECASE):
            attr += f";Name={match.group(1)}"
        elif match := re.search(r"\bseqid=([^;]+)", attr, flags=re.IGNORECASE):
            attr += f";Name={match.group(1)}"
        return attr


class GFFDataFrame:
    """Class for reading and processing GFF files."""

    # This class does all its work in the constructor, so maybe it could be a function,
    # but I like the idea of having a class that represents a GFF file
    # pylint: disable=too-few-public-methods

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
        self.splicing_info: list[SplicingInfo] | None = None

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
        # compression must be in dict format
        # https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
        compression = {"method": "gzip"} if FileUtils.is_gzipped(self.file_path) else None
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
        self.df["attributes"] = self.df["attributes"].apply(AttributeParser.normalize_name_attribute)

        # parsing and expanding make significant changes to the attributes column,
        # this makes it unable to be written back to the original GFF format
        self._expand_attributes()

    def _expand_attributes(self) -> None:
        """Expand attributes into separate columns, must be called after parsing attributes"""
        if self.df is None:
            return

        self.df["parsed_attributes"] = self.df["attributes"].apply(AttributeParser.parse_attributes)

        existing_cols = set(self.df.columns)
        attr_df = pd.DataFrame(self.df["parsed_attributes"].tolist())

        # Only add new columns
        new_cols = [col for col in attr_df.columns if col not in existing_cols]
        if new_cols:
            self.df = pd.concat([self.df, attr_df[new_cols]], axis=1)

        # I need Parent to be a column, even if it's empty
        # This is because I use it to groupby in the GFFFilter class
        if "Parent" not in self.df.columns:
            self.df["Parent"] = None

        self.df.drop(columns=["parsed_attributes"], inplace=True)

    def validate_dataframe(self, feature_type: str | None = None) -> bool:
        """Returns True if the dataframe is not empty and the feature type is not None"""
        if self.df is None or feature_type is None:
            self.logger.warning(
                "The GFF file is empty after filtering.\n"
                "This might mean that there are no records within the GFF that match the "
                "sequence ID(s) in the given Fasta file.\n"
                "This could also mean that there are no records within the GFF that match "
                f"the feature type '[cyan]{feature_type}[/cyan]'.\n"
                "Please check your inputs and try again."
            )
            return False
        return True

    def export_gff_to_file(self, path: str) -> None:
        """
        Writes out the df to a gff file according to the GFF3 spec.
        This means the combined attributes column is included, but the individual attribute columns are not.
        """
        assert self.df is not None, "DataFrame must be set before writing out GFF file"
        output_df = self.df[["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]

        with open(path, "w", encoding="utf-8") as f:
            assert self.header is not None, "Header must be set before writing out GFF file"
            f.write(self.header.export_header_text())
            f.write(output_df.to_csv(sep="\t", index=False, header=False))


class SequenceReader:
    """Handles reading sequence and gff files"""

    def __init__(self, logger: Logger, verbose: bool = False) -> None:
        self.logger = logger
        self.verbose = verbose

    def read_gff(self, file: Path) -> GFFDataFrame:
        """
        Reads a GFF file and returns a GffDataFrame object.

        Parameters
        ----------
        file : Path
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
        return GFFDataFrame(inputfile=file, logger=self.logger, verbose=self.verbose)

    def read_fasta(self, file: Path) -> list[SeqRecord]:
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
            self.logger.info(f"Parsing FASTA input file: '[green]{file.name}[/green]'")
        return list(SeqIO.parse(file, "fasta"))
