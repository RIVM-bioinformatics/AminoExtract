"""
Copyright (C) 2022  RIVM

https://github.com/RIVM-bioinformatics/AminoExtract
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from AminoExtract.args import validate_args
from AminoExtract.filter import GFFRecordFilter, SequenceFilter
from AminoExtract.logging import log
from AminoExtract.reader import SequenceReader
from AminoExtract.sequences import GFFDataFrame, SequenceExtractor
from AminoExtract.writer import FastaWriter


@dataclass
class Config:
    """Configuration class for AminoExtract."""

    # For a config class I dont care how many attributes it has,
    # so longs as they are all together in one place.
    # pylint: disable=too-many-instance-attributes

    input_fasta: Path
    input_gff: Path
    output: Path
    feature_type: str
    keep_gaps: bool
    verbose: bool
    name: str
    outtype_int: int
    outtype: str

    def __post_init__(self) -> None:
        """
        Add a string representation of the output type.
        This makes it easier to use the output type in the FastaWriter.
        """
        # 1 is dir, 0 is file
        if self.outtype_int == 0:
            self.outtype = "single_file"
        elif self.outtype_int == 1:
            self.outtype = "multiple_files"
        else:
            raise ValueError(f"Invalid output type: {self.outtype_int}")


class AminoAcidExtractor:
    """
    Extract amino acid sequences from GFF annotations and FASTA files.
    """

    # as this is a __main__ module, I just want it to run the main function
    # pylint: disable=too-few-public-methods

    def __init__(self, config: Config | None = None) -> None:
        self.log = log
        if not config:
            print("No config provided. Only reader access available.")
            self.reader = SequenceReader(logger=self.log, verbose=False)
        else:
            self.reader = SequenceReader(logger=self.log, verbose=config.verbose)
            self.config = config
            self.seq_records = self.reader.read_fasta(config.input_fasta)

    def run(self) -> None:
        """
        First load and filter the GFF file, then extract the sequences and write the output.
        """
        gff_data = self._load_and_filter_gff()
        unique_column_name = self.get_unique_col_name(gff_data)
        sequences = self._extract_sequences(gff_data, unique_column_name)
        self._write_output(sequences)

    @staticmethod
    def get_unique_col_name(gff: GFFDataFrame) -> str:
        if "ID" in gff.df.columns:
            return "ID"
        elif "gene" in gff.df.columns:
            return "gene"
        elif "seqid" in gff.df.columns:
            return "seqid"
        else:
            raise ValueError("No valid unique column found in GFF DataFrame.")

    def _load_and_filter_gff(self) -> GFFDataFrame:
        gff_obj = self.reader.read_gff(self.config.input_gff)
        unique_column_name = self.get_unique_col_name(gff_obj)

        gff_filter = GFFRecordFilter(gff_records=gff_obj, logger=self.log, verbose=self.config.verbose)
        filtered_gff = gff_filter.apply_filters(
            seq_records=self.seq_records, feature_type=self.config.feature_type, unique_col_name=unique_column_name
        )

        if not filtered_gff.validate_dataframe(self.config.feature_type):
            raise ValueError("Validation failed, either the GFF file is empty or the feature type is None")
        return filtered_gff

    def _extract_sequences(self, gff_data: GFFDataFrame, unique_column_name: str) -> dict[str, dict[str, Seq]]:
        seq_filter = SequenceFilter(
            seq_records=self.seq_records,
            logger=self.log,
            verbose=self.config.verbose,
        )
        filtered_seq_records = seq_filter.filter_sequences(gff_data)

        extractor = SequenceExtractor(
            logger=self.log,
            verbose=self.config.verbose,
            keep_gaps=self.config.keep_gaps,
        )
        return extractor.extract_aminoacids(gff_obj=gff_data, seq_records=filtered_seq_records, unique_col_name=unique_column_name)

    def _write_output(self, sequences: dict[str, dict[str, Seq]]) -> None:
        writer = FastaWriter(output_path=self.config.output, logger=self.log)
        writer.write(sequences, self.config.name, self.config.outtype)


def get_feature_name_attribute(input_gff: str, input_seq: str, feature_type: str) -> dict[str, list[str]]:
    """
    This function takes a GFF file, a FASTA file and a feature type,
    and returns a dictionary of the feature names for each sequence in the FASTA file.
    Allows

    Parameters
    ----------
    input_gff : str
        the path to the gff file
    input_seq : str
        The path to the fasta file containing the sequences you want to extract features from.
    feature_type : str
        This is the type of feature you want to extract from the GFF file.
        For example, if you want to
    extract the CDS features, you would enter "CDS". Options are "CDS", "gene" and "all".

        Returns
        -------
            A dict with the sequence id as the key and a list of the feature names as the value.

    """
    reader = SequenceReader(log, verbose=False)
    gff = reader.read_gff(Path(input_gff))
    seq = reader.read_fasta(Path(input_seq))

    unique_col_name = AminoAcidExtractor.get_unique_col_name(gff)

    gff_filter = GFFRecordFilter(gff_records=gff, logger=log, verbose=False)
    gff_records = gff_filter.apply_filters(seq_records=seq, feature_type=feature_type, unique_col_name=unique_col_name)
    assert gff_records.df is not None, "The GFF file is empty"

    seq_filter = SequenceFilter(seq_records=seq, logger=log, verbose=False)
    filtered_seqs = seq_filter.filter_sequences(gff_records)

    seq_attributes: dict[str, list[str]] = {record.id: [] for record in filtered_seqs if record.id is not None}

    valid_name_type_columns = ["Name", "gene", "ID"]
    name_type_col = next((col for col in valid_name_type_columns if col in gff_records.df.columns), None)
    if name_type_col is None:
        gff_records.df["Name"] = pd.NA
        name_type_col = "Name"
    for row in gff_records.df.itertuples():
        value = getattr(row, name_type_col)

        assert isinstance(row.seqid, str) and (
            isinstance(value, str) or pd.isna(value)
        ), f"Invalid row({row}): seqid needs to be a str and is {type(row.seqid)}, Name needs to be a str or NaN and is {type(value)}"
        seq_attributes[row.seqid].append(value)

    new_seq_attr = {}
    for k, v in seq_attributes.items():
        new_seq_attr[k] = [x for x in v if not pd.isna(x)]

    return new_seq_attr


def main(provided_args: list[str] | None = None) -> None:
    """
    Extract amino acid sequences from GFF annotations and FASTA files.

    Reads GFF and FASTA files, filters features based on type, extracts amino acid
    sequences and writes results to file.

    Parameters
    ----------
    provided_args : list[str] | None, optional
        Command line arguments as list of strings. If None, sys.argv[1:] is used.

    Returns
    -------
    None
        Writes output file with amino acid sequences.

    Examples
    --------
    >>> main(['--input', 'sequences.fa', '--features', 'annot.gff'])
    >>> main(['--input', 'seqs.fa', '--features', 'genes.gff', '--feature-type', 'CDS'])

    Notes
    -----
    Exit codes:
        0: Success
        1: Something went wrong
    """
    args = validate_args(provided_args if provided_args else sys.argv[1:])

    config = Config(
        input_fasta=args.input,
        input_gff=args.features,
        output=args.output,
        feature_type=args.feature_type,
        keep_gaps=args.keep_gaps,
        verbose=args.verbose,
        name=args.name,
        outtype_int=args.outtype,
        outtype="",
    )

    extractor = AminoAcidExtractor(config)
    extractor.run()
