"""
Module for filtering GFF records and SeqRecords based on sequence IDs and feature types.
"""

from logging import Logger
from typing import cast

import pandas as pd
from Bio.SeqRecord import SeqRecord

from AminoExtract.gff_data import SplicingInfo
from AminoExtract.reader import GFFDataFrame


class GFFFilter:
    """
    Class for filtering GFF records based on sequence IDs and feature types.
    """

    def __init__(self, df: pd.DataFrame, logger: Logger, verbose: bool = False) -> None:
        self.df = df
        self.logger = logger
        self.verbose = verbose

    def filter_by_seqids(self, seq_ids: list[str], df: pd.DataFrame) -> pd.DataFrame:
        """Filter GFF records of the GFFFilter class objext by some sequence IDs"""
        return df[df["seqid"].isin(seq_ids)]

    def filter_by_feature_type(self, feature_type: str, df: pd.DataFrame) -> pd.DataFrame:
        """Filter GFF records of the GFFFilter class object by some feature type"""
        if feature_type == "all":
            return df
        return df[df["type"] == feature_type]

    def get_splicing_info(self, filtered_df: pd.DataFrame, unique_col_name: str) -> list[SplicingInfo]:
        """
        Get splicing information from a filtered dataframe.
        Attributes are assumed to be parsed.
        """
        splicing = filtered_df.groupby(unique_col_name).agg({"start": tuple, "end": tuple, "Parent": set})

        start_lengths = splicing["start"].apply(len)
        end_lengths = splicing["end"].apply(len)
        if not start_lengths.eq(end_lengths).all():  # sanity check
            raise ValueError("There should be an equal amount coding start locations as coding end locations")

        splicing["CDSes"] = splicing.apply(lambda x: list(zip(x["start"], x["end"])), axis=1)

        # Im using cast here to tell mypy that the types are correct
        # This is necessary because if you access a row in a DataFrame, it returns as an Any object
        return [
            SplicingInfo(
                cds_locations=cast(list[tuple[int, int]], row.CDSes),
                parent_ids=cast(set[str], row.Parent),
                gene_id=cast(str, row.Index),
            )
            for row in splicing[["CDSes", "Parent"]].itertuples()
        ]


class GFFRecordFilter:
    """Wrapper class for GFFFilter that applies those filters to GFF records."""

    # This class is made to bundle the whole process of filtering GFF records,
    # so it only has one method
    # pylint: disable=too-few-public-methods

    def __init__(self, gff_records: GFFDataFrame, logger: Logger, verbose: bool = False) -> None:
        self.gff_records = gff_records
        assert gff_records.df is not None, "GFF records are empty"
        self.logger = logger
        self.filter = GFFFilter(gff_records.df, self.logger, verbose)
        self.verbose = verbose

    def apply_filters(self, seq_records: list[SeqRecord], feature_type: str, unique_col_name: str) -> GFFDataFrame:
        """
        Applies filters to sequences based on sequence IDs and feature types
        and stores them as a GFFDataFrame object.
        """
        seq_ids = [record.id for record in seq_records if record.id is not None]

        if self.verbose:
            self._log_filtering_info(seq_ids, feature_type)

        assert seq_ids, "No sequence IDs found in the given SeqRecords"
        filtered_df = self.filter.filter_by_seqids(seq_ids, self.filter.df)
        filtered_df = self.filter.filter_by_feature_type(feature_type, filtered_df)

        self.gff_records.df = filtered_df
        self.gff_records.splicing_info = self.filter.get_splicing_info(filtered_df, unique_col_name)

        return self.gff_records

    def _log_filtering_info(self, seq_ids: list[str], feature_type: str) -> None:
        """Logging for verbose mode"""
        self.logger.info("Filtering GFF records:")
        self.logger.info(f"Feature type: {feature_type}")
        self.logger.info(f"Sequence IDs: {', '.join(seq_ids)}")


class SequenceFilter:
    """
    Filter SeqRecord objects based on sequence IDs.
    """

    # This could maybe be a function instead of a class,
    # but all the other classes are classes and I wanted to keep it consistent
    # pylint: disable=too-few-public-methods

    def __init__(self, seq_records: list[SeqRecord], logger: Logger, verbose: bool = False) -> None:
        self.seq_records = seq_records
        self.logger = logger
        self.verbose = verbose

    def filter_sequences(self, gff: GFFDataFrame) -> list[SeqRecord]:
        """
        Takes a GffDataFrame object and a list of SeqRecord objects,
        and returns a list of SeqRecord objects that only contain
        the sequences that are specified in the GffDataFrame object.

        Parameters
        ----------
        gff : GffDataFrame
            The GffDataFrame object to filter
        SeqRecords : list
            A list of SeqRecord objects

        Returns
        -------
            A list of SeqRecord objects that only contain the sequences that
            are specified in the GffDataFrame object.

        """

        if gff.df is None:
            raise ValueError("The GFF file is empty")
        valid_ids = set(gff.df["seqid"])
        return [record for record in self.seq_records if record.id in valid_ids]
