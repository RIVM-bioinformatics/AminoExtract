from dataclasses import dataclass

import pandas as pd
from Bio.SeqRecord import SeqRecord

from AminoExtract.functions import log
from AminoExtract.gff_data import SplicingInfo
from AminoExtract.reader import GffDataFrame


class GFFFilter:
    def __init__(self, df: pd.DataFrame, verbose: bool = False):
        self.df = df
        self.verbose = verbose

    def filter_by_seqids(self, seq_ids: list[str]) -> pd.DataFrame:
        return self.df[self.df["seqid"].isin(seq_ids)]

    def filter_by_feature_type(self, feature_type: str) -> pd.DataFrame:
        if feature_type == "all":
            return self.df
        return self.df[self.df["type"] == feature_type]

    def get_splicing_info(self, filtered_df: pd.DataFrame) -> list[SplicingInfo | None]:
        """Get splicing information from a filtered dataframe. Attributes are assumed to be parsed."""

        splicing = filtered_df.groupby("ID").agg(
            {"start": tuple, "end": tuple, "Parent": set}
        )

        start_lengths = splicing["start"].apply(len)
        end_lengths = splicing["end"].apply(len)
        if not start_lengths.eq(end_lengths).all():  # sanity check
            raise ValueError(
                "There should be an equal amount coding start locations as coding end locations"
            )

        splicing["CDSes"] = splicing.apply(
            lambda x: list(zip(x["start"], x["end"])), axis=1
        )

        return [
            SplicingInfo(
                cds_locations=row.CDSes, parent_ids=row.Parent, gene_id=row.Index
            )
            for row in splicing[["CDSes", "Parent"]].itertuples()
        ]


class GFFRecordFilter:
    def __init__(self, gff_records: GffDataFrame, verbose: bool = False):
        self.gff_records = gff_records
        self.filter = GFFFilter(gff_records.df, verbose)
        self.verbose = verbose

    def apply_filters(
        self, seq_records: list[SeqRecord], feature_type: str
    ) -> GffDataFrame:
        seq_ids = [record.id for record in seq_records if record.id is not None]

        if self.verbose:
            self._log_filtering_info(seq_ids, feature_type)

        assert seq_ids, "No sequence IDs found in the given SeqRecords"
        filtered_df = self.filter.filter_by_seqids(seq_ids)
        filtered_df = self.filter.filter_by_feature_type(feature_type)

        self.gff_records.df = filtered_df
        self.gff_records.splicing_info = self.filter.get_splicing_info(filtered_df)

        return self.gff_records

    def _log_filtering_info(self, seq_ids: list[str], feature_type: str) -> None:
        """Logging for verbose mode"""
        log.info("Filtering GFF records:")
        log.info(f"Feature type: {feature_type}")
        log.info(f"Sequence IDs: {', '.join(seq_ids)}")


class SequenceFilter:
    def __init__(self, seq_records: list[SeqRecord], verbose: bool = False):
        self.seq_records = seq_records
        self.verbose = verbose

    def filter_sequences(self, gff: GffDataFrame) -> list[SeqRecord]:
        """Takes a GffDataFrame object and a list of SeqRecord objects, and returns a list of SeqRecord objects
        that only contain the sequences that are specified in the GffDataFrame object.

        Parameters
        ----------
        gff : GffDataFrame
            The GffDataFrame object to filter
        SeqRecords : list
            A list of SeqRecord objects

        Returns
        -------
            A list of SeqRecord objects that only contain the sequences that are specified in the GffDataFrame object.

        """

        if gff.df is None:
            raise ValueError("The GFF file is empty")
        valid_ids = set(gff.df["seqid"])
        return [record for record in self.seq_records if record.id in valid_ids]
