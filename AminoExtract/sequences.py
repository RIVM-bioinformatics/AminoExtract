"""Sequence extraction module for AminoExtract."""

from dataclasses import dataclass
from logging import Logger

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import Series

from AminoExtract.gff_data import SplicingInfo
from AminoExtract.logging import log
from AminoExtract.reader import GFFDataFrame


@dataclass
class ExonData:
    """Dataclass for exon information"""

    start: int
    end: int
    strand: str
    phase: int


@dataclass
class FeatureData:
    """Dataclass for feature information"""

    name: str
    exons: list[ExonData]
    sequence_id: str


class SequenceExtractor:
    """Extract amino acid sequences from sequence records based on GFF annotations."""

    def __init__(self, logger: Logger = log, verbose: bool = False, keep_gaps: bool = False) -> None:
        self.logger = logger
        self.verbose = verbose
        self.keep_gaps = keep_gaps

    def _process_sequence(self, sequence: Seq) -> Seq:
        """Process sequence based on gap handling preference"""
        if self.keep_gaps:
            return sequence.replace("-", "N")
        return sequence.replace("-", "")

    def _extract_single_exon(
        self,
        seq_dict: dict[str, Seq],
        exon: ExonData,
        feat: FeatureData,
    ) -> Seq:
        """Extract sequence for a single exon"""

        # GFF start and end values are 1-based indexed and inclusive
        # https://www.ensembl.org/info/website/upload/gff.html
        start = exon.start - 1  # Convert to 0-based
        sequence = seq_dict[feat.sequence_id][start : exon.end]
        return self._process_sequence(sequence)

    def _combine_exons(self, sequences: list[Seq], strand: str) -> Seq:
        """Combine exon sequences and handle strand orientation"""
        combined = sum(sequences, Seq(""))
        return combined.reverse_complement() if strand == "-" else combined

    def extract_feature(self, feature: FeatureData, seq_dict: dict[str, Seq]) -> Seq:
        """Extract and translate sequence for a single feature"""
        exon_sequences = [self._extract_single_exon(seq_dict, exon, feature) for exon in feature.exons]

        full_seq = self._combine_exons(exon_sequences, feature.exons[0].strand)
        return full_seq.translate(to_stop=True)

    def _get_splicing_detail(self, gff_obj: GFFDataFrame, row: Series) -> SplicingInfo:
        """
        Retrieve splicing details for a GFF row, using available identifiers.
        It needs to find a unique identifier, so it tries ID and gene first.
        seqid is a fallback because most GFF files have them, but it is not unique so it wont split features.
        """
        gene_id = getattr(row, "ID", None) or getattr(row, "gene", None) or getattr(row, "seqid", None)
        if gene_id is None:
            raise ValueError("Row must have an 'ID', 'gene', or 'seqid' attribute for splicing details.")

        assert gff_obj.splicing_info is not None, "No splicing information loaded"
        splicing_details = [x for x in gff_obj.splicing_info if gene_id == x.gene_id]
        if len(splicing_details) != 1:
            raise ValueError(f"Expected one splicing detail for gene_id '{gene_id}', found {len(splicing_details)}")
        return splicing_details[0]

    def _get_exons(self, row, splicing_info: SplicingInfo) -> list[ExonData]:
        """Extract exon information from GFF row"""
        return [ExonData(cds[0], cds[1], row.strand, getattr(row, "phase", 0)) for cds in splicing_info.cds_locations]

    def extract_aminoacids(self, gff_obj: GFFDataFrame, seq_records: list[SeqRecord], unique_col_name: str) -> dict[str, dict[str, Seq]]:
        """Extract amino acid sequences from sequence records based on GFF annotations.

        Parameters
        ----------
        gff_obj : GffDataFrame
            GFF annotation data
        seq_records : list[SeqRecord]
            Sequence records to process
        keep_gaps : bool
            Whether to preserve gaps in sequences
        verbose : bool
            Enable verbose logging

        Returns
        -------
        dict[str, dict[str, Seq]]
            Nested dictionary of sequence IDs and their features' amino acid sequences
        """

        # This next block is a bit weird, because it would be easier to use 2 dict comprehensions
        # But SeqRecord objects are not guaranteed to have IDs, so we need to check for that
        # Besides, now we only need to iterate over the sequence records once
        seq_dict: dict[str, Seq] = {}
        for record in seq_records:
            if record.id is None:
                raise ValueError("Sequence record ID cannot be None")
            seq_dict[record.id] = record.seq
            result: dict[str, dict[str, Seq]] = {}
            result[record.id] = {}

        assert gff_obj.df is not None, "No dataframe of the GFF loaded"  # sanity check

        # I do not see a way to avoid this for loop, I tried vectorizing it but it became unreadable
        for _, row in gff_obj.df.iterrows():
            name = getattr(row, "Name", f"ID-{row.seqid}-{getattr(row, unique_col_name, 'unknown')}")

            splicing_detail = self._get_splicing_detail(gff_obj, row)

            feature = FeatureData(
                name=name,
                sequence_id=row.seqid,
                exons=self._get_exons(row, splicing_detail),
            )

            result[feature.sequence_id][feature.name] = self.extract_feature(feature, seq_dict)

        return result
