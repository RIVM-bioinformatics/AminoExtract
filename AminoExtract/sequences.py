import logging
from dataclasses import dataclass

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas import Series

from AminoExtract.gff_data import SplicingInfo
from AminoExtract.reader import GffDataFrame


@dataclass
class ExonData:
    start: int
    end: int
    strand: str
    phase: int


@dataclass
class FeatureData:
    name: str
    exons: list[ExonData]
    sequence_id: str


class SequenceExtractor:
    def __init__(self, keep_gaps: bool = False, verbose: bool = False):
        self.keep_gaps = keep_gaps
        self.logger = logging.getLogger(__name__)
        self.verbose = verbose

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
        exon_sequences = [
            self._extract_single_exon(seq_dict, exon, feature) for exon in feature.exons
        ]

        full_seq = self._combine_exons(exon_sequences, feature.exons[0].strand)
        return full_seq.translate(to_stop=True)

    def _get_splicing_detail(
        self, gff_obj: GffDataFrame, row: Series
    ) -> SplicingInfo | None:
        if gff_obj.splicing_info:
            if not hasattr(row, "ID"):
                raise ValueError(
                    f"If there are splicing details, the GFF must have an 'ID' column"
                )
            splicing_details = [x for x in gff_obj.splicing_info if row.ID == x.gene_id]
            assert (
                len(splicing_details) == 1  # sanity check
            ), f"CDSes can only have one gene, found {len(splicing_details)}"

            splicing_detail = splicing_details[0]
        else:
            splicing_detail = None
        return splicing_detail

    def _get_exons(self, row, splicing_info: SplicingInfo) -> list[ExonData]:
        """Extract exon information from GFF row"""
        return [
            ExonData(cds[0], cds[1], row.strand, getattr(row, "phase", 0))
            for cds in splicing_info.cds_locations
        ]

    def extract_aminoacids(
        self,
        gff_obj: GffDataFrame,
        seq_records: list[SeqRecord],
    ) -> dict[str, dict[str, Seq]]:
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
        seq_dict: dict[str, Seq] = {record.id: record.seq for record in seq_records}
        result: dict[str, dict[str, Seq]] = {record.id: {} for record in seq_records}

        for row in gff_obj.df.itertuples():
            name = getattr(row, "Name", f"ID-{row.seqid}")

            splicing_detail = self._get_splicing_detail(gff_obj, row)

            feature = FeatureData(
                name=name,
                sequence_id=row.seqid,
                exons=self._get_exons(row, splicing_detail),
            )

            result[feature.sequence_id][feature.name] = self.extract_feature(
                feature, seq_dict
            )

        return result
