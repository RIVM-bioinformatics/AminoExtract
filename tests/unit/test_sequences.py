"""
Tests for the SequenceExtractor class and related functions.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from AminoExtract.__main__ import AminoAcidExtractor

# see tests/unit/test_reader.py for explanation
# pylint: disable=import-error
from AminoExtract.gff_data import SplicingInfo
from AminoExtract.logging import log
from AminoExtract.reader import GFFDataFrame
from AminoExtract.sequences import ExonData, FeatureData, SequenceExtractor

# This is a test file, so we can ignore protected-access
# pylint: disable=protected-access


# If a fixture is used in the same module in which it is defined,
# the function name of the fixture will be shadowed by
# the function arg that requests the fixture
# https://docs.pytest.org/en/stable/reference/reference.html#pytest-fixture
@pytest.fixture(name="sequence_extractor")
def fixture_sequence_extractor() -> SequenceExtractor:
    """Return a SequenceExtractor instance"""
    return SequenceExtractor(logger=log, verbose=False, keep_gaps=False)


@pytest.fixture(name="sample_sequence")
def fixture_sample_sequence() -> Seq:
    """Return a sample sequence with gaps"""
    return Seq("ATG-CC-TAG")


@pytest.fixture(name="sample_seq_dict")
def fixture_sample_seq_dict() -> dict[str, Seq]:
    """Return a sample sequence dictionary"""
    return {"seq1": Seq("ATGCCCTAGGGG")}


@pytest.fixture(name="sample_exon")
def fixture_sample_exon() -> ExonData:
    """Return a sample exon"""
    return ExonData(start=1, end=6, strand="+", phase=0)


@pytest.fixture(name="sample_feature")
def fixture_sample_feature() -> FeatureData:
    """Return a sample feature"""
    exon = ExonData(start=1, end=6, strand="+", phase=0)
    return FeatureData(name="gene1", exons=[exon], sequence_id="seq1")


@pytest.fixture(name="sample_splicing_info")
def fixture_sample_splicing_info() -> list[SplicingInfo]:
    """Return a sample splicing info list"""
    return [SplicingInfo(cds_locations=[(1, 6)], parent_ids={"parent1"}, gene_id="gene1")]


@pytest.fixture(name="sample_gff_data")
def fixture_sample_gff_data() -> GFFDataFrame:
    """Return a sample GFFDataFrame"""
    data = {
        "seqid": ["seq1"],
        "type": ["gene"],
        "start": [1],
        "end": [6],
        "strand": ["+"],
        "phase": [0],
        "ID": ["gene1"],
        "Name": ["test_gene"],
    }
    df = pd.DataFrame(data)
    gff = GFFDataFrame.__new__(GFFDataFrame)  # Create without __init__
    gff.df = df
    return gff


def test_process_sequence_no_gaps(sequence_extractor: SequenceExtractor, sample_sequence: Seq) -> None:
    """Test sequence processing with gaps disabled"""
    result = sequence_extractor._process_sequence(sample_sequence)
    assert str(result) == "ATGCCTAG"


def test_process_sequence_with_gaps(sample_sequence: Seq) -> None:
    """Test sequence processing with gaps enabled"""
    extractor = SequenceExtractor(keep_gaps=True)
    result = extractor._process_sequence(sample_sequence)
    assert str(result) == "ATGNCCNTAG"


def test_extract_single_exon(
    sequence_extractor: SequenceExtractor,
    sample_seq_dict: dict[str, Seq],
    sample_exon: ExonData,
    sample_feature: FeatureData,
) -> None:
    """Test extracting a single exon sequence"""
    result = sequence_extractor._extract_single_exon(sample_seq_dict, sample_exon, sample_feature)
    assert str(result) == "ATGCCC"  # 1-based to 0-based conversion


def test_combine_exons_forward(sequence_extractor: SequenceExtractor) -> None:
    """Test combining exons on forward strand"""
    sequences = [Seq("ATG"), Seq("CCC"), Seq("TAG")]
    result = sequence_extractor._combine_exons(sequences, "+")
    assert str(result) == "ATGCCCTAG"


def test_combine_exons_reverse(sequence_extractor: SequenceExtractor) -> None:
    """Test combining exons on reverse strand"""
    sequences = [Seq("ATG"), Seq("CCC"), Seq("TAG")]
    result = sequence_extractor._combine_exons(sequences, "-")
    assert str(result) == "CTAGGGCAT"  # Reverse complement


def test_extract_feature(sequence_extractor: SequenceExtractor, sample_seq_dict: dict[str, Seq]) -> None:
    """
    Test complete feature extraction.

    Notes:
    -------
    GFF start and end values are 1-based and inclusive
    https://www.ensembl.org/info/website/upload/gff.html
    """
    feature = FeatureData(
        name="gene1",
        exons=[
            ExonData(start=2, end=4, strand="+", phase=0),
            ExonData(start=8, end=10, strand="+", phase=0),
        ],
        sequence_id="seq1",
    )
    result = sequence_extractor.extract_feature(feature, sample_seq_dict)
    assert isinstance(result, Seq)
    assert str(result) == "CR"  # TGC + AGG


def test_get_splicing_detail(
    sequence_extractor: SequenceExtractor,
    sample_gff_data: GFFDataFrame,
    sample_splicing_info: list[SplicingInfo],
) -> None:
    """Test splicing detail extraction"""

    sample_gff_data.splicing_info = sample_splicing_info
    assert sample_gff_data.df is not None  # Ensure data is set before test
    row = sample_gff_data.df.iloc[0]
    result = sequence_extractor._get_splicing_detail(sample_gff_data, row)
    assert isinstance(result, SplicingInfo)
    assert result.gene_id == "gene1"


def test_extract_aminoacids_integration(
    sequence_extractor: SequenceExtractor,
    sample_gff_data: GFFDataFrame,
    sample_splicing_info: list[SplicingInfo],
) -> None:
    """Integration test for complete amino acid extraction"""
    seq_records = [SeqRecord(Seq("ATGCCCTAG"), id="seq1", name="test_seq")]

    sample_gff_data.splicing_info = sample_splicing_info

    unique_column_name = AminoAcidExtractor.get_unique_col_name(sample_gff_data)

    result = sequence_extractor.extract_aminoacids(sample_gff_data, seq_records, unique_col_name=unique_column_name)
    assert isinstance(result, dict)
    assert "seq1" in result
    assert isinstance(result["seq1"], dict)
