import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from AminoExtract.gff_data import SplicingInfo
from AminoExtract.reader import GffDataFrame
from AminoExtract.sequences import ExonData, FeatureData, SequenceExtractor


@pytest.fixture
def sequence_extractor():
    return SequenceExtractor(keep_gaps=False)


@pytest.fixture
def sample_sequence():
    return Seq("ATG-CC-TAG")


@pytest.fixture
def sample_seq_dict():
    return {"seq1": Seq("ATGCCCTAGGGG")}


@pytest.fixture
def sample_exon():
    return ExonData(start=1, end=6, strand="+", phase=0)


@pytest.fixture
def sample_feature():
    exon = ExonData(start=1, end=6, strand="+", phase=0)
    return FeatureData(name="gene1", exons=[exon], sequence_id="seq1")


@pytest.fixture
def sample_splicing_info():
    return [
        SplicingInfo(cds_locations=[(1, 6)], parent_ids={"parent1"}, gene_id="gene1")
    ]


@pytest.fixture
def sample_gff_data():
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
    gff = GffDataFrame.__new__(GffDataFrame)  # Create without __init__
    gff.df = df
    return gff


def test_process_sequence_no_gaps(sequence_extractor, sample_sequence):
    """Test sequence processing with gaps disabled"""
    result = sequence_extractor._process_sequence(sample_sequence)
    assert str(result) == "ATGCCTAG"


def test_process_sequence_with_gaps(sample_sequence):
    """Test sequence processing with gaps enabled"""
    extractor = SequenceExtractor(keep_gaps=True)
    result = extractor._process_sequence(sample_sequence)
    assert str(result) == "ATGNCCNTAG"


def test_extract_single_exon(
    sequence_extractor, sample_seq_dict, sample_exon, sample_feature
):
    """Test extracting a single exon sequence"""
    result = sequence_extractor._extract_single_exon(
        sample_seq_dict, sample_exon, sample_feature, 0
    )
    assert str(result) == "ATGCCC"  # 1-based to 0-based conversion


def test_combine_exons_forward(sequence_extractor):
    """Test combining exons on forward strand"""
    sequences = [Seq("ATG"), Seq("CCC"), Seq("TAG")]
    result = sequence_extractor._combine_exons(sequences, "+")
    assert str(result) == "ATGCCCTAG"


def test_combine_exons_reverse(sequence_extractor):
    """Test combining exons on reverse strand"""
    sequences = [Seq("ATG"), Seq("CCC"), Seq("TAG")]
    result = sequence_extractor._combine_exons(sequences, "-")
    assert str(result) == "CTAGGGCAT"  # Reverse complement


def test_extract_feature(sequence_extractor, sample_seq_dict):
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


def test_get_splicing_detail(sequence_extractor, sample_gff_data, sample_splicing_info):
    """Test splicing detail extraction"""

    sample_gff_data.splicing_info = sample_splicing_info

    row = sample_gff_data.df.iloc[0]
    result = sequence_extractor._get_splicing_detail(sample_gff_data, row)
    assert isinstance(result, SplicingInfo)
    assert result.gene_id == "gene1"


def test_extract_aminoacids_integration(
    sequence_extractor, sample_gff_data, sample_splicing_info
):
    """Integration test for complete amino acid extraction"""
    seq_records = [SeqRecord(Seq("ATGCCCTAG"), id="seq1", name="test_seq")]

    sample_gff_data.splicing_info = sample_splicing_info

    result = sequence_extractor.extract_aminoacids(sample_gff_data, seq_records)
    assert isinstance(result, dict)
    assert "seq1" in result
    assert isinstance(result["seq1"], dict)
