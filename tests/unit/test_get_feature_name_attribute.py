import pytest

from AminoExtract import get_feature_name_attribute


def test_get_feature_name_attribute_happy():
    feature_path = "tests/data/unit/test_features1.gff"
    seq_path = "tests/data/unit/test_reference1.fasta"
    feature_type = "CDS"

    a = get_feature_name_attribute(
        input_gff=feature_path, input_seq=seq_path, feature_type=feature_type
    )


def test_get_feature_name_attribute_unhappy_different_id_gff_fasta():
    feature_path = "tests/data/unit/test_features2.gff"
    seq_path = "tests/data/unit/test_reference2.fasta"
    feature_type = "CDS"

    with pytest.raises(ValueError):
        a = get_feature_name_attribute(
            input_gff=feature_path, input_seq=seq_path, feature_type=feature_type
        )
