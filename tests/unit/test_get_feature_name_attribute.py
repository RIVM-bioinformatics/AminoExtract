from AminoExtract import get_feature_name_attribute


def test_get_feature_name_attribute_viro():
    feature_path = "tests/data/unit/viro_input.gff"
    seq_path = "tests/data/unit/viro_input.fasta"
    feature_type = "CDS"

    a = get_feature_name_attribute(
        input_gff=feature_path, input_seq=seq_path, feature_type=feature_type
    )
