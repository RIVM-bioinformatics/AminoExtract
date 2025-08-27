from AminoExtract import get_feature_name_attribute


def test_get_feat_name_attr():
    result = get_feature_name_attribute(
        input_gff="tests/data/unit/test_eqa_reader.gff", input_seq="tests/data/unit/reference_genome.fasta", feature_type="all"
    )
    assert result == {"NC_045512.2": ["1_1", "1_2", "1_3", "1_4", "1_5", "1_6", "1_7", "1_8", "1_9"]}
