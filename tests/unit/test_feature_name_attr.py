from AminoExtract import get_feature_name_attribute


def test_get_feat_name_attr():
    result = get_feature_name_attribute(
        input_gff="tests/data/unit/test_eqa_reader.gff", input_seq="tests/data/unit/reference_genome.fasta", feature_type="all"
    )
    assert result == {"NC_045512.2": ["1_1", "1_2", "1_3", "1_4", "1_5", "1_6", "1_7", "1_8", "1_9"]}


def test_get_feat_name_attr_complex():
    result = get_feature_name_attribute(
        input_gff="tests/data/unit/test_eqa_reader_complex.gff", input_seq="tests/data/unit/reference_genome.fasta", feature_type="all"
    )
    assert result == {
        "NC_045512.2": [
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "ORF1ab",
            "S",
            "S",
            "ORF3a",
            "ORF3a",
            "E",
            "E",
            "M",
            "M",
            "ORF6",
            "ORF6",
            "ORF7a",
            "ORF7a",
            "ORF7b",
            "ORF7b",
            "ORF8",
            "ORF8",
            "N",
            "N",
            "ORF10",
            "ORF10",
            "ORF10",
            "ORF10",
        ]
    }
