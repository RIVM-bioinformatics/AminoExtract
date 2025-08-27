from AminoExtract.reader import AttributeParser


def test_parse_attributes():
    attr_string = "a=1;b=2;c=3"
    expected = {"a": "1", "b": "2", "c": "3"}
    assert AttributeParser.parse_attributes(attr_string) == expected


def test_normalize_name_attribute():
    attr_test_cases = [
        ("gene_name=ABC123;id=1_1;seqid=NC_045512.2", "gene_name=ABC123;id=1_1;seqid=NC_045512.2;Name=ABC123"),
        ("id=1_1;seqid=NC_045512.2", "id=1_1;seqid=NC_045512.2;Name=1_1"),
        ("seqid=NC_045512.2", "seqid=NC_045512.2;Name=NC_045512.2"),
    ]

    for attr, expected in attr_test_cases:
        assert AttributeParser.normalize_name_attribute(attr) == expected
