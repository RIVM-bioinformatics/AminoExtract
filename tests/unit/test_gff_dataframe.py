from AminoExtract import SequenceReader


def test_exporting_gff_files():
    reader = SequenceReader(None, False)
    gff_df_obj = reader.read_gff("tests/data/unit/test_eqa_reader.gff")
    gff_df_obj.export_gff_to_file("tests/data/unit/test_eqa_reader_exported.gff")
