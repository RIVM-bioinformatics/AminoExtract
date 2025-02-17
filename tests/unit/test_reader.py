""""Unit tests for the reader module of AminoExtract"""

from pathlib import Path

# pylint cannot find modules in subdirectories of a repo
# https://shorturl.at/AyTL0
# pylint: disable=import-error
from AminoExtract import reader


class TestReader:
    """Test the reader module"""

    def test_basic_gff(self):
        """Test basic GFF file"""
        gff = reader.GFFDataFrame(inputfile="tests/data/test_reader.gff")
        assert gff.file_path == Path("tests/data/test_reader.gff")
        assert gff.verbose is False
        assert (
            gff.header.raw_text
            == "##gff-version 3\n##sequence-region CY121680.1 1 1752\n"
        )
        assert gff.df.shape == (2, 10)
        assert (
            gff.df.Name[1] == "test2"
        )  # [1] has gene_name instead of name, which needs to be normalized

    def test_gzipped_gff(self):
        """Test gzipped GFF file"""
        gff = reader.GFFDataFrame(inputfile="tests/data/test_reader.gff.gz")
        assert gff.file_path == Path("tests/data/test_reader.gff.gz")
        assert gff.verbose is False
        assert (
            gff.header.raw_text
            == "##gff-version 3\n##sequence-region CY121680.1 1 1752\n"
        )
        assert gff.df.shape == (2, 10)
        assert gff.df.Name[1] == "test2"
