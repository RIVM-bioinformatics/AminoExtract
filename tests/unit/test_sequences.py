"""
Unit tests for the AminoExtract.sequences module, 
specifically the extract_aminoacids function.

Classes
-------
TestSequences
    Contains unit tests for the AminoExtract module.
"""

import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from AminoExtract import GffDataFrame, extract_aminoacids


class TestSequences:
    """
    Unit tests for the AminoExtract.sequences module.

    Attributes
    ----------
    mock_sequences : list of SeqRecord
        Mock sequences used for testing.
    gff_obj : GffDataFrame
        GFF data frame object used for testing.
    """

    mock_sequences = [
        SeqRecord(Seq("ATCCGAATCGGAATC"), id="seq1"),
        SeqRecord(
            Seq("GATTCCGATTCGGAT"), id="seq2"
        ),  # reverse complement of seq1, - strand in GFF
        SeqRecord(Seq("AUCGAGAUCGAGAUC"), id="seq3"),  # RNA sequence
        SeqRecord(Seq("ATC--GGAATCGGAATC-"), id="seq4"),  # with gaps
    ]
    gff_obj = GffDataFrame(inputfile="tests/data/test_sequences_module.gff")

    def test_extract_aminoacids_basic(self) -> None:
        """
        Test the extract_aminoacids function with basic parameters.

        Asserts
        -------
        bool
            True if the extracted amino acids are as expected.
        """
        assert extract_aminoacids(self.gff_obj, self.mock_sequences) == {
            "seq1": {"test1": Seq("IRIGI")},
            "seq2": {"test2": Seq("IRIGI")},
            "seq3": {"test3": Seq("IEIEI")},
            "seq4": {"test4": Seq("IGIG")},
        }

    def test_extract_aminoacids_keep_gaps(self) -> None:
        """
        Test the extract_aminoacids function with the keep_gaps parameter.

        Asserts
        -------
        bool
            True if the extracted amino acids with gaps are as expected.
        """
        assert extract_aminoacids(
            self.gff_obj, self.mock_sequences, keep_gaps=True
        ) == {
            "seq1": {"test1": Seq("IRIGI")},
            "seq2": {"test2": Seq("IRIGI")},
            "seq3": {"test3": Seq("IEIEI")},
            "seq4": {"test4": Seq("IXESE")},  # How does this work?
        }

    def test_extract_aminoacids_verbose(self, caplog):
        """
        Test the extract_aminoacids function with the verbose parameter.

        Parameters
        ----------
        caplog : pytest.LogCaptureFixture
            Fixture to capture log output.

        Asserts
        -------
        bool
            True if the expected log message is found.
        """
        with caplog.at_level(logging.INFO, logger="rich"):
            extract_aminoacids(self.gff_obj, self.mock_sequences, verbose=True)
            assert (
                "Extracting and translating the amino acid sequence(s) from the nucleotide sequence(s)"
                in caplog.text
            )
