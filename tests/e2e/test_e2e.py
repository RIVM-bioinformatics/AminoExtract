"""End-to-end tests for the AminoExtract package."""

from pathlib import Path

import pytest
from Bio import BiopythonWarning

from AminoExtract import ROOT_DIR

# see tests/unit/test_reader.py for explanation
# pylint: disable=import-error
from AminoExtract.__main__ import main


def _order_fasta_by_name(fasta_lines: list[str]) -> list[str]:
    """Function to order a fasta file by the read name."""
    chunks = [fasta_lines[i : i + 2] for i in range(0, len(fasta_lines), 2)]
    sorted_chunks = sorted(chunks, key=lambda x: x[0])
    sorted_fastq = [line for chunk in sorted_chunks for line in chunk]
    return sorted_fastq


def _compare_outputs(output_file: Path, expected_output_file: Path) -> None:
    """Compare the output fasta file with the expected output fasta file."""
    with open(output_file, "r", encoding="utf-8") as output, open(expected_output_file, "r", encoding="utf-8") as expected_output:
        output_lines = output.readlines()
        expected_output_lines = expected_output.readlines()
        # the fasta file can be unordered, so we need to sort it by the read name
        output_lines = _order_fasta_by_name(output_lines)
        expected_output_lines = _order_fasta_by_name(expected_output_lines)
        for line1, line2 in zip(output_lines, expected_output_lines):
            line1 = line1.strip()
            line2 = line2.strip()
            if line1 != line2:
                raise AssertionError(f"Output: {line1}, Expected Output: {line2}")


class TestE2E:
    """Test class for the end-to-end tests."""

    data_path = ROOT_DIR / "tests" / "data" / "e2e"

    def test_e2e_simple(self):
        """A simple end-to-end test for the AminoExtract package."""
        output_path = self.data_path / "simple_output.fa"
        if output_path.exists():
            output_path.unlink()
        args = [
            "-i",
            str(self.data_path / "simple_input.fa"),
            "-gff",
            str(self.data_path / "simple_input.gff"),
            "-o",
            str(output_path),
            "-n",
            "simple_output",
        ]
        main(args)
        assert output_path.exists()
        _compare_outputs(output_path, self.data_path / "simple_output.faa")

    def test_e2e_complex(self):
        """
        A complex end-to-end test for the AminoExtract package.

        Notes:
        ------
        Based on the example for here:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        """
        output_path = self.data_path / "complex_output.fa"
        if output_path.exists():
            output_path.unlink()

        args = [
            "-i",
            str(self.data_path / "complex_input.fa"),
            "-gff",
            str(self.data_path / "complex_input.gff"),
            "-o",
            str(output_path),
            "-n",
            "complex_input",
        ]

        with pytest.warns(BiopythonWarning):
            # The warning is that the sequence is not divisible by 3, so not a full codon
            main(args)

        assert output_path.exists()
        _compare_outputs(output_path, self.data_path / "complex_output.faa")

    def test_e2e_sars_wuhan(self):
        """End-to-end test for the SARS-CoV-2 Wuhan reference genome."""
        output_path = self.data_path / "sars_wuhan_output.fa"
        if output_path.exists():
            output_path.unlink()

        args = [
            "-i",
            str(self.data_path / "sars_wuhan_input.fa"),
            "-gff",
            str(self.data_path / "sars_wuhan_input.gff"),
            "-o",
            str(output_path),
            "-n",
            "sars_wuhan",
        ]
        main(args)
        assert output_path.exists()
        _compare_outputs(output_path, self.data_path / "sars_wuhan_output.faa")

    def test_e2e_error2(self):
        """End-to-end test for the SARS-CoV-2 error from 2025-09-03."""
        output_path = self.data_path / "error2.fa"
        if output_path.exists():
            output_path.unlink()

        args = [
            "-i",
            str(self.data_path / "error2_input.fa"),
            "-gff",
            str(self.data_path / "error2_input.gff"),
            "-o",
            str(output_path),
            "-n",
            "error2",
        ]
        main(args)
        assert output_path.exists()
        _compare_outputs(output_path, self.data_path / "error2.faa")
