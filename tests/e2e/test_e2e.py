"""End-to-end tests for the AminoExtract package."""

import os

# see tests/unit/test_reader.py for explanation
# pylint: disable=import-error
from AminoExtract.__main__ import main


def _order_fasta_by_name(fasta_lines: list[str]) -> list[str]:
    """Function to order a fasta file by the read name."""
    chunks = [fasta_lines[i : i + 2] for i in range(0, len(fasta_lines), 2)]
    sorted_chunks = sorted(chunks, key=lambda x: x[0])
    sorted_fastq = [line for chunk in sorted_chunks for line in chunk]
    return sorted_fastq


def _compare_outputs(output_file: str, expected_output_file: str) -> None:
    """Compare the output fasta file with the expected output fasta file."""
    with open(output_file, "r", encoding="utf-8") as output, open(
        expected_output_file, "r", encoding="utf-8"
    ) as expected_output:
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

    def test_e2e_simple(self):
        """A simple end-to-end test for the AminoExtract package."""
        if os.path.exists("tests/data/e2e/simple_output.fa"):
            os.remove("tests/data/e2e/simple_output.fa")
        args = [
            "-i",
            "tests/data/e2e/simple_input.fa",
            "-gff",
            "tests/data/e2e/simple_input.gff",
            "-o",
            "tests/data/e2e/simple_output.fa",
            "-n",
            "simple_output",
        ]
        main(args)
        assert os.path.exists("tests/data/e2e/simple_output.fa")
        _compare_outputs(
            "tests/data/e2e/simple_output.fa", "tests/data/e2e/simple_output.faa"
        )

    def test_e2e_complex(self):
        """
        A complex end-to-end test for the AminoExtract package.

        Notes:
        ------
        Based on the example for here:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        """
        if os.path.exists("tests/data/e2e/complex_output.fa"):
            os.remove("tests/data/e2e/complex_output.fa")
        args = [
            "-i",
            "tests/data/e2e/complex_input.fa",
            "-gff",
            "tests/data/e2e/complex_input.gff",
            "-o",
            "tests/data/e2e/complex_output.fa",
            "-n",
            "complex_input",
        ]
        main(args)
        assert os.path.exists("tests/data/e2e/complex_output.fa")
        _compare_outputs(
            "tests/data/e2e/complex_output.fa", "tests/data/e2e/complex_output.faa"
        )
