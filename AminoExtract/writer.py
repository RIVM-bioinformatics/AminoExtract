"""Writer module for AminoExtract."""

import pathlib
from logging import Logger
from pathlib import Path

from Bio.Seq import Seq

from AminoExtract.logging import log

# see comment in AminoExtract/args.py
# pylint: disable=logging-fstring-interpolation


class FastaWriter:
    """Class for writing sequences to FASTA files"""

    # I like having a writer class instead of a long function.
    # pylint: disable=too-few-public-methods

    def __init__(self, output_path: Path, logger: Logger) -> None:
        self.output_path = output_path
        self.logger = logger

    def write(self, sequences: dict[str, dict[str, Seq]], name: str, single_file: str) -> None:
        """Write sequences to file(s)"""
        if single_file == "single_file":
            self._write_single_file(sequences, name)
        else:
            self._write_multiple_files(sequences, name)

    def _write_single_file(self, sequences: dict[str, dict[str, Seq]], name: str) -> None:
        """Write all sequences to a single FASTA file"""
        with self.output_path.open("w") as f:
            for seq_id, features in sequences.items():
                for feature, aa in features.items():
                    self._write_sequence(f, feature=feature, sequence=aa, name=name)
                    self._log_write(seq_id, feature, self.output_path)

    def _write_multiple_files(self, sequences: dict[str, dict[str, Seq]], name: str) -> None:
        """Write each sequence to a separate FASTA file"""
        self.output_path.mkdir(exist_ok=True)

        for seq_id, features in sequences.items():
            for feature, aa in features.items():
                output_file = self.output_path / f"{name}_{feature}.faa"
                with output_file.open("a") as f:  # 'a' is append mode
                    self._write_sequence(f, feature=feature, sequence=aa, name=name)
                    self._log_write(seq_id, feature, output_file)

    def _write_sequence(self, file_handle, feature: str, sequence: Seq, name: str) -> None:
        """Write single sequence in FASTA format"""
        file_handle.write(f">{name}.{feature}\n{sequence}\n")

    def _log_write(self, seq_id: str, feature: str, filepath: Path) -> None:
        """Log sequence writing operation"""
        self.logger.info(f"Writing '[cyan]{seq_id} - {feature}[/cyan]' to " f"'[green]{filepath}[/green]'")


def write_aa_file(
    aa_dict: dict[str, dict[str, Seq]],
    output: pathlib.Path,
    name: str,
    outtype: int,
) -> None:
    """
    Takes a dictionary of dictionaries of Seq.Seq objects,
    a pathlib.Path object, a string, and an integer,
    and writes the Seq.Seq objects to either a single file or multiple files

    Parameters
    ----------
    AA_dict : dict[str, dict[str, Seq.Seq]]
        a dictionary of dictionaries. The first dictionary is keyed by the sequence ID,
        and the second dictionary is keyed by the feature name.
        The value of the second dictionary is the amino acid
        sequence as a biopython Seq() object.
    output : pathlib.Path
        The path to the output file or directory as a (POSIX)Path object
    name : str
        The name of the output file.
    outtype : int
        0 = single file, 1 = multiple files

    """

    log.info(f"{'='*20} Writing output(s) {'='*20}")
    if outtype == 0:
        # this should probably be changed so a user can
        # specify the fasta header per record instead of assuming {name}.{feature}
        with open(output, "w", encoding="utf-8") as out:
            for seq_id, features in aa_dict.items():
                for feature, aa in features.items():
                    log.info(f"Writing '[cyan]{seq_id} - {feature}[/cyan]' to file" f" '[green]{output}[/green]'")
                    out.write(f">{name}.{feature}\n{aa}\n")
    elif outtype == 1:
        if not output.exists():
            output.mkdir()
        for seq_id, features in aa_dict.items():
            for feature, aa in features.items():
                log.info(f"Writing '[cyan]{seq_id} - {feature}[/cyan]' to file" f" \"[green]{output / f'{name}_{feature}.faa'}[/green]\"")
                with open(output / f"{name}_{feature}.faa", "a", encoding="utf-8") as out:
                    out.write(f">{seq_id}.{feature}\n{aa}\n")
