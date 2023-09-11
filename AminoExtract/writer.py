import pathlib

from Bio import Seq

from AminoExtract.functions import log


def write_aa_file(
    AA_dict: dict[str, dict[str, Seq.Seq]],
    output: pathlib.Path,
    name: str,
    outtype: int,
) -> None:
    """Takes a dictionary of dictionaries of Seq.Seq objects, a pathlib.Path object, a string, and an
    integer, and writes the Seq.Seq objects to either a single file or multiple files

    Parameters
    ----------
    AA_dict : dict[str, dict[str, Seq.Seq]]
        a dictionary of dictionaries. The first dictionary is keyed by the sequence ID, and the second
    dictionary is keyed by the feature name. The value of the second dictionary is the amino acid
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
        # this should probably be changed so a user can specify the fasta header per record instead of assuming {name}.{feature}
        with open(output, "w") as out:
            for SeqID, features in AA_dict.items():
                for feature, aa in features.items():
                    log.info(
                        f"Writing '[cyan]{SeqID} - {feature}[/cyan]' to file '[green]{output}[/green]'"
                    )
                    out.write(f">{name}.{feature}\n{aa}\n")
    elif outtype == 1:
        if not output.exists():
            output.mkdir()
        for SeqID, features in AA_dict.items():
            for feature, aa in features.items():
                log.info(
                    f"Writing '[cyan]{SeqID} - {feature}[/cyan]' to file \"[green]{output / f'{name}_{feature}.faa'}[/green]\""
                )
                with open(output / f"{name}_{feature}.faa", "a") as out:
                    out.write(f">{SeqID}.{feature}\n{aa}\n")
