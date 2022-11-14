"""
Copyright (C) 2022  Florian Zwagemaker

https://github.com/florianzwagemaker/AminoExtract
"""

import sys

import pandas as pd
from rich import print

from AminoExtract.args import validate_args
from AminoExtract.filter import filter_gff, filter_sequences
from AminoExtract.functions import log
from AminoExtract.reader import GffDataFrame, read_fasta, read_gff
from AminoExtract.sequences import Extract_AminoAcids
from AminoExtract.writer import write_aa_file


def empty_dataframe(
    frame: pd.DataFrame = pd.DataFrame(), feature_type: str | None = None
) -> bool:
    if frame.empty or feature_type is None:
        log.warn(
            f"The GFF file is empty after filtering.\nThis might mean that there are no records within the GFF that match the sequence ID(s) in the given Fasta file.\nThis could also mean that there are no records within the GFF that match the feature type '[cyan]{feature_type}[/cyan]'.\nPlease check your inputs and try again."
        )
        return True
    return False


def get_feature_name_attribute(
    input_gff: str, input_seq: str, feature_type: str
) -> dict[str, list[str]]:
    gff = read_gff(file=input_gff, verbose=False)
    seq = read_fasta(input_seq)
    gff_records = filter_gff(GffRecords=gff, SeqRecords=seq, feature_type=feature_type)
    FiltSequences = filter_sequences(gff, seq)
    SeqAttributes = {record.id: [] for record in FiltSequences}
    for row in gff_records.df.itertuples():
        SeqAttributes[row.seqid].append(row.Name)
    return SeqAttributes


def main() -> None:
    args = validate_args(sys.argv[1:])

    GFF_Obj = read_gff(file=args.features, verbose=True)
    fasta_records = read_fasta(args.input, verbose=True)

    GFF_Obj = filter_gff(
        GffRecords=GFF_Obj,
        SeqRecords=fasta_records,
        feature_type=args.feature_type,
        verbose=True,
    )
    sys.exit(1) if empty_dataframe(GFF_Obj.df, args.feature_type) else None
    SeqRecords = filter_sequences(GFF_Obj, fasta_records)

    log.info(
        "Extracting and translating the amino acid sequence(s) from the nucleotide sequence(s)"
    )
    AA_dict = Extract_AminoAcids(GFFobj=GFF_Obj, SeqRecords=SeqRecords)

    write_aa_file(AA_dict, args.output, args.name, args.outtype)
