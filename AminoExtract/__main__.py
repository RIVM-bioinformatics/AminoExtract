"""
Copyright (C) 2022  RIVM

https://github.com/RIVM-bioinformatics/AminoExtract
"""

import sys

from AminoExtract.args import validate_args
from AminoExtract.filter import empty_dataframe, filter_gff, filter_sequences
from AminoExtract.reader import read_fasta, read_gff
from AminoExtract.sequences import extract_aminoacids
from AminoExtract.writer import write_aa_file


def get_feature_name_attribute(
    input_gff: str, input_seq: str, feature_type: str
) -> dict[str, list[str]]:
    """This function takes a GFF file, a FASTA file, and a feature type, and returns a dictionary of the
        feature names for each sequence in the FASTA file

        Parameters
        ----------
        input_gff : str
            the path to the gff file
        input_seq : str
            The path to the fasta file containing the sequences you want to extract features from.
        feature_type : str
            This is the type of feature you want to extract from the GFF file. For example, if you want to
    extract the CDS features, you would enter "CDS". Options are "CDS", "gene" and "all".

        Returns
        -------
            A dictionary with the sequence id as the key and a list of the feature names as the value.

    """
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

    AA_dict = extract_aminoacids(
        GFFobj=GFF_Obj, SeqRecords=SeqRecords, keep_gaps=args.keep_gaps, verbose=True
    )

    write_aa_file(AA_dict, args.output, args.name, args.outtype)
