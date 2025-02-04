"""
Copyright (C) 2022  RIVM

https://github.com/RIVM-bioinformatics/AminoExtract
"""

import sys

from AminoExtract.args import validate_args
from AminoExtract.filter import GFFRecordFilter, SequenceFilter
from AminoExtract.reader import SequenceReader
from AminoExtract.sequences import extract_aminoacids
from AminoExtract.writer import write_aa_file


def get_feature_name_attribute(
    input_gff: str, input_seq: str, feature_type: str
) -> dict[str, list[str]]:
    """
    This function takes a GFF file, a FASTA file and a feature type,
    and returns a dictionary of the feature names for each sequence in the FASTA file

    Parameters
    ----------
    input_gff : str
        the path to the gff file
    input_seq : str
        The path to the fasta file containing the sequences you want to extract features from.
    feature_type : str
        This is the type of feature you want to extract from the GFF file.
        For example, if you want to
    extract the CDS features, you would enter "CDS". Options are "CDS", "gene" and "all".

        Returns
        -------
            A dict with the sequence id as the key and a list of the feature names as the value.

    """
    reader = SequenceReader(verbose=False)
    gff = reader.read_gff(input_gff)
    seq = reader.read_fasta(input_seq)

    gff_filter = GFFRecordFilter(gff_records=gff, verbose=False)
    gff_records = gff_filter.apply_filters(seq_records=seq, feature_type=feature_type)

    seq_filter = SequenceFilter(seq_records=seq, verbose=False)
    filtered_seqs = seq_filter.filter_sequences(gff_records)

    seq_attributes: dict[str, list[str]] = {record.id: [] for record in filtered_seqs}
    for row in gff_records.df.itertuples():
        seq_attributes[row.seqid].append(row.Name)
    return seq_attributes


def main(provided_args: list[str] | None = None) -> None:
    """
    Extract amino acid sequences from GFF annotations and FASTA files.

    Reads GFF and FASTA files, filters features based on type, extracts amino acid
    sequences and writes results to file.

    Parameters
    ----------
    provided_args : list[str] | None, optional
        Command line arguments as list of strings. If None, sys.argv[1:] is used.

    Returns
    -------
    None
        Writes output file with amino acid sequences.

    Examples
    --------
    >>> main(['--input', 'sequences.fa', '--features', 'annot.gff'])
    >>> main(['--input', 'seqs.fa', '--features', 'genes.gff', '--feature-type', 'CDS'])

    Notes
    -----
    Exit codes:
        0: Success
        1: Something went wrong
    """
    if provided_args:
        args = validate_args(provided_args)
    else:
        args = validate_args(sys.argv[1:])

    reader = SequenceReader(verbose=args.verbose)
    gff_obj = reader.read_gff(args.features)
    fasta_records = reader.read_fasta(args.input)

    filter = GFFRecordFilter(gff_records=gff_obj, verbose=args.verbose)
    gff_obj = filter.apply_filters(
        seq_records=fasta_records, feature_type=args.feature_type
    )
    assert gff_obj.validate_dataframe(
        args.feature_type
    ), "Validation failed, either the GFF file is empty or the feature type is None"

    seq_filter = SequenceFilter(seq_records=fasta_records, verbose=args.verbose)
    seq_records = seq_filter.filter_sequences(gff_obj)

    aa_dict = extract_aminoacids(
        gff_obj=gff_obj, seq_records=seq_records, keep_gaps=args.keep_gaps, verbose=True
    )

    write_aa_file(aa_dict, args.output, args.name, args.outtype)
