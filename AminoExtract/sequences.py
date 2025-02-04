from Bio.Seq import Seq

from AminoExtract.functions import log
from AminoExtract.reader import GffDataFrame


def extract_aminoacids(
    gff_obj: GffDataFrame,
    seq_records: list,
    keep_gaps: bool = False,
    verbose: bool = False,
) -> dict[str, dict[str, Seq]]:
    """
    Extract amino acids from the SeqRecord objects based on the start and end positions of the GFFobj.df dataframe

    Parameters
    ----------
    GFFobj : GffDataFrame object
        GffDataFrame
    SeqRecords : list
        list of SeqRecord objects
    keep_gaps : bool, optional
        If True, gaps ('-') in the nucleotide sequence will not be removed before AA translation.
        If False, gaps will be removed from the nucleotide sequence before translation.
        (default is False)
    verbose : bool, optional
        bool = False

    Returns
    -------
        A dictionary with the sequence ID as the key and a dictionary as the value. The dictionary has the
    name of the feature as the key and the amino acid sequence as the value.

    """

    if verbose:
        log.info(
            "Extracting and translating the amino acid sequence(s) from the nucleotide sequence(s)"
        )

    # create a dictionary with the SeqRecord.id as the key and the SeqRecord.seq as the value
    SeqDict = {record.id: record.seq for record in seq_records}
    # create an empty dictionary
    aa_dict: dict[str, dict[str, Seq]] = {record.id: {} for record in seq_records}

    tempdf = gff_obj.df.copy()

    # iterate through the dataframe
    for row in tempdf.itertuples():
        try:
            name = row.Name
        except AttributeError:
            if verbose:
                log.warning(
                    "No '[green]Name[/green]' attribute found in GFF records. Using '[cyan]ID[/cyan]' instead"
                )
            name = f"ID-{row.ID}"

        # get the sequence ID from the row
        seq_id = row.seqid

        # get the start and end positions from the row
        # subtract 1 from the start position to account for 0-based indexing
        # end position is also 0-based but we don't need to subtract 1 because of the way python slices function
        if gff_obj.splicing_table is None:
            start, end = row.start - 1, row.end

            # get the nucleotide sequence from the sequence dictionary
            NucSequence = SeqDict[seq_id]

            # get the sequence slice from the start to the end position
            full_seq = (
                NucSequence[start:end].replace("-", "N")
                if keep_gaps
                else NucSequence[start:end].replace("-", "")
            )
        else:
            splicing_info = gff_obj.splicing_table.loc[row.ID]
            exon_count = len(splicing_info.CDSes)
            full_seq_str = ""
            for i in range(exon_count):
                start, end = splicing_info.CDSes[i]
                if i == 0:
                    start = start - 1

                NucSequence = SeqDict[seq_id]
                seq_part = (
                    NucSequence[start:end].replace("-", "N")
                    if keep_gaps
                    else NucSequence[start:end].replace("-", "")
                )
                full_seq_str += seq_part
            full_seq = Seq(full_seq_str)

        if row.strand == "-":
            full_seq = full_seq.reverse_complement()

        aa_sequence = full_seq.translate(to_stop=True)

        aa_dict[seq_id][name] = aa_sequence
    return aa_dict
