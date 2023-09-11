from Bio.Seq import Seq

from AminoExtract.functions import log
from AminoExtract.reader import GffDataFrame, _split_attributes_column


def Reverse_complement(seq: str) -> Seq:
    """Reverse complement a sequence

    Parameters
    ----------
    seq : str
        The sequence to reverse complement

    Returns
    -------
    Seq
        The reverse complement of the input sequence
    """
    seq_obj = Seq(seq)  # create a Seq object from the sequence
    return seq_obj.reverse_complement()  # type: ignore as BioPythons seq object is weird with typehints


def extract_aminoacids(
    GFFobj: GffDataFrame,
    SeqRecords: list,
    keep_gaps: bool = False,
    verbose: bool = False,
) -> dict:
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

    log.info(
        "Extracting and translating the amino acid sequence(s) from the nucleotide sequence(s)"
    ) if verbose else None

    # create a dictionary with the SeqRecord.id as the key and the SeqRecord.seq as the value
    SeqDict = {record.id: record.seq for record in SeqRecords}
    # create an empty dictionary
    AA_dict = {record.id: {} for record in SeqRecords}

    tempdf = GFFobj.df.copy()
    if "attributes" in tempdf.columns:
        tempdf = _split_attributes_column(df=tempdf)

    # iterate through the dataframe
    for row in tempdf.itertuples():
        try:
            name = row.Name
        except AttributeError:
            log.warn(
                "No '[green]Name[/green]' attribute found in GFF records. Using '[cyan]ID[/cyan]' instead"
            ) if verbose else None
            name = f"ID-{row.ID}"

        # get the sequence ID from the row
        seq_id = row.seqid

        # get the start and end positions from the row
        # subtract 1 from the start position to account for 0-based indexing
        # end position is also 0-based but we don't need to subtract 1 because of the way python slices function
        start, end = row.start - 1, row.end

        # get the nucleotide sequence from the sequence dictionary
        NucSequence = SeqDict[seq_id]

        # get the sequence slice from the start to the end position
        seq_slice = (
            NucSequence[start:end].replace("-", "N")
            if keep_gaps
            else NucSequence[start:end].replace("-", "")
        )

        # convert the sequence slice to a string
        seq_slice_str = str(seq_slice)

        # reverse complement the sequence slice if the strand is negative
        if row.strand == "-":
            seq_slice_str = Reverse_complement(seq_slice_str)

        # create a Seq object from the sequence slice
        seq_slice_obj = Seq(seq_slice_str)

        # translate the sequence slice to amino acids
        AASequence = seq_slice_obj.translate(to_stop=True)

        # add the amino acid sequence to the amino_acid_dict dictionary
        AA_dict[seq_id][name] = AASequence
    return AA_dict
