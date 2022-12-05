import pandas as pd

from AminoExtract.functions import log
from AminoExtract.reader import GffDataFrame


def filter_name(frame: pd.DataFrame, Sequence_IDs: list) -> pd.DataFrame:
    """Takes a dataframe and a list of Biopythons SeqRecord objects, and returns a dataframe with only the rows that
    match the sequence IDs within the SeqRecord objects

    Parameters
    ----------
    frame : pd.DataFrame
        the dataframe you want to filter
    seq_ids : list
        a list of SeqRecord objects

    Returns
    -------
        A dataframe with the rows that have the seqids in the fasta_ids list.

    """
    return frame[frame["seqid"].isin(Sequence_IDs)]


def filter_feature_type(frame: pd.DataFrame, feature_type: str) -> pd.DataFrame:
    """Return a new dataframe containing only the rows of the input dataframe that have the specified
    feature type.

    Parameters
    ----------
    frame : pd.DataFrame
        the dataframe to filter
    feature_type : str
        The type of feature you want to filter for.

    Returns
    -------
        A dataframe with only the rows that containing the specified feature type.

    """
    # keep all rows if feature_type is "all"
    return frame if feature_type == "all" else frame[frame["type"] == feature_type]


def filter_gff(
    GffRecords: GffDataFrame, SeqRecords: list, feature_type: str, verbose: bool = False
) -> GffDataFrame:
    """Filter the GFF dataframe by feature type and sequence name

    Parameters
    ----------
    GffRecords : GffDataFrame object
        GffDataFrame
    SeqRecords : list
        list of SeqRecord objects
    feature_type : str
        str

    Returns
    -------
        A GffDataFrame object.

    """
    Sequence_IDs = [record.id for record in SeqRecords]
    log.info(
        f"Filtering GFF records to only contain the following information:\n * Feature type: '[green]{feature_type}[/green]'\n * Sequence IDs: '[green]{', '.join(Sequence_IDs)}[/green]'"
    ) if verbose else None
    GffRecords.df = filter_feature_type(
        filter_name(GffRecords.df, Sequence_IDs), feature_type
    )
    return GffRecords


def filter_sequences(gff: GffDataFrame, SeqRecords: list):
    """Takes a GffDataFrame object and a list of SeqRecord objects, and returns a list of SeqRecord objects
    that only contain the sequences that are specified in the GffDataFrame object.

    Parameters
    ----------
    gff : GffDataFrame
        The GffDataFrame object to filter
    SeqRecords : list
        A list of SeqRecord objects

    Returns
    -------
        A list of SeqRecord objects that only contain the sequences that are specified in the GffDataFrame object.

    """
    return [record for record in SeqRecords if record.id in gff.df["seqid"].tolist()]


def empty_dataframe(
    frame: pd.DataFrame = pd.DataFrame(), feature_type: str | None = None
) -> bool:
    if frame.empty or feature_type is None:
        log.warn(
            f"The GFF file is empty after filtering.\nThis might mean that there are no records within the GFF that match the sequence ID(s) in the given Fasta file.\nThis could also mean that there are no records within the GFF that match the feature type '[cyan]{feature_type}[/cyan]'.\nPlease check your inputs and try again."
        )
        return True
    return False
