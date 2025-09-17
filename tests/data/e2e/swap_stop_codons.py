"""Quick script to swap the stop codons of the edenprotein.1 CDS in the complex_input.fa file"""

import pandas as pd

with open("tests/data/e2e/complex_input.fa", "r", encoding="utf-8") as f:
    complex_output = f.read()

column_names = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
table = pd.read_csv(
    "tests/data/e2e/complex_input.gff",
    sep="\t",
    comment="#",
    header=None,
    names=column_names,
)


cds_table = table[table["type"] == "CDS"]
eden1_table = cds_table[cds_table["attributes"].str.contains("Name=edenprotein.1")]
eden1_table["start"] = eden1_table["start"].astype(int) - 1  # Convert to 0-based
eden1_table["end"] = eden1_table["end"].astype(int) - 1

just_nucl = complex_output.strip().split(">")[1:][0].split("\n", 1)[1:][0].replace("\n", "")
eden1_table["nucl"] = eden1_table.apply(lambda row: just_nucl[row["start"] : row["end"]], axis=1)

# combine all the CDS sequences into one
OLD_EDEN1_SEQ = "".join(eden1_table["nucl"].tolist())

NEW_EDEN1_SEQ = "ATG" + OLD_EDEN1_SEQ[3:-3] + "TAA"

# edit the nucleotide sequence
lowest = min(eden1_table["start"])
highest = max(eden1_table["end"])
just_nucl_2 = just_nucl[:lowest] + "ATG" + just_nucl[lowest + 3 : highest + 1] + "TAA" + just_nucl[highest + 4 :]
