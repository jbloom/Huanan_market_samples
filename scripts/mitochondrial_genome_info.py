"""Get CSV with information on mitochondrial genomes."""


import re

import Bio.SeqIO

import pandas as pd


regex = re.compile(
    "(?P<id>[^\s]+) (?P<species>.+) "
    + "(?:mitochondrion|mitochondria|mitochondrial|chloroplast|kinetoplast|genome assembly).+"
)


records = []
for seq in Bio.SeqIO.parse(snakemake.input.fasta, "fasta"):
    desc = seq.description
    m = regex.fullmatch(desc)
    if m:
        records.append((m.group("id"), m.group("species")))
    elif desc == "NC_036144.1 Trichoderma hamatum, complete genome":
        records.append(desc.split(",")[0].split(maxsplit=1))
    else:
        raise ValueError(f"Cannot recognize {desc}")

pd.DataFrame(records, columns=["id", "species"]).to_csv(
    snakemake.output.csv, index=False,
)
