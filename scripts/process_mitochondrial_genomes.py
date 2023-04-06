"""Get CSV with information on mitochondrial genomes."""


import os
import re

import Bio.Entrez
import Bio.SeqIO

import pandas as pd


regex = re.compile(
    "(?P<species>.+) "
    + "(?:mitochondrion|mitochondria|mitochondrial|chloroplast|kinetoplast|genome assembly).+"
)


os.makedirs(snakemake.output.per_genome_seqs, exist_ok=True)

extra_seqs = [
    Bio.SeqIO.read(Bio.Entrez.efetch(db="nuccore", id=acc, rettype="gb"), "genbank")
    for acc in snakemake.params.extra_genomes
]

records = []
outfiles = []
seqs = []
for seq in list(Bio.SeqIO.parse(snakemake.input.gb, "genbank")) + extra_seqs:
    desc = seq.description
    m = regex.fullmatch(desc)
    if m:
        species = m.group("species")
    elif desc == "Trichoderma hamatum, complete genome":
        species = desc.split(",")[0]
    else:
        raise ValueError(f"Cannot recognize\n{desc}\nfor {seq}")
    try:
        outfasta = os.path.join(snakemake.output.per_genome_seqs, f"{seq.id}.fa")
        Bio.SeqIO.write([seq], outfasta, "fasta")
        outgb = os.path.join(snakemake.output.per_genome_seqs, f"{seq.id}.gb")
        Bio.SeqIO.write([seq], outgb, "genbank")
        records.append((seq.id, species))
        outfiles.append(outfasta)
        seqs.append(seq)
    except Bio.Seq.UndefinedSequenceError:
        print(f"No sequence for {seq.id} {desc}")

pd.DataFrame(records, columns=["id", "species"]).to_csv(
    snakemake.output.csv, index=False,
)

with open(snakemake.output.per_genome_fasta_list, "w") as f:
    f.write("\n".join(outfiles))

Bio.SeqIO.write(seqs, snakemake.output.fasta, "fasta")
