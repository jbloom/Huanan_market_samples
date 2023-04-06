"""Trim 3' polyA."""


import Bio.SeqIO


seq = Bio.SeqIO.read(snakemake.input.fasta, "fasta")

if snakemake.params.trim_3_polyA:
    while seq[-1] == "A":
        seq = seq[: -1]

Bio.SeqIO.write([seq], snakemake.output.fasta, "fasta")
