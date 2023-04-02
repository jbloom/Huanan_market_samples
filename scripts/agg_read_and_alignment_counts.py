"""Aggregate read and alignment counts."""


import pandas as pd


records = []
for fname in snakemake.input:
    with open(fname) as f:
        run, n, n_preprocessed, n_aligned = [line.strip() for line in f]
    records.append((run, int(n), int(n_preprocessed), int(n_aligned)))

(
    pd.DataFrame(
        records,
        columns=["Run accession", "total_reads", "preprocessed_reads", "aligned_reads"],
    )
    .to_csv(snakemake.output.csv, index=False)
)
