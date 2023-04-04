"""Aggregate read counts."""


import pandas as pd


records = []
for fname in snakemake.input:
    with open(fname) as f:
        run, n, n_preprocessed = [line.strip() for line in f]
    records.append((run, int(n), int(n_preprocessed)))

(
    pd.DataFrame(
        records,
        columns=["Run accession", "total_reads", "preprocessed_reads"],
    )
    .to_csv(snakemake.output.csv, index=False)
)
