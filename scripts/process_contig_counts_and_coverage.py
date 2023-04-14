"""Process ``coverm`` counts and coverage at species level."""


import pandas as pd

dfs = []

for tsv, (accession, d) in zip(snakemake.input.tsvs, snakemake.params.d.items()):

    sample = d["sample"]
    species = list(d["species"])

    assert accession in tsv

    df = (
        pd.read_csv(
            tsv,
            sep="\t",
            header=0,
            names=["reference", "n_reads", "covered_bases"],
        )
        .assign(species=lambda x: x["reference"].str.split("-").str[0])
        .groupby("species", as_index=False)
        .aggregate(
            aligned_contigs=pd.NamedAgg("n_reads", "sum"),
            covered_bases=pd.NamedAgg("covered_bases", "sum"),
        )
        .assign(accession=accession, sample=sample)
    )

    assert set(df["species"]) == set(species), f"{df['species'].unique()=}, {species=}"

    dfs.append(df)

pd.concat(dfs).to_csv(snakemake.output.csv, index=False)
