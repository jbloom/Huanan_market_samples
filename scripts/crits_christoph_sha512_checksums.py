"""Merge sample information and SHA-512 check sums from Crits-Christoph et al data."""


import pandas as pd


file_info = pd.read_csv(snakemake.input.file_info)

sample_names = pd.read_csv(snakemake.input.sample_names)


merged = (
    file_info
    .rename(columns={"GISAID Accession ID": "GISAID accession"})
    .merge(
        (
            sample_names
            [["GISAID accession", "Sample Name", "Sample Date"]]
            .rename(
                columns={"Sample Name": "Sample title", "Sample Date": "Collection date"},
            )
            .drop_duplicates()
        ),
        on="GISAID accession",
        validate="many_to_one",
        how="left",
    )
)

assert merged.notnull().all().all(), merged.notnull().all()

merged.to_csv(snakemake.output.checksums, index=False)
