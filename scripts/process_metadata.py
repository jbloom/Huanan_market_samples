"""Process metadata from Excel."""


import os
import itertools

import pandas as pd


excel_metadata = pd.read_excel(snakemake.input.excel, sheet_name=None)

assert set(excel_metadata) == {"Sample", "Experiment", "Run"}

samples = excel_metadata["Sample"]
samples.to_csv(snakemake.output.samples, index=False)

runs = excel_metadata["Run"]
runs.to_csv(snakemake.output.runs, index=False)

experiments = excel_metadata["Experiment"]
experiments.to_csv(snakemake.output.experiments, index=False)

assert set(samples["Sample name"]) == set(experiments["BioSample name"])
assert set(samples.columns).intersection(experiments.columns) == {"Accession", "ID"}

metadata = (
    samples
    .drop(columns=["ID", "Public description"])
    .rename(columns={"Accession": "BioSample accession"})
    .merge(
        (
            experiments
            .drop(columns="ID")
            .rename(
                columns={
                    "Accession": "Experiment accession",
                    "Library Construction / Experimental Design": "description",
                },
            )
        ),
        on="BioSample accession",
        validate="one_to_many",
    )
    .drop(
        columns=[
            "accession_in_other_db",
            "other_db",
            "other_db_url",
            "BioProject accession",
        ]
    )
    .merge(
        runs.drop(columns="ID").rename(columns={"Accession": "Run accession"}),
        on="Experiment accession",
        validate="one_to_many",
    )
)

# drop uninformative columns
uninformative_columns = [
    c for c in metadata.columns if metadata[c].nunique(dropna=False) <= 1
]
metadata = metadata.drop(columns=uninformative_columns)

# drop duplicated columns with different names
dup_cols = [
    cols[1] for cols in itertools.combinations(metadata.columns, 2)
    if (metadata[cols[0]] == metadata[cols[1]]).all()
]
metadata = metadata.drop(columns=dup_cols)

# drop other irrelevant columns
metadata = metadata.drop(columns=["Planned number of cycles", "Run title"])

# group download paths and MD5 checksums
metadata = (
    metadata
    .drop(columns=[col for col in metadata.columns if col.startswith("File name")])
    .assign(
        download_urls=lambda x: (
            x[[c for c in metadata.columns if c.startswith("DownLoad")]].apply(
                lambda vals: [val for val in vals if not pd.isnull(val)],
                axis=1,
            )
        ),
        md5_checksums=lambda x: (
            x[[c for c in metadata.columns if c.startswith("MD5")]].apply(
                lambda vals: [val for val in vals if not pd.isnull(val)],
                axis=1,
            )
        ),
        fastqs=lambda x: x["download_urls"].map(lambda fs: [os.path.basename(f) for f in fs]),
    )
    .drop(columns=[col for col in metadata.columns if col.startswith("DownLoad")])
    .drop(columns=[col for col in metadata.columns if col.startswith("MD5")])
)

# get most relevant columns first
first_cols = [
    "Sample name",
    "Sample title",
    "description",
    "Collection date",
    "Run accession",
    "BioSample accession",
    "fastqs",
]

metadata = metadata[first_cols + [c for c in metadata.columns if c not in first_cols]]

print(f"Writing merged metadata to {snakemake.output.metadata}")
metadata.to_csv(snakemake.output.metadata, index=False)

fastqs_to_get = (
    metadata
    [["fastqs", "download_urls", "md5_checksums"]]
    .explode(column=["fastqs", "download_urls", "md5_checksums"])
)
assert fastqs_to_get["fastqs"].nunique() == len(fastqs_to_get)
print(f"Writing FASTQs to get to {snakemake.output.fastqs}")
fastqs_to_get.to_csv(snakemake.output.fastqs, index=False)
