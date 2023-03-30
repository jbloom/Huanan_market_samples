"""Check downloaded FASTQ MD5 checksums versus metadata."""


import os

import pandas as pd


metadata = (
    pd.read_csv(snakemake.input.fastq_metadata)
    .drop(columns="download_urls")
    .rename(
        columns={
            "fastqs": "fastq",
            "md5_checksums": "expected_md5_checksum",
        }
    )
)

records = []
for fname in snakemake.input.checksums:
    with open(fname) as f:
        md5, fullpath = f.read().strip().split()
        records.append((os.path.basename(fullpath), md5))
actual = pd.DataFrame(records, columns=["fastq", "actual_md5_checksum"])

assert set(actual["fastq"]) == set(metadata["fastq"])

merged = (
    metadata
    .merge(actual, on="fastq", validate="one_to_one", how="outer")
    .assign(
        md5_matches=lambda x: x["actual_md5_checksum"] == x["expected_md5_checksum"],
    )
)

if (not merged["md5_matches"].all()) or merged["md5_matches"].isnull().any():
    raise ValueError(
        "Some FASTQs have mismatched MD5 checksums:\n"
        + str(merged.query("md5_matches != True"))
    )

merged.to_csv(snakemake.output.csv, index=False)
