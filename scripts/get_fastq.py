"""Download the FASTQs."""


import subprocess

import pandas as pd


fastq = snakemake.wildcards.fastq

fastq_to_url = (
    pd.read_csv(snakemake.input.csv)
    # note we have to fix the URL specified in the metadata by substituting
    # "/gsa2/" for "/gsa/" as the medata page gives wrong URL
    .assign(download_urls=lambda x: x["download_urls"].str.replace("/gsa/", "/gsa2/", regex=False))
    .set_index("fastqs")
    ["download_urls"]
    .to_dict()
)

url = fastq_to_url[fastq]

print(f"Downloading {fastq} from {url} to {snakemake.output.fastq}")

cmds = ["curl", "-s", url, "-o", snakemake.output.fastq]

res = subprocess.run(cmds, capture_output=True, check=False)

if res.returncode != 0:
    raise ValueError(f"Command failed:\n{res}")
