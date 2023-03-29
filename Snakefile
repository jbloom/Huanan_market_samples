"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import pandas as pd


configfile: "config.yaml"


rule all:
    input:
        "_temp"


checkpoint process_metadata:
    """Process metadata sample sheet to aggregate metadata and get list of FASTQs."""
    input:
        excel="data/CRA010170.xlsx",
    output:
        samples="results/metadata/samples.csv",
        experiments="results/metadata/experiments.csv",
        runs="results/metadata/runs.csv",
        metadata="results/metadata/merged_metadata.csv",
        fastqs="results/metadata/fastqs.csv",
    conda:
        "environment.yml"
    script:
        "scripts/process_metadata.py"


def fastq_info(wildcards):
    """Return list of tuples of FASTQs, their URLs, and MD5 checksums."""
    fname = checkpoints.process_metadata.get().output.fastqs
    return list(pd.read_csv(fname).itertuples(index=False))


rule get_fastq:
    """Download a FASTQ file."""
    input:
        csv=rules.process_metadata.output.fastqs,
    output:
        fastq=protected("results/fastqs/{fastq}"),
    script:
        "scripts/get_fastq.py"


rule agg_fastqs:
    input:
        lambda wc: [f"results/fastqs/{fastq}" for (fastq, _, _) in fastq_info(wc)],
    output:
        "_temp"
    shell:
        "echo not_implemented"
