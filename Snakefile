"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import pandas as pd


configfile: "config.yaml"


rule all:
    input:
        "results/fastqs_md5/check_vs_metadata.csv",
#        "_temp",


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
    conda:
        "environment.yml"
    script:
        "scripts/get_fastq.py"


rule fastq_checksum:
    """Get checksum for FASTQ file."""
    input:
        fastq=rules.get_fastq.output.fastq
    output:
        checksum="results/fastqs_{checksumtype}/{fastq}.{checksumtype}",
    conda:
        "environment.yml"
    shell:
        "{wildcards.checksumtype}sum {input.fastq} > {output.checksum}"


rule check_fastq_md5s:
    """Check MD5s for all downloaded FASTQs versus metadata, raise error if mismatch."""
    input:
        checksums=lambda wc: [
            f"results/fastqs_md5/{fastq}.md5" for (fastq, _, _) in fastq_info(wc)
        ],
        fastq_metadata=rules.process_metadata.output.fastqs,
    output:
        csv="results/fastqs_md5/check_vs_metadata.csv",
    conda:
        "environment.yml"
    script:
        "scripts/check_fastq_md5s.py"


rule agg_fastq_sha512:
    input:
        checksums=lambda wc: [
            f"results/fastqs_sha512/{fastq}.sha512" for (fastq, _, _) in fastq_info(wc)
        ],
    output:
        "_temp",
    conda:
        "environment.yml"
    shell:
        "echo not_implemented"
