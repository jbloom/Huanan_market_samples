"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import pandas as pd


configfile: "config.yaml"


rule all:
    input:
        "results/fastqs_md5/check_vs_metadata.csv",
        "results/crits_christoph_data/check_sha512_vs_crits_christoph.csv",


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
    """Get checksum for FASTQ file, both gzipped and not (input assumed gzipped)."""
    input:
        fastq=rules.get_fastq.output.fastq
    output:
        checksum="results/fastqs_{checksumtype}/{fastq}.{checksumtype}",
        checksum_nogz="results/fastqs_{checksumtype}/{fastq}_unzipped.{checksumtype}",
    params:
        fastq_nogz=lambda wc: os.path.splitext(wc.fastq)[0],
    conda:
        "environment.yml"
    shell:
        """
        {wildcards.checksumtype}sum {input.fastq} > {output.checksum}
        gzip -cd {input.fastq} | {wildcards.checksumtype}sum > {output.checksum_nogz}
        sed -i 's/-/{params.fastq_nogz}/g' {output.checksum_nogz}
        """


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


rule get_crits_christoph_data:
    """Get the data from the Crits-Christoph et al report."""
    params:
        url=lambda wc: config["crits_christoph_data"][wc.cc_data],
    output:
        csv="results/crits_christoph_data/{cc_data}.csv",
    conda:
        "environment.yml"
    shell:
        "curl -s {params.url} -o {output.csv}"


rule crits_christoph_sha512_checksums:
    """Aggregate Crits-Christoph et al SHA-512 checksums and sample names/dates."""
    input:
        file_info="results/crits_christoph_data/file_info.csv",
        sample_names="results/crits_christoph_data/sample_names.csv",
    output:
        checksums="results/crits_christoph_data/sha512_checksums.csv",
    conda:
        "environment.yml"
    script:
        "scripts/crits_christoph_sha512_checksums.py"


rule check_sha512_vs_crits_christoph:
    """Check the SHA-512 checksums of downloaded FASTQs vs those from Crits-Christoph."""
    input:
        checksums=lambda wc: [
            f"results/fastqs_sha512/{fastq}.sha512" for (fastq, _, _) in fastq_info(wc)
        ],
        checksums_nogz=lambda wc: [
            f"results/fastqs_sha512/{fastq}_unzipped.sha512"
            for (fastq, _, _) in fastq_info(wc)
        ],
        cc_checksums=rules.crits_christoph_sha512_checksums.output.checksums,
        metadata=rules.process_metadata.output.metadata,
    output:
        csv="results/crits_christoph_data/check_sha512_vs_crits_christoph.csv",
    log:
        notebook="results/crits_christoph_data/check_sha512_vs_crits_christoph.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/check_sha512_vs_crits_christoph.py.ipynb"
