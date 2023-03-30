"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import ast

import pandas as pd


configfile: "config.yaml"


rule all:
    input:
        "results/fastqs_md5/check_vs_metadata.csv",
        "results/crits_christoph_data/check_sha512_vs_crits_christoph.csv",
        "_temp",


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


def fastqs(wildcards):
    """Return list of FASTQs."""
    fname = checkpoints.process_metadata.get().output.fastqs
    fastqs = pd.read_csv(fname)["fastqs"].tolist()
    assert len(fastqs) == len(set(fastqs))
    return fastqs


def accessions(wildcards):
    """Return list of run accessions."""
    fname = checkpoints.process_metadata.get().output.metadata
    accs = pd.read_csv(fname)["Run accession"].tolist()
    assert len(accs) == len(set(accs))
    return accs


def accession_fastqs(wildcards):
    """Given {accession} returns FASTQs as dict keyed by "r1" and (optionally) "r2"."""
    fname = checkpoints.process_metadata.get().output.metadata
    acc_fastqs = (
        pd.read_csv(fname, converters={"fastqs": ast.literal_eval})
        .set_index("Run accession")
        ["fastqs"]
        .to_dict()
        [wildcards.accession]
    )
    if 1 <= len(acc_fastqs) <= 2:
        return dict(zip(["r1", "r2"], [f"results/fastqs/{f}" for f in acc_fastqs]))
    else:
        raise ValueError(f"Not 1 or 2 FASTQs\n{acc_fastqs=}\n{wildcards.accession=}")


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
        checksums=lambda wc: [f"results/fastqs_md5/{fastq}.md5" for fastq in fastqs(wc)],
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
            f"results/fastqs_sha512/{fastq}.sha512" for fastq in fastqs(wc)
        ],
        checksums_nogz=lambda wc: [
            f"results/fastqs_sha512/{fastq}_unzipped.sha512" for fastq in fastqs(wc)
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


rule preprocess_single_fastq:
    """Pre-process the FASTQ files."""
    input:
        unpack(lambda wc: accession_fastqs(wc)),
    output:
        fastq="results/fastqs_preprocessed/{accession}.fq.gz",
    threads: 2
    conda:
        "environment.yml"
    shell:
        "fastp -i {input.r1} -o {output.fastq} -w {threads}"


rule preprocess_paired_fastq:
    """Pre-process the FASTQ files."""
    input:
        unpack(lambda wc: accession_fastqs(wc)),
    output:
        r1="results/fastqs_preprocessed/{accession}_R1.fq.gz",
        r2="results/fastqs_preprocessed/{accession}_R2.fq.gz",
        json="results/fastqs_preprocessed/{accession}.json",
        html="results/fastqs_preprocessed/{accession}.html",
    threads: 2
    conda:
        "environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            --json {output.json} \
            --html {output.html}
        """


rule align_fastq:
    """Align FASTQs for an accession, aligning single or paired end as depending on data."""
    input:
        unpack(
            lambda wc: (
                {"r1": f"results/fastqs_preprocessed/{wc.accession}.fq.gz"}
                if len(accession_fastqs(wc)) == 1
                else {
                    "r1": f"results/fastqs_preprocessed/{wc.accession}_R1.fq.gz",
                    "r2": f"results/fastqs_preprocessed/{wc.accession}_R2.fq.gz",
                }
            )
        ),
    output:
        "results/alignment/{accession}",
    conda:
        "environment.yml"
    shell:
        "echo not_implemented"


rule agg_alignments:
    input:
        lambda wc: [f"results/alignment/{accession}" for accession in accessions(wc)],
    output:
        "_temp"
    conda:
        "environment.yml"
    shell:
        "echo not_implemented"
