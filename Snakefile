"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import ast

import pandas as pd


configfile: "config.yaml"


# output files with aggregated counts
aggregated_counts_csvs = [
    "results/aggregated_counts/sars2_mito_aligned_by_run.csv",
    "results/aggregated_counts/sars2_mito_aligned_by_sample.csv",
]

rule all:
    input:
        "results/metadata/merged_metadata.csv",
        "results/fastqs_md5/check_vs_metadata.csv",
        "results/crits_christoph_data/check_sha512_vs_crits_christoph.csv",
        "results/mitochondrial_genomes/retained.csv",
        aggregated_counts_csvs,


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


def fastqs(_):
    """Return list of FASTQs."""
    fname = checkpoints.process_metadata.get().output.fastqs
    fastqs = pd.read_csv(fname)["fastqs"].tolist()
    assert len(fastqs) == len(set(fastqs))
    return fastqs


def accessions(_):
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
        fastq=protected("results/fastqs_preprocessed/{accession}.fq.gz"),
        json="results/fastqs_preprocessed/{accession}.json",
        html="results/fastqs_preprocessed/{accession}.html",
    threads: 3
    conda:
        "environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.fastq} \
            -w {threads} \
            --json {output.json} \
            --html {output.html}
        """


rule preprocess_paired_fastq:
    """Pre-process the FASTQ files."""
    input:
        unpack(lambda wc: accession_fastqs(wc)),
    output:
        r1=protected("results/fastqs_preprocessed/{accession}_R1.fq.gz"),
        r2=protected("results/fastqs_preprocessed/{accession}_R2.fq.gz"),
        json="results/fastqs_preprocessed/{accession}.json",
        html="results/fastqs_preprocessed/{accession}.html",
    threads: 3
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


rule get_mitochondrial_genomes:
    """Get the mitochondrial genomes."""
    params:
        url=config["mitochondrial_genomes"],
    output:
        fasta="results/mitochondrial_genomes/all.fasta",
    conda:
        "environment.yml"
    shell:
        "curl -s {params.url} | gzip -cd > {output.fasta}"


checkpoint process_mitochondrial_genomes:
    """Extract info for mitochondrial genomes to CSV and write a per-genome FASTA."""
    input:
        fasta=rules.get_mitochondrial_genomes.output.fasta,
    output:
        csv="results/mitochondrial_genomes/all.csv",
        per_genome_fastas=directory("results/mitochondrial_genomes/per_genome_fastas"),
        per_genome_fasta_list="results/mitochondrial_genomes/per_genome_fasta_list.csv",
    conda:
        "environment.yml"
    script:
        "scripts/process_mitochondrial_genomes.py"


def mitochondrial_genome_ids(_):
    """Get list of all mitochondrial genome IDs."""
    fname = checkpoints.process_mitochondrial_genomes.get().output.csv
    return pd.read_csv(fname)["id"].tolist()


rule mash_dist_mitochondrial_genome:
    """Compute Mash distance of one mitochondrial genome to all others."""
    input:
        genome_list=rules.process_mitochondrial_genomes.output.per_genome_fasta_list,
        genome="results/mitochondrial_genomes/per_genome_fastas/{mito_id}.fa",
    output:
        tsv="results/mitochondrial_genomes/per_genome_mash/{mito_id}.tsv",
    conda:
        "environment.yml"
    shell:
        "mash dist {input.genome} {input.genome_list} -l -t > {output.tsv} -s 5000"


rule mitochondrial_genomes_to_retain:
    """Retain sufficiently unique mitochondrial genomes."""
    input:
        mashes=lambda wc: [
            f"results/mitochondrial_genomes/per_genome_mash/{mito_id}.tsv"
            for mito_id in mitochondrial_genome_ids(wc)
        ],
        fasta=rules.get_mitochondrial_genomes.output.fasta,
        info_csv=rules.process_mitochondrial_genomes.output.csv,
    output:
        csv="results/mitochondrial_genomes/retained.csv",
        fasta="results/mitochondrial_genomes/retained.fasta",
    params:
        min_mash_dist=config["mitochondrial_genome_min_mash_dist"],
        to_keep=config["mitochondrial_genomes_to_keep"],
    log:
        notebook="results/mitochondrial_genomes/mitochondrial_genomes_to_retain.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/mitochondrial_genomes_to_retain.py.ipynb"


rule get_sars2_ref:
    """Get SARS-CoV-2 reference genome."""
    params:
        url=config["sars2_ref"],
    output:
        fasta="results/sars2_ref/ref.fa",
    conda:
        "environment.yml"
    shell:
        "curl -s {params.url} | gzip -cd > {output.fasta}"


rule minimap2_ref:
    """Build ``minimap2`` reference genome."""
    input:
        rules.get_sars2_ref.output.fasta,
        rules.mitochondrial_genomes_to_retain.output.fasta,
    output:
        fasta="results/minimap2_ref/ref.fa",
        mmi="results/minimap2_ref/ref.mmi",
    conda:
        "environment.yml",
    shell:
        """
        cat {input} > {output.fasta}
        minimap2 -x sr -d {output.mmi} {output.fasta}
        """


rule minimap2_alignments:
    """Align FASTQs for accession, aligning single or paired end as depending on data."""
    input:
        fastqs=lambda wc: (
            [f"results/fastqs_preprocessed/{wc.accession}.fq.gz"]
            if len(accession_fastqs(wc)) == 1
            else [
                f"results/fastqs_preprocessed/{wc.accession}_R1.fq.gz",
                f"results/fastqs_preprocessed/{wc.accession}_R2.fq.gz",
            ]
        ),
        ref=rules.minimap2_ref.output.mmi
    output:
        sam=temp("results/minimap2_alignments/{accession}.sam"),
        unsorted_bam=temp("results/minimap2_alignments/{accession}.bam"),
        bam=protected("results/minimap2_alignments/{accession}_sorted.bam"),
    threads: 3
    params:
        extra_threads=lambda _, threads: threads - 1
    conda:
        "environment.yml"
    shell:
        """
        minimap2 \
            -t {threads} \
            -x sr \
            --secondary=yes \
            --sam-hit-only \
            -a {input.ref} \
            {input.fastqs} \
            > {output.sam}
        samtools view -@ {params.extra_threads} -b -o {output.unsorted_bam} {output.sam}
        samtools sort -@ {params.extra_threads} -o {output.bam} {output.unsorted_bam}
        """


rule read_counts:
    """Count total and pre-processed reads."""
    input:
        fastqs=lambda wc: accession_fastqs(wc).values(),
        preprocessed_fastqs=lambda wc: (
            [f"results/fastqs_preprocessed/{wc.accession}.fq.gz"]
            if len(accession_fastqs(wc)) == 1
            else [
                f"results/fastqs_preprocessed/{wc.accession}_R1.fq.gz",
                f"results/fastqs_preprocessed/{wc.accession}_R2.fq.gz",
            ]
        ),
    output:
        counts="results/read_counts/{accession}.txt",
    conda:
        "environment.yml"
    shell:
        """
        echo "{wildcards.accession}" > {output.counts}
        zcat {input.fastqs} | echo $((`wc -l`/4)) >> {output.counts}
        if [ -s {input.preprocessed_fastqs[0]} ]; then
            zcat {input.preprocessed_fastqs} | echo $((`wc -l`/4)) >> {output.counts}
        else
            echo 0 >> {output.counts}
        fi
        """


rule agg_read_counts:
    """Aggregate total and pre-processed reads."""    
    input:
        lambda wc: [
            f"results/read_counts/{accession}.txt"
            for accession in accessions(wc)
        ],
    output:
        csv="results/read_counts/aggregated_counts.csv",
    conda:
        "environment.yml"
    script:
        "scripts/agg_read_counts.py"


rule tally_alignment_counts:
    """Tally alignment counts to each reference."""
    input:
        bamfile=rules.minimap2_alignments.output.bam,
    output:
        csv="results/tallied_alignment_counts/{accession}.csv",
    conda:
        "environment.yml"
    script:
        "scripts/tally_alignment_counts.py"


rule aggregate_all_counts:
    """Aggregate all counts."""
    input:
        tallied_counts=lambda wc: [
            f"results/tallied_alignment_counts/{accession}.csv"
            for accession in accessions(wc)
        ],
        mito_ref_info=rules.mitochondrial_genomes_to_retain.output.csv,
        metadata=rules.process_metadata.output.metadata,
        read_counts=rules.agg_read_counts.output.csv,
    output:
        **{os.path.splitext(os.path.basename(f))[0]: f for f in aggregated_counts_csvs}
    params:
        sars2_ref_id=config["sars2_ref_id"],
        mito_genomes_to_keep=config["mitochondrial_genomes_to_keep"],
    log:
        notebook="results/aggregated_counts/aggregate_all_counts.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/aggregate_all_counts.py.ipynb"
