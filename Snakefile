"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import ast

import pandas as pd

import yaml


wildcard_constraints:
    accession="CRR\d+",

configfile: "config.yaml"

with open(config["docs_plot_annotations"]) as f:
    docs_plot_annotations = yaml.safe_load(f)

# output files with aggregated counts
aggregated_counts_csvs = {
    name: f"results/aggregated_counts/{name}.csv"
    for name in [
        "sars2_aligned_by_run",
        "sars2_aligned_by_sample",
        "mito_composition_by_run",
        "mito_composition_by_sample",
    ]
}

# plot files
plot_htmls = {
    name: f"results/plots/{name}.html"
    for name in [
        "crits_christoph_vs_current_run_corr",
        "mito_composition",
        "genomic_contig_composition",
        "sars2_aligned",
        "sars2_aligned_vertical",
        "per_species_corr_faceted",
        "per_species_corr_single",
        "overall_corr",
    ]
}

rule all:
    input:
        "results/metadata/merged_metadata.csv",
        "results/fastqs_md5/check_vs_metadata.csv",
        "results/crits_christoph_data/check_sha512_vs_crits_christoph.csv",
        "results/mitochondrial_genomes/retained.csv",
        aggregated_counts_csvs.values(),
        "results/contigs/counts_and_coverage/processed_counts.csv",
        plot_htmls.values(),
        "results/plots/susceptible_table.csv",
        "results/plots/susceptible_mammal_table.csv",
        "results/plots/raccoon_dog_long.csv",
        "results/plots/susceptible_table.tex",
        "results/plots/susceptible_mammal_table.tex",
        "results/rt_qpcr/rt_qpcr.csv",
        "results/rt_qpcr/ct_vs_reads.html",
        expand("docs/{plot}.html", plot=list(plot_htmls) + ["ct_vs_reads"]),
        "docs/index.html",


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
            --trim_poly_g \
            --trim_poly_x \
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
            --trim_poly_g \
            --trim_poly_x \
            --json {output.json} \
            --html {output.html}
        """


rule get_mitochondrial_genomes:
    """Get the mitochondrial genomes."""
    params:
        url=config["mitochondrial_genomes"],
    output:
        gb="results/mitochondrial_genomes/all.gb",
    conda:
        "environment.yml"
    shell:
        "curl -s {params.url} | gzip -cd > {output.gb}"


checkpoint process_mitochondrial_genomes:
    """Extract info for mitochondrial genomes to CSV and write a per-genome FASTA."""
    input:
        gb=rules.get_mitochondrial_genomes.output.gb,
    output:
        csv="results/mitochondrial_genomes/all.csv",
        fasta="results/mitochondrial_genomes/all.fa",
        per_genome_seqs=directory("results/mitochondrial_genomes/per_genome_seqs"),
        per_genome_fasta_list="results/mitochondrial_genomes/per_genome_fasta_list.csv",
    params:
        extra_genomes=config["extra_mitochondrial_genomes"],
    conda:
        "environment.yml"
    script:
        "scripts/process_mitochondrial_genomes.py"


def mitochondrial_genome_ids(_):
    """Get list of all mitochondrial genome IDs."""
    fname = checkpoints.process_mitochondrial_genomes.get().output.csv
    return pd.read_csv(fname)["id"].tolist()


rule mitochondrial_genome_taxa:
    """Get the taxa (phylum and subphylum) for a mitochondrial genome."""
    input:
        gb="results/mitochondrial_genomes/per_genome_seqs/{mito_id}.gb",
    output:
        csv="results/mitochondrial_genomes/taxa/{mito_id}.csv",
    conda:
        "environment.yml"
    script:
        "scripts/mitochondrial_genome_taxa.py"


rule mash_dist_mitochondrial_genome:
    """Compute Mash distance of one mitochondrial genome to all others."""
    input:
        genome_list=rules.process_mitochondrial_genomes.output.per_genome_fasta_list,
        genome="results/mitochondrial_genomes/per_genome_seqs/{mito_id}.fa",
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
        taxa=lambda wc: [
            f"results/mitochondrial_genomes/taxa/{mito_id}.csv"
            for mito_id in mitochondrial_genome_ids(wc)
        ],
        fasta=rules.process_mitochondrial_genomes.output.fasta,
        info_csv=rules.process_mitochondrial_genomes.output.csv,
    output:
        csv="results/mitochondrial_genomes/retained.csv",
        fasta="results/mitochondrial_genomes/retained.fasta",
    params:
        taxa_to_keep=config["mitochondrial_taxa_to_keep"],
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
        fasta="results/sars2_ref/untrimmed_ref.fa",
    conda:
        "environment.yml"
    shell:
        "curl -s {params.url} | gzip -cd > {output.fasta}"


rule trim_sars2_ref:
    input:
        fasta=rules.get_sars2_ref.output.fasta,
    output:
        fasta="results/sars2_ref/trimmed_ref.fa",
    params:
        trim_3_polyA=config["sars2_ref_trim_3_polyA"],
    conda:
        "environment.yml"
    script:
        "scripts/trim_sars2_ref.py"


rule minimap2_ref:
    """Build ``minimap2`` reference genome."""
    input:
        rules.trim_sars2_ref.output.fasta,
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
        sam=temp("results/minimap2_alignments/not_mapq_filtered/{accession}.sam"),
        unsorted_bam=temp("results/minimap2_alignments/not_mapq_filtered/{accession}.bam"),
        bam=protected("results/minimap2_alignments/not_mapq_filtered/{accession}_sorted.bam"),
    threads: 3
    params:
        extra_threads=lambda _, threads: threads - 1,
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
        samtools view \
            -@ {params.extra_threads} \
            -b \
            -o {output.unsorted_bam} \
            {output.sam}
        samtools sort -@ {params.extra_threads} -o {output.bam} {output.unsorted_bam}
        """


rule mapq_filter_bam:
    """Filter BAM for only reads above a certain mapping quality."""
    input:
        bam=rules.minimap2_alignments.output.bam,
    output:
        bam="results/minimap2_alignments/mapq_filtered/{accession}_sorted.bam",
    params:
        min_mapq=config["min_mapq"],
    conda:
        "environment.yml"
    shell:
        "samtools view -q {params.min_mapq} -b -o {output.bam} {input.bam}"


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


rule alignment_counts_and_coverage:
    """Get alignment counts to each reference."""
    input:
        bam=rules.mapq_filter_bam.output.bam,
    output:
        tsv="results/alignment_counts_and_coverage/{accession}.tsv",
    params:
        **config["coverm_flags"],
    conda:
        "environment.yml"
    shell:
        """
        coverm contig \
            -b {input.bam} \
            -m count covered_bases \
            --min-read-aligned-length {params.min_read_aligned_length} \
            --contig-end-exclusion {params.contig_end_exclusion} \
            --min-read-percent-identity {params.min_read_percent_identity} \
            > {output.tsv}
        """


rule aggregate_all_counts:
    """Aggregate all counts and compute composition of reads."""
    input:
        counts_and_coverage=lambda wc: [
            f"results/alignment_counts_and_coverage/{accession}.tsv"
            for accession in accessions(wc)
        ],
        mito_ref_info=rules.mitochondrial_genomes_to_retain.output.csv,
        metadata=rules.process_metadata.output.metadata,
        read_counts=rules.agg_read_counts.output.csv,
    output:
        **aggregated_counts_csvs,
    params:
        sars2_ref_id=config["sars2_ref_id"],
        mito_genomes_to_keep=config["mitochondrial_genomes_to_keep"],
        mito_composition_filters=config["mito_composition_filters"],
        metagenomic_descriptions=config["metagenomic_descriptions"],
    log:
        notebook="results/aggregated_counts/aggregate_all_counts.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/aggregate_all_counts.py.ipynb"


rule build_contigs:
    """Build contigs using ``Trinity``."""
    input:
        fastqs=lambda wc: (
            [f"results/fastqs_preprocessed/{wc.accession}.fq.gz"]
            if len(accession_fastqs(wc)) == 1
            else [
                f"results/fastqs_preprocessed/{wc.accession}_R1.fq.gz",
                f"results/fastqs_preprocessed/{wc.accession}_R2.fq.gz",
            ]
        ),
    output:
        outdir=directory("results/contigs/trinity_{accession}"),
        contigs=protected("results/contigs/trinity_{accession}/Trinity.fasta"),
    params: 
        readin=lambda wc, input: (
            f"--single {input.fastqs}"
            if len(input) == 1
            else f"--left {input.fastqs[0]} --right {input.fastq[1]}"
        ),
    threads: 12
    envmodules:
        "Trinity/2.12.0-foss-2020b"
    shell:
        """
        Trinity \
            --CPU {threads} \
            --seqType fq \
            {params.readin} \
            --output {output.outdir} \
            --max_memory 50G
        """


rule get_contig_species_ref:
    """Get reference for a particular species for contig alignment, update header."""
    output:
        fasta="results/contigs/refs/{accession}_{species}.fa",
    params:
        url=lambda wc: config["metagenomic_contigs"][wc.accession]["species"][wc.species],
    conda:
        "environment.yml",
    shell:
        # since we will make concatenated genome, add species to beginning of chromosome name
        """
        curl -s {params.url} \
        | gzip -cd \
        | sed -E 's/>/>{wildcards.species}-/g' \
        > {output.fasta}
        """


rule build_contig_ref:
    """Build reference genome for each contig alignment."""
    input:
        fastas=lambda wc: [
            f"results/contigs/refs/{{accession}}_{species}.fa"
            for species in config["metagenomic_contigs"][wc.accession]["species"]
        ],
    output:
        fasta="results/contigs/refs/{accession}.fa",
        mmi="results/contigs/refs/{accession}.mmi",
    conda:
        "environment.yml",
    shell:
        """
        cat {input.fastas} > {output.fasta}
        minimap2 -I 40G -x splice -d {output.mmi} {output.fasta}
        """


rule align_contigs:
    """Align contigs to specified reference genomes."""
    input:
        contigs=rules.build_contigs.output.contigs,
        ref=rules.build_contig_ref.output.mmi,
    output:
        sam=temp("results/contigs/alignments/{accession}.sam"),
        unsorted_bam=temp("results/contigs/alignments/{accession}.bam"),
        bam=protected("results/contigs/alignments/{accession}_sorted.bam"),
    threads: 4
    params:
        extra_threads=lambda _, threads: threads - 1,
    conda:
        "environment.yml",
    shell:
        """
        minimap2 \
            -I 40G \
            -t {threads} \
            -x splice \
            --secondary=no \
            --sam-hit-only \
            -a {input.ref} \
            {input.contigs} \
            > {output.sam}
        samtools view \
            -@ {params.extra_threads} \
            -b \
            -o {output.unsorted_bam} \
            {output.sam}
        samtools sort -@ {params.extra_threads} -o {output.bam} {output.unsorted_bam}
        """


rule mapq_filter_contigs_bam:
    """Filter aligned contig BAM for only reads above a certain mapping quality."""
    input:
        bam=rules.align_contigs.output.bam,
    output:
        bam="results/contigs/alignments/{accession}_mapq_filtered_sorted.bam",
    params:
        min_mapq=config["metagenomic_contigs_min_mapq"],
    conda:
        "environment.yml"
    shell:
        "samtools view -q {params.min_mapq} -b -o {output.bam} {input.bam}"


rule contig_counts_and_coverage:
    """Get contig alignment counts to each reference."""
    input:
        bam=rules.mapq_filter_contigs_bam.output.bam,
    output:
        tsv="results/contigs/counts_and_coverage/{accession}.tsv",
    params:
        **config["metagenomic_contigs_coverm_flags"],
    conda:
        "environment.yml"
    shell:
        """
        coverm contig \
            -b {input.bam} \
            -m count covered_bases \
            --min-read-aligned-length {params.min_read_aligned_length} \
            --contig-end-exclusion {params.contig_end_exclusion} \
            --min-read-percent-identity {params.min_read_percent_identity} \
            > {output.tsv}
        """


rule process_contig_counts_and_coverage:
    """Process contig counts and coverage at species level."""
    input:
        tsvs=expand(
            rules.contig_counts_and_coverage.output.tsv,
            accession=config["metagenomic_contigs"],
        ),
    output:
        csv="results/contigs/counts_and_coverage/processed_counts.csv",
    params:
        d=lambda wc: config["metagenomic_contigs"],
    conda:
        "environment.yml"
    script:
        "scripts/process_contig_counts_and_coverage.py"


rule rt_qpcr:
    """Compare deep sequencing to China CDC table of RT-qPCR for metagenomic samples."""
    input:
        positive_table="data/positive_table.csv",
        sars2_aligned_by_sample=aggregated_counts_csvs["sars2_aligned_by_sample"],
    output:
        csv="results/rt_qpcr/rt_qpcr.csv",
        html="results/rt_qpcr/ct_vs_reads.html",
    params:
        metagenomic_descriptions=config["metagenomic_descriptions"],
    log:
        notebook="results/rt_qpcr/rt_qpcr.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/rt_qpcr.py.ipynb"


rule make_plots:
    """Make final plots."""
    input:
        **aggregated_counts_csvs,
        crits_christoph_read_counts="results/crits_christoph_data/read_counts.csv",
        ngdc_to_crits_christoph=rules.check_sha512_vs_crits_christoph.output.csv,
        contig_counts=rules.process_contig_counts_and_coverage.output.csv,
    params:
        metagenomic_descriptions=config["metagenomic_descriptions"],
        crits_christoph_plotted_species=config["crits_christoph_plotted_species"],
        susceptible_table=config["susceptible_table"],
    output:
        **plot_htmls,
        susceptible_csv="results/plots/susceptible_table.csv",
        susceptible_tex="results/plots/susceptible_table.tex",
        susceptible_mammal_csv="results/plots/susceptible_mammal_table.csv",
        susceptible_mammal_tex="results/plots/susceptible_mammal_table.tex",
        raccoon_dog_long="results/plots/raccoon_dog_long.csv",
    log:
        notebook="results/plots/make_plots.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/make_plots.py.ipynb"


rule format_plot_for_docs:
    """Format a specific plot for the GitHub pages docs."""
    input:
        plot=lambda wc: (
            plot_htmls[wc.plot]
            if wc.plot != "ct_vs_reads"
            else rules.rt_qpcr.output.html
        ),
        script="scripts/format_altair_html.py",
    output:
        plot="docs/{plot}.html",
        markdown=temp("docs/{plot}.md"),
        google_analytics_tag=temp("docs/google_analytics_tag_{plot}.txt"),
    params:
        annotations=lambda wc: docs_plot_annotations["plots"][wc.plot],
        url="https://jbloom.github.io/Huanan_market_samples",
        legend_suffix=docs_plot_annotations["legend_suffix"],
        google_analytics_tag=docs_plot_annotations["google_analytics_tag"].replace('"', '\"'),
    conda:
        "environment.yml"
    shell:
        """
        echo "## {params.annotations[title]}\n" > {output.markdown}
        echo "{params.annotations[legend]}\n\n" >> {output.markdown}
        echo "{params.legend_suffix}" >> {output.markdown}
        echo "{params.google_analytics_tag}" > {output.google_analytics_tag}
        python {input.script} \
            --chart {input.plot} \
            --markdown {output.markdown} \
            --site {params.url} \
            --title "{params.annotations[title]}" \
            --description "{params.annotations[title]}" \
            --google_analytics_tag {output.google_analytics_tag} \
            --output {output.plot}
        """


rule docs_index:
    """Write index for GitHub pages docs."""
    output:
        html="docs/index.html",
    params:
        plot_annotations=docs_plot_annotations,
    conda:
        "environment.yml"
    script:
        "scripts/docs_index.py"
