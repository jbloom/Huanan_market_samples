# Association between SARS-CoV-2 and metagenomic content of samples from the Huanan Seafood Market

This repository contains a fully reproducible computational pipeline for analyzing the SARS-CoV-2 and chordate mitochondrial metagenomic content of the deep sequencing of environmental samples taken from the Huanan Seafood market by [Liu et al (2023)](https://www.nature.com/articles/s41586-023-06043-2).

The pipeline was created by Jesse Bloom.

The analysis and results are described in **XXX**.

## Results and plots
Key results are in the [./results/](results) subdirectory.
These include:

 - [results/merged_metadata.csv](results/merged_metadata.csv): metadata about the samples extracted by processing files provided on the NGDC by [Liu et al (2023)](https://www.nature.com/articles/s41586-023-06043-2).
 - [results/crits_christoph_data/check_sha512_vs_crits_christoph.csv](results/crits_christoph_data/check_sha512_vs_crits_christoph.csv): comparison of SHA-512 hashes for the FASTQ files downloaded from the NGDC to those reported in the earlier analysis by [Crits-Christoph et al (2023)](https://zenodo.org/record/7754299#.ZEghB-zMKX0).
 - [results/mitochondrial_genomes/retained.csv](results/mitochondrial_genomes/retained.csv): the set of chordate mitochondrial genomes to which reads were aligned for the metagenomic analysis.
 - [results/aggregated_counts/sars2_aligned_by_run.csv](results/aggregated_counts/sars2_aligned_by_run.csv): number of aligned SARS-CoV-2 reads for each sequencing run.
 - [results/aggregated_counts/sars2_aligned_by_sample.csv](results/aggregated_counts/sars2_aligned_by_sample.csv): number of aligned SARS-CoV-2 reads for each sample.
 - [results/aggregated_counts/mito_composition_by_run.csv](results/aggregated_counts/mito_composition_by_run.csv): chordate mitochondrial composition for each sequencing run.
 - [results/aggregated_counts/mito_composition_by_sample.csv](results/aggregated_counts/mito_composition_by_sample.csv): chordate mitochondrial composition for each sample.
 - [results/rt_qpcr/rt_qpcr.csv](results/rt_qpcr/rt_qpcr.csv): SARS-CoV-2 content of samples determined in current sequencing and Ct values from RT-qPCR reported by [Liu et al (2023)](https://www.nature.com/articles/s41586-023-06043-2).
 - [results/plots/susceptible_table.csv](results/plots/susceptible_table.csv): SARS-CoV-2 content of samples with high mitochondrial composition from susceptible species sold live at the market.
 - [results/contigs/counts_and_coverage/processed_counts.csv](results/contigs/counts_and_coverage/processed_counts.csv): results for aligning assembled contigs to full genomes for selected samples and genomes.

Note that the pipeline also produces many other results files (some of which are very large) that are not tracked in this repo.

Interactive plots of the results created using [Altair](https://altair-viz.github.io/) are rendered via GitHub Pages at [https://jbloom.github.io/Huanan_market_samples/](https://jbloom.github.io/Huanan_market_samples/)

## Understanding and running the pipeline
The entire analysis can be run in automated fashion using [snakemake](https://snakemake.readthedocs.io/).

The pipeline itself is in [Snakefile](Snakefile).
The configuration for the pipeline is specified in [config.yaml](config.yaml).
The pipeline uses the [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml), which specifies the precise versions of all software used.
The one exception is that the rule `build_contigs` in [Snakefile](Snakefile) uses an [environment module](https://modules.readthedocs.io/en/latest/) that is pre-built on the Fred Hutch computing cluster to run the [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to build contigs---to run this rule, you will need to specify a comparable module for whatever computing system you are using, or skip the contig building by commenting out the file `results/contigs/counts_and_coverage/processed_counts.csv` as an input to the `all` rule in [Snakefile](Snakefile).
The scripts and Jupyter notebooks used by the pipeline are in [./scripts/](scripts) and [./notebooks/](notebooks), respectively.

Most data used by the pipeline is downloaded by the pipeline, but it takes the following to input files, both found in [./data/](data):

  - [data/CRA010170.xlsx](data/CRA010170.xlsx) is the GSA BioProject metadata sheet downloaded from the NGDC GSA page [https://ngdc.cncb.ac.cn/gsa/browse/CRA010170](https://ngdc.cncb.ac.cn/gsa/browse/CRA010170) on March-29-2023.

  - [data/positive_table.csv](data/positive_table.csv) is a version of Supplementary Table 2 from [Liu et al (2023)](https://www.nature.com/articles/s41586-023-06043-2), taken from [this link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06043-2/MediaObjects/41586_2023_6043_MOESM4_ESM.docx) (archived [here](https://web.archive.org/web/20230405155400/https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06043-2/MediaObjects/41586_2023_6043_MOESM4_ESM.docx)).

To run the pipeline on the Fred Hutch computing cluster, use the commands in [run_Hutch_cluster.bash](run_Hutch_cluster.bash).

Below is a rulegraph of the pipeline built with:

    snakemake --forceall --rulegraph | dot -Tpng > rulegraph.png

![Rulegraph of `snakemake` pipeline](results/rulegraph.png)
