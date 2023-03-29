#!/bin/bash
#
#SBATCH -c 1

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

mkdir -p slurm_logs

printf "Running snakemake...\n"
# set ntasks-per-node=1 because pysradb download issues if too many downloads from an IP
snakemake \
    -j 500 \
    --cluster "sbatch -c 1 -t 4-0 -J Huanan_market_samples -o slurm_logs/slurm-%j.out -e slurm_logs/slurm-%j.out --ntasks-per-node=1" \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --retries 3 \
    --scheduler greedy \
    --rerun-incomplete
printf "Run of snakemake complete.\n"
