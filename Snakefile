"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


configfile: "config.yaml"


checkpoint process_metadata:
    """Process metadata sample sheet."""
    input:
        excel="data/CRA010170.xlsx",
    output:
        samples="results/metadata/samples.csv",
        experiments="results/metadata/experiments.csv",
        runs="results/metadata/runs.csv",
    conda:
        "environment.yml"
    script:
        "scripts/process_metadata.py"

