"""Tally counts of alignments to each target.

We tally counts of:

 - primary alignments that are also highest DP score of all alignments.

 - primary alignments regardless of whether there are equivalent DP score secondaries.

 - all alignments with DP scores >= score of primary alignment.


"""


import collections

import pandas as pd

import pysam


n_unique_best = n_multi_best = 0

reads = collections.defaultdict(list)
    
hits_by_ref = collections.defaultdict(
    lambda: {
        "n_primary": 0,
        "n_primary_and_unique_best": 0,
        "n_primary_or_equally_good_secondary": 0,
    }
)

with pysam.AlignmentFile(snakemake.input.bamfile, "rb") as bam:
    for read in bam:
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue  # ignore supplementary alignments
        reads[read.query_name].append(
            (
                not read.is_secondary,
                read.get_tag("AS"),  # DP alignment score
                read.reference_id,
            )
        )

    for read_name, read_list in reads.items():
        read_list.sort(reverse=True)
        assert read_list[0][0] is True
        hits_by_ref[read_list[0][2]]["n_primary"] += 1
        if (len(read_list) == 1) or (read_list[0][1] > read_list[1][1]):
            hits_by_ref[read_list[0][2]]["n_primary_and_unique_best"] += 1
        i = 0
        while (i < len(read_list)) and (read_list[i][1] >= read_list[0][1]):
            hits_by_ref[read_list[i][2]]["n_primary_or_equally_good_secondary"] += 1
            i += 1

    # switch from reference IDs to reference names
    hits_by_ref = {
        bam.get_reference_name(ref_id): ref_dict
        for ref_id, ref_dict in hits_by_ref.items()
    }

df = (
    pd.DataFrame.from_dict(hits_by_ref, orient="index")
    .assign(accession=snakemake.wildcards.accession)
    .rename(columns={"accession": "Run accession"})
    .rename_axis("alignment_reference")
    .reset_index()
)

df.to_csv(snakemake.output.csv, index=False)
