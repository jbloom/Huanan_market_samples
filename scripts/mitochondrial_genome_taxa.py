"""Get taxa for a mitochondrial genome."""


import re
import xml

import Bio.Entrez
import Bio.SeqIO

import pandas as pd


seq = Bio.SeqIO.read(snakemake.input.gb, format="genbank")

source_feature = [f for f in seq.features if f.type == "source"]
source_feature = source_feature[0]
taxon_qualifier = [
    q for q in source_feature.qualifiers["db_xref"] if q.startswith("taxon:")
]
assert len(taxon_qualifier) == 1
taxon_key, taxonid = taxon_qualifier[0].split(":")
assert taxon_key == "taxon", source_feature.qualifiers
print(f"Getting taxa for {taxonid=}")
tax_xml = Bio.Entrez.efetch(db="taxonomy", id=taxonid, retmode="xml").read().decode("utf-8")

taxa_groups = ["phylum", "subphylum"]
taxa = []
for taxa_group in taxa_groups:
    m = re.search(
        "\s+<TaxId>\d+</TaxId>\n"
        + "\s+<ScientificName>(?P<taxa>\w+)</ScientificName>\n"
        + f"\s+<Rank>{taxa_group}</Rank>\n",
        tax_xml
    )
    if not m:
        taxa.append("")
        assert taxa_group not in tax_xml.lower(), f"{taxa_group=}\n{tax_xml}"
    else:
        taxa.append(m.group("taxa"))

with open(snakemake.output.csv, "w") as f:
    f.write("id,taxonid," + ",".join(taxa_groups) + "\n")
    f.write(f"{seq.id},{taxonid}," + ",".join(taxa) + "\n")
