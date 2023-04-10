"""Build index for GitHub pages docs."""


import markdown


html_file = snakemake.output.html
plot_annotations = snakemake.params.plot_annotations

text = [
    f"## {plot_annotations['index_title']}",
    plot_annotations["index_abstract"],
]

for plot, plot_d in plot_annotations["plots"].items():
    text.append(
        f"  - [{plot_d['title']}]({plot}.html)"
    )

text.append("\n" + plot_annotations["legend_suffix"])

html = markdown.markdown("\n".join(text))

with open(html_file, "w") as f:
    f.write(html)

