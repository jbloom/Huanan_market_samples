"""Process metadata from Excel."""


import pandas as pd


metadata = pd.read_excel(snakemake.input.excel, sheet_name=None)

for sheetname, sheet in metadata.items():
    output_name = sheetname.lower() + "s"  # name we give output file in snakemake
    output_csv = snakemake.output[output_name]
    print(f"Writing {sheetname=} to {output_csv=}")
    sheet.to_csv(output_csv, index=False)
