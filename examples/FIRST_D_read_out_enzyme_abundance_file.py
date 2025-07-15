import z_add_path  # noqa: F401

from cobrak.io import json_write

with open(
    "examples/common_needed_external_resources/schmidt_2016_mass_abundance_mapped_NNLS.csv",
    encoding="utf-8",
) as f:
    csvlines = f.read().split("\n")

mg_column = csvlines[0].split(",").index("glucose_MG1655_g_gDW")
enzid_column = -1

enzjson: dict[str, float] = {}
for csvline in csvlines[1:]:
    if len(csvline) == 0:
        continue
    enzvalue = csvline.split(",")[mg_column]
    enzid = csvline.split(",")[enzid_column]
    if enzvalue == "":
        continue
    enzjson[enzid] = float(enzvalue)
json_write(
    "examples/common_needed_external_resources/Schmidt_2016_full_data.json", enzjson
)
