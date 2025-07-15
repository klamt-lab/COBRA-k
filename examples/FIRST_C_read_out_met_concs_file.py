"""This script contains the transformation of raw data from the supplementary data of...
'Bennett, Bryson D., et al. "Absolute metabolite concentrations and implied enzyme active site occupancy in Escherichia coli."
Nature chemical biology 5.8 (2009): 593-599.'
...into a easily machine-readable JSON, i.e., a list of BiGG (iCH360 model) metabolite IDs and the measured minimal
and maximal concentrations (as given in the mentioned publication.)
"""

import cobra

# IMPORT SECTION #
import z_add_path  # noqa: F401

from cobrak.io import json_write

# ACTUAL LOGIC SECTION #

# Load (and clean) raw Bennett data (taken from the supplement of the mentioned Bennett et al., 2009)
with open(
    "examples/common_needed_external_resources/Bennett_2009_full_raw_data.txt",
    encoding="utf-8",
) as f:
    lines = f.readlines()
lines = [x.replace("\n", "") for x in lines]
names = [x.split(";")[0].split(" (")[0].split("  ")[:-1][0].lower() for x in lines]
print(names)

# Get BiGG IDs
model = cobra.io.read_sbml_model(
    "examples/iCH360/external_resources/EC_iCH360_fitted_kapps.xml"
)
name_to_bigg_ids = {}
for name in names:
    for metabolite in model.metabolites:
        lower_metname = metabolite.name.lower()
        if (
            (name == lower_metname)
            or (f"l-{name}" == lower_metname)
            or (f"(s)-{name}" == lower_metname)
            or (f"d-{name}" == lower_metname)
            or (name.split(" ")[0] == lower_metname.split(" ")[0])
        ):
            name_to_bigg_ids[name] = metabolite.id
            break
    if name not in name_to_bigg_ids:
        name_to_bigg_ids[name] = "N/A"

# Manual setting for found other metabolites where the automatic BiGG ID identification fails
name_to_bigg_ids["2,3-dihydroxybenzoic acid"] = "23dhb_c"  # KEGG ID C00196
name_to_bigg_ids["3-phosphoglycerate"] = "3pg_c"  # KEGG ID C00197
name_to_bigg_ids["6-phosphogluconate"] = "6pgc_c"  # KEGG ID C00345
name_to_bigg_ids["acetylphosphate"] = "actp_c"  # KEGG ID C00227
name_to_bigg_ids["adp-glucose"] = "adpglc_c"  # KEGG ID C00498
name_to_bigg_ids["alpha-ketoglutarate"] = "akg_c"  # KEGG ID C00026
name_to_bigg_ids["carbamylaspartate"] = "cbasp_c"  # KEGG ID C00438
name_to_bigg_ids["deoxyribose-5-p"] = "2dr5p_c"  # KEGG ID C00673
name_to_bigg_ids["fad"] = "fad_c"  # KEGG ID C00016
name_to_bigg_ids["fructose-1,6-bisphosphate"] = "fdp_c"  # MetaNetX MNXM417
name_to_bigg_ids["glucosamine-6-phosphate"] = "gam6p_c"  # MetaNetX MNXM370
name_to_bigg_ids["glycerate"] = "glyc__R"
name_to_bigg_ids["glycerol-3-phosphate"] = "glyc3p_c"
name_to_bigg_ids["n-acetyl-glucosamine-1p"] = "acgam1p_c"
name_to_bigg_ids["n-acetyl-ornithine"] = "acorn_c"
name_to_bigg_ids["propionyl-coa"] = "ppcoa_c"
name_to_bigg_ids["prpp"] = "prpp_c"
name_to_bigg_ids["udp-glucuronate"] = "udpglcur_c"
name_to_bigg_ids["udp-glucose"] = "udpg_c"
name_to_bigg_ids["nad+"] = "nad_c"
name_to_bigg_ids["nadh"] = "nadh_c"
name_to_bigg_ids["nadp+"] = "nadp_c"
name_to_bigg_ids["nadph"] = "nadph_c"
# gluconolactone does not seem to exist in iCH260
# hexose-pa is a general term
# pentose-pd is a general term

for name in name_to_bigg_ids:
    if name_to_bigg_ids[name] == "N/A":
        print("MISSING:", name)

with open(
    "examples/common_needed_external_resources/Bennett_2009_full_raw_data.txt",
    encoding="utf-8",
) as f:
    lines = f.readlines()
lines = [x.replace("\n", "") for x in lines]

lbs = []
ubs = []
full_data_dict = {}
for line in lines:
    linesplit = line.split("  ")
    name = linesplit[0]
    bigg_id = name_to_bigg_ids[name.lower()]

    bigg_id = (bigg_id + "\b").replace("_e\b", "_c")
    bigg_id = (bigg_id + "\b").replace("_p\b", "_c")
    bigg_id = bigg_id.replace("\b", "")

    concentration_raw = linesplit[1]
    mean_conc = float(concentration_raw.split(" ")[0])
    lb_conc = float(concentration_raw.split(" ")[1].replace("(", ""))
    ub_conc = float(concentration_raw.split(" ")[3].replace(")", ""))
    full_data_dict[bigg_id] = {
        "lb": lb_conc,
        "ub": ub_conc,
        "mean": mean_conc,
        "name": name,
    }
    lbs.append(lb_conc)
    ubs.append(ub_conc)
    if ub_conc > 0.02:
        print("OVER", full_data_dict[bigg_id])
    if lb_conc < 1e-6:
        print("UNDER", full_data_dict[bigg_id])

print("MIN", min(lbs))
print("MAX", max(ubs))

# Glutathione and glutathione disulfide special treatment
# Gluthathione disulfide seems to be mixed together with oxidized/reduced glutathione in iML1515
# Therefore, only an upper bound for all affected metabolites is set at the sum of the measured
# upper bound of gluthathione and glutathione sulfate concentrations.
# Since the distribution of oxidized and reduced glutathione remains unclear, the default lower
# bound is set.
ub_glutathione = 1.79e-2
ub_glutathione_disulfide = 2.90e-3
sum_ubs_glutathiones = ub_glutathione + ub_glutathione_disulfide
# -> For oxidized glutathione (gthox_c)
full_data_dict["gthox_c"] = {
    "lb": 1e-6,
    "ub": sum_ubs_glutathiones,
    "mean": (1e-6 + sum_ubs_glutathiones) / 2,
    "name": name,
}
# -> For reduced glutathione (gthrd_c)
full_data_dict["gthrd_c"] = {
    "lb": 1e-6,
    "ub": sum_ubs_glutathiones,
    "mean": (1e-6 + sum_ubs_glutathiones) / 2,
    "name": name,
}

json_write(
    "./examples/common_needed_external_resources/Bennett_2009_full_data.json",
    full_data_dict,
)
