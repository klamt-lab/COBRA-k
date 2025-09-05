# IMPORTS SECTION #  # noqa: D100
from math import log

import z_add_path  # noqa: F401

from cobrak.dataclasses import EnzymeReactionData
from cobrak.io import json_load, json_write
from cobrak.model_instantiation import get_cobrak_model_from_sbml_and_thermokinetic_data
from cobrak.utilities import delete_orphaned_metabolites_and_enzymes
from examples.iCH360.A_general_model_constants import (
    conc_ranges,
    kinetic_ignored_metabolites,
)

output_folder = "examples/iCH360/prepared_external_resources/"

dG0s = json_load(f"{output_folder}dG0s.json", dict[str, float])
enzyme_reaction_data: dict[str, EnzymeReactionData | None] = json_load(
    f"{output_folder}enzyme_reaction_data_COMBINED.json",
    dict[str, EnzymeReactionData | None],
)

mws = json_load(f"{output_folder}mws.json", dict[str, float])
cobrak_model = get_cobrak_model_from_sbml_and_thermokinetic_data(
    sbml_path=f"{output_folder}EC_iCH360_unadjusted_kcats_cleaned.xml",
    extra_linear_constraints=[],
    dG0s=dG0s,
    dG0_uncertainties={},
    conc_ranges=conc_ranges,
    enzyme_molecular_weights=mws,
    enzyme_reaction_data=enzyme_reaction_data,
    max_prot_pool=0.39,
    kinetic_ignored_metabolites=kinetic_ignored_metabolites,
    fwd_suffix="_fw",
    rev_suffix="_bw",
)

# Delete NAD-dependent fatty acid EAR reactions (if not essential)
reac_ids = list(cobrak_model.reactions.keys())
for reac_id in reac_ids:
    if not reac_id.startswith("EAR"):
        continue
    if reac_id.endswith(("y_fw", "y_bw")):
        continue
    if reac_id.replace("x_", "y_") in reac_ids:
        print("DELETED:", reac_id)
        del cobrak_model.reactions[reac_id]
    else:
        print("KEPT:", reac_id)

cobrak_model = delete_orphaned_metabolites_and_enzymes(cobrak_model)

# Widen concentration range of metabolites
for met_id, metabolite in cobrak_model.metabolites.items():
    if met_id.startswith(("h_", "h2o_")):
        continue
    if not met_id.endswith("_e"):
        continue
    metabolite.log_max_conc = log(0.2)

cobrak_model.max_conc_sum = 0.4
cobrak_model.conc_sum_ignore_prefixes = ["h2o_", "h_"]
cobrak_model.conc_sum_include_suffixes = ["_c"]
cobrak_model.conc_sum_max_rel_error = 0.1
cobrak_model.conc_sum_min_abs_error = 1e-6

json_write(
    f"{output_folder}iCH360_cobrak_prepstepA_bottlenecked_unfilled_uncalibrated.json",
    cobrak_model,
)
