# IMPORTS SECTION #

import z_add_path  # noqa: F401

from cobrak.constants import (
    ALL_OK_KEY,
    OBJECTIVE_VAR_NAME,
)
from cobrak.dataclasses import Model
from cobrak.io import json_load, json_write
from cobrak.lps import (
    perform_lp_optimization,
    perform_lp_thermodynamic_bottleneck_analysis,
)
from cobrak.printing import print_optimization_result
from cobrak.utilities import (
    get_reverse_reac_id_if_existing,
)

# RUNNING SCRIPT SECTION #
working_folder = "examples/iCH360/prepared_external_resources/"

biomass_reac_id = "Biomass_fw"
cobrak_model: Model = json_load(
    f"{working_folder}iCH360_cobrak_prepstepA_bottlenecked_unfilled_uncalibrated.json",
    Model,
)
original_conc_sum = cobrak_model.max_conc_sum
cobrak_model.max_conc_sum = float("inf")

with cobrak_model as original_prot_pool_model:
    print("Running FBA with maximal growth rate and original protein pool as test...")
    original_prot_pool_model.max_prot_pool = 0.315157359265934
    fba_result_original = perform_lp_optimization(
        original_prot_pool_model,
        objective_target=biomass_reac_id,
        objective_sense=+1,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=False,
        with_loop_constraints=False,
        variability_dict={},
        ignored_reacs=[],
        verbose=False,
    )
    assert fba_result_original[ALL_OK_KEY]
    print(round(fba_result_original[biomass_reac_id], 3))
    assert round(fba_result_original[biomass_reac_id], 3) == 1.733
    assert round(fba_result_original[biomass_reac_id], 5) == round(
        fba_result_original[OBJECTIVE_VAR_NAME], 5
    )

print()
print("Running FBA for maximal growth rate...")
fba_result = perform_lp_optimization(
    cobrak_model,
    objective_target=biomass_reac_id,
    objective_sense=+1,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=False,
    with_loop_constraints=False,
    variability_dict={},
    ignored_reacs=[],
    verbose=False,
)
print("FBA result with maximal growth rate:")
print_optimization_result(cobrak_model, fba_result, ignore_unused=True)
assert fba_result[ALL_OK_KEY]

print()
print("Running thermodynamic bottleneck analysis...")
with cobrak_model as enforced_biomass_model:
    enforced_biomass_model.reactions[biomass_reac_id].min_flux = (
        fba_result[biomass_reac_id] - 1e-6
    )
    bottleneck_reactions, _ = perform_lp_thermodynamic_bottleneck_analysis(
        cobrak_model=enforced_biomass_model,
        with_enzyme_constraints=True,
        verbose=False,
    )
print(
    f"Identified thermodynamic bottlenecks @ growth rate {fba_result[biomass_reac_id] - 1e-6} 1/h:",
    bottleneck_reactions,
)
assert len(bottleneck_reactions) == 2
assert set(bottleneck_reactions) == {"AIRC3_bw", "Htex_bw"}

print(
    "Add H2O transport reactions (as H2O has a fixed concentration) and biomass and exchange reactions (as they are pseudo-reactions)"
)
bottleneck_reactions += ["H2Otex_fw", "H2Otpp_fw", "Biomass_fw"] + [
    reac_id for reac_id in cobrak_model.reactions if reac_id.startswith("EX_")
]
bottleneck_reactions_plus_reverse = []
for bottleneck_reaction in bottleneck_reactions:
    reverse_id = get_reverse_reac_id_if_existing(
        bottleneck_reaction, cobrak_model.fwd_suffix, cobrak_model.rev_suffix
    )
    for reac_id in [bottleneck_reaction, reverse_id]:
        if reac_id not in cobrak_model.reactions:
            continue
        bottleneck_reactions_plus_reverse.append(reac_id)
        cobrak_model.reactions[reac_id].dG0 = None
        cobrak_model.reactions[reac_id].enzyme_reaction_data = None
        print(f"->Removing ΔG'° and enzyme data of {reac_id}")

print("Remove enzyme data of remaining unenzymatic diffusion reactions...")
diffusion_reactions = [
    reac_id
    for reac_id, reaction in cobrak_model.reactions.items()
    if (
        "diffusion" in reaction.name
        and reac_id not in bottleneck_reactions_plus_reverse
    )
    or (reac_id.startswith("ATPM_"))
]
for diffusion_reaction in diffusion_reactions:
    cobrak_model.reactions[diffusion_reaction].enzyme_reaction_data = None
    print(f"->Remove enzyme data of diffusion reaction {diffusion_reaction}")

cobrak_model.max_conc_sum = original_conc_sum
json_write(
    f"{working_folder}iCH360_cobrak_prestepB_unfilled_uncalibrated.json",
    cobrak_model,
)
json_write(
    f"{working_folder}bottleneck_and_pseudo_and_diffusion_reactions.json",
    diffusion_reactions + bottleneck_reactions_plus_reverse,
)
