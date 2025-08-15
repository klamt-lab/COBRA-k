# IMPORTS SECTION
import time

import z_add_path  # noqa: F401

from cobrak.dataclasses import Model
from cobrak.io import json_load, json_write
from cobrak.lps import (
    perform_lp_variability_analysis,
)
from cobrak.standard_solvers import CPLEX_FOR_VARIABILITY_ANALYSIS
from cobrak.utilities import get_model_with_filled_missing_parameters

# RUNNING SCRIPT SECTION #
biomass_reac_id = "Biomass_fw"
cobrak_model: Model = json_load(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak_prestepB_unfilled_uncalibrated.json",
    Model,
)

print("Create filled model and run variability analysis with it...")
filled_cobrak_model = get_model_with_filled_missing_parameters(
    cobrak_model,
    ignore_prefixes=json_load(
        "examples/iCH360/prepared_external_resources/bottleneck_and_pseudo_and_diffusion_reactions.json",
        list[str],
    ),
    add_dG0_extra_constraints=True,
    use_median_for_kms=True,
    use_median_for_kcats=True,
)

json_write(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak_prestepC_uncalibrated.json",
    filled_cobrak_model,
)
with filled_cobrak_model as unconcrestricted_filled_cobrak_model:
    t3 = time.time()
    variability_dict = perform_lp_variability_analysis(
        unconcrestricted_filled_cobrak_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        min_flux_cutoff=1e-6,
        solver=CPLEX_FOR_VARIABILITY_ANALYSIS,
        active_reactions=["Biomass_fw"],
        min_active_flux=1e-5,
    )
    t4 = time.time()
    json_write(
        "examples/iCH360/prepared_external_resources/variability_dict_enforced_growth_uncalibrated.json",
        variability_dict,
    )
    print("...done!")
    print("TIME [s]:", t4 - t3)

with filled_cobrak_model as unconcrestricted_filled_cobrak_model:
    t3 = time.time()
    variability_dict = perform_lp_variability_analysis(
        unconcrestricted_filled_cobrak_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        min_flux_cutoff=1e-6,
        solver=CPLEX_FOR_VARIABILITY_ANALYSIS,
        active_reactions=["Biomass_fw", "EX_ac_e_fw"],
        min_active_flux=1e-5,
    )
    t4 = time.time()
    json_write(
        "examples/iCH360/prepared_external_resources/variability_dict_enforced_growth_and_acetate_uncalibrated.json",
        variability_dict,
    )
    print("...done!")
    print("TIME [s]:", t4 - t3)
