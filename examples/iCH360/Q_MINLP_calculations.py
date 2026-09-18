from copy import deepcopy

import z_add_path  # noqa: D100, F401

from cobrak.dataclasses import Model, Solver, VarResult
from cobrak.io import json_load
from cobrak.nlps import perform_nlp_reversible_optimization
from cobrak.standard_solvers import BARON, IPOPT_MA57, SCIP  # noqa: F401

cobrak_model: Model = json_load(
    "examples/iCH360/RESULTS_GLCUPTAKE/used_cobrak_model__1_maxglc9.65.json",
    # "examples/iCH360/RESULTS_MAXAC/used_cobrak_model__1__maxac.json",
    Model,
)
variability_dict: VarResult = json_load(
    "examples/iCH360/RESULTS_GLCUPTAKE/variability_dict__1_maxglc9.65.json"
    # "examples/iCH360/RESULTS_MAXAC/variability_dict__1__maxac.json"
)
for varkey in deepcopy(list(variability_dict.keys())):
    if varkey in cobrak_model.reactions:
        variability_dict[varkey] = (
            0.0 if variability_dict[varkey][0] < 1e-3 else variability_dict[varkey][0],
            max(variability_dict[varkey][1], 0.001),
        )
    # if varkey.startswith(ENZYME_VAR_PREFIX):
    #     variability_dict[varkey] = (
    #         0.0 if variability_dict[varkey][0] < 1e-3 else variability_dict[varkey][0],
    #         1e-3 if variability_dict[varkey][1] < 1e-3 else variability_dict[varkey][1],
    #     )
# variability_dict["Biomass_fw"] = (0.3, 0.7)
nlp_result = perform_nlp_reversible_optimization(
    approximation_value=0.01,
    cobrak_model=cobrak_model,
    objective_target="Biomass_fw",  # 0.6879331919321442 is th efound optimum
    objective_sense=+1,
    variability_dict=variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_iota=False,
    with_alpha=False,
    strict_mode=False,
    verbose=True,
    solver=Solver(
        name="pounce",
        solver_options={
            "nlp_scaling_method": "equilibration-based",
            "hsllib": IPOPT_MA57.solver_options["hsllib"],
            "linear_solver": "ma97",
            "max_iter": 400_000,
        },
    ),  # , solver_options={"mu_strategy": "adaptive", "acceptable_tol": 1e-4, "acceptable_constr_viol_tol": 0.01, "corrector_type": "primal-dual",}),  # Change with BARON to test this solver
    with_flux_sum_var=False,
    show_variable_count=True,
)
