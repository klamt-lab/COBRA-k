import z_add_path  # noqa: D100, F401

from cobrak.dataclasses import Model, VarResult
from cobrak.io import json_load
from cobrak.nlps import perform_nlp_reversible_optimization
from cobrak.standard_solvers import BARON, SCIP  # noqa: F401

cobrak_model: Model = json_load(
    "examples/iCH360/RESULTS_GLCUPTAKE/used_cobrak_model__1_maxglc1000.json",
    # "examples/iCH360/RESULTS_MAXAC/used_cobrak_model__1__maxac.json",
    Model,
)
variability_dict: VarResult = json_load(
    "examples/iCH360/RESULTS_GLCUPTAKE/variability_dict__1_maxglc1000.json"
    # "examples/iCH360/RESULTS_MAXAC/variability_dict__1__maxac.json"
)

nlp_result = perform_nlp_reversible_optimization(
    cobrak_model=cobrak_model,
    objective_target="EX_ac_e_fw",
    objective_sense=+1,
    variability_dict=variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_iota=False,
    with_alpha=False,
    strict_mode=False,
    verbose=True,
    solver=SCIP,  # Change with BARON to test this solver
    with_flux_sum_var=False,
    show_variable_count=True,
)
