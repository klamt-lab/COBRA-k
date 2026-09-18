"""Runs all analyses for the toymodel as shown in COBRA-k's initial publication"""

# from cobrak.io import json_write
import time

from cobrak.constants import OBJECTIVE_VAR_NAME

try:  # noqa: SIM105
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

from math import log

from cobrak._evolution2 import (
    EvolutionSettings,
    perform_nlp_evolutionary_optimization,
)
from cobrak.dataclasses import (
    ExtraLinearConstraint,
)
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_variability_analysis
from cobrak.nlps import (  # noqa: F401
    perform_nlp_irreversible_optimization_with_active_reacs_only,
    perform_nlp_reversible_optimization,
)
from cobrak.standard_solvers import BARON, GUROBI, IPOPT, IPOPT_MA57, SCIP  # noqa: F401

IPOPT.solver_options["acceptable_tol"] = 1e-10
IPOPT.solver_options["max_iter"] = 100_000
IPOPT.solver_options["mu_strategy"] = "adaptive"
IPOPT.solver_options["corrector_type"] = "primal-dual"

side_reac_id = "Glycolysis"
main_reac_ids = ["Respiration", "Overflow"]

toy_model.extra_linear_constraints = [
    ExtraLinearConstraint(
        stoichiometries={
            "x_ATP": 1.0,
            "x_ADP": -1.0,
        },
        lower_value=log(3.0),
    )
]

# ecTFVA #
print("[b]Run variability analysis with thermodynamic and enzyme constraints...[/b]")
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
    min_flux_cutoff=1e-7,
)

print("----------------------------------------")

variability_dict["EX_S"] = (0.0, 14.0)
t0 = time.time()

for metabolite in toy_model.metabolites.values():
    metabolite.molar_mass = 1.0
toy_model.max_met_mass_sum = 1000.0
toy_model.max_conc_sum = 0.8
toy_model.include_mets_in_prot_pool = True
toy_model.conc_sum_include_suffixes = list(toy_model.metabolites)

result = perform_nlp_evolutionary_optimization(
    EvolutionSettings(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        variability_dict=variability_dict,
        num_gens=10,
        population_size=10,
        with_kappa=True,
        with_gamma=True,
        with_alpha=False,
        with_iota=False,
        lp_solver=GUROBI,
        nlp_solver=IPOPT,
        use_original_ectfba_nlp_binaries_as_template=(True, True, False),
    )
)
t1 = time.time()
print(
    f"max(ATP_Consumption) from evolutionary algorithm under EX_S <= 14: {result.best_nlps[0][OBJECTIVE_VAR_NAME]}"
)
assert result.best_nlps[0][OBJECTIVE_VAR_NAME] > 32.718
assert result.best_nlps[0][OBJECTIVE_VAR_NAME] < 32.719
print("TIME FOR COBRA-k evolutionary algorithm:", t1 - t0)
