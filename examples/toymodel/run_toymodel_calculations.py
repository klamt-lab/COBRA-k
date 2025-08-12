"""Runs all analyses for the toymodel as shown in COBRA-k's initial publication"""

import os
import time
from copy import deepcopy

try:  # noqa: SIM105
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

from math import log

from cobrak.constants import (
    OBJECTIVE_VAR_NAME,
    Z_VAR_PREFIX,
)
from cobrak.dataclasses import (
    ExtraLinearConstraint,
)
from cobrak.evolution import (
    perform_nlp_evolutionary_optimization,
    perform_nlp_irreversible_optimization_with_active_reacs_only,
)
from cobrak.example_models import toy_model
from cobrak.io import (
    json_write,
    load_annotated_sbml_model_as_cobrak_model,
    save_cobrak_model_as_annotated_sbml_model,
)
from cobrak.lps import perform_lp_optimization, perform_lp_variability_analysis
from cobrak.nlps import perform_nlp_reversible_optimization  # noqa: F401
from cobrak.printing import (
    print_variability_result,
)
from cobrak.standard_solvers import BARON, IPOPT, SCIP  # noqa: F401

side_reac_id = "Glycolysis"
main_reac_ids = ["Respiration", "Overflow"]

save_cobrak_model_as_annotated_sbml_model(
    toy_model, filepath="examples/toymodel/sbml_model.xml"
)
toy_model = load_annotated_sbml_model_as_cobrak_model(
    filepath="examples/toymodel/sbml_model.xml"
)

toy_model.extra_linear_constraints = [
    ExtraLinearConstraint(
        stoichiometries={
            "x_ATP": 1.0,
            "x_ADP": -1.0,
        },
        lower_value=log(3.0),
    )
]

os.environ["PATH"] += os.pathsep + "/usr/local/net/GAMS/38.1/"

# ecTFVA #
print("[b]Run variability analysis with thermodynamic and enzyme constraints...[/b]")
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
    min_flux_cutoff=1e-7,
)
json_write(
    "examples/toymodel/variability_dict.json",
    variability_dict,
)
print("[b]...done! Variability result:[/b]")
print_variability_result(toy_model, variability_dict)
print()

# MINLP #
print("Run MINLP...")
test_variability_dict = deepcopy(variability_dict)

"""
# uncomment to run MINLP test runs
t0 = time.time()
nlp_result_rev = perform_nlp_reversible_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    variability_dict=test_variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_iota=False,
    with_alpha=False,
    verbose=True,
    solver=SCIP,
    # solver=BARON,
    show_variable_count=True,
)
t1 = time.time()
print("...done! First MINLP result time:", t1 - t0)
print()
"""

print("----------------------------------------")

for max_glc_uptake in [50, 14.0]:
    print(f"~~~MAX GLC UPTAKE: {max_glc_uptake}~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    for exclusions in [[main_reac_ids[0]], [main_reac_ids[1]], []]:
        active_reacs = " & ".join(
            [
                reac_id
                for reac_id in toy_model.reactions
                if (not reac_id.startswith("EX_"))
                and (reac_id not in exclusions)
                and (reac_id in main_reac_ids)
            ]
        )

        test_variability_dict = deepcopy(variability_dict)
        for exclusion in exclusions:
            test_variability_dict[exclusion] = (0.0, 0.0)
        test_variability_dict["EX_S"] = (0.0, max_glc_uptake)

        run_cobrak_model = deepcopy(toy_model)
        if not exclusions:
            run_cobrak_model.extra_linear_constraints = [
                ExtraLinearConstraint(
                    stoichiometries={
                        f"{Z_VAR_PREFIX}{reac_id}": 1.0,
                    },
                    lower_value=1.0,
                    upper_value=1.0,
                )
                for reac_id in [side_reac_id, *main_reac_ids]
            ]

        print(run_cobrak_model.reactions.keys())
        lp_result = perform_lp_optimization(
            cobrak_model=run_cobrak_model,
            objective_target="ATP_Consumption",
            objective_sense=+1,
            variability_dict=test_variability_dict,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            with_loop_constraints=False,
            verbose=False,
        )
        print(
            f"ecTFBA | max(ATP_Consumption) with {active_reacs}:",
            round(lp_result[OBJECTIVE_VAR_NAME], 4),
            "uptake:",
            round(lp_result["EX_S"], 2),
        )

        if not exclusions:
            for reac_id in [
                side_reac_id,
                *main_reac_ids,
                "EX_S",
                "EX_M",
                "EX_C",
                "EX_D",
                "EX_P",
                "ATP_Consumption",
            ]:
                lp_result[reac_id] = 1.0
        try:
            nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
                deepcopy(toy_model),
                objective_target="ATP_Consumption",
                objective_sense=+1,
                optimization_dict=lp_result,
                variability_dict=test_variability_dict,
                with_kappa=True,
                with_gamma=True,
                with_iota=False,
                with_alpha=False,
                verbose=False,
                solver=IPOPT,
            )
        except ValueError:
            continue

        from math import exp

        print(
            f"NLP    | max(ATP_Consumption) with {active_reacs}:",
            round(nlp_result[OBJECTIVE_VAR_NAME], 4),
            "uptake:",
            round(nlp_result["EX_S"], 2),
            [
                (enzyme_id, round(nlp_result[enzyme_id], 100))
                for enzyme_id in nlp_result
                if enzyme_id.startswith(("gamma_", "f_"))
            ],
            [
                (enzyme_id, round(exp(nlp_result[enzyme_id]), 9))
                for enzyme_id in nlp_result
                if enzyme_id.startswith(("x_",))
            ],
        )
        print("----------------------------------------")

variability_dict["EX_S"] = (0.0, 14.0)
t0 = time.time()
result = perform_nlp_evolutionary_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    variability_dict=variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_alpha=False,
    with_iota=False,
    sampling_wished_num_feasible_starts=2,
    objvalue_json_path="examples/toymodel/evo_objvalues_genetic.json",
    evolution_num_gens=10,
)
t1 = time.time()
print(
    f"max(ATP_Consumption) from evolutionary algorithm under EX_S <= 14: {list(result.keys())[0]}"
)
print("TIME FOR COBRA-k evolutionary algorithm:", t1 - t0)
