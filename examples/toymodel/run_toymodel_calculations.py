import os
import time
from copy import deepcopy

try:
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

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
from cobrak.nlps import (
    perform_nlp_reversible_optimization,
)
from cobrak.plotting import plot_objvalue_evolution
from cobrak.printing import (
    print_model,
    print_optimization_result,
    print_variability_result,
)
from cobrak.spreadsheet_functionality import (
    OptimizationDataset,
    VariabilityDataset,
    create_cobrak_spreadsheet,
)
from cobrak.standard_solvers import IPOPT, SCIP

side_reac_id = "Glycolysis"
main_reac_ids = ["Respiration", "Overflow"]

save_cobrak_model_as_annotated_sbml_model(
    toy_model, filepath="examples/toymodel/sbml_model.xml"
)
toy_model = load_annotated_sbml_model_as_cobrak_model(
    filepath="examples/toymodel/sbml_model.xml"
)

print_model(toy_model)

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

nlp_result_rev = perform_nlp_reversible_optimization(
    cobrak_model=toy_model,
    objective_target="EX_ATP",
    objective_sense=+1,
    variability_dict=test_variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_iota=False,
    with_alpha=False,
    verbose=False,
    solver=SCIP,
    # solver=BARON,
    show_variable_count=True,
)
print("...done! First MINLP result:")
print_optimization_result(toy_model, nlp_result_rev)
print()

print("----------------------------------------")

for max_glc_uptake in [30, 10.0]:
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
        test_variability_dict["EX_A"] = (0.0, max_glc_uptake)

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

        lp_result = perform_lp_optimization(
            cobrak_model=run_cobrak_model,
            objective_target="EX_ATP",
            objective_sense=+1,
            variability_dict=test_variability_dict,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            with_loop_constraints=False,
            verbose=False,
        )
        print(
            f"ecTFBA | max(EX_ATP) with {active_reacs}:",
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
                "EX_ATP",
            ]:
                lp_result[reac_id] = 1.0
        nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
            deepcopy(toy_model),
            objective_target="EX_ATP",
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

        print(
            f"NLP    | max(EX_ATP) with {active_reacs}:",
            round(nlp_result[OBJECTIVE_VAR_NAME], 4),
            "uptake:",
            round(nlp_result["EX_S"], 2),
            [
                (enzyme_id, round(nlp_result[enzyme_id], 100))
                for enzyme_id in nlp_result
                if enzyme_id.startswith(("gamma_", "f_"))
            ],
        )
        print("----------------------------------------")

print("----------------------------------------")

create_cobrak_spreadsheet(
    path="examples/toymodel/result.xlsx",
    cobrak_model=toy_model,
    variability_datasets={
        "tfva": VariabilityDataset(variability_dict, with_df=True),
    },
    optimization_datasets={
        "nlp_irr": OptimizationDataset(
            nlp_result,
            with_df=True,
            with_vplus=True,
            with_kappa=True,
            with_gamma=True,
            with_iota=False,
            with_kinetic_differences=True,
        ),
        "nlp_rev": OptimizationDataset(
            nlp_result_rev,
            with_df=True,
            with_vplus=True,
            with_kappa=True,
            with_gamma=True,
            with_iota=False,
            with_kinetic_differences=True,
        ),
    },
    is_maximization=True,
)

variability_dict["EX_A"] = (0.0, 14.0)
for algorithm in ("genetic",):
    t0 = time.time()
    result = perform_nlp_evolutionary_optimization(
        cobrak_model=toy_model,
        objective_target="EX_ATP",
        objective_sense=+1,
        variability_dict=variability_dict,
        with_kappa=True,
        with_gamma=True,
        with_alpha=False,
        with_iota=False,
        sampling_wished_num_feasible_starts=2,
        algorithm=algorithm,  # type: ignore
        objvalue_json_path=f"examples/toymodel/evo_objvalues_{algorithm}.json",
        evolution_num_gens=10,
    )
    t1 = time.time()
    print("TIME FOR COBRA-k:", t1 - t0)
    plot_objvalue_evolution(
        f"examples/toymodel/evo_objvalues_{algorithm}.json",
        f"examples/toymodel/evo_objvalues_{algorithm}.png",
        algorithm=algorithm,
    )  # type: ignore

x1 = time.time()
nlp_result_rev = perform_nlp_reversible_optimization(
    cobrak_model=toy_model,
    objective_target="EX_ATP",
    objective_sense=+1,
    variability_dict=test_variability_dict,
    with_kappa=True,
    with_gamma=True,
    with_iota=False,
    with_alpha=False,
    verbose=False,
    solver=SCIP,
    # solver=BARON,
)
x2 = time.time()
print("MINLP TIME", x2 - x1)
json_write("examples/toymodel/minlp_result.json", nlp_result_rev)
print("...done! MINLP result:")
print("Objective", nlp_result_rev["EX_ATP"])
print()
