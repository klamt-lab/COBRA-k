# IMPORTS SECTION #
from os.path import exists
from sys import argv

import z_add_path  # noqa: F401

from cobrak.constants import PROT_POOL_REAC_NAME
from cobrak.dataclasses import (
    Model,
    VarResult,
)
from cobrak.evolution import (
    perform_nlp_evolutionary_optimization,
    postprocess,
)
from cobrak.io import (
    ensure_folder_existence,
    json_load,
    json_write,
    json_zip_write,
)
from cobrak.plotting import plot_objvalue_evolution
from cobrak.spreadsheet_functionality import (
    OptimizationDataset,
)
from cobrak.standard_solvers import CPLEX, IPOPT_MA57
from cobrak.utilities import (
    create_cnapy_scenario_out_of_optimization_dict,
)

try:
    corr_round = int(argv[1])
except (ValueError, IndexError):
    corr_round = None

# RUNNING SCRIPT SECTION #
biomass_reac_id = "Biomass_fw"
cobrak_model: Model = json_load(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak_prestepC_uncalibrated.json",
    Model,
)

variability_dict: VarResult = json_load(
    "examples/iCH360/prepared_external_resources/variability_dict_enforced_growth_uncalibrated.json",
    VarResult,
)

variability_dict["Biomass_fw"] = (0.7, variability_dict["Biomass_fw"][1])
mustring = "minmu0.7_X"
for arg in argv:
    if arg.startswith("minmu"):
        variability_dict["Biomass_fw"] = (
            float(arg[len("minmu") :]),
            variability_dict["Biomass_fw"][1],
        )
        mustring = arg
        print("Min µ is", variability_dict["Biomass_fw"][1])
    if arg.startswith("protpool"):
        cobrak_model.max_prot_pool = float(arg[len("protpool") :])


results_folder = f"examples/iCH360/results_CALIBRATION_RUN_{mustring}/"
ensure_folder_existence(results_folder)

optimization_datasets: dict[str, OptimizationDataset] = {}
corr_rounds = [0] if corr_round is None else [corr_round]

for corr_round in corr_rounds:
    file_suffix: str = f"_{corr_round}"
    objvalue_json_path = f"{results_folder}objvalues{file_suffix}.json"
    final_best_result_json_path = (
        f"{results_folder}final_best_result_{file_suffix}.json"
    )
    evolution_best_result_json_path = (
        f"{results_folder}best_evolution_result{file_suffix}.json"
    )

    if not exists(final_best_result_json_path):
        objvalue_png_path = f"{results_folder}objvalues{file_suffix}.png"

        if not exists(evolution_best_result_json_path):
            results = perform_nlp_evolutionary_optimization(
                cobrak_model=cobrak_model,
                objective_target=PROT_POOL_REAC_NAME,
                objective_sense=-1,
                variability_dict=variability_dict,
                sampling_always_deactivated_reactions=[],
                with_kappa=True,
                with_gamma=True,
                with_iota=False,
                with_alpha=False,
                sampling_wished_num_feasible_starts=5,
                evolution_num_gens=100,
                lp_solver=CPLEX,
                nlp_solver=IPOPT_MA57,
                objvalue_json_path=objvalue_json_path,
                pop_size=64,
            )
            evolution_best_result = results[list(results.keys())[0]][0]
            json_write(evolution_best_result_json_path, evolution_best_result)
        else:
            evolution_best_result = json_load(evolution_best_result_json_path)

        create_cnapy_scenario_out_of_optimization_dict(
            f"{results_folder}best_evolution_result{file_suffix}.scen",
            cobrak_model,
            evolution_best_result,
        )
        json_zip_write(
            f"{results_folder}full_evolution_results{file_suffix}.json",
            results,
        )

        plot_objvalue_evolution(
            json_path=objvalue_json_path,
            output_path=objvalue_png_path,
            ylabel="correction error",
        )

        last_result = evolution_best_result
        postprocess_round = 0
        while True:
            postprocess_results, best_postprocess_result = postprocess(
                cobrak_model=cobrak_model,
                opt_dict=last_result,
                objective_target=PROT_POOL_REAC_NAME,
                objective_sense=-1,
                variability_data=variability_dict,
                lp_solver=CPLEX,
                nlp_solver=IPOPT_MA57,
            )

            if len(postprocess_results) == 0:
                break
            json_write(
                f"{results_folder}postprocess_round{postprocess_round}_full_result_{file_suffix}.json",
                postprocess_results,
            )
            json_write(
                f"{results_folder}postprocess_round{postprocess_round}_best_result_{file_suffix}.json",
                best_postprocess_result,
            )
            last_result = best_postprocess_result
            postprocess_round += 1
        json_write(
            f"{results_folder}final_best_result_{file_suffix}.json",
            last_result,
        )
