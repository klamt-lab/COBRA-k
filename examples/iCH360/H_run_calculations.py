# IMPORTS SECTION #
from copy import deepcopy
from dataclasses import dataclass
from os.path import exists
from sys import argv
from time import time

import z_add_path  # noqa: F401
from pydantic import Field

from cobrak.constants import OBJECTIVE_VAR_NAME
from cobrak.dataclasses import (
    ExtraLinearConstraint,
    Model,
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
    standardize_folder,
)
from cobrak.lps import perform_lp_variability_analysis
from cobrak.plotting import plot_objvalue_evolution
from cobrak.standard_solvers import CPLEX, CPLEX_FOR_VARIABILITY_ANALYSIS, IPOPT_MA57
from cobrak.utilities import (
    create_cnapy_scenario_out_of_optimization_dict,
    delete_orphaned_metabolites_and_enzymes,
    get_model_with_varied_parameters,
    is_objsense_maximization,
)

try:
    corr_round = int(argv[1])
except (ValueError, IndexError):
    corr_round = None


@dataclass
class RunConfig:
    # Model changes
    manually_changed_kms: dict[str, dict[str, float]]
    manually_changed_kcats: dict[str, float]
    manually_changed_dG0s: dict[str, float]
    # Folder settings
    results_folder: str
    # ecTFVA settings
    ectfva_active_reacs: list[str]
    # Evolutionary algorithm settings
    round_num: int
    objective_target: str | dict[str, float]
    objective_sense: int
    deactivated_reacs: list[str]
    set_bounds: dict[str, tuple[float, float]]
    sampling_rounds_per_metaround: int = 2
    sampling_wished_num_feasible_starts: int = 5
    sampling_max_metarounds: int = 3
    evolution_num_gens: int = 150
    pop_size: int = 32
    protein_pool: float | None = None
    uses_bennett_concs: bool = False
    working_results: list[dict[str, float]] = Field(default_factory=list)
    max_conc_sum: float | None = None
    nameaddition: str | None = None
    kicked_reacs: list[str] = Field(default_factory=list)
    do_parameter_variation: bool = False
    varied_reacs: list[str] = Field(default_factory=list)
    max_km_variation: float | None = None
    max_kcat_variation: float | None = None
    max_ki_variation: float | None = None
    max_ka_variation: float | None = None
    max_dG0_variation: float | None = None
    with_iota: bool = False
    with_alpha: bool = False


# LOAD RUN CONFIGURATION #
run_config: RunConfig = json_load(argv[1], RunConfig)
file_suffix: str = (
    f"_{run_config.round_num}"
    if run_config.nameaddition is None
    else f"_{run_config.round_num}_{run_config.nameaddition}"
)

# CREATE RESULTS FOLDER #
run_config.results_folder = (
    f"{'' if run_config.results_folder.startswith('.') else '.'}"
    + standardize_folder(run_config.results_folder)
)
ensure_folder_existence(run_config.results_folder)

# LOAD AND CHANGE COBRAK MODEL #
used_cobrak_model_path = (
    f"{run_config.results_folder}used_cobrak_model_{file_suffix}.json"
)
if not exists(used_cobrak_model_path):
    cobrak_model: Model = json_load(
        "examples/iCH360/prepared_external_resources/iCH360_cobrak.json",
        Model,
    )
    if run_config.do_parameter_variation:
        cobrak_model = get_model_with_varied_parameters(
            cobrak_model,
            run_config.max_km_variation,
            run_config.max_kcat_variation,
            run_config.max_ki_variation,
            run_config.max_ka_variation,
            run_config.max_dG0_variation,
            run_config.varied_reacs,
        )

    for kicked_reac in run_config.kicked_reacs:
        if kicked_reac in cobrak_model.reactions:
            del cobrak_model.reactions[kicked_reac]
    if len(run_config.kicked_reacs) > 0:
        cobrak_model = delete_orphaned_metabolites_and_enzymes(cobrak_model)

    if run_config.protein_pool is not None:
        assert run_config.protein_pool > 0
        cobrak_model.max_prot_pool = run_config.protein_pool
    for reac_id, changed_kms_dict in run_config.manually_changed_kms.items():
        if reac_id in run_config.kicked_reacs:
            continue
        for met_id, changed_km in changed_kms_dict.items():
            assert changed_km > 0
            cobrak_model.reactions[reac_id].enzyme_reaction_data.k_ms[met_id] = (
                changed_km
            )
    for reac_id, changed_kcat in run_config.manually_changed_kcats.items():
        if reac_id in run_config.kicked_reacs:
            continue
        assert changed_kcat > 0
        cobrak_model.reactions[reac_id].enzyme_reaction_data.k_cat = changed_kcat
    for reac_id, changed_dG0 in run_config.manually_changed_dG0s.items():
        if reac_id in run_config.kicked_reacs:
            continue
        # assert float(changed_dG0)
        cobrak_model.reactions[reac_id].dG0 = changed_dG0

    for deactivated_reac in run_config.deactivated_reacs:
        cobrak_model.reactions[deactivated_reac].min_flux = 0.0
        cobrak_model.reactions[deactivated_reac].max_flux = 0.0

    for var_id, bound_tuple in run_config.set_bounds.items():
        if not var_id:
            continue
        cobrak_model.extra_linear_constraints.append(
            ExtraLinearConstraint(
                stoichiometries={
                    var_id: 1.0,
                },
                lower_value=bound_tuple[0],
                upper_value=bound_tuple[1],
            )
        )

    if run_config.max_conc_sum is not None:
        assert run_config.max_conc_sum > 0
        cobrak_model.max_conc_sum = run_config.max_conc_sum

    json_write(used_cobrak_model_path, cobrak_model)

cobrak_model: Model = json_load(
    used_cobrak_model_path,
    Model,
)


var_dict_filepath = f"{run_config.results_folder}variability_dict_{file_suffix}.json"
if not exists(var_dict_filepath):
    cobrak_model_with_slightly_more_protein = deepcopy(cobrak_model)
    cobrak_model_with_slightly_more_protein.max_prot_pool *= 1.02
    t0 = time()
    variability_dict = perform_lp_variability_analysis(
        cobrak_model_with_slightly_more_protein,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        min_flux_cutoff=1e-6,
        solver=CPLEX_FOR_VARIABILITY_ANALYSIS,
        active_reactions=run_config.ectfva_active_reacs,
        min_active_flux=1e-5,
    )
    json_write(var_dict_filepath, variability_dict)
    t1 = time()
    json_write(
        f"{run_config.results_folder}variability_dict_{file_suffix}_used_time.json",
        [t1 - t0],
    )
else:
    variability_dict = json_load(var_dict_filepath, dict[str, tuple[float, float]])

objvalue_json_path = f"{run_config.results_folder}objvalues_over_time{file_suffix}.json"
final_best_result_json_path = (
    f"{run_config.results_folder}final_best_result_{file_suffix}.json"
)
evolution_best_result_json_path = (
    f"{run_config.results_folder}best_evolution_result{file_suffix}.json"
)

if exists(final_best_result_json_path):
    print("INFO: Final best result (after genetic algorithm and postprocessing)")
    print(f"already exists in {final_best_result_json_path}")
    print("Hence, no calculation is needed.")
    print(
        "If you want to recalculate the postprocesing, delete the mentioned .json file."
    )
    print(
        "If you also want to recalculate the genetic algorithm (done before the postprocessing),"
    )
    print(f"also delete {evolution_best_result_json_path}")
else:
    evolution_best_result_time_json_path = (
        f"{run_config.results_folder}best_evolution_result{file_suffix}_used_time.json"
    )
    objvalue_png_path = f"{run_config.results_folder}objvalues{file_suffix}.png"

    if exists(evolution_best_result_json_path):
        print(
            f"INFO: Best genetic algorithm result already exists as {evolution_best_result_json_path}."
        )
        print(
            "Hence, only the postprocessing is performed. If you want to recalculate everything,"
        )
        print("delete the mentioned .json file.")
        evolution_best_result = json_load(evolution_best_result_json_path)
    else:
        t0 = time()
        results = perform_nlp_evolutionary_optimization(
            cobrak_model=cobrak_model,
            objective_target=run_config.objective_target,
            objective_sense=run_config.objective_sense,
            variability_dict=variability_dict,
            sampling_always_deactivated_reactions=[],
            with_kappa=True,
            with_gamma=True,
            with_iota=run_config.with_iota,
            with_alpha=run_config.with_alpha,
            evolution_num_gens=run_config.evolution_num_gens,
            lp_solver=CPLEX,
            nlp_solver=IPOPT_MA57,
            objvalue_json_path=objvalue_json_path,
            pop_size=run_config.pop_size,
            sampling_rounds_per_metaround=run_config.sampling_rounds_per_metaround,
            sampling_wished_num_feasible_starts=run_config.sampling_wished_num_feasible_starts,
            sampling_max_metarounds=run_config.sampling_max_metarounds,
            working_results=run_config.working_results,
        )
        evolution_best_result = results[list(results.keys())[0]][0]
        json_write(evolution_best_result_json_path, evolution_best_result)
        json_zip_write(
            f"{run_config.results_folder}full_evolution_results{file_suffix}.json",
            results,
        )
        t1 = time()
        json_write(evolution_best_result_time_json_path, [t1 - t0])

    plot_objvalue_evolution(
        json_path=objvalue_json_path,
        output_path=objvalue_png_path,
        ylabel="correction error",
    )

    t0 = time()
    last_working_result = deepcopy(evolution_best_result)
    postprocess_round = 0
    while True:
        postprocess_results, best_postprocess_result = postprocess(
            cobrak_model=cobrak_model,
            opt_dict=last_working_result,
            objective_target=run_config.objective_target,
            objective_sense=run_config.objective_sense,
            variability_data=variability_dict,
            lp_solver=CPLEX,
            nlp_solver=IPOPT_MA57,
            with_iota=run_config.with_iota,
            with_alpha=run_config.with_alpha,
        )

        if (
            len(postprocess_results) == 0
            or best_postprocess_result == {}
            or OBJECTIVE_VAR_NAME not in best_postprocess_result
        ):
            break
        set_new_best = False
        if (
            is_objsense_maximization(run_config.objective_sense)
            and best_postprocess_result[OBJECTIVE_VAR_NAME]
            > last_working_result[OBJECTIVE_VAR_NAME] + 1e-6
        ) or best_postprocess_result[OBJECTIVE_VAR_NAME] < last_working_result[
            OBJECTIVE_VAR_NAME
        ] - 1e-6:
            set_new_best = True
        if not set_new_best:
            break
        # json_write(
        #    f"{run_config.results_folder}postprocess_round{postprocess_round}_{file_suffix}_full_result.json",
        #    postprocess_results,
        # )
        json_write(
            f"{run_config.results_folder}postprocess_round{postprocess_round}_{file_suffix}_best_result.json",
            best_postprocess_result,
        )
        last_working_result: list[float | int] = deepcopy(best_postprocess_result)
        postprocess_round += 1

    json_write(
        f"{run_config.results_folder}final_best_result_{file_suffix}.json",
        last_working_result,
    )
    t1 = time()
    json_write(
        f"{run_config.results_folder}final_best_result_{file_suffix}_time.json",
        [t1 - t0],
    )

    create_cnapy_scenario_out_of_optimization_dict(
        f"{run_config.results_folder}final_best_result_cnapy_{file_suffix}.scen",
        cobrak_model,
        last_working_result,
    )
