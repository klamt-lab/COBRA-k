""""""
from copy import deepcopy
from dataclasses import dataclass
from os import cpu_count
from pyomo.common.errors import ApplicationError
from cobrak.constants import ALL_OK_KEY, OBJECTIVE_VAR_NAME, Z_VAR_PREFIX
from cobrak.dataclasses import ExtraLinearConstraint, Model, Solver, CorrectionConfig
from cobrak.evolution import is_objsense_maximization
from cobrak.lps import perform_lp_optimization
from cobrak.nlps import delete_unused_reactions_in_optimization_dict, perform_nlp_irreversible_optimization
from cobrak.standard_solvers import SCIP, IPOPT
from pydantic import validate_call
from joblib import Parallel, delayed
from cobrak.utilities import get_stoichiometrically_coupled_reactions
from random import randint, shuffle, choices, choice

@dataclass
class LpNlpBlockResult:
    """"""
    original_binaries: tuple[int, ...]
    lp_result: dict[str, float] = {}
    lp_binaries: tuple[int, ...]
    nlp_result: dict[str, float] = {}
    nlp_binaries: tuple[int, ...]


######## PRIVATE FUNCTIONS ########
@validate_call(validate_return=True)
def _get_binaries_from_opt_result(
    opt_result: dict[str, float],
    reac_couples_list: list[tuple[str, ...]],
) -> tuple[int, ...]:
    binaries = [0 for _ in range(len(reac_couples_list))]
    for i, reac_couple in enumerate(reac_couples_list):
        if reac_couple[0] in opt_result:
            binaries[i] = 1
    return tuple(binaries)


@validate_call(validate_return=True)
def _ectfba_block(
    cobrak_model_with_deletions: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    lp_solver: Solver,
    correction_config: CorrectionConfig,
    ignore_nonlinear_extra_terms_in_ectfbas: bool,
    binaries: tuple[int],
) -> tuple[dict[str, float], tuple[int, ...]]:
    try:
        ectfba_dict = perform_lp_optimization(
            cobrak_model=cobrak_model_with_deletions,
            objective_target=objective_target,
            objective_sense=objective_sense,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            with_loop_constraints=True,
            variability_dict=variability_dict,
            solver=lp_solver,
            correction_config=correction_config,
            ignore_nonlinear_terms=ignore_nonlinear_extra_terms_in_ectfbas,
        )
    except (ApplicationError, AttributeError, ValueError):
        return {}, binaries
    if not ectfba_dict[ALL_OK_KEY]:
        return {}, binaries
    return ectfba_dict, binaries


@validate_call(validate_return=True)
def _nlp_block(
    cobrak_model_with_deletions: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    nlp_solver: Solver,
    nlp_strict_mode: bool,
    nlp_single_strict_reacs: list[str],
    correction_config: CorrectionConfig,
) -> dict[str, float]:
    try:
        nlp_result: dict[str, float] = perform_nlp_irreversible_optimization(
            cobrak_model=cobrak_model_with_deletions,
            objective_target=objective_target,
            objective_sense=objective_sense,
            variability_dict=variability_dict,
            with_kappa=with_kappa,
            with_gamma=with_gamma,
            with_iota=with_iota,
            with_alpha=with_alpha,
            solver=nlp_solver,
            correction_config=correction_config,
            strict_mode=nlp_strict_mode,
            single_strict_reacs=nlp_single_strict_reacs,
        )
    except (ApplicationError, AttributeError, ValueError):
        return {}
    if not nlp_result[ALL_OK_KEY]:
        return {}

    return nlp_result

@validate_call(validate_return=True)
def _ectfba_nlp_block(
    cobrak_model_with_deletions: Model,
    binaries: tuple[int],
    reac_couples_list: list[tuple[str, ...]],
    lp_objective_target: str | dict[str, float],
    lp_objective_sense: int,
    lp_extra_linear_constraints: list[ExtraLinearConstraint],
    nlp_objective_target: str | dict[str, float],
    nlp_objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    lp_solver: Solver,
    nlp_solver: Solver,
    nlp_strict_mode: bool,
    nlp_single_strict_reacs: list[str],
    correction_config: CorrectionConfig,
    ignore_nonlinear_extra_terms_in_ectfbas: bool,
) -> LpNlpBlockResult:
    with cobrak_model_with_deletions as cobrak_model_with_deletions_and_extra_constraints:
        cobrak_model_with_deletions_and_extra_constraints.extra_linear_constraints += lp_extra_linear_constraints
        if lp_objective_target == "MAXZ":
            lp_objective_target = {
                f"{Z_VAR_PREFIX}{reac_id}": 1.0
                for (reac_id, reac_data) in cobrak_model_with_deletions_and_extra_constraints.reactions.items()
                if (reac_data.dG0 is not None)
                and (variability_dict[reac_id][1] > 0.0)
            }
        ectfba_dict: dict[str, float] = _ectfba_block(
            cobrak_model_with_deletions=cobrak_model_with_deletions_and_extra_constraints,
            objective_target=lp_objective_target,
            objective_sense=lp_objective_sense,
            variability_dict=variability_dict,
            lp_solver=lp_solver,
            correction_config=correction_config,
            ignore_nonlinear_extra_terms_in_ectfbas=ignore_nonlinear_extra_terms_in_ectfbas,
            binaries=binaries,
        )[1]
    error_target_missing: bool = any(errortarget not in ectfba_dict for errortarget in correction_config.error_scenario)
    if not ectfba_dict or not ectfba_dict[ALL_OK_KEY] or error_target_missing:
        return LpNlpBlockResult(original_binaries=binaries, lp_binaries=(), nlp_binaries=())

    nlp_result: dict[str, float] = _nlp_block(
        cobrak_model_with_deletions=delete_unused_reactions_in_optimization_dict(cobrak_model_with_deletions, ectfba_dict),
        objective_target=nlp_objective_target,
        objective_sense=nlp_objective_sense,
        variability_dict=variability_dict,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        nlp_solver=nlp_solver,
        nlp_strict_mode=nlp_strict_mode,
        nlp_single_strict_reacs=nlp_single_strict_reacs,
        correction_config=correction_config,
    )
    if not nlp_result or not nlp_result[ALL_OK_KEY]:
        return LpNlpBlockResult(
            original_binaries=binaries,
            lp_result=ectfba_dict,
            lp_binaries=_get_binaries_from_opt_result(ectfba_dict, reac_couples_list),
            nlp_binaries=(),
        )

    return LpNlpBlockResult(
        original_binaries=binaries,
        lp_result=ectfba_dict,
        lp_binaries=_get_binaries_from_opt_result(ectfba_dict, reac_couples_list),
        nlp_result=nlp_result,
        nlp_binaries=_get_binaries_from_opt_result(nlp_result, reac_couples_list)
    )

@validate_call(validate_return=True)
def _get_cobrak_model_with_deleted_binary_zero_reacs(
    cobrak_model: Model,
    binaries: tuple[int],
    reac_couples_list: list[tuple[str, ...]],
) -> Model:
    model_with_deletions: Model = deepcopy(cobrak_model)
    for couple_idx, binary in enumerate(binaries):
        if binary == 0:
            del model_with_deletions.reactions[reac_couples_list[couple_idx][0]]
    return model_with_deletions


@validate_call(validate_return=True)
def _add_eligible_binaries_and_get_best_nlp_solution(
    best_nlp_solution: dict[str, float],
    lpnlpblock_results: list[LpNlpBlockResult],
    binary_results: dict[tuple[int, ...], float | None],
    is_maximization: bool,
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    for result in lpnlpblock_results:
        if not result.lp_result:
            binary_results[result.original_binaries] = None
            continue
        if not result.nlp_result:
            binary_results[result.original_binaries] = None
            binary_results[result.lp_binaries] = None
            continue
        binary_results[result.original_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
        binary_results[result.lp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
        binary_results[result.nlp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]

        if not best_nlp_solution or\
           (is_maximization and result.nlp_result[OBJECTIVE_VAR_NAME] > best_nlp_solution[OBJECTIVE_VAR_NAME]) or\
           (not is_maximization and result.nlp_result[OBJECTIVE_VAR_NAME] < best_nlp_solution[OBJECTIVE_VAR_NAME]):
            best_nlp_solution = deepcopy(result.nlp_result)

    return best_nlp_solution, binary_results

@validate_call(validate_return=True)
def _sampling(
    cobrak_model: Model,
    reac_couples_list: list[tuple[str, ...]],
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    correction_config: CorrectionConfig,
    lp_solver: Solver,
    nlp_solver: Solver,
    nlp_strict_mode: bool,
    nlp_single_strict_reacs: list[str],
    ignore_nonlinear_extra_terms_in_lps: bool,
    num_cpus: int = 1,
    max_targeted_reacs: int = 5,
    max_calculations: int = 50,
    num_wished_starts: int = 3,
) -> dict[tuple[int, ...], float | None]:
    num_reacs = len(cobrak_model.reactions)
    samples: dict[tuple[int, ...], float | None] = {}
    best_nlp_solution: dict[str, float] = {}
    num_calculations = 0
    while num_calculations < max_calculations:
        binary_scenarios: list[tuple[int, ...]] = []
        for _ in range(num_cpus):
            num_zeroes: int = randint(0, min(num_reacs, max_targeted_reacs))
            zero_one_mix: list[int] = [1] * (num_reacs - num_zeroes) + [0] * num_zeroes
            shuffle(zero_one_mix)
            binary_scenarios.append(
                tuple(zero_one_mix)
            )
        binary_scenarios = [scenario for scenario in binary_scenarios if scenario not in samples]
        results: list[LpNlpBlockResult] = Parallel(n_jobs=-1, verbose=10)(
            delayed(_ectfba_nlp_block)(
                _get_cobrak_model_with_deleted_binary_zero_reacs(
                    cobrak_model,
                    binary_scenario,
                    reac_couples_list
                ),
                binary_scenario,
                reac_couples_list,
                objective_target,
                objective_sense,
                [],
                objective_target,
                objective_sense,
                variability_dict,
                with_kappa,
                with_gamma,
                with_iota,
                with_alpha,
                lp_solver,
                nlp_solver,
                nlp_strict_mode,
                nlp_single_strict_reacs,
                correction_config,
                ignore_nonlinear_extra_terms_in_lps,
            )
            for binary_scenario in binary_scenarios
        )
        best_nlp_solution, samples = _add_eligible_binaries_and_get_best_nlp_solution(
            best_nlp_solution=best_nlp_solution,
            lpnlpblock_results=results,
            binary_results=samples,
            is_maximization=is_objsense_maximization(objective_sense),
        )

        if len(samples) >= num_wished_starts:
            break
        num_calculations += num_cpus
    return samples


@validate_call(validate_return=True)
def _get_binaries_according_to_selection(
    sorted_results: dict[tuple[int, ...], float],
    selection_method: str
) -> tuple[int, ...]:
    binaries: tuple[int, ...]
    match selection_method:
        case "weighted":
            binaries = deepcopy(choices(
                population=list(sorted_results.keys()),
                weights=list(sorted_results.values()),
                k=1,
            )[0])
        case "elite":
            binaries = deepcopy(choice(list(sorted_results.keys())[:9]))
        case "random":
            binaries = deepcopy(choice(list(sorted_results.keys())))
        case _:
            raise ValueError
    return binaries


@validate_call(validate_return=True)
def _get_evolution_binaries(
    population_size: int,
    num_reac_couples: int,
    evolution_results: dict[tuple[int, ...], float],
    fractions_genetic_method: dict[str, float],
    fractions_population_selection: dict[str, float],
    is_maximization: bool,
) -> list[tuple[int, ...]]:
    non_na_evolution_results: dict[tuple[int, ...], float] = {
        key: value
        for key, value in evolution_results.items()
        if value is not None
    }
    sorted_results: dict[tuple[int, ...], float] = dict(
        sorted(
            non_na_evolution_results.items(),
            key=lambda item: item[1],
            reverse=is_maximization,
        )
    )
    binaries: list[tuple[int, ...]] = []
    for _ in range(population_size):
        selection_method: str = choices(
            population=list(fractions_population_selection.keys()),
            weights=list(fractions_population_selection.values()),
            k=1,
        )[0]
        first_binaries: tuple[int, ...] = _get_binaries_according_to_selection(
            sorted_results=sorted_results,
            selection_method=selection_method,
        )

        genetic_method: str  = choices(
            population=list(fractions_genetic_method.keys()),
            weights=list(fractions_genetic_method.values()),
            k=1,
        )[0]
        match genetic_method:
            case "neighborhood":
                num_tries = 0
                while first_binaries in sorted_results:
                    flip_location = randint(0, num_reac_couples-1)
                    first_binaries[flip_location] = tuple(
                        list(first_binaries[:flip_location]) + [int(not binaries[-1][flip_location])] + list(first_binaries[flip_location+1:])
                    )
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(first_binaries)
            case "random":
                binaries.append(
                    tuple([randint(0, 1) for _ in range(num_reac_couples)])
                )
            case "multimutation":
                num_tries = 0
                while first_binaries in sorted_results:
                    flip_locations: list[int] = [randint(0, num_reac_couples-1) for _ in range(3)]
                    for flip_location in flip_locations:
                        first_binaries[flip_location] = tuple(
                            list(first_binaries[:flip_location]) + [int(not binaries[-1][flip_location])] + list(first_binaries[flip_location+1:])
                        )
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(first_binaries)
            case "crossover":
                second_binaries: tuple[int, ...] = _get_binaries_according_to_selection(
                    sorted_results=sorted_results,
                    selection_method=selection_method,
                )
                num_tries = 0
                crossed_over_binaries: list[int]
                while first_binaries in sorted_results:
                    crossover_point = randint(0, num_reac_couples-1)
                    crossed_over_binaries = first_binaries[:crossover_point] + second_binaries[crossover_point:]
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(crossed_over_binaries)
            case _:
                raise ValueError
    return binaries


@validate_call(validate_return=True)
def _evolution(
    cobrak_model: Model,
    reac_couples_list: list[tuple[str, ...]],
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    pre_results: dict[tuple[int, ...], float | None],
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    correction_config: CorrectionConfig,
    lp_solver: Solver,
    nlp_solver: Solver,
    nlp_strict_mode: bool,
    nlp_single_strict_reacs: list[str],
    num_gens: int,
    population_size: int,
    ignore_nonlinear_extra_terms_in_lps: bool,
    max_rounds_same_objvalue: float = float("inf"),
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    opt_selector = max if is_objsense_maximization(objective_sense) else min
    evolution_results: dict[tuple[int, ...], float | None] = deepcopy(pre_results)
    if type(objective_target) is str:
        objective_target_as_dict: dict[str, int | float] = {objective_target: 1.0}
    elif type(objective_target) is dict:
        objective_target_as_dict: dict[str, int | float] = objective_target

    best_nlp_solution: dict[str, float] = {}
    current_best_objvalue: float = opt_selector([value for value in evolution_results.values() if value is not None])
    num_rounds_with_same_objvalue = 0
    for _ in range(num_gens):
        tested_binaries = _get_evolution_binaries(
            population_size=population_size,
            num_reac_couples=len(reac_couples_list),
            evolution_results=evolution_results,
            fractions_genetic_method={
                "neighborhood": 1/8,
                "random": 1/8,
                "multimutation": 1/4,
                "crossover": 1/2,
            },
            fractions_population_selection={
                "weighted": 1/2,
                "elite": 1/4,
                "random": 1/4,
            },
            is_maximization=is_objsense_maximization(objective_sense),
        )

        ectfba_results: list[dict[str, float]] = Parallel(n_jobs=-1, verbose=10)(
            delayed(_ectfba_block)(
                _get_cobrak_model_with_deleted_binary_zero_reacs(
                    cobrak_model=cobrak_model,
                    binaries=tested_binary,
                    reac_couples_list=reac_couples_list,
                ),
                objective_target,
                objective_sense,
                variability_dict,
                lp_solver,
                correction_config,
                ignore_nonlinear_extra_terms_in_lps,
                tested_binary,
            )
            for tested_binary in tested_binaries
        )
        eligible_binaries_with_objvalue: list[tuple[tuple[int, ...], float]] = []
        for (original_binaries, ectfba_result) in ectfba_results:
            if not ectfba_result or not ectfba_result[ALL_OK_KEY]:
                evolution_results[original_binaries] = None
            else:
                eligible_binaries_with_objvalue.append((_get_binaries_from_opt_result(ectfba_result), ectfba_results[OBJECTIVE_VAR_NAME]))

        results: list[LpNlpBlockResult] = Parallel(n_jobs=-1, verbose=10)(
            delayed(_ectfba_nlp_block)(
                _get_cobrak_model_with_deleted_binary_zero_reacs(
                    cobrak_model=cobrak_model,
                    binaries=eligible_binary,
                    reac_couples_list=reac_couples_list,
                ),
                eligible_binary,
                "MAXZ",
                objective_sense,
                objective_target,
                [
                    ExtraLinearConstraint(
                        stoichiometries=objective_target_as_dict,
                        lower_value=objvalue - 1e-9,
                        upper_value=objvalue + 1e-9,
                    )
                ],
                objective_sense,
                variability_dict,
                with_kappa,
                with_gamma,
                with_iota,
                with_alpha,
                lp_solver,
                nlp_solver,
                nlp_strict_mode,
                nlp_single_strict_reacs,
                correction_config,
                ignore_nonlinear_extra_terms_in_lps,
            )
            for (eligible_binary, objvalue) in eligible_binaries_with_objvalue
        )
        best_nlp_solution, evolution_results = _add_eligible_binaries_and_get_best_nlp_solution(
            best_nlp_solution=best_nlp_solution,
            lpnlpblock_results=results,
            binary_results=evolution_results,
            is_maximization=is_objsense_maximization(objective_sense),
        )

        gen_best_objvalue = opt_selector(evolution_results.values())
        if gen_best_objvalue != current_best_objvalue:
            current_best_objvalue = gen_best_objvalue
            num_rounds_with_same_objvalue = 0
        else:
            num_rounds_with_same_objvalue += 1
        if num_rounds_with_same_objvalue >= max_rounds_same_objvalue:
            break

    return best_nlp_solution, evolution_results


@validate_call(validate_return=True)
def _get_idx_to_reac_ids(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    variability_dict: dict[str, tuple[float, float]],
    error_scenario: dict[str, tuple[float, float]],
) -> tuple[tuple[str, ...], ...]:
    reac_couples: list[list[str]] = get_stoichiometrically_coupled_reactions(
        cobrak_model=cobrak_model,
        rounding=10,
    )
    reac_couples_list: list[tuple[str, ...]] = []
    for reac_ids in reac_couples:
        # Discard couple withn blocked, essential, non-kinetic reactions and ones with error targets
        if any(
            variability_dict[reac_id][1] <= 0.0 or
            variability_dict[reac_id][0] > 0.0 or
            (not cobrak_model.reactions[reac_id].dG0 and not cobrak_model.reactions[reac_id].enzyme_reaction_data) or
            reac_id in error_scenario
            for reac_id in reac_ids
        ):
            continue
        # Discard couple with objective target(s)
        objective_target_strlist: list[str]
        if type(objective_target) is str:
            objective_target_strlist = [objective_target]
        elif type(objective_target) is dict:
            objective_target_strlist = list(objective_target.keys())
        if any(objective_target in reac_ids for objective_target in objective_target_strlist):
            continue
        reac_couples_list.append(tuple(reac_ids))

    return tuple(reac_couples_list)


######## PUBLIC FUNCTIONS ########
@validate_call(validate_return=True)
def perform_nlp_evolutionary_optimization(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    num_gens: int,
    population_size: int,
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
    lp_solver: Solver = SCIP,
    nlp_solver: Solver = IPOPT,
    nlp_strict_mode: bool = False,
    nlp_single_strict_reacs: list[str] = [],
    max_rounds_same_objvalue: float = float("inf"),
    ignore_nonlinear_extra_terms_in_lps: bool = True,
    existing_evolution_result: dict[float, list[dict[str, float]]] = {},
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    """"""
    # PHASE 1: BUILD INDEX TO REAC COUPLES DATA, AND CPU DATA
    reac_couples_list = _get_idx_to_reac_ids(
        cobrak_model=cobrak_model,
        objective_target=objective_target,
        variability_dict=variability_dict,
        error_scenario=correction_config.error_scenario,
    )
    num_cpus_raw: int | None = cpu_count()
    num_cpus: int = 1 if num_cpus_raw is None else num_cpus_raw

    # PHASE 2: SAMPLING (IF NO EXISTING SOLUTION GIVEN)
    if not existing_evolution_result:
        sampling_solutions = _sampling(
            cobrak_model=cobrak_model,
            reac_couples_list=reac_couples_list,
            objective_target=objective_target,
            objective_sense=objective_sense,
            variability_dict=variability_dict,
            with_kappa=with_kappa,
            with_gamma=with_gamma,
            with_iota=with_iota,
            with_alpha=with_alpha,
            correction_config=correction_config,
            lp_solver=lp_solver,
            nlp_solver=nlp_solver,
            nlp_strict_mode=nlp_strict_mode,
            nlp_single_strict_reacs=nlp_single_strict_reacs,
            ignore_nonlinear_extra_terms_in_lps=ignore_nonlinear_extra_terms_in_lps,
            num_cpus=num_cpus,
            max_targeted_reacs=5,
            max_calculations=50,
            num_wished_starts=3,
        )
        if all(result is None for result in sampling_solutions.values()):
            print("ERROR: No working result given or found in sampling or given solutions")
            raise ValueError
        existing_evolution_result = sampling_solutions

    # PHASE 3: ACTUAL EVOLUTION ALGORITHM (USING GIVEN OR SAMPLED RESULTS AS STARTING POINTS)
    return _evolution(
        cobrak_model=cobrak_model,
        reac_couples_list=reac_couples_list,
        objective_target=objective_target,
        objective_sense=objective_sense,
        variability_dict=variability_dict,
        pre_results=existing_evolution_result,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        num_gens=num_gens,
        population_size=population_size,
        correction_config=correction_config,
        lp_solver=lp_solver,
        nlp_solver=nlp_solver,
        nlp_strict_mode=nlp_strict_mode,
        nlp_single_strict_reacs=nlp_single_strict_reacs,
        ignore_nonlinear_extra_terms_in_lps=ignore_nonlinear_extra_terms_in_lps,
        max_rounds_same_objvalue=max_rounds_same_objvalue,
    )
