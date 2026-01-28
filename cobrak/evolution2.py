""""""
from audioop import mul
from copy import deepcopy
from dataclasses import dataclass
from os import cpu_count
from numpy.distutils.npy_pkg_config import parse_sections
from pyomo.common.errors import ApplicationError
from cobrak.constants import ALL_OK_KEY, OBJECTIVE_VAR_NAME
from cobrak.dataclasses import Model, Solver, CorrectionConfig
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
    original_binaries: tuple[int]
    lp_result: dict[str, float] = {}
    lp_binaries: tuple[int] = []
    nlp_result: dict[str, float] = {}
    nlp_binaries: tuple[int] = []


######## PRIVATE FUNCTIONS ########
@validate_call(validate_return=True)
def _get_binaries_from_opt_result(
    opt_result: dict[str, float],
    reac_couples_list: list[tuple[str, ...]],
) -> tuple[int]:
    binaries = [0 for _ in range(len(reac_couples_list))]
    for i, reac_couple in enumerate(reac_couples_list):
        if reac_couple[0] in opt_result and opt_result > reac_couple[0]:
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
) -> tuple[tuple[int], dict[str, float]]:
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
    binaries: list[int],
    reac_couples_list: list[tuple[str, ...]],
    lp_objective_target: str | dict[str, float],
    lp_objective_sense: int,
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
    ectfba_dict: dict[str, float] = _ectfba_block(
        cobrak_model_with_deletions=cobrak_model_with_deletions,
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
        return LpNlpBlockResult(original_binaries=binaries)

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
            ectfba_dict=ectfba_dict,
            lp_binaries=_get_binaries_from_opt_result(ectfba_dict, reac_couples_list),
        )

    return LpNlpBlockResult(
        original_binaries=binaries,
        ectfba_dict=ectfba_dict,
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
    num_calculations = 0
    while num_calculations < max_calculations:
        binary_scenarios: list[tuple[int]] = []
        for _ in range(num_cpus):
            num_zeroes = randint(0, min(num_reacs, max_targeted_reacs))
            binary_scenarios.append(
                tuple(shuffle([1] * (num_reacs - num_zeroes) + [0] * num_zeroes))
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
        for result in results:
            if not result.lp_result:
                samples[result.original_binaries] = None
                continue
            if not result.nlp_result:
                samples[result.original_binaries] = None
                samples[result.lp_binaries] = None
                continue
            samples[result.original_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
            samples[result.lp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
            samples[result.nlp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]

        if len(samples) >= num_wished_starts:
            break
        num_calculations += num_cpus
    return samples


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
        first_binaries: tuple[int, ...]
        second_binaries: tuple[int, ...]
        match selection_method:
            case "weighted":
                first_binaries = deepcopy(choices(
                    population=list(sorted_results.keys()),
                    weights=list(sorted_results.values()),
                    k=1,
                )[0])
                second_binaries = deepcopy(choices(
                    population=list(sorted_results.keys()),
                    weights=list(sorted_results.values()),
                    k=1,
                )[0])
            case "elite":
                first_binaries = deepcopy(choice(list(sorted_results.keys())[:9]))
                second_binaries = deepcopy(choice(list(sorted_results.keys())[:9]))
            case "random":
                first_binaries = deepcopy(choice(list(sorted_results.keys())))
                second_binaries = deepcopy(choice(list(sorted_results.keys())))
            case _:
                raise ValueError
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
                    first_binaries[flip_location] = tuple(int(not binaries[-1][flip_location]))
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(first_binaries)
            case "random":
                binaries.append(
                    [randint(0, 1) for _ in range(num_reac_couples)]
                )
            case "multimutation":
                num_tries = 0
                while first_binaries in sorted_results:
                    flip_locations: list[int] = [randint(0, num_reac_couples-1) for _ in range(3)]
                    for flip_location in flip_locations:
                        first_binaries[flip_location] = tuple(int(not binaries[-1][flip_location]))
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(first_binaries)
            case "crossover":
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
    ignore_nonlinear_extra_terms_in_lps: bool,
    num_gens: int,
    population_size: int,
    ignore_nonlinear_extra_terms_in_ectfbas: bool,
    max_rounds_same_objvalue: float = float("inf"),
    existing_evolution_result: dict[tuple[int, ...], float] = {},
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    opt_selector = max if is_objsense_maximization(objective_sense) else min
    evolution_results: dict[tuple[int, ...], float | None] = deepcopy(pre_results)

    current_best_objvalue = opt_selector(evolution_results.values())
    num_rounds_with_same_objvalue = 0
    for _ in range(len(num_gens)):
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
                ignore_nonlinear_extra_terms_in_ectfbas,
                tested_binary,
            )
            for tested_binary in tested_binaries
        )
        eligible_binaries: list[tuple[int, ...]] = []
        for (original_binaries, ectfba_result) in ectfba_results:
            if not ectfba_result or not ectfba_result[ALL_OK_KEY]:
                evolution_results[original_binaries] = None
            else:
                eligible_binaries.append(_get_binaries_from_opt_result(ectfba_result))

        results: list[LpNlpBlockResult] = Parallel(n_jobs=-1, verbose=10)(
            delayed(_ectfba_nlp_block)(
                _get_cobrak_model_with_deleted_binary_zero_reacs(
                    cobrak_model=cobrak_model,
                    binaries=eligible_binary,
                    reac_couples_list=reac_couples_list,
                ),
                eligible_binary,
                objective_target,
                objective_sense,
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
            for eligible_binary in eligible_binaries
        )
        for result in results:
            if not result.lp_result:
                evolution_results[result.original_binaries] = None
                continue
            if not result.nlp_result:
                evolution_results[result.original_binaries] = None
                evolution_results[result.lp_binaries] = None
                continue
            evolution_results[result.original_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
            evolution_results[result.lp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]
            evolution_results[result.nlp_binaries] = result.nlp_result[OBJECTIVE_VAR_NAME]

        gen_best_objvalue = opt_selector(evolution_results.values())
        if gen_best_objvalue != current_best_objvalue:
            current_best_objvalue = gen_best_objvalue
            num_rounds_with_same_objvalue = 0
        else:
            num_rounds_with_same_objvalue += 1
        if num_rounds_with_same_objvalue >= max_rounds_same_objvalue:
            break

    return evolution_results


@validate_call(validate_return=True)
def _get_idx_to_reac_ids(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    variability_dict: dict[str, tuple[float, float]],
    error_scenario: dict[str, tuple[float, float]],
) -> tuple[tuple[tuple[str, ...]], dict[str, int]]:
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
        else:
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
        print("ERROR: No working result given or found in sampling")
        raise ValueError

    # PHASE 3: ACTUAL ALGORITHM
