""""""
from cobrak.io import json_write
import random
from ast import literal_eval
from dataclasses import dataclass, field
from pyomo.common.errors import ApplicationError
from cobrak.constants import ALL_OK_KEY, OBJECTIVE_VAR_NAME, Z_VAR_PREFIX
from cobrak.dataclasses import ExtraLinearConstraint, Model, Solver, CorrectionConfig
from cobrak.evolution import is_objsense_maximization
from cobrak.lps import perform_lp_optimization
from cobrak.nlps import perform_nlp_irreversible_optimization
from cobrak.standard_solvers import SCIP, IPOPT
from pydantic import validate_call, PositiveInt, NonNegativeFloat, PositiveFloat
from joblib import Parallel, delayed
from cobrak.utilities import get_stoichiometrically_coupled_reactions, delete_orphaned_metabolites_and_enzymes, delete_unused_reactions_in_optimization_dict
from random import randint, choices, choice, sample, uniform


@dataclass
class LpNlpBlockResult:
    """"""
    original_binaries: tuple[int, ...]
    lp_binaries: tuple[int, ...]
    nlp_binaries: tuple[int, ...]
    lp_result: dict[str, float] = field(default_factory=dict)
    nlp_result: dict[str, float] = field(default_factory=dict)


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
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    lp_solver: Solver,
    correction_config: CorrectionConfig,
    ignore_nonlinear_extra_terms_in_ectfbas: bool,
    binaries: tuple[int, ...],
    reac_couples_list: list[tuple[str, ...]],
    do_reac_deletions: bool = True,
) -> tuple[dict[str, float | None], tuple[int, ...]]:
    with cobrak_model as cobrak_model_with_deletions:
        if do_reac_deletions:
            for couple_idx, binary in enumerate(binaries):
                if binary == 0:
                    for reac_id in reac_couples_list[couple_idx]:
                        del cobrak_model_with_deletions.reactions[reac_id]
        try:
            ectfba_dict: dict[str, int | float] = perform_lp_optimization(
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
    if not ectfba_dict[ALL_OK_KEY] or None in ectfba_dict.values():
        return {}, binaries
    return ectfba_dict, binaries  # ty:ignore[invalid-return-type]


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
) -> dict[str, float | None]:
    try:
        nlp_result: dict[str, float | None] = perform_nlp_irreversible_optimization(
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
    cobrak_model: Model,
    binaries: tuple[int, ...],
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
    delete_nonthermodynamic_reacs_for_nlp: bool,
) -> LpNlpBlockResult:
    with cobrak_model as cobrak_model_with_deletions_and_extra_constraints:
        for couple_idx, binary in enumerate(binaries):
            if binary == 0:
                for reac_id in reac_couples_list[couple_idx]:
                    del cobrak_model_with_deletions_and_extra_constraints.reactions[reac_id]
        cobrak_model_with_deletions_and_extra_constraints.extra_linear_constraints += lp_extra_linear_constraints
        if lp_objective_target == "MAXZ" or lp_objective_target == "MINZ":
            lp_objective_sense = +1 if lp_objective_target == "MAXZ" else -1
            lp_objective_target = {
                f"{Z_VAR_PREFIX}{reac_id}": 1.0
                for (reac_id, reac_data) in cobrak_model_with_deletions_and_extra_constraints.reactions.items()
                if (reac_data.dG0 is not None)
                and (variability_dict[reac_id][1] > 0.0)
            }
        ectfba_dict: dict[str, float | None] = _ectfba_block(
            cobrak_model=cobrak_model_with_deletions_and_extra_constraints,
            objective_target=lp_objective_target,
            objective_sense=lp_objective_sense,
            variability_dict=variability_dict,
            lp_solver=lp_solver,
            correction_config=correction_config,
            ignore_nonlinear_extra_terms_in_ectfbas=ignore_nonlinear_extra_terms_in_ectfbas,
            binaries=binaries,
            reac_couples_list=reac_couples_list,
            do_reac_deletions=False,
        )[0]
    error_target_missing: bool = any(errortarget not in ectfba_dict for errortarget in correction_config.error_scenario)
    if not ectfba_dict or not ectfba_dict[ALL_OK_KEY] or error_target_missing or None in ectfba_dict.values():
        return LpNlpBlockResult(original_binaries=binaries, lp_binaries=(), nlp_binaries=())

    nlp_result: dict[str, float | None] = _nlp_block(
        cobrak_model_with_deletions=delete_unused_reactions_in_optimization_dict(cobrak_model, ectfba_dict, delete_nonthermodynamic_reacs=delete_nonthermodynamic_reacs_for_nlp),
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
    print("++++", lp_objective_sense, nlp_result[OBJECTIVE_VAR_NAME])
    print("++++++")
    return LpNlpBlockResult(
        original_binaries=binaries,
        lp_result=ectfba_dict,
        lp_binaries=_get_binaries_from_opt_result(ectfba_dict, reac_couples_list),
        nlp_result=nlp_result,
        nlp_binaries=_get_binaries_from_opt_result(nlp_result, reac_couples_list)
    )


@validate_call(validate_return=True)
def _add_eligible_binaries_and_get_best_nlp_solution(
    best_nlp_solution: dict[str, float],
    lpnlpblock_results: list[LpNlpBlockResult],
    binary_results: dict[tuple[int, ...], float | None],
    is_maximization: bool,
    add_nlp_result_only: bool,
    min_abs_objvalue: float,
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    for result in lpnlpblock_results:
        if not result.lp_result:
            binary_results[result.original_binaries] = None
            continue
        if not result.nlp_result or (None in result.nlp_result.values()):
            binary_results[result.original_binaries] = None
            binary_results[result.lp_binaries] = None
            continue
        if abs(result.nlp_result[OBJECTIVE_VAR_NAME]) < min_abs_objvalue:
            binary_results[result.original_binaries] = None
            binary_results[result.lp_binaries] = None
            binary_results[result.nlp_binaries] = None
            continue
        keys_to_update = [result.nlp_binaries]
        if not add_nlp_result_only:
            keys_to_update.extend([result.lp_binaries])
        compare = max if is_maximization else min
        nlp_objvalue = result.nlp_result[OBJECTIVE_VAR_NAME]
        for key in keys_to_update:
            if key in binary_results and binary_results[key] is None:
                binary_results[key] = nlp_objvalue
            else:
                binary_results[key] = compare(binary_results[key], nlp_objvalue) if key in binary_results else nlp_objvalue  # ty:ignore[invalid-argument-type]

        if not best_nlp_solution or\
           (is_maximization and result.nlp_result[OBJECTIVE_VAR_NAME] > best_nlp_solution[OBJECTIVE_VAR_NAME]) or\
           (not is_maximization and result.nlp_result[OBJECTIVE_VAR_NAME] < best_nlp_solution[OBJECTIVE_VAR_NAME]):
            best_nlp_solution = result.nlp_result

    return best_nlp_solution, binary_results


@validate_call(validate_return=True)
def _get_binaries_according_to_selection(
    sorted_results: dict[tuple[int, ...], float],
    selection_method: str
) -> tuple[int, ...]:
    keylist: list[tuple[int, ...]] = list(sorted_results.keys())
    match selection_method:
        case "weighted":
            return choices(
                population=list(sorted_results.keys()),
                weights=list(sorted_results.values()),
                k=1,
            )[0]
        case "worst_75_pct":
            return choice(keylist[int(len(keylist) * .25):])
        case "top_25_pct":
            return choice(keylist[:int(len(keylist) * .25)])
        case "top_3":
            return choice(keylist[:6])
        case "random":
            return choice(keylist)
        case _:
            raise ValueError


@validate_call(validate_return=True)
def _get_evolution_binaries(
    population_size: int,
    num_reac_couples: int,
    evolution_results: dict[tuple[int, ...], float | None],
    fractions_genetic_method: dict[str, float],
    fractions_population_selection: dict[str, float],
    sampling_p_random: float,
    sampling_max_knockouts: int,
    sampling_start_solutions: int,
    is_maximization: bool,
    num_rounds_with_same_objvalue: int,
) -> list[tuple[int, ...]]:
    if not evolution_results or len([value for value in evolution_results.values() if value is not None]) < sampling_start_solutions:
        sampling_binaries = []
        for _ in range(population_size):
            if not random.random() < sampling_p_random:
                # Completely random
                sampling_binaries.append(tuple([randint(0, 1) for _ in range(num_reac_couples)]))
            else:
                # max. 5 deletions
                sampling_binary: list[int] = [1] * num_reac_couples
                zero_indices: list[int] = sample(
                    range(num_reac_couples),
                    randint(0, min(num_reac_couples, sampling_max_knockouts))
                )
                for zero_index in zero_indices:
                    sampling_binary[zero_index] = 0
                sampling_binaries.append(tuple(sampling_binary))
        return sampling_binaries

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

        min_change_p = 0.1 * 0.95**num_rounds_with_same_objvalue
        max_change_p = 0.1 * 1.05**num_rounds_with_same_objvalue
        change_p = uniform(min_change_p, max_change_p)
        change_p = max(0.001, change_p)
        change_p = min(0.999, change_p)
        genetic_method: str  = choices(
            population=list(fractions_genetic_method.keys()),
            weights=list(fractions_genetic_method.values()),
            k=1,
        )[0]
        match genetic_method:
            case "extend":
                mutated_x = []
                for x in first_binaries:
                    if x == 1:
                        mutated_x.append(1)
                        continue
                    if uniform(0.0, 1.0) < change_p:
                        mutated_x.append(1)
                    else:
                        mutated_x.append(x)
                binaries.append(tuple(mutated_x))
            case "decrease":
                mutated_x = []
                for x in first_binaries:
                    if x == 0:
                        mutated_x.append(0)
                        continue
                    if uniform(0.0, 1.0) < change_p:
                        mutated_x.append(0)
                    else:
                        mutated_x.append(x)
                binaries.append(tuple(mutated_x))
            case "extend_and_decrease":
                mutated_x = []
                for x in first_binaries:
                    if x == 1:
                        if uniform(0.0, 1.0) < change_p:
                            mutated_x.append(0)
                        else:
                            mutated_x.append(x)
                    else:
                        if uniform(0.0, 1.0) < change_p:
                            mutated_x.append(1)
                        else:
                            mutated_x.append(x)
                binaries.append(tuple(mutated_x))
            case "neighborhood":
                num_tries = 0
                while first_binaries in sorted_results:
                    flip_location: int = randint(0, num_reac_couples-1)
                    first_binaries: tuple[int, ...] = tuple(
                        list(first_binaries[:flip_location]) + [int(not first_binaries[flip_location])] + list(first_binaries[flip_location+1:])
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
                        first_binaries = tuple(
                            list(first_binaries[:flip_location]) + [int(not first_binaries[flip_location])] + list(first_binaries[flip_location+1:])
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
                crossed_over_binaries: tuple[int, ...]
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
    reac_couples_list: tuple[tuple[str, ...], ...],
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    evolution_results: dict[tuple[int, ...], float | None],
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
    fractions_genetic_method: dict[str, float],
    fractions_population_selection: dict[str, float],
    sampling_p_random: float,
    sampling_max_knockouts: int,
    sampling_start_solutions: int,
    min_abs_objvalue: float,
    inner_lp_objectives: tuple[str, ...],
    max_rounds_same_objvalue: int,
    verbose: bool,
    round_result_json_path: str,
    do_sampling_only: bool,
) -> tuple[dict[str, float], dict[tuple[int, ...], float | None]]:
    opt_selector = max if is_objsense_maximization(objective_sense) else min
    if type(objective_target) is str:
        objective_target_as_dict: dict[str, int | float] = {objective_target: 1.0}
    elif type(objective_target) is dict:
        objective_target_as_dict: dict[str, int | float] = objective_target

    best_nlp_solution: dict[str, float] = {}
    if evolution_results:
        current_best_objvalue: float = opt_selector([value for value in evolution_results.values() if value is not None])
    else:
        current_best_objvalue: float = -float("inf") if is_objsense_maximization(objective_sense) else float("inf")
    num_rounds_with_same_objvalue = 0
    for current_round in range(num_gens):
        if do_sampling_only and len([value for value in evolution_results.values() if value is not None]) >= sampling_start_solutions:
            print("ENDING SAMPLING (do_sampling_only argument is set to True)")
            break
        tested_binaries = _get_evolution_binaries(
            population_size=population_size,
            num_reac_couples=len(reac_couples_list),
            evolution_results=evolution_results,
            fractions_genetic_method=fractions_genetic_method,
            fractions_population_selection=fractions_population_selection,
            sampling_p_random=sampling_p_random,
            sampling_max_knockouts=sampling_max_knockouts,
            sampling_start_solutions=sampling_start_solutions,
            is_maximization=is_objsense_maximization(objective_sense),
            num_rounds_with_same_objvalue=num_rounds_with_same_objvalue,
        )

        ectfba_results: list[dict[str, float]] = Parallel(n_jobs=-1, verbose=0)(
            delayed(_ectfba_block)(
                cobrak_model,
                objective_target,
                objective_sense,
                variability_dict,
                lp_solver,
                correction_config,
                ignore_nonlinear_extra_terms_in_lps,
                tested_binary,
                reac_couples_list,
            )
            for tested_binary in tested_binaries
        )
        eligible_binaries_with_objvalue: dict[tuple[int, ...], float] = {}
        for (ectfba_result, original_binaries) in ectfba_results:
            if not ectfba_result or not ectfba_result[ALL_OK_KEY]:
                evolution_results[original_binaries] = None
            elif original_binaries not in evolution_results:
                eligible_binaries_with_objvalue[original_binaries] = ectfba_result[OBJECTIVE_VAR_NAME]
        results: list[LpNlpBlockResult] = Parallel(n_jobs=-1, verbose=0)(
            delayed(_ectfba_nlp_block)(
                cobrak_model,
                eligible_binary,
                reac_couples_list,
                inner_lp_objective,
                +1,
                [
                    ExtraLinearConstraint(
                        stoichiometries=objective_target_as_dict,
                        lower_value=objvalue - 1e-8,
                        upper_value=objvalue + 1e-8,
                    )
                ],
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
                False,
            )
            for eligible_binary, objvalue in eligible_binaries_with_objvalue.items()
            for inner_lp_objective in inner_lp_objectives
        )
        best_nlp_solution, evolution_results = _add_eligible_binaries_and_get_best_nlp_solution(
            best_nlp_solution=best_nlp_solution,
            lpnlpblock_results=results,
            binary_results=evolution_results,
            is_maximization=is_objsense_maximization(objective_sense),
            add_nlp_result_only=False,
            min_abs_objvalue=min_abs_objvalue,
        )
        if verbose:
            print(f"ROUND {current_round} NON-NONE OBJECTIVE VALUES: {sorted(x for x in set(evolution_results.values()) if x is not None)}")
        if round_result_json_path:
            json_write(round_result_json_path, [best_nlp_solution, {str(key): value for key, value in evolution_results.items()}])

        non_none_objvalues = [value for value in evolution_results.values() if value is not None]
        if len(non_none_objvalues) > 0:
            gen_best_objvalue = opt_selector(non_none_objvalues)
            if gen_best_objvalue != current_best_objvalue:
                current_best_objvalue = gen_best_objvalue
                num_rounds_with_same_objvalue = 0
            else:
                num_rounds_with_same_objvalue += 1
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
            reac_id in error_scenario
            for reac_id in reac_ids
        ):
            continue
        if all(not cobrak_model.reactions[reac_id].dG0 and not cobrak_model.reactions[reac_id].enzyme_reaction_data for reac_id in reac_ids):
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
    ignore_nonlinear_extra_terms_in_lps: bool = True,
    existing_evolution_results: dict[str, float | None] = {},
    fractions_genetic_method: dict[str, NonNegativeFloat] = {
        # "neighborhood": 0.0,
        # "random": 1/2,
        # "multimutation": 1/2,
        "crossover": 1/4,
        "extend": 1/4,
        "decrease": 1/4,
        "extend_and_decrease": 1/4,
    },
    fractions_population_selection: dict[str, NonNegativeFloat] = {
        # "weighted": 1/4,
        # "random": 1/4,
        "top_3": 1/4,
        "top_25_pct": 1/2,
        "worst_75_pct": 1/4,
    },
    sampling_p_random: NonNegativeFloat = 0.33,
    sampling_max_knockouts: PositiveInt = 5,
    sampling_start_solutions: PositiveInt = 2,
    min_abs_objvalue: PositiveFloat = 1e-8,
    inner_lp_objectives: tuple[str, ...] = ("MAXZ",),
    max_rounds_same_objvalue: PositiveInt = 1_000_000,
    verbose: bool = False,
    round_result_json_path: str = "",
    do_sampling_only: bool = False,
) -> tuple[dict[str, float], dict[str, float | None]]:
    """"""
    # PHASE 1: BUILD INDEX TO REAC COUPLES DATA, AND CPU DATA
    reac_couples_list: tuple[tuple[str, ...], ...] = _get_idx_to_reac_ids(
        cobrak_model=cobrak_model,
        objective_target=objective_target,
        variability_dict=variability_dict,
        error_scenario=correction_config.error_scenario,
    )

    # PHASE 2: ACTUAL EVOLUTION ALGORITHM (USING GIVEN OR SAMPLED RESULTS AS STARTING POINTS)
    best_result, binary_results = _evolution(
        cobrak_model=cobrak_model,
        reac_couples_list=reac_couples_list,
        objective_target=objective_target,
        objective_sense=objective_sense,
        variability_dict=variability_dict,
        evolution_results={
            literal_eval(key): value
            for key, value in existing_evolution_results.items()
        },
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
        fractions_genetic_method=fractions_genetic_method,
        fractions_population_selection=fractions_population_selection,
        sampling_p_random=sampling_p_random,
        sampling_max_knockouts=sampling_max_knockouts,
        sampling_start_solutions=sampling_start_solutions,
        min_abs_objvalue=min_abs_objvalue,
        inner_lp_objectives=inner_lp_objectives,
        max_rounds_same_objvalue=max_rounds_same_objvalue,
        verbose=verbose,
        round_result_json_path=round_result_json_path,
        do_sampling_only=do_sampling_only,
    )

    return best_result, {str(key): value for key, value in binary_results.items()}
