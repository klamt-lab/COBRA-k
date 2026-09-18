"""Experimental! Not for usage currently"""

import random
from ast import literal_eval
from copy import deepcopy
from math import ceil, floor
from random import choice, choices, randint, sample, uniform

from joblib import Parallel, delayed, parallel_config
from pydantic import NonNegativeFloat, PositiveFloat, PositiveInt
from pydantic.dataclasses import Field, dataclass
from pyomo.common.errors import ApplicationError

from cobrak.constants import ALL_OK_KEY, BIG_M, OBJECTIVE_VAR_NAME, Z_VAR_PREFIX
from cobrak.dataclasses import CorrectionConfig, ExtraLinearConstraint, Model, Solver
from cobrak.evolution import is_objsense_maximization
from cobrak.io import json_write
from cobrak.lps import perform_lp_optimization
from cobrak.nlps import (
    perform_nlp_irreversible_optimization,
    perform_nlp_irreversible_optimization_with_active_reacs_only,
)
from cobrak.utilities import (
    delete_orphaned_metabolites_and_enzymes,
    get_stoichiometrically_coupled_reactions,
)


# DATACLASSES #
@dataclass
class EvolutionSettings:
    cobrak_model: Model
    objective_target: str | dict[str, float]
    objective_sense: int
    variability_dict: dict[str, tuple[float, float]]
    num_gens: int
    population_size: int
    with_kappa: bool = True
    with_gamma: bool = True
    with_iota: bool = False
    with_alpha: bool = False
    correction_config: CorrectionConfig = Field(default_factory=CorrectionConfig)
    lp_solver: Solver = Field(default_factory=Solver)
    nlp_solver: Solver = Field(default_factory=Solver)
    nlp_strict_mode: bool = False
    nlp_single_strict_reacs: list[str] = Field(default_factory=list)
    ignore_nonlinear_extra_terms_in_lps: bool = True
    fractions_genetic_method: dict[str, NonNegativeFloat] = Field(
        default_factory=lambda: {
            # "neighborhood": 0.0,
            # "random": 1/2,
            # "multimutation": 1/2,
            "crossover": 1 / 4,
            "extend": 1 / 4,
            "decrease": 1 / 4,
            "extend_and_decrease": 1 / 4,
        }
    )
    fractions_population_selection: dict[str, NonNegativeFloat] = Field(
        default_factory=lambda: {
            # "weighted": 1/4,
            # "random": 1/4,
            "top_3": 1 / 4,
            "top_25_pct": 1 / 2,
            "worst_75_pct": 1 / 4,
        }
    )
    use_original_ectfba_nlp_binaries_as_template: tuple[bool, bool, bool] = Field(
        default_factory=lambda: (True, False, False)
    )
    sampling_p_random: NonNegativeFloat = 0.33
    sampling_max_knockouts: PositiveInt = 5
    sampling_start_solutions: PositiveInt = 2
    min_abs_objvalue: PositiveFloat = 1e-8
    inner_lp_objectives: tuple[str, ...] = ("MAXZ",)
    max_rounds_same_objvalue: PositiveInt = 1_000_000
    verbose: bool = False
    round_result_json_path: str = ""
    do_sampling_only: bool = False
    num_used_cpu_cores: float = -1
    approximation_value: float = 0.0001


@dataclass
class EvolutionResult:
    best_nlps: list[dict[str, float]] = Field(default_factory=list)
    original_binary_results: dict[str, list[float | None]] = Field(default_factory=dict)
    ectfba_binary_results: dict[str, float | None] = Field(default_factory=dict)
    nlp_binary_results: dict[str, float | None] = Field(default_factory=dict)
    reac_couples_list: tuple[tuple[str, ...], ...] = Field(default_factory=tuple)


@dataclass
class _InternalEvolutionResult:
    best_nlps: list[dict[str, float]] = Field(default_factory=list)
    original_binary_results: dict[tuple[int, ...], list[float | None]] = Field(
        default_factory=dict
    )
    ectfba_binary_results: dict[tuple[int, ...], float | None] = Field(
        default_factory=dict
    )
    nlp_binary_results: dict[tuple[int, ...], float | None] = Field(
        default_factory=dict
    )
    reac_couples_list: tuple[tuple[str, ...], ...] = Field(default_factory=tuple)


@dataclass
class LpNlpBlockResult:
    """"""

    original_binaries: tuple[int, ...]
    lp_binaries: tuple[int, ...]
    nlp_binaries: tuple[int, ...]
    lp_result: dict[str, float] | None = Field(default_factory=dict)
    nlp_result: dict[str, float] | None = Field(default_factory=dict)


@dataclass
class PostprocessingSetting:
    target_type: str
    target_reac: str
    lp_objective: str
    num_active_deviation: int
    num_inactive_deviation: int
    extra_zvars: tuple[str, ...]
    extra_linear_constraints: list[ExtraLinearConstraint]


@dataclass
class PostprocessingOutcome:
    target_reac: str
    target_type: str
    lp_objective: str
    num_active_deviation: int
    num_inactive_deviation: int
    ectfba_result: float | None
    nlp_result: float | None


# DATACLASS CONVERTORS #
def _internal_to_external_evolution_result(
    internal_evolution_result: _InternalEvolutionResult,
) -> EvolutionResult:
    return EvolutionResult(
        best_nlps=list(
            {
                frozenset(d.items()): d for d in internal_evolution_result.best_nlps
            }.values()
        ),  # return only unique solutions
        original_binary_results={
            str(key): value
            for key, value in internal_evolution_result.original_binary_results.items()
        },
        ectfba_binary_results={
            str(key): value
            for key, value in internal_evolution_result.ectfba_binary_results.items()
        },
        nlp_binary_results={
            str(key): value
            for key, value in internal_evolution_result.nlp_binary_results.items()
        },
        reac_couples_list=internal_evolution_result.reac_couples_list,
    )


def _external_to_internal_evolution_result(
    external_evolution_result: EvolutionResult,
) -> _InternalEvolutionResult:
    return _InternalEvolutionResult(
        best_nlps=external_evolution_result.best_nlps,
        original_binary_results={
            literal_eval(key): value
            for key, value in external_evolution_result.original_binary_results.items()
        },
        ectfba_binary_results={
            literal_eval(key): value
            for key, value in external_evolution_result.ectfba_binary_results.items()
        },
        nlp_binary_results={
            literal_eval(key): value
            for key, value in external_evolution_result.nlp_binary_results.items()
        },
        reac_couples_list=external_evolution_result.reac_couples_list,
    )


# PUBLIC FUNCTIONS #
# @validate_call(validate_return=True)
def delete_unused_reactions_in_optimization_dict_2(
    cobrak_model: Model,
    optimization_dict: dict[str, float],
    exception_prefix: str = "",
    delete_missing_reactions: bool = True,
    min_abs_flux: NonNegativeFloat = 1e-15,
    do_not_delete_with_z_var_one: bool = True,
    delete_nonthermodynamic_reacs: bool = True,
) -> Model:
    """Delete unused reactions in a COBRAk model based on an optimization dictionary.

    This function creates a deep copy of the provided COBRAk model and removes reactions that are either not present
    in the optimization dictionary or have flux values below a specified threshold. Optionally,
    reactions with a specific prefix can be excluded from deletion.
    Additionally, orphaned metabolites (those not used in any remaining reactions) are also removed.

    Args:
        cobrak_model (Model): COBRAk model containing reactions and metabolites.
        optimization_dict (dict[str, float]): Dictionary mapping reaction IDs to their optimized flux values.
        exception_prefix (str, optional): A prefix for reaction IDs that should not be deleted. Defaults to "".
        delete_missing_reactions (bool, optional): Whether to delete reactions not present in the optimization dictionary. Defaults to True.
        min_abs_flux (float, optional): The minimum absolute flux value below which reactions are considered unused. Defaults to 1e-10.

    Returns:
        Model: A new COBRAk model with unused reactions and orphaned metabolites removed.
    """
    cobrak_model = deepcopy(cobrak_model)
    reacs_to_delete: list[str] = []
    for reac_id in cobrak_model.reactions:
        to_delete = False
        if (reac_id not in optimization_dict) and delete_missing_reactions:
            to_delete = True
        elif (reac_id in optimization_dict) and abs(
            optimization_dict[reac_id]
        ) <= min_abs_flux:
            z_var_id = f"{Z_VAR_PREFIX}{reac_id}"
            if z_var_id in optimization_dict:
                if do_not_delete_with_z_var_one and (
                    optimization_dict[z_var_id] <= 1e-6
                ):
                    to_delete = True
                else:
                    to_delete = False
            else:
                if (
                    not delete_nonthermodynamic_reacs
                    and cobrak_model.reactions[reac_id].dG0 is None
                ):
                    to_delete = False
                else:
                    to_delete = True
        if to_delete:
            reacs_to_delete.append(reac_id)
    for reac_to_delete in reacs_to_delete:
        if (exception_prefix) and (reac_to_delete.startswith(exception_prefix)):
            continue
        del cobrak_model.reactions[reac_to_delete]
    return delete_orphaned_metabolites_and_enzymes(cobrak_model)


######## PRIVATE FUNCTIONS ########
# @validate_call(validate_return=True)
def _get_binaries_from_opt_result(
    opt_result: dict[str, float],
    reac_couples_list: list[tuple[str, ...]] | tuple[tuple[str, ...], ...],
) -> tuple[int, ...]:
    binaries = [0 for _ in range(len(reac_couples_list))]
    for i, reac_couple in enumerate(reac_couples_list):
        if reac_couple[0] in opt_result:
            binaries[i] = 1
    return tuple(binaries)


# @validate_call(validate_return=True)
def _ectfba_block(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    lp_solver: Solver,
    correction_config: CorrectionConfig,
    ignore_nonlinear_extra_terms_in_ectfbas: bool,
    binaries: tuple[int, ...],
    reac_couples_list: tuple[tuple[str, ...], ...],
    do_reac_deletions: bool = True,
    verbose: bool = False,
) -> tuple[dict[str, float], tuple[int, ...]]:
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
                verbose=verbose,
            )
        except (ApplicationError, AttributeError, ValueError):
            return {}, binaries
    if not ectfba_dict[ALL_OK_KEY] or None in ectfba_dict.values():
        return {}, binaries
    return ectfba_dict, binaries


# @validate_call(validate_return=True)
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
    approximation_value: float,
    verbose: bool = False,
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
            verbose=verbose,
            approximation_value=approximation_value,
        )
    except (ApplicationError, AttributeError, ValueError):
        return {}
    if not nlp_result[ALL_OK_KEY]:
        return {}
    if None in nlp_result.values():
        return {}

    return nlp_result


# @validate_call(validate_return=True)
def _ectfba_nlp_block(
    cobrak_model: Model,
    binaries: tuple[int, ...],
    reac_couples_list: tuple[tuple[str, ...], ...],
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
    approximation_value: float,
    verbose: bool = False,
    lp_extra_binary_vars: list[str] = [],
) -> LpNlpBlockResult:
    with cobrak_model as cobrak_model_with_deletions_and_extra_constraints:
        for couple_idx, binary in enumerate(binaries):
            if binary == 0:
                for reac_id in reac_couples_list[couple_idx]:
                    del cobrak_model_with_deletions_and_extra_constraints.reactions[
                        reac_id
                    ]
        cobrak_model_with_deletions_and_extra_constraints.extra_linear_constraints += (
            lp_extra_linear_constraints
        )
        cobrak_model_with_deletions_and_extra_constraints.extra_binary_vars += (
            lp_extra_binary_vars
        )
        if lp_objective_target == "MAXZ" or lp_objective_target == "MINZ":
            lp_objective_sense = +1 if lp_objective_target == "MAXZ" else -1
            lp_objective_target = {
                f"{Z_VAR_PREFIX}{reac_id}": 1.0
                for (
                    reac_id,
                    reac_data,
                ) in cobrak_model_with_deletions_and_extra_constraints.reactions.items()
                if (reac_data.dG0 is not None) and (variability_dict[reac_id][1] > 0.0)
            } | {
                binvar: 1.0
                for binvar in lp_extra_binary_vars
                if binvar.startswith("PPBIN_")
            }
        ectfba_dict: dict[str, float] = _ectfba_block(
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
            verbose=verbose,
        )[0]

    error_target_missing: bool = any(
        errortarget not in ectfba_dict
        for errortarget in correction_config.error_scenario
    )
    if (
        not ectfba_dict
        or not ectfba_dict[ALL_OK_KEY]
        or error_target_missing
        or None in ectfba_dict.values()
    ):
        return LpNlpBlockResult(
            original_binaries=binaries, lp_binaries=(), nlp_binaries=()
        )

    nlp_result: dict[str, float] = _nlp_block(
        cobrak_model_with_deletions=delete_unused_reactions_in_optimization_dict_2(
            cobrak_model,
            ectfba_dict,
            delete_nonthermodynamic_reacs=delete_nonthermodynamic_reacs_for_nlp,
        ),
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
        approximation_value=approximation_value,
        verbose=verbose,
    )
    if not nlp_result or not nlp_result[ALL_OK_KEY] or None in nlp_result.values():
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
        nlp_binaries=_get_binaries_from_opt_result(nlp_result, reac_couples_list),
    )


# @validate_call(validate_return=True)
def _add_eligible_binaries_and_get_best_nlp_solution(
    evolution_result: _InternalEvolutionResult,
    lpnlpblock_results: list[LpNlpBlockResult],
    is_maximization: bool,
    min_abs_objvalue: float,
) -> _InternalEvolutionResult:
    for result in lpnlpblock_results:
        if not result.lp_result:
            if result.original_binaries not in evolution_result.original_binary_results:
                evolution_result.original_binary_results.setdefault(
                    result.original_binaries, []
                ).append(None)
            continue
        if not result.nlp_result or (None in result.nlp_result.values()):
            if result.original_binaries not in evolution_result.original_binary_results:
                evolution_result.original_binary_results.setdefault(
                    result.original_binaries, []
                ).append(None)
            evolution_result.ectfba_binary_results[result.lp_binaries] = None
            continue
        if abs(result.nlp_result[OBJECTIVE_VAR_NAME]) < min_abs_objvalue:
            if result.original_binaries not in evolution_result.original_binary_results:
                evolution_result.original_binary_results.setdefault(
                    result.original_binaries, []
                ).append(None)
            evolution_result.ectfba_binary_results[result.lp_binaries] = None
            evolution_result.nlp_binary_results[result.nlp_binaries] = None
            continue

        # All eligible ⇒ ecTFBA and NLP results exist :D
        nlp_objvalue = result.nlp_result[OBJECTIVE_VAR_NAME]
        evolution_result.ectfba_binary_results[result.lp_binaries] = nlp_objvalue
        evolution_result.nlp_binary_results[result.nlp_binaries] = nlp_objvalue
        evolution_result.original_binary_results.setdefault(
            result.original_binaries, []
        ).append(nlp_objvalue)

        no_nlp_solution = not evolution_result.best_nlps
        better_nlp_solution = (not no_nlp_solution) and (
            (
                is_maximization
                and result.nlp_result[OBJECTIVE_VAR_NAME]
                > evolution_result.best_nlps[0][OBJECTIVE_VAR_NAME]
            )
            or (
                not is_maximization
                and result.nlp_result[OBJECTIVE_VAR_NAME]
                < evolution_result.best_nlps[0][OBJECTIVE_VAR_NAME]
            )
        )
        if no_nlp_solution or better_nlp_solution:
            evolution_result.best_nlps = [deepcopy(result.nlp_result)]
        elif not no_nlp_solution and (
            result.nlp_result[OBJECTIVE_VAR_NAME]
            == evolution_result.best_nlps[0][OBJECTIVE_VAR_NAME]
        ):
            evolution_result.best_nlps.append(deepcopy(result.nlp_result))

    return evolution_result


# @validate_call(validate_return=True)
def _get_binaries_according_to_selection(
    sorted_results: dict[tuple[int, ...], float], selection_method: str
) -> tuple[int, ...]:
    keylist: list[tuple[int, ...]] = list(sorted_results.keys())
    if selection_method == "weighted":
        return choices(
            population=list(sorted_results.keys()),
            weights=list(sorted_results.values()),
            k=1,
        )[0]
    if selection_method.startswith("worst_") and selection_method.endswith("_pct"):
        cleaned_str_float = (
            float(selection_method.replace("worst_", "").replace("_pct", "")) / 100
        )
        return choice(keylist[floor(len(keylist) * cleaned_str_float) :])
    if selection_method.startswith("top_") and selection_method.endswith("_pct"):
        cleaned_str_float = (
            float(selection_method.replace("top_", "").replace("_pct", "")) / 100
        )
        return choice(keylist[: ceil(len(keylist) * cleaned_str_float)])
    if selection_method.startswith("top_") and not selection_method.endswith("_pct"):
        cleaned_str_int = int(selection_method.replace("top_", ""))
        return choice(keylist[:cleaned_str_int])
    if selection_method == "random":
        return choice(keylist)
    print(f"ERROR: Invalid selection method {selection_method}")
    raise ValueError


# @validate_call(validate_return=True)
def _get_binaries_to_test(
    population_size: int,
    num_reac_couples: int,
    evolution_results: _InternalEvolutionResult,
    fractions_genetic_method: dict[str, float],
    fractions_population_selection: dict[str, float],
    sampling_p_random: float,
    sampling_max_knockouts: int,
    sampling_start_solutions: int,
    is_maximization: bool,
    num_rounds_with_same_objvalue: int,
    use_original_ectfba_nlp_binaries_as_template: tuple[bool, bool, bool],
) -> list[tuple[int, ...]]:
    # SAMPLING (IF ELIGIBLE)
    if (
        not evolution_results
        or len(
            [
                value
                for valuelist in evolution_results.original_binary_results.values()
                for value in valuelist
                if value is not None
            ]
        )
        < sampling_start_solutions
    ):
        sampling_binaries = []
        for _ in range(population_size):
            if not random.random() < sampling_p_random:
                # Completely random
                sampling_binaries.append(
                    tuple([randint(0, 1) for _ in range(num_reac_couples)])
                )
            else:
                # max. 5 deletions
                sampling_binary: list[int] = [1] * num_reac_couples
                zero_indices: list[int] = sample(
                    range(num_reac_couples),
                    randint(0, min(num_reac_couples, sampling_max_knockouts)),
                )
                for zero_index in zero_indices:
                    sampling_binary[zero_index] = 0
                sampling_binaries.append(tuple(sampling_binary))
        return sampling_binaries

    # NO SAMPLING ⇒ NORMAL EVOLUTION MUTATIONS ;D
    comparison_function = max if is_maximization else min
    non_na_evolution_results: dict[tuple[int, ...], float] = {}
    if use_original_ectfba_nlp_binaries_as_template[0]:
        for binaries, values in evolution_results.original_binary_results.items():
            nonna_values = [value for value in values if value is not None]
            if nonna_values:
                non_na_evolution_results[binaries] = comparison_function(nonna_values)
    if use_original_ectfba_nlp_binaries_as_template[1]:
        for binaries, value in evolution_results.ectfba_binary_results.items():
            if value:
                non_na_evolution_results[binaries] = value
    if use_original_ectfba_nlp_binaries_as_template[2]:
        for binaries, value in evolution_results.nlp_binary_results.items():
            if value:
                non_na_evolution_results[binaries] = value
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
        genetic_method: str = choices(
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
                    flip_location: int = randint(0, num_reac_couples - 1)
                    first_binaries: tuple[int, ...] = tuple(
                        list(first_binaries[:flip_location])
                        + [int(not first_binaries[flip_location])]
                        + list(first_binaries[flip_location + 1 :])
                    )
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(first_binaries)
            case "random":
                binaries.append(tuple([randint(0, 1) for _ in range(num_reac_couples)]))
            case "multimutation":
                num_tries = 0
                while first_binaries in sorted_results:
                    flip_locations: list[int] = [
                        randint(0, num_reac_couples - 1) for _ in range(3)
                    ]
                    for flip_location in flip_locations:
                        first_binaries = tuple(
                            list(first_binaries[:flip_location])
                            + [int(not first_binaries[flip_location])]
                            + list(first_binaries[flip_location + 1 :])
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
                    crossover_point = randint(0, num_reac_couples - 1)
                    crossed_over_binaries = (
                        first_binaries[:crossover_point]
                        + second_binaries[crossover_point:]
                    )
                    num_tries += 1
                    if num_tries == 100:
                        break
                if num_tries < 100:
                    binaries.append(crossed_over_binaries)
            case _:
                raise ValueError
    return binaries


# @validate_call(validate_return=True)
def _evolution(
    settings: EvolutionSettings,
    existing_result: _InternalEvolutionResult,
) -> _InternalEvolutionResult:
    existing_result = deepcopy(existing_result)

    comparison_function = (
        max if is_objsense_maximization(settings.objective_sense) else min
    )
    if type(settings.objective_target) is str:
        objective_target_as_dict: dict[str, int | float] = {
            settings.objective_target: 1.0
        }
    elif type(settings.objective_target) is dict:
        objective_target_as_dict: dict[str, int | float] = settings.objective_target

    if existing_result.original_binary_results:
        current_best_objvalue: float = comparison_function(
            [
                comparison_function([x for x in value if x is not None])
                for value in existing_result.original_binary_results.values()
                if value
            ]
        )
    else:
        current_best_objvalue: float = (
            -float("inf")
            if is_objsense_maximization(settings.objective_sense)
            else float("inf")
        )

    num_rounds_with_same_objvalue = 0
    for current_round in range(settings.num_gens):
        if (
            settings.do_sampling_only
            and len(
                [
                    value
                    for value in existing_result.original_binary_results.values()
                    if value is not None
                ]
            )
            >= settings.sampling_start_solutions
        ):
            print("ENDING AFTER SAMPLING (do_sampling_only argument is set to True)")
            break
        tested_binaries = _get_binaries_to_test(
            population_size=settings.population_size,
            num_reac_couples=len(existing_result.reac_couples_list),
            evolution_results=existing_result,
            fractions_genetic_method=settings.fractions_genetic_method,
            fractions_population_selection=settings.fractions_population_selection,
            sampling_p_random=settings.sampling_p_random,
            sampling_max_knockouts=settings.sampling_max_knockouts,
            sampling_start_solutions=settings.sampling_start_solutions,
            is_maximization=is_objsense_maximization(settings.objective_sense),
            num_rounds_with_same_objvalue=num_rounds_with_same_objvalue,
            use_original_ectfba_nlp_binaries_as_template=settings.use_original_ectfba_nlp_binaries_as_template,
        )
        ectfba_results: list[tuple[dict[str, float], tuple[int, ...]]] = Parallel(
            n_jobs=settings.num_used_cpu_cores, verbose=0
        )(
            delayed(_ectfba_block)(
                settings.cobrak_model,
                settings.objective_target,
                settings.objective_sense,
                settings.variability_dict,
                settings.lp_solver,
                settings.correction_config,
                settings.ignore_nonlinear_extra_terms_in_lps,
                tested_binary,
                existing_result.reac_couples_list,
            )
            for tested_binary in tested_binaries
        )
        eligible_binaries_with_objvalue: dict[tuple[int, ...], float] = {}
        for ectfba_result, original_binaries in ectfba_results:
            if not ectfba_result or not ectfba_result.get(ALL_OK_KEY, False):
                existing_result.original_binary_results.setdefault(
                    original_binaries, []
                ).append(None)
            elif original_binaries not in existing_result.original_binary_results:
                eligible_binaries_with_objvalue[original_binaries] = ectfba_result[
                    OBJECTIVE_VAR_NAME
                ]
        with parallel_config(backend="loky", inner_max_num_threads=1):
            results: list[LpNlpBlockResult] = Parallel(
                n_jobs=settings.num_used_cpu_cores, verbose=0
            )(
                delayed(_ectfba_nlp_block)(
                    settings.cobrak_model,
                    eligible_binary,
                    existing_result.reac_couples_list,
                    inner_lp_objective,
                    +1,
                    [
                        ExtraLinearConstraint(
                            stoichiometries=objective_target_as_dict,
                            lower_value=objvalue - 1e-8,
                            upper_value=objvalue + 1e-8,
                        )
                    ],
                    settings.objective_target,
                    settings.objective_sense,
                    settings.variability_dict,
                    settings.with_kappa,
                    settings.with_gamma,
                    settings.with_iota,
                    settings.with_alpha,
                    settings.lp_solver,
                    settings.nlp_solver,
                    settings.nlp_strict_mode,
                    settings.nlp_single_strict_reacs,
                    settings.correction_config,
                    settings.ignore_nonlinear_extra_terms_in_lps,
                    False,
                    settings.approximation_value,
                )
                for eligible_binary, objvalue in eligible_binaries_with_objvalue.items()
                for inner_lp_objective in settings.inner_lp_objectives
            )
        existing_result = _add_eligible_binaries_and_get_best_nlp_solution(
            evolution_result=existing_result,
            lpnlpblock_results=results,
            is_maximization=is_objsense_maximization(settings.objective_sense),
            min_abs_objvalue=settings.min_abs_objvalue,
        )
        if settings.verbose:
            print(
                f"ROUND {current_round} OBJECTIVE VALUES: {existing_result.original_binary_results.values()}"
            )
        if settings.round_result_json_path:
            json_write(settings.round_result_json_path, existing_result)

        non_none_objvalues = [
            value
            for valuelist in existing_result.original_binary_results.values()
            for value in valuelist
            if value is not None
        ]
        if non_none_objvalues:
            gen_best_objvalue = comparison_function(non_none_objvalues)
            if gen_best_objvalue != current_best_objvalue:
                current_best_objvalue = gen_best_objvalue
                num_rounds_with_same_objvalue = 0
            else:
                num_rounds_with_same_objvalue += 1
        else:
            num_rounds_with_same_objvalue += 1
        if num_rounds_with_same_objvalue >= settings.max_rounds_same_objvalue:
            break

    return existing_result


##########################
def _get_postprocessing_calculation_settings(
    cobrak_model: Model,
    best_binary: tuple[int, ...],
    best_result: dict[str, float],
    reac_couples_list: tuple[tuple[str, ...], ...],
    num_active_deviations: tuple[int, ...],
    num_inactive_deviations: tuple[int, ...],
    lp_targets: tuple[str, ...],
    objective_target: dict[str, float],
    objective_sense: int,
    target_reaction_selection: list[str],
) -> list[PostprocessingSetting]:
    first_reacs_in_couple = tuple([reacs[0] for reacs in reac_couples_list])
    extra_binary_var_prefix = "PPBIN_"
    extra_binary_vars = {
        reac: f"{extra_binary_var_prefix}{reac}"
        for reac in first_reacs_in_couple
        if (cobrak_model.reactions[reac].dG0 is None)
    }

    if not target_reaction_selection:
        target_reacs = ["BASETEST"] + list(first_reacs_in_couple)
    else:
        target_reacs = target_reaction_selection
    postprocessing_calculation_settings: list[
        PostprocessingSetting
    ] = []  # list[tuple[str, str, int, int, tuple[str, ...]
    for num_active_deviation in num_active_deviations:
        for num_inactive_deviation in num_inactive_deviations:
            for target_reac in target_reacs:
                target_type = ""
                extra_inactive_zvars: list[str] = []
                existing_inactive_zvars: list[str] = []
                extra_active_zvars: list[str] = []
                existing_active_zvars: list[str] = []
                extra_linear_constraints: list[ExtraLinearConstraint] = []
                for i, zvalue in enumerate(best_binary):
                    reacname = first_reacs_in_couple[i]
                    if target_reac == reacname and reacname != "BASETEST":
                        if zvalue == 0:  # or best_result.get(reacname, 0.0) < 1e-8:
                            target_type = "ACT."
                            extra_linear_constraints.append(
                                ExtraLinearConstraint(
                                    stoichiometries={
                                        first_reacs_in_couple[i]: 1.0,
                                    },
                                    lower_value=1e-7,
                                    upper_value=None,
                                )
                            )
                        else:
                            target_type = "INACT."
                            extra_linear_constraints.append(
                                ExtraLinearConstraint(
                                    stoichiometries={
                                        first_reacs_in_couple[i]: 1.0,
                                    },
                                    lower_value=None,
                                    upper_value=0.0,
                                )
                            )
                        continue
                    if cobrak_model.reactions[first_reacs_in_couple[i]].dG0 is None:
                        binary_varname = (
                            extra_binary_vars[first_reacs_in_couple[i]]
                            if (
                                cobrak_model.reactions[first_reacs_in_couple[i]].dG0
                                is None
                            )
                            else f"{Z_VAR_PREFIX}{first_reacs_in_couple[i]}"
                        )
                        if zvalue == 0:  # or best_result.get(reacname, 0.0) < 1e-12:
                            extra_inactive_zvars.append(binary_varname)
                        else:
                            extra_active_zvars.append(binary_varname)
                        extra_linear_constraints.append(
                            ExtraLinearConstraint(
                                stoichiometries={
                                    first_reacs_in_couple[i]: 1.0,
                                    binary_varname: -BIG_M,
                                },
                                lower_value=None,
                                upper_value=0.0,
                            )
                        )
                    else:
                        z_varname = f"{Z_VAR_PREFIX}{first_reacs_in_couple[i]}"
                        if zvalue == 0:  # or best_result.get(reacname, 0.0) < 1e-12:
                            existing_inactive_zvars.append(z_varname)
                        else:
                            existing_active_zvars.append(z_varname)
                extra_linear_constraints.append(
                    ExtraLinearConstraint(
                        stoichiometries=dict.fromkeys(
                            extra_inactive_zvars + existing_inactive_zvars, 1.0
                        ),
                        lower_value=0.0,
                        upper_value=num_inactive_deviation,
                    )
                )
                extra_linear_constraints.append(
                    ExtraLinearConstraint(
                        stoichiometries=dict.fromkeys(
                            extra_active_zvars + existing_active_zvars, 1.0
                        ),
                        lower_value=max(
                            0,
                            len(existing_active_zvars)
                            + len(extra_active_zvars)
                            - num_active_deviation,
                        ),
                        upper_value=len(existing_active_zvars)
                        + len(extra_active_zvars)
                        + num_active_deviation,
                    )
                )
                if is_objsense_maximization(objective_sense):
                    extra_linear_constraints.append(
                        ExtraLinearConstraint(
                            stoichiometries=objective_target,
                            lower_value=best_result[OBJECTIVE_VAR_NAME],
                        )
                    )
                else:
                    extra_linear_constraints.append(
                        ExtraLinearConstraint(
                            stoichiometries=objective_target,
                            upper_value=best_result[OBJECTIVE_VAR_NAME],
                        )
                    )
                for lp_target in lp_targets:
                    postprocessing_calculation_settings.append(
                        PostprocessingSetting(
                            target_type,
                            target_reac,
                            lp_target,
                            num_active_deviation,
                            num_inactive_deviation,
                            tuple(extra_inactive_zvars + extra_active_zvars),
                            deepcopy(extra_linear_constraints),
                        ),
                    )

    return postprocessing_calculation_settings


# @validate_call
def _postprocessing_block(
    evo_settings: EvolutionSettings,
    pp_setting: PostprocessingSetting,
    reac_couples_list: tuple[tuple[str, ...], ...],
    objective_target_as_dict: dict[str, float],
    verbose: bool,
) -> tuple[LpNlpBlockResult, PostprocessingSetting]:
    setting_model = deepcopy(evo_settings.cobrak_model)
    setting_model.extra_binary_vars = list(pp_setting.extra_zvars)
    setting_model.extra_linear_constraints += pp_setting.extra_linear_constraints

    ectfba_result, binaries = _ectfba_block(
        cobrak_model=setting_model,
        objective_target=evo_settings.objective_target,
        objective_sense=evo_settings.objective_sense,
        variability_dict=evo_settings.variability_dict,
        lp_solver=evo_settings.lp_solver,
        correction_config=evo_settings.correction_config,
        ignore_nonlinear_extra_terms_in_ectfbas=evo_settings.ignore_nonlinear_extra_terms_in_lps,
        binaries=(),
        reac_couples_list=reac_couples_list,
        do_reac_deletions=False,
        verbose=False,
    )

    if OBJECTIVE_VAR_NAME not in ectfba_result or None in ectfba_result.values():
        if verbose:
            print(
                f"{pp_setting.target_type} {pp_setting.target_reac} {pp_setting.lp_objective} | AC {pp_setting.num_active_deviation} INAC {pp_setting.num_inactive_deviation} | ecTFBA infeasible"
            )
        return LpNlpBlockResult(
            original_binaries=binaries,
            lp_binaries=(),
            nlp_binaries=(),
        ), pp_setting
    ectfba_nlp_block_result = _ectfba_nlp_block(
        cobrak_model=evo_settings.cobrak_model,
        binaries=(),
        reac_couples_list=reac_couples_list,
        lp_objective_target=pp_setting.lp_objective,
        lp_objective_sense=evo_settings.objective_sense,
        lp_extra_linear_constraints=[
            ExtraLinearConstraint(
                stoichiometries=objective_target_as_dict,
                lower_value=ectfba_result[OBJECTIVE_VAR_NAME] - 1e-8,
                upper_value=ectfba_result[OBJECTIVE_VAR_NAME] + 1e-8,
            )
        ]
        + pp_setting.extra_linear_constraints,
        nlp_objective_target=evo_settings.objective_target,
        nlp_objective_sense=evo_settings.objective_sense,
        variability_dict=evo_settings.variability_dict,
        with_kappa=evo_settings.with_kappa,
        with_gamma=evo_settings.with_gamma,
        with_iota=evo_settings.with_iota,
        with_alpha=evo_settings.with_alpha,
        lp_solver=evo_settings.lp_solver,
        nlp_solver=evo_settings.nlp_solver,
        nlp_strict_mode=evo_settings.nlp_strict_mode,
        nlp_single_strict_reacs=evo_settings.nlp_single_strict_reacs,
        correction_config=evo_settings.correction_config,
        ignore_nonlinear_extra_terms_in_ectfbas=evo_settings.ignore_nonlinear_extra_terms_in_lps,
        delete_nonthermodynamic_reacs_for_nlp=False,
        approximation_value=evo_settings.approximation_value,
        verbose=False,
        lp_extra_binary_vars=list(pp_setting.extra_zvars),
    )

    if verbose:
        if (
            ectfba_nlp_block_result.lp_result
            and OBJECTIVE_VAR_NAME in ectfba_nlp_block_result.lp_result
        ):
            ectfba_value = round(
                ectfba_nlp_block_result.lp_result[OBJECTIVE_VAR_NAME], 8
            )
        else:
            ectfba_value = None
        if (
            ectfba_nlp_block_result.nlp_result
            and OBJECTIVE_VAR_NAME in ectfba_nlp_block_result.nlp_result
        ):
            nlp_value = round(ectfba_nlp_block_result.nlp_result[OBJECTIVE_VAR_NAME], 8)
        else:
            nlp_value = None
        print(
            f"{pp_setting.target_type} {pp_setting.target_reac} {pp_setting.lp_objective} | AC {pp_setting.num_active_deviation} INAC {pp_setting.num_inactive_deviation} | "
            f"ecTFBA {ectfba_value} "
            f"NLP: {nlp_value}"
        )
    return ectfba_nlp_block_result, pp_setting


# @validate_call(validate_return=True)
def postprocessing(
    settings: EvolutionSettings,
    evolution_results: EvolutionResult,
    num_active_deviations: tuple[int, ...],
    num_inactive_deviations: tuple[int, ...],
    lp_targets: tuple[str, ...] = (
        "MAXZ",
        "MINZ",
    ),
    idx_of_best_nlp: int = 0,
    verbose: bool = False,
    target_reaction_selection: list[str] = [],
) -> tuple[EvolutionResult, list[PostprocessingOutcome]]:
    evolution_results: _InternalEvolutionResult = (
        _external_to_internal_evolution_result(deepcopy(evolution_results))
    )

    best_nlp = evolution_results.best_nlps[idx_of_best_nlp]

    if type(settings.objective_target) is str:
        objective_target_as_dict: dict[str, int | float] = {
            settings.objective_target: 1.0
        }
    elif type(settings.objective_target) is dict:
        objective_target_as_dict: dict[str, int | float] = settings.objective_target

    best_binary = next(
        key
        for key, value in evolution_results.original_binary_results.items()
        if best_nlp[OBJECTIVE_VAR_NAME] in value
    )
    postprocessing_settings: list[PostprocessingSetting] = (
        _get_postprocessing_calculation_settings(
            cobrak_model=settings.cobrak_model,
            best_binary=best_binary,
            best_result=best_nlp,
            reac_couples_list=evolution_results.reac_couples_list,
            num_active_deviations=num_active_deviations,
            num_inactive_deviations=num_inactive_deviations,
            lp_targets=lp_targets,
            objective_target=objective_target_as_dict,
            objective_sense=settings.objective_sense,
            target_reaction_selection=target_reaction_selection,
        )
    )

    with parallel_config(backend="loky", inner_max_num_threads=1):
        results: list[tuple[LpNlpBlockResult, PostprocessingSetting]] = Parallel(
            n_jobs=settings.num_used_cpu_cores, verbose=10
        )(
            delayed(_postprocessing_block)(
                settings,
                pp_setting,
                evolution_results.reac_couples_list,
                objective_target_as_dict,
                verbose,
            )
            for pp_setting in postprocessing_settings
        )
    postprocessing_outcomes: list[PostprocessingOutcome] = []
    for result, pp_setting in results:
        if result.lp_result and not result.nlp_result:
            evolution_results.original_binary_results.setdefault(
                result.original_binaries, []
            ).append(None)
            evolution_results.ectfba_binary_results[result.lp_binaries] = None
        elif not result.lp_result:
            evolution_results.original_binary_results.setdefault(
                result.original_binaries, []
            ).append(None)
        if result.nlp_result:
            evolution_results.original_binary_results.setdefault(
                result.original_binaries, []
            ).append(result.nlp_result[OBJECTIVE_VAR_NAME])
            evolution_results.ectfba_binary_results[result.nlp_binaries] = (
                result.nlp_result[OBJECTIVE_VAR_NAME]
            )
            evolution_results.nlp_binary_results[result.nlp_binaries] = (
                result.nlp_result[OBJECTIVE_VAR_NAME]
            )
            postprocessing_outcomes.append(
                PostprocessingOutcome(
                    target_reac=pp_setting.target_reac,
                    target_type=pp_setting.target_type,
                    lp_objective=pp_setting.lp_objective,
                    num_active_deviation=pp_setting.num_active_deviation,
                    num_inactive_deviation=pp_setting.num_inactive_deviation,
                    ectfba_result=result.lp_result[OBJECTIVE_VAR_NAME]
                    if result.lp_result
                    else None,
                    nlp_result=result.nlp_result[OBJECTIVE_VAR_NAME]
                    if result.nlp_result
                    else None,
                )
            )
            if (
                (
                    is_objsense_maximization(settings.objective_sense)
                    and result.nlp_result[OBJECTIVE_VAR_NAME]
                    > evolution_results.best_nlps[0][OBJECTIVE_VAR_NAME]
                )
                or (not is_objsense_maximization(settings.objective_sense))
                and result.nlp_result[OBJECTIVE_VAR_NAME]
                < evolution_results.best_nlps[0][OBJECTIVE_VAR_NAME]
            ):
                evolution_results.best_nlps = [result.nlp_result]
            elif (
                result.nlp_result[OBJECTIVE_VAR_NAME]
                == evolution_results.best_nlps[0][OBJECTIVE_VAR_NAME]
            ):
                evolution_results.best_nlps.append(result.nlp_result)

    return _internal_to_external_evolution_result(
        evolution_results
    ), postprocessing_outcomes


##########################


# @validate_call(validate_return=True)
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
            variability_dict[reac_id][1] <= 0.0
            or variability_dict[reac_id][0] > 0.0
            or reac_id in error_scenario
            for reac_id in reac_ids
        ):
            continue
        if all(
            (cobrak_model.reactions[reac_id].dG0 is None)
            and not cobrak_model.reactions[reac_id].enzyme_reaction_data
            for reac_id in reac_ids
        ):
            continue
        # Discard couple with objective target(s)
        objective_target_strlist: list[str]
        if type(objective_target) is str:
            objective_target_strlist = [objective_target]
        elif type(objective_target) is dict:
            objective_target_strlist = list(objective_target.keys())
        if any(
            objective_target in reac_ids
            for objective_target in objective_target_strlist
        ):
            continue
        reac_couples_list.append(tuple(reac_ids))
    return tuple(reac_couples_list)


######## PUBLIC FUNCTIONS ########
# @validate_call(validate_return=True)
def perform_nlp_evolutionary_optimization(
    settings: EvolutionSettings,
    existing_result: EvolutionResult = EvolutionResult(),
) -> EvolutionResult:
    """"""
    # CONVERT JSON-ABLE EVOLUTION RESULT TO INTERNAL EVOLUTION RESULT
    internal_existing_result = _external_to_internal_evolution_result(existing_result)

    # PHASE 1: BUILD INDEX TO REAC COUPLES DATA, AND CPU DATA
    if not existing_result.reac_couples_list:
        internal_existing_result.reac_couples_list = _get_idx_to_reac_ids(
            cobrak_model=settings.cobrak_model,
            objective_target=settings.objective_target,
            variability_dict=settings.variability_dict,
            error_scenario=settings.correction_config.error_scenario,
        )
    else:
        internal_existing_result.reac_couples_list = deepcopy(
            existing_result.reac_couples_list
        )

    # PHASE 2: ACTUAL EVOLUTION ALGORITHM (USING GIVEN OR SAMPLED RESULTS AS STARTING POINTS)
    internal_evolution_result = _evolution(
        settings=settings,
        existing_result=internal_existing_result,
    )

    return _internal_to_external_evolution_result(internal_evolution_result)


###################
###################
def _simple_postprocessing_block(
    evo_settings: EvolutionSettings,
    target_reac: str,
    opt_dict: dict[str, float],
    objective_target_as_dict: dict[str, float],
    objective_sense: int,
    reac_couples_list: tuple[tuple[str, ...], ...],
    verbose: bool,
) -> tuple[str, dict[str, float]]:
    pp_lp_model = deepcopy(evo_settings.cobrak_model)
    pp_lp_vardict = deepcopy(evo_settings.variability_dict)
    pp_optdict = deepcopy(opt_dict)

    if target_reac != "BASETEST":
        is_ko = pp_optdict.get(target_reac, 0.0) > 5e-9

        pp_lp_vardict[target_reac] = (
            0.0 if is_ko else min(1e-5, pp_lp_vardict[target_reac][1]),
            0.0 if is_ko else pp_lp_vardict[target_reac][1],
        )
        pp_lp_model.extra_binary_vars.append("EXTRA_BINVAR_ONE")
        pp_lp_model.extra_linear_constraints.append(
            ExtraLinearConstraint(
                stoichiometries={
                    "EXTRA_BINVAR_ONE": 1.0,
                },
                lower_value=1.0 - 1e-6,
                upper_value=1.0 + 1e-6,
            )
        )
        for couple_reacs in reac_couples_list:
            reac_id = couple_reacs[0]
            if reac_id == target_reac:
                continue
            binary_varname = f"EXTRAZ_{reac_id}"
            pp_lp_model.extra_binary_vars.append(binary_varname)
            if pp_optdict.get(couple_reacs[0], 0.0) > 5e-9:
                # Originally active reac: vᵢ ≤ M(1-zᵢ) ⇒ vᵢ ≤ M-M⋅zᵢ ⇒ vᵢ - M + M⋅zᵢ ≤ 0
                # i.e. zᵢ=1→vᵢ≤0
                pp_lp_model.extra_linear_constraints.append(
                    ExtraLinearConstraint(
                        stoichiometries={
                            reac_id: 1.0,
                            "EXTRA_BINVAR_ONE": -BIG_M,
                            binary_varname: +BIG_M,
                        },
                        lower_value=None,
                        upper_value=0.0,
                    )
                )
            else:
                # Originally inactive reac: vᵢ ≤ M⋅zᵢ
                # i.e. zᵢ=1→vᵢ≥0
                pp_lp_model.extra_linear_constraints.append(
                    ExtraLinearConstraint(
                        stoichiometries={
                            reac_id: 1.0,
                            binary_varname: -BIG_M,
                        },
                        lower_value=None,
                        upper_value=0.0,
                    )
                )
        pp_lp_model.extra_linear_constraints.append(
            ExtraLinearConstraint(
                stoichiometries=objective_target_as_dict,
                lower_value=opt_dict[OBJECTIVE_VAR_NAME] - 1e-8
                if is_objsense_maximization(objective_sense)
                else None,
                upper_value=opt_dict[OBJECTIVE_VAR_NAME] + 1e-8
                if not is_objsense_maximization(objective_sense)
                else None,
            )
        )
        try:
            result = perform_lp_optimization(
                cobrak_model=pp_lp_model,
                objective_target={
                    binvar: 1.0
                    for binvar in pp_lp_model.extra_binary_vars
                    if binvar.startswith("EXTRAZ_")
                },
                objective_sense=-1,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=True,
                variability_dict=pp_lp_vardict,
                solver=evo_settings.lp_solver,
                ignore_nonlinear_terms=evo_settings.ignore_nonlinear_extra_terms_in_lps,
                correction_config=evo_settings.correction_config,
                verbose=False,
            )
        except (ApplicationError, AttributeError, ValueError):
            return target_reac, {}
        if not result[ALL_OK_KEY]:
            return target_reac, {}
        pp_optdict[target_reac] = 0.0 if is_ko else 1.0
        for key, value in result.items():
            if not key.startswith("EXTRAZ_"):
                continue
            if value > 0.99:
                pp_optdict[key[len("EXTRAZ_") :]] = 1.0
                print("AAA", key[len("EXTRAZ_") :])
    try:
        nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
            cobrak_model=evo_settings.cobrak_model,
            objective_target=evo_settings.objective_target,
            objective_sense=evo_settings.objective_sense,
            optimization_dict=pp_optdict,
            variability_dict=evo_settings.variability_dict,
            solver=evo_settings.nlp_solver,  # Solver(name="pounce"), # IPOPT_MA57,
            verbose=False,
            with_kappa=evo_settings.with_kappa,
            with_gamma=evo_settings.with_gamma,
            with_iota=evo_settings.with_iota,
            with_alpha=evo_settings.with_alpha,
            strict_mode=evo_settings.nlp_strict_mode,
            single_strict_reacs=evo_settings.nlp_single_strict_reacs,
            approximation_value=evo_settings.approximation_value,
        )
    except (ApplicationError, AttributeError, ValueError):
        return target_reac, {}
    if not nlp_result[ALL_OK_KEY]:
        return target_reac, {}
    if verbose:
        print(
            f"POSTPROCESS {'KO' if is_ko else 'ACT'} {target_reac}",
            nlp_result[OBJECTIVE_VAR_NAME],
        )

    return target_reac, nlp_result


def simple_postprocessing(
    settings: EvolutionSettings,
    evolution_results: EvolutionResult,
    idx_of_best_nlp: int = 0,
    verbose: bool = False,
    target_reaction_selection: list[str] = [],
) -> list[tuple[str, dict[str, float]]]:
    if type(settings.objective_target) is str:
        objective_target_as_dict: dict[str, int | float] = {
            settings.objective_target: 1.0
        }
    elif type(settings.objective_target) is dict:
        objective_target_as_dict: dict[str, int | float] = settings.objective_target

    best_nlp = evolution_results.best_nlps[idx_of_best_nlp]
    ko_targets = [
        reac_id
        for reac_id in settings.cobrak_model.reactions
        if best_nlp.get(reac_id, 0.0) > 5e-9
        and settings.cobrak_model.reactions[reac_id].min_flux == 0.0
    ]
    if target_reaction_selection:
        ko_targets = [
            ko_target
            for ko_target in ko_targets
            if ko_target in target_reaction_selection
        ]
    results: list[tuple[str, dict[str, float]]] = Parallel(
        n_jobs=settings.num_used_cpu_cores, verbose=10
    )(
        delayed(_simple_postprocessing_block)(
            settings,
            ko_target,
            best_nlp,
            objective_target_as_dict,
            settings.objective_sense,
            evolution_results.reac_couples_list,
            verbose,
        )
        for ko_target in ko_targets
    )
    working_results = [
        result for result in results if result[1].get(OBJECTIVE_VAR_NAME, None)
    ]
    better_results = [
        result
        for result in working_results
        if (
            is_objsense_maximization(settings.objective_sense)
            and (
                result[1][OBJECTIVE_VAR_NAME]
                > evolution_results.best_nlps[idx_of_best_nlp][OBJECTIVE_VAR_NAME]
            )
        )
        or (
            not is_objsense_maximization(settings.objective_sense)
            and (
                result[1][OBJECTIVE_VAR_NAME]
                < evolution_results.best_nlps[idx_of_best_nlp][OBJECTIVE_VAR_NAME]
            )
        )
    ]
    return better_results


###################
###################
