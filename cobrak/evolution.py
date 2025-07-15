from copy import deepcopy
from random import sample
from tempfile import TemporaryDirectory
from time import time
from typing import Any, Literal

from joblib import Parallel, cpu_count, delayed
from numpy.random import randint
from pyomo.common.errors import ApplicationError

from .constants import ALL_OK_KEY, OBJECTIVE_VAR_NAME, Z_VAR_PREFIX
from .dataclasses import CorrectionConfig, ExtraLinearConstraint, Model, Solver
from .genetic import COBRAKGENETIC
from .io import json_load, json_write
from .lps import perform_lp_optimization, perform_lp_variability_analysis
from .nlps import perform_nlp_irreversible_optimization_with_active_reacs_only
from .pso_parallel import COBRAKPAPSO
from .standard_solvers import IPOPT, SCIP
from .utilities import (
    get_active_reacs_from_optimization_dict,
    get_files,
    get_stoichiometrically_coupled_reactions,
    is_objsense_maximization,
    split_list,
    standardize_folder,
)


class COBRAKProblem:
    def __init__(
        self,
        cobrak_model: Model,
        objective_target: dict[str, float],
        objective_sense: int,
        variability_dict: dict[str, tuple[float, float]],
        nlp_dict_list: list[dict[str, float]],
        best_value: float,
        with_kappa: bool = True,
        with_gamma: bool = True,
        with_iota: bool = True,
        with_alpha: bool = True,
        num_gens: int = 5,
        algorithm: Literal["pso", "genetic"] = "pso",
        lp_solver: Solver = SCIP,
        nlp_solver: Solver = IPOPT,
        objvalue_json_path: str = "",
        max_rounds_same_objvalue: float = float("inf"),
        correction_config: CorrectionConfig = CorrectionConfig(),
        min_abs_objvalue: float = 1e-6,
        pop_size: int | None = None,
    ) -> None:
        self.original_cobrak_model: Model = deepcopy(cobrak_model)
        self.objective_target = objective_target
        self.objective_sense = objective_sense
        self.blocked_reacs: list[str] = []
        self.initial_xs_list: list[list[int | float]] = [
            [] for _ in range(len(nlp_dict_list))
        ]
        self.minimal_xs_dict: dict[float, list[float]] = {}
        self.variability_data = deepcopy(variability_dict)
        self.algorithm = algorithm

        reac_couples = get_stoichiometrically_coupled_reactions(
            self.original_cobrak_model
        )

        objective_target_ids = list(objective_target.keys())
        filtered_reac_couples: list[tuple[str, ...]] = []
        for reac_couple in reac_couples:
            filtered_reac_couple = [
                reac_id
                for reac_id in reac_couple
                if (abs(self.variability_data[reac_id][1]) > 0.0)
                and (abs(self.variability_data[reac_id][0]) <= 0.0)
                and not (
                    cobrak_model.reactions[reac_id].dG0 is None
                    and cobrak_model.reactions[reac_id].enzyme_reaction_data is None
                )
            ]

            found_invalid_id = False
            for objective_target_id in objective_target_ids:
                if objective_target_id in filtered_reac_couple:
                    found_invalid_id = True
            for var_id in correction_config.error_scenario:
                if var_id in filtered_reac_couple:
                    found_invalid_id = True

            if found_invalid_id:
                continue

            if len(filtered_reac_couple) > 0:
                filtered_reac_couples.append(tuple(filtered_reac_couple))

        self.idx_to_reac_ids: dict[int, tuple[str, ...]] = {}
        couple_idx = 0
        for filtered_reac_couplex in filtered_reac_couples:
            nlp_idx = 0
            for nlp_dict in nlp_dict_list:
                first_reac_id = filtered_reac_couplex[0]
                if (first_reac_id in nlp_dict) and (nlp_dict[first_reac_id] > 0.0):
                    self.initial_xs_list[nlp_idx].append(
                        1.0 if algorithm == "pso" else 1
                    )
                else:
                    self.initial_xs_list[nlp_idx].append(
                        0.0 if algorithm == "pso" else 0
                    )

                nlp_idx += 1  # noqa: SIM113

            self.idx_to_reac_ids[couple_idx] = filtered_reac_couplex
            couple_idx += 1

        self.with_kappa = with_kappa
        self.with_gamma = with_gamma
        self.with_iota = with_iota
        self.with_alpha = with_alpha
        self.dim = couple_idx
        self.num_gens = num_gens
        self.lp_solver = lp_solver
        self.nlp_solver = nlp_solver
        self.temp_directory_name = ""
        self.best_value = best_value
        self.objvalue_json_path = objvalue_json_path
        self.max_rounds_same_objvalue = max_rounds_same_objvalue
        self.correction_config = correction_config
        self.min_abs_objvalue = min_abs_objvalue
        self.pop_size = pop_size

    def fitness(
        self,
        x: list[float | int],
    ) -> list[tuple[float, list[float | int]]]:
        # Preliminary TFBA :3
        deactivated_reactions: list[str] = []
        for couple_idx, reac_ids in self.idx_to_reac_ids.items():
            if x[couple_idx] <= 0.02:
                deactivated_reactions.extend(reac_ids)

        try:
            first_ectfba_dict = perform_lp_optimization(
                cobrak_model=self.original_cobrak_model,
                objective_target=self.objective_target,
                objective_sense=self.objective_sense,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=True,
                variability_dict=deepcopy(self.variability_data),
                ignored_reacs=deactivated_reactions,
                solver=self.lp_solver,
                correction_config=self.correction_config,
            )
        except (ApplicationError, AttributeError, ValueError):
            first_ectfba_dict = {ALL_OK_KEY: False}
        if not first_ectfba_dict[ALL_OK_KEY]:
            return [(1_000_000.0, [])]

        nlp_results: list[dict[str, float]] = []

        if is_objsense_maximization(self.objective_sense):
            lower_value = first_ectfba_dict[OBJECTIVE_VAR_NAME] - 1e-12
            upper_value = None
        else:
            lower_value = None
            upper_value = first_ectfba_dict[OBJECTIVE_VAR_NAME] + 1e-12

        maxz_model = deepcopy(self.original_cobrak_model)
        maxz_model.extra_linear_constraints = [
            ExtraLinearConstraint(
                stoichiometries=self.objective_target,
                lower_value=lower_value,
                upper_value=upper_value,
            )
        ]
        maxz_model.extra_linear_constraints += [
            ExtraLinearConstraint(
                stoichiometries={f"{Z_VAR_PREFIX}{reac_id}": 1.0},
                upper_value=0.0,
            )
            for (reac_id, reac_data) in self.original_cobrak_model.reactions.items()
            if (reac_data.dG0 is not None) and (reac_id in deactivated_reactions)
        ]
        eligible_z_sum_objective = {
            f"{Z_VAR_PREFIX}{reac_id}": 1.0
            for (reac_id, reac_data) in self.original_cobrak_model.reactions.items()
            if (reac_data.dG0 is not None)
            and (self.variability_data[reac_id][1] > 0.0)
            and (reac_id not in deactivated_reactions)
        }
        try:
            maxz_ectfba_dict = perform_lp_optimization(
                cobrak_model=maxz_model,
                objective_target=eligible_z_sum_objective,
                objective_sense=+1,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=True,
                variability_dict=deepcopy(self.variability_data),
                ignored_reacs=deactivated_reactions,
                solver=self.lp_solver,
                correction_config=self.correction_config,
            )
        except (ApplicationError, AttributeError, ValueError):
            maxz_ectfba_dict = {ALL_OK_KEY: False}
        if maxz_ectfba_dict[ALL_OK_KEY]:
            used_maxz_tfba_dict: dict[str, float] = {}
            for var_id in maxz_ectfba_dict:
                if var_id not in self.original_cobrak_model.reactions:
                    continue
                reaction = self.original_cobrak_model.reactions[var_id]
                if (
                    (reaction.dG0 is None) and (var_id not in deactivated_reactions)
                ) or (
                    (reaction.dG0 is not None)
                    and (maxz_ectfba_dict[f"{Z_VAR_PREFIX}{var_id}"] > 0.0)
                ):
                    used_maxz_tfba_dict[var_id] = 1.0
                else:
                    used_maxz_tfba_dict[var_id] = 0.0

            used_maxz_tfba_dict[ALL_OK_KEY] = True

            if used_maxz_tfba_dict[ALL_OK_KEY]:
                try:
                    second_nlp_dict = (
                        perform_nlp_irreversible_optimization_with_active_reacs_only(
                            cobrak_model=self.original_cobrak_model,
                            objective_target=self.objective_target,
                            objective_sense=self.objective_sense,
                            optimization_dict=deepcopy(used_maxz_tfba_dict),
                            variability_dict=deepcopy(self.variability_data),
                            with_kappa=self.with_kappa,
                            with_gamma=self.with_gamma,
                            with_iota=self.with_iota,
                            with_alpha=self.with_alpha,
                            solver=self.nlp_solver,
                            correction_config=self.correction_config,
                        )
                    )
                    if second_nlp_dict[ALL_OK_KEY] and (
                        abs(second_nlp_dict[OBJECTIVE_VAR_NAME]) > self.min_abs_objvalue
                    ):
                        nlp_results.append(second_nlp_dict)
                except (ApplicationError, AttributeError, ValueError):
                    pass

        ####
        try:
            minz_ectfba_dict = perform_lp_optimization(
                cobrak_model=maxz_model,
                objective_target=eligible_z_sum_objective,
                objective_sense=-1,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=True,
                variability_dict=deepcopy(self.variability_data),
                ignored_reacs=deactivated_reactions,
                solver=self.lp_solver,
                correction_config=self.correction_config,
            )
        except (ApplicationError, AttributeError, ValueError):
            minz_ectfba_dict = {ALL_OK_KEY: False}
        if minz_ectfba_dict[ALL_OK_KEY]:
            used_minz_tfba_dict: dict[str, float] = {}
            for var_id in minz_ectfba_dict:
                if var_id not in self.original_cobrak_model.reactions:
                    continue
                reaction = self.original_cobrak_model.reactions[var_id]
                if (
                    (reaction.dG0 is None) and (var_id not in deactivated_reactions)
                ) or (
                    (reaction.dG0 is not None)
                    and (minz_ectfba_dict[f"{Z_VAR_PREFIX}{var_id}"] > 0.0)
                ):
                    used_minz_tfba_dict[var_id] = 1.0
                else:
                    used_minz_tfba_dict[var_id] = 0.0

            used_minz_tfba_dict[ALL_OK_KEY] = True

            if used_minz_tfba_dict[ALL_OK_KEY]:
                try:
                    third_nlp_dict = (
                        perform_nlp_irreversible_optimization_with_active_reacs_only(
                            cobrak_model=self.original_cobrak_model,
                            objective_target=self.objective_target,
                            objective_sense=self.objective_sense,
                            optimization_dict=deepcopy(used_minz_tfba_dict),
                            variability_dict=deepcopy(self.variability_data),
                            with_kappa=self.with_kappa,
                            with_gamma=self.with_gamma,
                            with_iota=self.with_iota,
                            with_alpha=self.with_alpha,
                            solver=self.nlp_solver,
                            correction_config=self.correction_config,
                        )
                    )
                    if third_nlp_dict[ALL_OK_KEY] and (
                        abs(third_nlp_dict[OBJECTIVE_VAR_NAME]) > self.min_abs_objvalue
                    ):
                        nlp_results.append(third_nlp_dict)
                except (ApplicationError, AttributeError, ValueError):
                    pass
        ####

        output: list[tuple[float, list[float | int]]] = [(1_000_000, [])]
        for nlp_result in nlp_results:
            objvalues = [nlp_result[OBJECTIVE_VAR_NAME] for nlp_result in nlp_results]
            if is_objsense_maximization(self.objective_sense):
                opt_idx = objvalues.index(max(objvalues))
            else:
                opt_idx = objvalues.index(min(objvalues))
            opt_nlp_dict = nlp_results[opt_idx]

            objective_value = opt_nlp_dict[OBJECTIVE_VAR_NAME]

            if self.temp_directory_name:
                filename = f"{self.temp_directory_name}{objective_value}{time()}{randint(0, 1_000_000_000)}.json"
                json_write(filename, opt_nlp_dict)

            if is_objsense_maximization(self.objective_sense):
                objective_value *= -1

            print("No error, objective value is:", objective_value)

            active_nlp_x: list[float | int] = [
                0 for _ in range(len(list(self.idx_to_reac_ids.keys())))
            ]
            for couple_idx, reac_ids in self.idx_to_reac_ids.items():
                reac_id = reac_ids[0]
                if reac_id not in opt_nlp_dict or opt_nlp_dict[reac_id] < 1e-11:
                    set_value = 0
                else:
                    set_value = 1
                active_nlp_x[couple_idx] = set_value
            output.append((objective_value, active_nlp_x))

        return output

    def optimize(self) -> dict[float, list[dict[str, float]]]:
        temp_directory = TemporaryDirectory()
        self.temp_directory_name = standardize_folder(temp_directory.name)

        match self.algorithm:
            case "pso":
                evolution = COBRAKPAPSO(
                    fitness_function=self.fitness,
                    xs_dim=self.dim,
                    extra_xs=self.initial_xs_list,
                    gen=self.num_gens,
                    objvalue_json_path=self.objvalue_json_path,
                    max_rounds_same_objvalue=self.max_rounds_same_objvalue,
                    pop_size=self.pop_size,
                )
            case "genetic":
                evolution = COBRAKGENETIC(
                    fitness_function=self.fitness,
                    xs_dim=self.dim,
                    extra_xs=self.initial_xs_list,
                    gen=self.num_gens,
                    objvalue_json_path=self.objvalue_json_path,
                    max_rounds_same_objvalue=self.max_rounds_same_objvalue,
                    pop_size=self.pop_size,
                )
            case _:
                print(
                    f"ERROR: Evolution algorithm {self.algorithm} does not exist! Use either 'pso' or 'genetic'."
                )
                raise ValueError
        evolution.run()

        result_dict: dict[float, list[dict[str, float]]] = {}
        for json_filename in get_files(self.temp_directory_name):
            json_data = json_load(f"{self.temp_directory_name}{json_filename}", Any)
            objective_value = json_data[OBJECTIVE_VAR_NAME]
            if objective_value not in result_dict:
                result_dict[objective_value] = []
            result_dict[objective_value].append(deepcopy(json_data))

        temp_directory.cleanup()

        return {
            key: result_dict[key] for key in sorted(result_dict.keys(), reverse=True)
        }


from pyomo.environ import Binary, Constraint, Reals, Var

from .constants import BIG_M
from .lps import add_statuses_to_optimziation_dict, get_lp_from_cobrak_model
from .pyomo_functionality import (
    add_objective_to_model,
    get_solver,
)
from .utilities import (
    add_objective_value_as_extra_linear_constraint,
    apply_variability_dict,
    delete_orphaned_metabolites_and_enzymes,
    get_pyomo_solution_as_dict,
)


def _postprocess_batch(
    reac_couples: list[str],
    target_couples: list[tuple[str, list[str]]],
    active_reacs: list[str],
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_data: dict[str, tuple[float, float]],
    pyomo_lp_solver,
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    lp_solver: Solver = SCIP,
    nlp_solver: Solver = IPOPT,
    verbose: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
    onlytested: str = "",
):
    original_var_data = deepcopy(variability_data)
    feasible_switches = []
    for target_type, target_couple, max_allowed_changes in target_couples:
        if onlytested and (not any(onlytested in reac_id for reac_id in target_couple)):
            continue

        for at_maximum in (False,):
            target_cobrak_model = deepcopy(cobrak_model)
            variability_data = deepcopy(original_var_data)
            if target_type == "deac":
                print(
                    f"DEACTIVATING {target_couple} {at_maximum}, max changes: {max_allowed_changes}"
                )
                for reac_id in target_couple:
                    del target_cobrak_model.reactions[reac_id]
                target_cobrak_model = delete_orphaned_metabolites_and_enzymes(
                    target_cobrak_model
                )
            else:
                print(
                    f"ACTIVATING {target_couple} {at_maximum}, max changes: {max_allowed_changes}"
                )
                variability_data[target_couple[0]] = (
                    1e-5,
                    variability_data[target_couple[0]][1],
                )

            baselp = get_lp_from_cobrak_model(
                cobrak_model=target_cobrak_model,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=False,
                with_flux_sum_var=False,
                add_extra_linear_constraints=True,
                correction_config=correction_config,
            )
            baselp = apply_variability_dict(
                baselp,
                target_cobrak_model,
                variability_data,
                error_scenario=correction_config.error_scenario,
            )

            active_z_var_changes_sum = 0.0
            for active_reac in active_reacs:
                if (
                    active_reac in target_couple
                    or target_cobrak_model.reactions[active_reac].dG0 is None
                ):
                    continue
                z_var_id = f"{Z_VAR_PREFIX}{active_reac}"
                z_var_change_id = f"{z_var_id}_change"
                setattr(baselp, z_var_change_id, Var(within=Binary))
                setattr(
                    baselp,
                    f"fix_{z_var_id}",
                    Constraint(
                        expr=getattr(baselp, z_var_id)
                        == 1.0 - getattr(baselp, z_var_change_id)
                    ),
                )
                active_z_var_changes_sum += getattr(baselp, z_var_change_id)
            max_z_var_changes_id = "max_active_z_var_changes"
            setattr(
                baselp,
                max_z_var_changes_id,
                Var(within=Reals),
            )
            getattr(baselp, max_z_var_changes_id).lb = 0.0
            getattr(baselp, max_z_var_changes_id).ub = max_allowed_changes
            setattr(
                baselp,
                "active_original_z_vars_sum",
                Constraint(
                    expr=active_z_var_changes_sum
                    <= getattr(baselp, max_z_var_changes_id)
                ),
            )

            baselp = add_objective_to_model(
                baselp,
                objective_target,
                objective_sense,
                "ORIGINAL_OBJ",
                "ORIGINAL_OBJ_VAR",
            )
            results = pyomo_lp_solver.solve(
                baselp, tee=verbose, **lp_solver.solve_extra_options
            )
            lp_resultdict = add_statuses_to_optimziation_dict(
                get_pyomo_solution_as_dict(baselp), results
            )
            if not lp_resultdict[ALL_OK_KEY]:
                continue
            new_active_reacs = get_active_reacs_from_optimization_dict(
                target_cobrak_model, lp_resultdict
            )
            getattr(baselp, "ORIGINAL_OBJ").deactivate()
            if at_maximum:
                getattr(baselp, "ORIGINAL_OBJ_VAR").lb = (
                    getattr(baselp, "ORIGINAL_OBJ_VAR").value - 1e-6
                )

            for active_reac in new_active_reacs:
                if (
                    active_reac in target_couple
                    or target_cobrak_model.reactions[active_reac].dG0 is None
                ):
                    continue
                if active_reac in active_reacs:
                    continue
                getattr(baselp, f"{Z_VAR_PREFIX}{active_reac}").lb = 1.0

            EXTRAZ_PREFIX = "extra_z_var_"
            extraz_names: list[str] = []
            for reac_couple in reac_couples:
                if target_couple == reac_couple:
                    continue
                if reac_couple[0] in active_reacs:
                    continue
                if all(
                    target_cobrak_model.reactions[reac_id].enzyme_reaction_data is None
                    and target_cobrak_model.reactions[reac_id].dG0 is None
                    for reac_id in reac_couple
                ):
                    continue
                if variability_data[reac_couple[0]][1] == 0.0:
                    continue

                extraz_names.append(f"{EXTRAZ_PREFIX}_of_couple_of_{reac_couple[0]}")
                setattr(
                    baselp,
                    extraz_names[-1],
                    Var(
                        within=Binary,
                    ),
                )
                setattr(
                    baselp,
                    f"{EXTRAZ_PREFIX}_constraint_of_couple_of_{reac_couple[0]}",
                    Constraint(
                        expr=getattr(baselp, reac_couple[0])
                        <= BIG_M * getattr(baselp, extraz_names[-1])
                    ),
                )
                couple_reacs_with_dG0 = [
                    reac_id
                    for reac_id in target_cobrak_model.reactions
                    if target_cobrak_model.reactions[reac_id].dG0 is not None
                ]
                if len(couple_reacs_with_dG0) == 0:
                    continue
                setattr(
                    baselp,
                    f"{EXTRAZ_PREFIX}_z_var_constraint_of_couple_of_{reac_couple[0]}",
                    Constraint(
                        expr=getattr(
                            baselp, f"{Z_VAR_PREFIX}{couple_reacs_with_dG0[0]}"
                        )
                        >= getattr(baselp, extraz_names[-1])
                    ),
                )

            extraz_objective = dict.fromkeys(extraz_names, 1.0)

            for extraz_objective_sense in (-1, +1):
                baselp = add_objective_to_model(
                    baselp,
                    extraz_objective,
                    extraz_objective_sense,
                    f"postprocess_obj_{extraz_objective_sense}",
                    f"postprocess_obj_var_{extraz_objective_sense}",
                )
                getattr(
                    baselp, f"postprocess_obj_{extraz_objective_sense}"
                ).deactivate()

            for extraz_objective_sense in (-1, +1):
                getattr(baselp, f"postprocess_obj_{extraz_objective_sense}").activate()
                try:
                    results = pyomo_lp_solver.solve(
                        baselp, tee=verbose, **lp_solver.solve_extra_options
                    )
                except (ApplicationError, AttributeError, ValueError):
                    lp_resultdict = {ALL_OK_KEY: False}
                try:
                    lp_resultdict = add_statuses_to_optimziation_dict(
                        get_pyomo_solution_as_dict(baselp), results
                    )
                except (ApplicationError, AttributeError, ValueError):
                    lp_resultdict = {ALL_OK_KEY: False}
                getattr(
                    baselp, f"postprocess_obj_{extraz_objective_sense}"
                ).deactivate()
                if not lp_resultdict[ALL_OK_KEY]:
                    continue

                try:
                    nlp_resultdict = (
                        perform_nlp_irreversible_optimization_with_active_reacs_only(
                            cobrak_model=target_cobrak_model,
                            objective_target=objective_target,
                            objective_sense=objective_sense,
                            optimization_dict=deepcopy(lp_resultdict),
                            variability_dict=deepcopy(variability_data),
                            with_kappa=with_kappa,
                            with_gamma=with_gamma,
                            with_iota=with_iota,
                            with_alpha=with_alpha,
                            solver=nlp_solver,
                            correction_config=correction_config,
                        )
                    )
                    if nlp_resultdict[ALL_OK_KEY]:
                        feasible_switches.append(
                            [
                                target_couple,
                                tuple(
                                    extraz_var
                                    for extraz_var, value in lp_resultdict.items()
                                    if extraz_var.startswith(EXTRAZ_PREFIX)
                                    and value > 0.1
                                ),
                                nlp_resultdict,
                                f"{target_type}{at_maximum}",
                                max_allowed_changes,
                            ],
                        )
                except (ApplicationError, AttributeError, ValueError):
                    pass
    return feasible_switches


def postprocess(
    cobrak_model: Model,
    opt_dict: dict[str, float],
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_data: dict[str, tuple[float, float]],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    lp_solver: Solver = SCIP,
    nlp_solver: Solver = IPOPT,
    verbose: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
    onlytested: str = "",
) -> tuple[float, list[float | int]]:
    if variability_data == {}:
        variability_data = perform_lp_variability_analysis(
            cobrak_model=cobrak_model,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            active_reactions=[],
            solver=lp_solver,
        )
    else:
        variability_data = deepcopy(variability_data)

    pyomo_lp_solver = get_solver(
        lp_solver.name,
        lp_solver.solver_options,
        lp_solver.solver_attrs,
    )

    cobrak_model = deepcopy(cobrak_model)
    if type(objective_target) is str:
        objective_target = {objective_target: 1.0}
    obj_value = 0.0
    for obj_target_id, obj_target_multiplier in objective_target.items():
        obj_value += opt_dict[obj_target_id] * obj_target_multiplier
    epsilon = 1e-6 if is_objsense_maximization(objective_sense) else -1e-6
    cobrak_model = add_objective_value_as_extra_linear_constraint(
        cobrak_model,
        obj_value + epsilon,
        objective_target,
        objective_sense,
    )

    reac_couples = get_stoichiometrically_coupled_reactions(cobrak_model)
    active_reacs = [
        active_reac
        for active_reac in get_active_reacs_from_optimization_dict(
            cobrak_model, opt_dict
        )
        if (active_reac in [reac_ids[0] for reac_ids in reac_couples])
        and (active_reac not in objective_target)  # and opt_dict[active_reac] > 1e-8
    ]
    active_reac_couples: list[list[str]] = [
        reac_couple
        for reac_couple in reac_couples
        if reac_couple[0] in active_reacs and variability_data[reac_couple[0]][0] == 0.0
    ]
    inactive_reac_couples: list[list[str]] = [
        reac_couple
        for reac_couple in reac_couples
        if reac_couple[0] not in active_reacs
        and variability_data[reac_couple[0]][1] > 0.0
    ]
    targets = []
    for max_target_num in (0, 5):
        targets += [("deac", x, max_target_num) for x in active_reac_couples]
        targets += [("ac", x, max_target_num) for x in inactive_reac_couples]
    all_feasible_switches_metalist = Parallel(n_jobs=-1, verbose=10)(
        delayed(_postprocess_batch)(
            reac_couples,
            targets_batch,
            active_reacs,
            cobrak_model,
            objective_target,
            objective_sense,
            variability_data,
            pyomo_lp_solver,
            with_kappa,
            with_gamma,
            with_iota,
            with_alpha,
            lp_solver,
            nlp_solver,
            verbose,
            correction_config,
            onlytested,
        )
        for targets_batch in split_list(targets, cpu_count())
    )
    all_feasible_switches = []
    for sublist in all_feasible_switches_metalist:
        all_feasible_switches.extend(sublist)

    if len(all_feasible_switches) > 0:
        best_result = all_feasible_switches[0][2]
        for result in [x[2] for x in all_feasible_switches[1:]]:
            if is_objsense_maximization(objective_sense):
                if result[OBJECTIVE_VAR_NAME] > best_result[OBJECTIVE_VAR_NAME]:
                    best_result = result
            else:
                if result[OBJECTIVE_VAR_NAME] < best_result[OBJECTIVE_VAR_NAME]:
                    best_result = result
    else:
        best_result = {}

    return all_feasible_switches, best_result
    ####


def _sampling_routine(
    cobrak_model: Model,
    objective_target: str,
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    deactivated_reaction_lists: list[list[str]],
    lp_solver: Solver,
    nlp_solver: Solver,
    correction_config: CorrectionConfig,
    min_abs_objvalue: float,
) -> list[dict[str, float]]:
    working_dicts: list[dict[str, float]] = []
    for deactivated_reaction_set in deactivated_reaction_lists:
        try:
            ectfba_dict = perform_lp_optimization(
                cobrak_model=cobrak_model,
                objective_target=objective_target,
                objective_sense=objective_sense,
                with_enzyme_constraints=True,
                with_thermodynamic_constraints=True,
                with_loop_constraints=False,
                variability_dict=variability_dict,
                ignored_reacs=deactivated_reaction_set,
                solver=lp_solver,
                correction_config=correction_config,
            )
        except (ApplicationError, AttributeError, ValueError):
            continue

        if not ectfba_dict[ALL_OK_KEY]:
            continue

        active_ectfba_reacs = get_active_reacs_from_optimization_dict(
            cobrak_model, ectfba_dict
        )
        missing_errortarget = False
        for errortarget in correction_config.error_scenario:
            if errortarget not in cobrak_model.reactions:
                continue
            if errortarget not in active_ectfba_reacs:
                missing_errortarget = True
        if missing_errortarget:
            continue

        ###################################################
        # 3. NLP for initial set of rmax guesses
        try:
            nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
                cobrak_model=cobrak_model,
                objective_target=objective_target,
                objective_sense=objective_sense,
                optimization_dict=ectfba_dict,
                variability_dict=variability_dict,
                with_kappa=with_kappa,
                with_gamma=with_gamma,
                with_iota=with_iota,
                with_alpha=with_alpha,
                verbose=False,
                solver=nlp_solver,
                correction_config=correction_config,
            )
        except (ApplicationError, AttributeError, ValueError):
            continue
        if not nlp_result[ALL_OK_KEY]:
            continue
        if abs(nlp_result[OBJECTIVE_VAR_NAME]) < min_abs_objvalue:
            continue
        if nlp_result[ALL_OK_KEY]:
            working_dicts.append(deepcopy(nlp_result))
            print("Working sampling result:", nlp_result[OBJECTIVE_VAR_NAME])
    return working_dicts


def perform_nlp_evolutionary_optimization(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]] = {},
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    sampling_wished_num_feasible_starts: int = 3,
    sampling_max_metarounds: int = 3,
    sampling_rounds_per_metaround: int = 2,
    sampling_max_deactivated_reactions: int = 5,
    sampling_always_deactivated_reactions: list[str] = [],
    evolution_num_gens: int = 5,
    algorithm: Literal["pso", "genetic"] = "genetic",
    lp_solver: Solver = SCIP,
    nlp_solver: Solver = IPOPT,
    objvalue_json_path: str = "",
    max_rounds_same_objvalue: float = float("inf"),
    correction_config: CorrectionConfig = CorrectionConfig(),
    min_abs_objvalue: float = 1e-13,
    pop_size: int | None = None,
    working_results: list[dict[str, float]] = [],
) -> dict[float, list[dict[str, float]]]:
    if variability_dict == {}:
        variability_dict = perform_lp_variability_analysis(
            cobrak_model=cobrak_model,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            active_reactions=[],
            solver=lp_solver,
        )
    else:
        variability_dict = deepcopy(variability_dict)

    # Initial sampling
    if isinstance(objective_target, str):
        objective_target = {objective_target: 1.0}
    objective_target_ids = list(objective_target.keys())  # type: ignore

    deactivatable_reactions = [
        var_id
        for var_id in variability_dict
        if (var_id in cobrak_model.reactions)
        and (variability_dict[var_id][0] == 0.0)
        and (var_id not in objective_target_ids)
        and (var_id not in sampling_always_deactivated_reactions)
        and (var_id not in correction_config.error_scenario)
    ]
    distinct_feasible_start_solutions: dict[tuple[str, ...], dict[str, float]] = {}
    for current_round in range(sampling_max_metarounds):
        # Get deactivated reaction lists
        all_deactivated_reaction_lists = [
            [
                sample(
                    deactivatable_reactions,
                    randint(1, sampling_max_deactivated_reactions + 1),
                )
                + sampling_always_deactivated_reactions
                for _ in range(sampling_rounds_per_metaround)
            ]
            for _ in range(cpu_count())
        ]
        if current_round == 0:
            all_deactivated_reaction_lists[0][0] = deepcopy(
                sampling_always_deactivated_reactions
            )

        # run sampling
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(_sampling_routine)(
                cobrak_model,
                objective_target,
                objective_sense,
                variability_dict,
                with_kappa,
                with_gamma,
                with_iota,
                with_alpha,
                deactivated_reaction_lists,
                lp_solver,
                nlp_solver,
                correction_config,
                min_abs_objvalue,
            )
            for deactivated_reaction_lists in all_deactivated_reaction_lists
        )
        if len(working_results) > 0:
            results.append(working_results)
        best_result = (
            -float("inf") if is_objsense_maximization(objective_sense) else float("inf")
        )
        for result in results:
            for tfba_dict in result:
                active_reacs_tuple = tuple(
                    sorted(
                        get_active_reacs_from_optimization_dict(cobrak_model, tfba_dict)
                    )
                )
                distinct_feasible_start_solutions[active_reacs_tuple] = deepcopy(
                    tfba_dict
                )
                if is_objsense_maximization(objective_sense):
                    best_result = max(tfba_dict[OBJECTIVE_VAR_NAME], best_result)
                else:
                    best_result = min(tfba_dict[OBJECTIVE_VAR_NAME], best_result)

        if (
            len(distinct_feasible_start_solutions.keys())
            >= sampling_wished_num_feasible_starts
        ):
            break

    print(best_result)

    if len(distinct_feasible_start_solutions.keys()) == 0:
        print(
            "ERROR: No feasible sampling solution found! Check feasibility of problem and/or adjust sampling settings."
        )
        raise ValueError
    if (
        len(distinct_feasible_start_solutions.keys())
        < sampling_wished_num_feasible_starts
    ):
        print("INFO: Fewer feasible sampling solutions found than wished.")

    # Evolutionary algorithm
    problem = COBRAKProblem(
        cobrak_model=cobrak_model,
        objective_target=objective_target,  # type: ignore
        objective_sense=objective_sense,
        variability_dict=variability_dict,
        nlp_dict_list=list(distinct_feasible_start_solutions.values()),
        best_value=best_result,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        num_gens=evolution_num_gens,
        algorithm=algorithm,
        lp_solver=lp_solver,
        nlp_solver=nlp_solver,
        objvalue_json_path=objvalue_json_path,
        max_rounds_same_objvalue=max_rounds_same_objvalue,
        correction_config=correction_config,
        min_abs_objvalue=min_abs_objvalue,
        pop_size=pop_size,
    )

    return problem.optimize()
