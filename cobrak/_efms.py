from copy import deepcopy
from random import choice

from joblib import Parallel, cpu_count, delayed
from pyomo.environ import (
    Binary,
    ConcreteModel,
    Constraint,
    Objective,
    Reals,
    SolverFactoryClass,
    Var,
    maximize,
)

from .constants import BIG_M, STANDARD_MIN_MDF
from .dataclasses import Model, Solver
from .lps import (
    get_lp_from_cobrak_model,
    get_objective,
)
from .pyomo_functionality import get_solver
from .standard_solvers import SCIP
from .utilities import (
    delete_orphaned_metabolites_and_enzymes,
    delete_unused_reactions_in_variability_dict,
    get_stoichiometrically_coupled_reactions,
)

STANDARD_EFM_OBJ_NAME = "EFM_STANDARD_OBJECTIVE"


def _efm_batch_routine(
    pyomo_model: ConcreteModel,
    pyomo_solver: SolverFactoryClass,
    objective_name: str,
    cobrak_solver: Solver,
    objective_sense: int,
) -> None:
    getattr(pyomo_model, objective_name).activate()

    # OBJ ACTIVATE
    pyomo_solver.solve(
        pyomo_model, tee=True, warmstart=True, **cobrak_solver.solve_extra_options
    )

    if results.solver.termination_condition == TerminationCondition.maxTimeLimit:
        result = None
    else:
        result = value(getattr(model, target_id))

    getattr(pyomo_model, objective_name).deactivate()


def calculate_efms(
    cobrak_model: Model,
    with_inhomogenous_constraints: bool = False,
    with_enzyme_constraints: bool = False,
    with_thermodynamic_constraints: bool = False,
    solver: Solver = SCIP,
    variability_dict: dict[str, tuple[float, float]] = {},
    max_num_efms: int = 1_000_000_000,
    min_mdf: float = STANDARD_MIN_MDF,
    min_efm_flux_bound: float = 1e-5,
) -> set[tuple[str, ...]]:
    enforced_reacs = [
        reac_id
        for reac_id in variability_dict
        if variability_dict[reac_id][0] > 0.0 and reac_id in cobrak_model.reactions
    ]
    if not with_inhomogenous_constraints:
        for enforced_reac in enforced_reacs:
            cobrak_model.reactions[enforced_reac].min_flux = 0.0
            variability_dict[enforced_reac] = (0.0, variability_dict[enforced_reac][1])

    reduced_model = deepcopy(cobrak_model)
    if variability_dict != {}:
        reduced_model = delete_unused_reactions_in_variability_dict(
            reduced_model, variability_dict
        )
    reduced_model = delete_orphaned_metabolites_and_enzymes(reduced_model)

    all_couples = get_stoichiometrically_coupled_reactions(
        reduced_model,
    )
    eligible_couples: list[list[str]] = []
    for couple in all_couples:
        eligible = True
        for enforced_reac in enforced_reacs:
            if enforced_reac in couple:
                eligible = False
                break
        if eligible:
            eligible_couples.append(couple)

    pyomo_model = get_lp_from_cobrak_model(
        cobrak_model=reduced_model,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=with_thermodynamic_constraints,
        with_loop_constraints=True,
        min_mdf=min_mdf,
    )
    zc_var_sum = 0.0
    zc_var_ids: list[str] = []
    objective_names: list[str] = []
    for eligible_couple in eligible_couples:
        zc_var_id = "zC_var_" + "_".join(eligible_couple)
        zc_var_ids.append(zc_var_id)

        setattr(pyomo_model, zc_var_id, Var(within=Binary))
        setattr(
            pyomo_model,
            zc_var_id + "_constraint",
            Constraint(
                rule=getattr(pyomo_model, eligible_couple[0])
                <= BIG_M * getattr(pyomo_model, zc_var_id)
            ),
        )
        setattr(
            pyomo_model,
            zc_var_id + "_constraint_2",
            Constraint(
                rule=getattr(pyomo_model, eligible_couple[0])
                >= min_efm_flux_bound * getattr(pyomo_model, zc_var_id)
            ),
        )
        zc_var_sum += getattr(pyomo_model, zc_var_id)

        objective_names.append(f"objective_{eligible_couple[0]}")
        setattr(
            pyomo_model,
            objective_names[-1],
            Objective(expr=getattr(pyomo_model, eligible_couple[0]), sense=maximize),
        )
        getattr(pyomo_model, objective_names[-1]).deactivate()

    objective_names.append(f"objective_{eligible_couple[0]}")
    setattr(
        pyomo_model,
        STANDARD_EFM_OBJ_NAME,
        Objective(expr=0.0, sense=maximize),
    )
    getattr(pyomo_model, STANDARD_EFM_OBJ_NAME).deactivate()

    setattr(pyomo_model, "zc_var_sum", Var(within=Reals, bounds=(1.0, 10_000.0)))
    setattr(
        pyomo_model,
        "zc_var_sum_constraint_1",
        Constraint(rule=getattr(pyomo_model, "zc_var_sum") == zc_var_sum),
    )
    pyomo_model.obj = get_objective(pyomo_model, "zc_var_sum", -1)

    integer_cut_i = 0
    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)
    efms: set[tuple[str, ...]] = set()
    for i in range(max_num_efms):
        opt_targets = [STANDARD_EFM_OBJ_NAME] + [
            choice(objective_names) for _ in range(cpu_count() - 1)
        ]
        min_z_results = Parallel(n_jobs=-1, verbose=10)(
            delayed(_efm_batch_routine)(
                pyomo_model,
                pyomo_solver,
                opt_target,
                +1,
            )
            for opt_target in opt_targets
        )

        new_efms = min_z_results.difference(efms)
        if len(new_efms) == 0:
            break
        for new_efm in new_efms:
            new_integer_cut = 0.0
            new_integer_cut_sum = 0.0
            for zc_var_id in zc_var_ids:
                if getattr(pyomo_model, zc_var_id).value <= 0.1:
                    continue
                new_integer_cut += getattr(pyomo_model, zc_var_id)
                new_integer_cut_sum += 1
            setattr(
                pyomo_model,
                f"integer_cut_{integer_cut_i}",
                Constraint(expr=new_integer_cut <= new_integer_cut_sum - 1),
            )
            integer_cut_i += 1
        efms = efms.union(new_efms)

    return efms
