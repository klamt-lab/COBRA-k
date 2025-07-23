"""
This file contains all non-linear programs (NLP) functions, including the evolutionary NLP optimization algorithm,
that can be used with COBRAk models.
With NLPs, all types of constraints (stoichiomnetric, enzymatic, κ, γ, ι, ...) can be integrated.
However, NLPs can be very slow.
For linear-programs (LP) and mixed-integer linear programs (MILP),
see lps.py in the same folder.
"""

# IMPORTS SECTION #
from copy import deepcopy
from itertools import chain

from joblib import Parallel, delayed
from pyomo.environ import (
    Binary,
    ConcreteModel,
    Constraint,
    Objective,
    Reals,
    Var,
    exp,
    maximize,
    minimize,
)

from .constants import (
    ALL_OK_KEY,
    ALPHA_VAR_PREFIX,
    BIG_M,
    ENZYME_VAR_PREFIX,
    ERROR_VAR_PREFIX,
    GAMMA_VAR_PREFIX,
    IOTA_VAR_PREFIX,
    KAPPA_VAR_PREFIX,
    LNCONC_VAR_PREFIX,
    MDF_VAR_ID,
    OBJECTIVE_VAR_NAME,
    SOLVER_STATUS_KEY,
    STANDARD_MIN_MDF,
    TERMINATION_CONDITION_KEY,
    Z_VAR_PREFIX,
)
from .dataclasses import CorrectionConfig, Model, Solver
from .lps import (
    _add_concentration_vars_and_constraints,
    _add_df_and_dG0_var_for_reaction,
    _add_error_sum_to_model,
    _add_extra_linear_constraints_to_lp,
    _add_kappa_substrates_and_products_vars,
    _apply_error_scenario,
    _get_dG0_highbound,
    _get_km_bounds,
    get_lp_from_cobrak_model,
)
from .pyomo_functionality import get_model_var_names, get_objective, get_solver
from .standard_solvers import IPOPT, SCIP
from .utilities import (
    add_statuses_to_optimziation_dict,
    apply_variability_dict,
    delete_unused_reactions_in_optimization_dict,
    get_full_enzyme_id,
    get_pyomo_solution_as_dict,
    get_reaction_enzyme_var_id,
    get_stoichiometrically_coupled_reactions,
    have_all_unignored_km,
    is_any_error_term_active,
    split_list,
)


# FUNCTIONS SECTION #
def add_loop_constraints_to_nlp(
    model: ConcreteModel,
    cobrak_model: Model,
) -> ConcreteModel:
    """Adds loop constraints to a non-linear program (NLP) model.

    The loop constraints are of the nonlinear form v_fwd * v_rev = 0.0
    for any forward/reverse pair of split reversible reactions.

    Parameters
    * `model` (`ConcreteModel`): The NLP model to add constraints to.
    * `cobrak_model` (`Model`): The COBRAk model associated with the NLP model.

    Returns
    * `ConcreteModel`: The NLP model with added loop constraints.
    """
    model_var_names = [v.name for v in model.component_objects(Var)]
    for reac_id, reaction in cobrak_model.reactions.items():
        if reaction.dG0 is not None:
            continue
        if not reac_id.endswith(cobrak_model.rev_suffix):
            continue
        other_reac_id = reac_id.replace(
            cobrak_model.rev_suffix, cobrak_model.fwd_suffix
        )
        if other_reac_id not in model_var_names:
            continue

        setattr(
            model,
            f"loop_constraint_{reac_id}",
            Constraint(
                rule=getattr(model, reac_id) * getattr(model, other_reac_id) == 0.0
            ),
        )

    return model


def get_nlp_from_cobrak_model(
    cobrak_model: Model,
    ignored_reacs: list[str] = [],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    approximation_value: float = 0.0001,
    irreversible_mode: bool = False,
    variability_data: dict[str, tuple[float, float]] = {},
    strict_mode: bool = False,
    irreversible_mode_min_mdf: float = STANDARD_MIN_MDF,
    with_flux_sum_var: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
) -> ConcreteModel:
    """Creates a pyomo non-linear program (NLP) model instance from a COBRAk Model.

    For more, see COBRAk's NLP documentation chapter.

    # Parameters
    * `cobrak_model` (`Model`): The COBRAk model to create the NLP model from.
    * `ignored_reacs` (`list[str]`, optional): List of reaction IDs to ignore. Defaults to `[]`.
    * `with_kappa` (`bool`, optional): Whether to include κ saturation term terms. Defaults to `True`.
    * `with_gamma` (`bool`, optional): Whether to include γ thermodynamic terms. Defaults to `True`.
    * `with_iota` (`bool`, optional): Whether to include ι inhibition terms. Defaults to `False` and untested!
    * `with_alpha` (`bool`, optional): Whether to include α activation terms. Defaults to `False` and untested!
    * `approximation_value` (`float`, optional): Approximation value for κ, γ, ι, and α terms. Defaults to `0.0001`. This value is the
       minimal value for κ, γ, ι, and α terms, and can lead to an overapproximation in this regard.
    * `irreversible_mode` (`bool`, optional): Whether to use irreversible mode. Defaults to `False`.
    * `variability_data` (`dict[str, tuple[float, float]]`, optional): Variability data for reactions. Defaults to `{}`.
    * `strict_mode` (`bool`, optional): Whether to use strict mode. Defaults to `False`.
    * `irreversible_mode_min_mdf` (`float`, optional): Minimum MDF value for irreversible mode. Defaults to `STANDARD_MIN_MDF`.
    * `with_flux_sum_var` (`bool`, optional): Whether to include a flux sum variable of name ```cobrak.constants.FLUX_SUM_VAR```. Defaults to `False`.
    * `correction_config` (`CorrectionConfig`, optional): Parameter correction configuration. Defaults to `CorrectionConfig()`.

    # Returns
    * `ConcreteModel`: The created NLP model.
    """
    cobrak_model = deepcopy(cobrak_model)

    reac_ids = list(cobrak_model.reactions.keys())
    enforced_reacs: list[str] = []
    ignored_reacs = deepcopy(ignored_reacs)
    for reac_id in variability_data:
        if reac_id not in reac_ids:
            continue
        min_flux = variability_data[reac_id][0]
        if min_flux < 1e-6:
            continue
        enforced_reacs.append(reac_id)

        if reac_id.endswith("_REV"):
            other_id = reac_id.replace("_REV", "_FWD")
        elif reac_id.endswith("_FWD"):
            other_id = reac_id.replace("_FWD", "_REV")
        else:
            continue
        if other_id in reac_ids:
            ignored_reacs.append(other_id)

    model = get_lp_from_cobrak_model(
        cobrak_model=cobrak_model,
        ignored_reacs=ignored_reacs,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=False,
        with_loop_constraints=False,
        add_extra_linear_constraints=False,
        with_flux_sum_var=with_flux_sum_var,
        correction_config=CorrectionConfig(
            add_kcat_times_e_error_term=correction_config.add_kcat_times_e_error_term,
            kcat_times_e_error_cutoff=correction_config.kcat_times_e_error_cutoff,
            max_rel_kcat_times_e_correction=correction_config.max_rel_kcat_times_e_correction,
            add_error_sum_term=False,
        ),
    )
    model = _add_concentration_vars_and_constraints(model, cobrak_model)

    if correction_config.add_kcat_times_e_error_term:
        model_vars = get_model_var_names(model)

    if correction_config.add_km_error_term:
        kms_lowbound, kms_highbound = _get_km_bounds(
            cobrak_model, correction_config.km_error_cutoff
        )
    else:
        kms_lowbound, kms_highbound = 0.0, 0.0

    if correction_config.add_dG0_error_term:
        dG0_highbound = _get_dG0_highbound(
            cobrak_model, correction_config.dG0_error_cutoff
        )
    else:
        dG0_highbound = 0.0

    setattr(
        model,
        MDF_VAR_ID,
        Var(within=Reals, bounds=(irreversible_mode_min_mdf, 1_000_000)),
    )
    # Set "MM" constraints
    if not irreversible_mode:
        reaction_couples = get_stoichiometrically_coupled_reactions(
            cobrak_model=cobrak_model,
        )
        reac_id_to_reac_couple_id: dict[str, str] = {}
        for couple in reaction_couples:
            for reac_id in couple:
                reac_id_to_reac_couple_id[reac_id] = "".join(couple)
        created_z_vars = []

    for reac_id, reaction in cobrak_model.reactions.items():
        if reac_id in ignored_reacs:
            continue

        if reaction.dG0 is not None:
            model, f_var_name = _add_df_and_dG0_var_for_reaction(
                model,
                reac_id,
                reaction,
                cobrak_model,
                strict_df_equality=strict_mode,
                add_error_term=correction_config.add_dG0_error_term
                and (reaction.dG0 >= dG0_highbound),
                max_abs_dG0_correction=correction_config.max_abs_dG0_correction,
            )

            if (
                not irreversible_mode
                and variability_data[reac_id][0] == 0.0
                and variability_data[reac_id][1] != 0.0
            ):
                z_varname = f"{Z_VAR_PREFIX}{reac_id_to_reac_couple_id[reac_id]}"
                if z_varname not in created_z_vars:
                    setattr(model, z_varname, Var(within=Binary))
                    created_z_vars.append(z_varname)

                # Big-M 0: r_i <= lb * z_i
                bigm_optmdfpathway_0_constraint = getattr(
                    model, reac_id
                ) <= reaction.max_flux * getattr(model, z_varname)
                setattr(
                    model,
                    f"bigm_optmdfpathway_0_{reac_id}",
                    Constraint(rule=bigm_optmdfpathway_0_constraint),
                )

                # Big-M 1: f_i + (1-z_i) * M_i >= var_B
                bigm_optmdfpathway_1_constraint = getattr(model, f_var_name) + (
                    1 - getattr(model, z_varname)
                ) * BIG_M >= getattr(model, MDF_VAR_ID)

                setattr(
                    model,
                    f"bigm_optmdfpathway_1_{reac_id}",
                    Constraint(rule=bigm_optmdfpathway_1_constraint),
                )
            elif variability_data[reac_id][1] != 0.0:
                mdf_constraint = getattr(model, f_var_name) >= getattr(
                    model, MDF_VAR_ID
                )

                setattr(
                    model,
                    f"mdf_constraint_{reac_id}",
                    Constraint(rule=mdf_constraint),
                )

        if reaction.enzyme_reaction_data is None:
            continue

        # Determine whether or not κ, γ, ι and α are possible to add to the reaction
        # given its current kinetic and thermodynamic data.
        has_gamma = True
        has_kappa = True
        if not have_all_unignored_km(
            reaction, cobrak_model.kinetic_ignored_metabolites
        ):
            has_kappa = False
        if reaction.dG0 is None:
            has_gamma = False
        if (not has_kappa) and (not has_gamma):
            continue
        has_iota = reaction.enzyme_reaction_data.k_is != {}
        has_alpha = reaction.enzyme_reaction_data.k_as != {}

        reac_full_enzyme_id = get_full_enzyme_id(
            reaction.enzyme_reaction_data.identifiers
        )
        if not reac_full_enzyme_id:  # E.g., in ATPM
            continue
        enzyme_var_id = get_reaction_enzyme_var_id(reac_id, reaction)

        # V+
        k_cat = reaction.enzyme_reaction_data.k_cat

        if correction_config.add_kcat_times_e_error_term:
            kcat_times_e_error_var_id = f"{ERROR_VAR_PREFIX}_kcat_times_e_{reac_id}"
            if kcat_times_e_error_var_id in model_vars:
                v_plus = getattr(model, enzyme_var_id) * k_cat + getattr(
                    model, kcat_times_e_error_var_id
                )
            else:
                v_plus = getattr(model, enzyme_var_id) * k_cat
        else:
            v_plus = getattr(model, enzyme_var_id) * k_cat

        # κ (for solver stability, with a minimal value of 0.0001)
        if has_kappa and with_kappa:
            model, kappa_substrates_var_id, kappa_products_var_id = (
                _add_kappa_substrates_and_products_vars(
                    model,
                    reac_id,
                    reaction,
                    cobrak_model,
                    strict_kappa_products_equality=strict_mode,
                    add_error_term=correction_config.add_km_error_term,
                    max_rel_km_correction=correction_config.max_rel_km_correction,
                    kms_lowbound=kms_lowbound,
                    kms_highbound=kms_highbound,
                )
            )

            kappa_var_id = f"{KAPPA_VAR_PREFIX}{reac_id}"
            setattr(
                model,
                kappa_var_id,
                Var(within=Reals, bounds=(approximation_value, 1.0)),
            )
            kappa_rhs = approximation_value + exp(
                getattr(model, kappa_substrates_var_id)
            ) / (
                1
                + exp(getattr(model, kappa_substrates_var_id))
                + exp(getattr(model, kappa_products_var_id))
            )
            if strict_mode:
                kappa_constraint = getattr(model, kappa_var_id) == kappa_rhs
            else:
                kappa_constraint = getattr(model, kappa_var_id) <= kappa_rhs
            setattr(
                model, f"kappa_constraint_{reac_id}", Constraint(rule=kappa_constraint)
            )

        # γ (for solver stability, with a minimal value of 0.0001)
        if has_gamma and with_gamma:
            gamma_var_name = f"{GAMMA_VAR_PREFIX}{reac_id}"

            min_gamma_value = (
                approximation_value if irreversible_mode else -float("inf")
            )
            setattr(
                model,
                gamma_var_name,
                Var(within=Reals, bounds=(min_gamma_value, 1.0)),
            )
            f_by_RT = getattr(model, f_var_name) / (cobrak_model.R * cobrak_model.T)

            if irreversible_mode:
                gamma_rhs = approximation_value + (1 - exp(-f_by_RT))
            else:
                gamma_rhs = (
                    approximation_value
                    + (
                        1
                        - exp(
                            -f_by_RT
                        )  # * getattr(model, f"{Z_VAR_PREFIX}{reac_id_to_reac_couple_id[reac_id]}")
                    )
                )  # (f_by_RT**2) / (1 + (f_by_RT**2)) would be a rough approximation

            if strict_mode:
                gamma_var_constraint_0 = getattr(model, gamma_var_name) == gamma_rhs
            else:
                gamma_var_constraint_0 = getattr(model, gamma_var_name) <= gamma_rhs
            setattr(
                model,
                f"gamma_var_constraint_{reac_id}_0",
                Constraint(rule=gamma_var_constraint_0),
            )

        # ι (for solver stability, with a minimal value of 0.0001)
        if with_iota and has_iota:
            iota_product = 1.0
            for met_id, k_i in reaction.enzyme_reaction_data.k_is.items():
                if met_id in cobrak_model.kinetic_ignored_metabolites:
                    continue
                iota_product *= 1 / (
                    1 + exp(getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}")) / k_i
                )

            iota_var_name = f"{IOTA_VAR_PREFIX}{reac_id}"
            setattr(
                model,
                iota_var_name,
                Var(within=Reals, bounds=(approximation_value, 1.0)),
            )
            if strict_mode:
                iota_var_constraint_0 = (
                    getattr(model, iota_var_name) == approximation_value + iota_product
                )
            else:
                iota_var_constraint_0 = (
                    getattr(model, iota_var_name) <= approximation_value + iota_product
                )
            setattr(
                model,
                f"iota_var_constraint_{reac_id}_0",
                Constraint(rule=iota_var_constraint_0),
            )

        # α (for solver stability, with a minimal value of 0.0001)
        if with_alpha and has_alpha:
            alpha_product = 1.0
            for met_id, k_a in reaction.enzyme_reaction_data.k_as.items():
                if met_id in cobrak_model.kinetic_ignored_metabolites:
                    continue
                # alpha_product *= 1 + k_a * exp(
                #    getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}")
                # )
                # alpha_product *= ((exp(getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}"))) / (k_a + exp(getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}"))))
                alpha_product *= 1 / (
                    1 + (k_a / exp(getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}")))
                )
                # alpha_product *= (1 - (k_a / (k_a + exp(getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}")))))

            alpha_var_name = f"{ALPHA_VAR_PREFIX}{reac_id}"
            setattr(
                model,
                alpha_var_name,
                Var(within=Reals, bounds=(approximation_value, 1.0)),
            )

            if strict_mode:
                alpha_var_constraint_0 = (
                    getattr(model, alpha_var_name)
                    == approximation_value + alpha_product
                )
            else:
                alpha_var_constraint_0 = (
                    getattr(model, alpha_var_name)
                    <= approximation_value + alpha_product
                )
            setattr(
                model,
                f"alpha_var_constraint_{reac_id}_0",
                Constraint(rule=alpha_var_constraint_0),
            )
        # Build kinetic term for reaction according to included parts
        kinetic_rhs = v_plus
        if has_kappa and with_kappa:
            kinetic_rhs *= getattr(model, kappa_var_id)
        if has_gamma and with_gamma:
            kinetic_rhs *= getattr(model, gamma_var_name)
        if has_iota and with_iota:
            kinetic_rhs *= getattr(model, iota_var_name)
        if has_alpha and with_alpha:
            kinetic_rhs *= getattr(model, alpha_var_name)

        # Apply strict mode
        if strict_mode:
            setattr(
                model,
                f"full_reac_constraint_{reac_id}",
                Constraint(rule=getattr(model, reac_id) == kinetic_rhs),
            )
        else:
            setattr(
                model,
                f"full_reac_constraint_{reac_id}",
                Constraint(rule=getattr(model, reac_id) <= kinetic_rhs),
            )

    model = _add_extra_linear_constraints_to_lp(model, cobrak_model)
    if is_any_error_term_active(correction_config):
        if correction_config.error_scenario != {}:
            _apply_error_scenario(
                model,
                cobrak_model,
                correction_config,
            )
        if correction_config.add_error_sum_term:
            model = _add_error_sum_to_model(
                model,
                cobrak_model,
                correction_config,
            )

    """
    if cobrak_model.max_conc_sum < float("inf"):
        model = _add_conc_sum_constraints(cobrak_model, model)
    """
    ########################
    if cobrak_model.max_conc_sum < float("inf"):
        met_sum_ids: list[str] = []
        for var_id in get_model_var_names(model):
            if not var_id.startswith(LNCONC_VAR_PREFIX):
                continue
            if not any(
                var_id.endswith(suffix)
                for suffix in cobrak_model.conc_sum_include_suffixes
            ):
                continue
            if any(
                var_id.replace(LNCONC_VAR_PREFIX, "").startswith(prefix)
                for prefix in cobrak_model.conc_sum_ignore_prefixes
            ):
                continue
            met_sum_ids.append(var_id)

        conc_sum_expr = 0.0
        for met_sum_id in met_sum_ids:
            met_id = met_sum_id[len(LNCONC_VAR_PREFIX) :]
            # exp_var_id = f"expvar_{met_sum_id}"
            # setattr(model, exp_var_id, Var(within=Reals, bounds=(1e-5, 1e6)),)
            # setattr(
            #     model,
            #     f"expvarconstraint_{met_sum_id}",
            #     Constraint(rule=getattr(model, exp_var_id) >= getattr(model, met_sum_id)),
            # )
            conc_sum_expr += exp(getattr(model, met_sum_id))
        setattr(
            model,
            "met_sum_var",
            Var(within=Reals, bounds=(1e-5, cobrak_model.max_conc_sum)),
        )
        setattr(
            model,
            "met_sum_constraint",
            Constraint(rule=conc_sum_expr <= getattr(model, "met_sum_var")),
        )
    ################

    return model


def perform_nlp_reversible_optimization(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    approximation_value: float = 0.0001,
    strict_mode: bool = False,
    verbose: bool = False,
    solver: Solver = SCIP,
    with_flux_sum_var: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
    show_variable_count: bool = False,
) -> dict[str, float]:
    """Performs a reversible MILP-based non-linear program (NLP) optimization on a COBRAk model.

    For more on the MINLP, see the COBRAk documentation's NLP chapter.

    #### Parameters
    * `cobrak_model` (`Model`): The COBRAk model to optimize.
    * `objective_target` (`str | dict[str, float]`): The objective target (reaction ID or dictionary of reaction IDs and coefficients).
    * `objective_sense` (`int`): The objective sense (1 for maximization, -1 for minimization).
    * `variability_dict` (`dict[str, tuple[float, float]]`): Dictionary of reaction IDs and their variability (lower and upper bounds).
    * `with_kappa` (`bool`, optional): Whether to include κ saturation terms. Defaults to `True`.
    * `with_gamma` (`bool`, optional): Whether to include γ thermodynamic terms. Defaults to `True`.
    * `with_iota` (`bool`, optional): Whether to include ι inhibition terms. Defaults to `False` and untested!
    * `with_alpha` (`bool`, optional): Whether to include α activation terms. Defaults to `False` and untested!
    * `approximation_value` (`float`, optional): Approximation value for κ, γ, ι, and α terms. Defaults to `0.0001`. This value is the
       minimal value for κ, γ, ι, and α terms, and can lead to an overapproximation in this regard.
    * `strict_mode` (`bool`, optional): Whether to use strict mode. Defaults to `False`.
    * `verbose` (`bool`, optional): Whether to print solver output. Defaults to `False`.
    * `solver_name` (`str`, optional): Used MINLP solver. Defaults to SCIP,
    * `with_flux_sum_var` (`bool`, optional): Whether to include a reaction flux sum variable of name ```cobrak.constants.FLUX_SUM_VAR```. Defaults to `False`.
    * `correction_config` (`CorrectionConfig`, optional): Parameter correction configuration. Defaults to `CorrectionConfig()`.

    #### Returns
    * `dict[str, float]`: The optimization results.
    """
    nlp_model = get_nlp_from_cobrak_model(
        cobrak_model,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        approximation_value=approximation_value,
        irreversible_mode=False,
        variability_data=variability_dict,
        strict_mode=strict_mode,
        with_flux_sum_var=with_flux_sum_var,
        correction_config=correction_config,
    )

    nlp_model = apply_variability_dict(
        nlp_model, cobrak_model, variability_dict, correction_config.error_scenario
    )
    nlp_model.obj = get_objective(nlp_model, objective_target, objective_sense)
    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)

    if show_variable_count:
        float_vars = [v for v in nlp_model.component_objects(Var) if v.domain == Reals]
        num_float_vars = sum(1 for v in float_vars for i in v)
        binary_vars = [
            v for v in nlp_model.component_objects(Var) if v.domain == Binary
        ]
        num_binary_vars = sum(1 for v in binary_vars for i in v)
        print("# FLOAT VARS:", num_float_vars)
        print("# BINARY VARS:", num_binary_vars)

    results = pyomo_solver.solve(nlp_model, tee=verbose, **solver.solve_extra_options)

    nlp_result = get_pyomo_solution_as_dict(nlp_model)
    return add_statuses_to_optimziation_dict(nlp_result, results)


def perform_nlp_irreversible_optimization(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    approximation_value: float = 0.0001,
    verbose: bool = False,
    strict_mode: bool = False,
    min_mdf: float = STANDARD_MIN_MDF,
    solver: Solver = IPOPT,
    multistart_num_iterations: int = 10,
    multistart_solver_name: str = "ipopt",
    min_flux: float = 0.0,
    with_flux_sum_var: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
) -> dict[str, float]:
    """Performs an irreversible non-linear program (NLP) optimization on a COBRAk model.

    For more about the NLP, see the COBRAk documentation's NLP chapter.

    # Parameters
    * `cobrak_model` (`Model`): The COBRAk model to optimize.
    * `objective_target` (`str | dict[str, float]`): The objective target (reaction ID or dictionary of reaction IDs and coefficients).
    * `objective_sense` (`int`): The objective sense (1 for maximization, -1 for minimization).
    * `variability_dict` (`dict[str, tuple[float, float]]`): Dictionary of reaction IDs and their variability (lower and upper bounds).
    * `with_kappa` (`bool`, optional): Whether to include κ saturation terms. Defaults to `True`.
    * `with_gamma` (`bool`, optional): Whether to include γ thermodynamic terms. Defaults to `True`.
    * `with_iota` (`bool`, optional): Whether to include ι inhibition terms. Defaults to `False` and untested!
    * `with_alpha` (`bool`, optional): Whether to include α activation terms. Defaults to `False` and untested!
    * `approximation_value` (`float`, optional): Approximation value for κ, γ, ι, and α terms. Defaults to `0.0001`. This value is the
       minimal value for κ, γ, ι, and α terms, and can lead to an overapproximation in this regard.
    * `verbose` (`bool`, optional): Whether to print solver output. Defaults to `False`.
    * `strict_mode` (`bool`, optional): Whether to use strict mode. Defaults to `False`.
    * `min_mdf` (`float`, optional): Minimum MDF value. Defaults to `STANDARD_MIN_MDF`.
    * `solver_name` (Solver, optional): Used NLP solver. Defaults to IPOPT.
    * `multistart_num_iterations` (`int`, optional): Number of iterations for multistart solver. Defaults to `10`.
    * `multistart_solver_name` (`str`, optional): Solver name for multistart. Defaults to `"ipopt"`.
    * `min_flux` (`float`, optional): Minimum flux value. Defaults to `0.0`.
    * `with_flux_sum_var` (`bool`, optional): Whether to include a reaction flux sum variable of name ```cobrak.constants.FLUX_SUM_VAR```. Defaults to `False`.
    * `correction_config` (`CorrectionConfig`, optional): Parameter correction configuration. Defaults to `CorrectionConfig()`.

    # Returns
    * `dict[str, float]`: The optimization results.
    """
    nlp_model = get_nlp_from_cobrak_model(
        cobrak_model,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        approximation_value=approximation_value,
        irreversible_mode=True,
        variability_data=variability_dict,
        strict_mode=strict_mode,
        irreversible_mode_min_mdf=min_mdf,
        with_flux_sum_var=with_flux_sum_var,
        correction_config=correction_config,
    )
    variability_dict = deepcopy(variability_dict)
    if min_flux != 0.0:
        for reac_id in cobrak_model.reactions:
            if (reac_id in variability_dict) and (
                (variability_dict[reac_id][0] == 0.0)
                and (variability_dict[reac_id][1] >= min_flux)
            ):
                variability_dict[reac_id] = (min_flux, variability_dict[reac_id][1])

    nlp_model = apply_variability_dict(
        nlp_model, cobrak_model, variability_dict, correction_config.error_scenario
    )
    nlp_model.obj = get_objective(nlp_model, objective_target, objective_sense)
    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)
    if solver.name.startswith("multistart"):
        all_results = pyomo_solver.solve(
            nlp_model,
            solver=multistart_solver_name,
            solver_args={"options": solver.solver_options},
            iterations=multistart_num_iterations,
        )
        mmtfba_dict = get_pyomo_solution_as_dict(nlp_model)
        if isinstance(all_results, tuple):  # cobrak_multistart
            results = all_results[0]
            objs = all_results[1]
            solver_statuses = all_results[2]
            termination_conditions = all_results[3]
            mmtfba_dict = add_statuses_to_optimziation_dict(mmtfba_dict, results)
            mmtfba_dict[OBJECTIVE_VAR_NAME + "_MULTISTART"] = objs
            mmtfba_dict[SOLVER_STATUS_KEY + "_MULTISTART"] = solver_statuses
            mmtfba_dict[TERMINATION_CONDITION_KEY + "_MULTISTART"] = (
                termination_conditions
            )
        else:
            results = all_results
            mmtfba_dict = add_statuses_to_optimziation_dict(mmtfba_dict, results)
        return mmtfba_dict
    results = pyomo_solver.solve(nlp_model, tee=verbose, **solver.solve_extra_options)
    mmtfba_dict = get_pyomo_solution_as_dict(nlp_model)
    return add_statuses_to_optimziation_dict(mmtfba_dict, results)


def perform_nlp_irreversible_optimization_with_active_reacs_only(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    optimization_dict: dict[str, float],
    variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    approximation_value: float = 0.0001,
    verbose: bool = False,
    strict_mode: bool = False,
    min_mdf: float = STANDARD_MIN_MDF,
    solver: Solver = IPOPT,
    multistart_num_iterations: int = 10,
    do_not_delete_with_z_var_one: bool = False,
    correction_config: CorrectionConfig = CorrectionConfig(),
) -> dict[str, float]:
    """Performs an irreversible non-linear program (NLP) optimization on a COBRAk model, considering only active reactions of the optimization dict.

    For more about the NLP, see the COBRAk documentation's NLP chapter.

    # Parameters
    * `cobrak_model` (`Model`): The COBRAk model to optimize.
    * `objective_target` (`str | dict[str, float]`): The objective target (reaction ID or dictionary of reaction IDs and coefficients).
    * `objective_sense` (`int`): The objective sense (1 for maximization, -1 for minimization).
    * `optimization_dict` (`dict[str, float]`): Dictionary of reaction IDs and their optimization values.
    * `variability_dict` (`dict[str, tuple[float, float]]`): Dictionary of reaction IDs and their variability (lower and upper bounds).
    * `with_kappa` (`bool`, optional): Whether to include κ terms. Defaults to `True`.
    * `with_gamma` (`bool`, optional): Whether to include γ terms. Defaults to `True`.
    * `with_iota` (`bool`, optional): Whether to include ι inhibition terms. Defaults to `False` and untested!
    * `with_alpha` (`bool`, optional): Whether to include α activation terms. Defaults to `False` and untested!
    * `approximation_value` (`float`, optional): Approximation value for κ, γ, ι, and α terms. Defaults to `0.0001`. This value is the
       minimal value for κ, γ, ι, and α terms, and can lead to an overapproximation in this regard.
    * `verbose` (`bool`, optional): Whether to print solver output. Defaults to `False`.
    * `strict_mode` (`bool`, optional): Whether to use strict mode. Defaults to `False`.
    * `min_mdf` (`float`, optional): Minimum MDF value. Defaults to `STANDARD_MIN_MDF`.
    * `solver` (Solver, optional): Used NLP solver. Defaults to IPOPT.
    * `multistart_num_iterations` (`int`, optional): Number of iterations for multistart solver. Defaults to `10`.
    * `do_not_delete_with_z_var_one` (`bool`, optional): Whether to delete reactions with associated Z variables (in the optimization dics) equal to one.
      Defaults to `False`.
    * `correction_config` (`CorrectionConfig`, optional): Paramter correction configuration. Defaults to `CorrectionConfig()`.

    # Returns
    * `dict[str, float]`: The optimization results.
    """
    nlp_cobrak_model = delete_unused_reactions_in_optimization_dict(
        cobrak_model=cobrak_model,
        optimization_dict=optimization_dict,
        do_not_delete_with_z_var_one=do_not_delete_with_z_var_one,
    )
    return perform_nlp_irreversible_optimization(
        cobrak_model=nlp_cobrak_model,
        objective_target=objective_target,
        objective_sense=objective_sense,
        variability_dict=variability_dict,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        approximation_value=approximation_value,
        verbose=verbose,
        strict_mode=strict_mode,
        min_mdf=min_mdf,
        solver=solver,
        multistart_num_iterations=multistart_num_iterations,
        correction_config=correction_config,
    )


def _batch_nlp_variability_optimization(
    batch: list[tuple[str, str]],
    cobrak_model: Model,
    with_kappa: bool,
    with_gamma: bool,
    with_iota: bool,
    with_alpha: bool,
    approximation_value: float,
    tfba_variability_dict: dict[str, float],
    strict_mode: bool,
    min_mdf: float,
    solver: Solver,
) -> list[tuple[bool, str, float | None]]:
    """Performs a batch of non-linear program (NLP) variability optimizations on a COBRAk model.

    # Parameters
    batch (list[tuple[str, str]]): List of tuples containing objective names and target IDs.
    cobrak_model (Model): The COBRAk model to optimize.
    with_kappa (bool): Whether to include κ terms.
    with_gamma (bool): Whether to include γ terms.
    with_iota (bool): Whether to include ι terms.
    with_alpha (bool): Whether to include α terms.
    approximation_value (float): Approximation value for κ, γ, ι, and α terms.
    tfba_variability_dict (dict[str, float]): Dictionary of reaction IDs and their TFBA variability.
    strict_mode (bool): Whether to use strict mode.
    min_mdf (float): Minimum MDF value.
    solver (Solver): Used NLP name.
    """
    i = 0
    resultslist: list[tuple[bool, str, float | None]] = []
    for objective_name, target_id in batch:
        i += 1
        objective_sense = +1 if objective_name.startswith("MAX_OBJ_") else -1
        try:
            results = perform_nlp_irreversible_optimization(
                deepcopy(cobrak_model),
                objective_target=target_id,
                objective_sense=objective_sense,
                with_kappa=with_kappa,
                with_gamma=with_gamma,
                with_iota=with_iota,
                with_alpha=with_alpha,
                approximation_value=approximation_value,
                variability_dict=deepcopy(tfba_variability_dict),
                strict_mode=strict_mode,
                min_mdf=min_mdf,
                solver=solver,
                verbose=False,
            )
        except Exception:
            print("EXCEPTION", objective_name)
            resultslist.append((objective_name.startswith("MIN_OBJ_"), target_id, None))
            continue
        result = None if not results[ALL_OK_KEY] else results[OBJECTIVE_VAR_NAME]
        resultslist.append((objective_name.startswith("MIN_OBJ_"), target_id, result))
    return resultslist


def perform_nlp_irreversible_variability_analysis_with_active_reacs_only(
    cobrak_model: Model,
    optimization_dict: dict[str, float],
    tfba_variability_dict: dict[str, tuple[float, float]],
    with_kappa: bool = True,
    with_gamma: bool = True,
    with_iota: bool = False,
    with_alpha: bool = False,
    active_reactions: list[str] = [],
    min_active_flux: float = 1e-5,
    calculate_reacs: bool = True,
    calculate_concs: bool = True,
    calculate_rest: bool = True,
    extra_tested_vars_max: list[str] = [],
    extra_tested_vars_min: list[str] = [],
    strict_mode: bool = False,
    min_mdf: float = STANDARD_MIN_MDF,
    min_flux_cutoff: float = 1e-8,
    solver: Solver = IPOPT,
    do_not_delete_with_z_var_one: bool = False,
    parallel_verbosity_level: int = 0,
    approximation_value: float = 0.0001,
) -> dict[str, tuple[float, float]]:
    """Performs an irreversible non-linear program (NLP) variability analysis on a COBRAk model, considering only active reactions.

    This function calculates the minimum and maximum values of reaction fluxes, metabolite concentrations, and other variables in the model,
    given a set of active reactions and a variability dictionary.
    It uses a combination of NLP optimizations and parallel processing to efficiently compute the variability of the model.

    # Parameters
    * `cobrak_model` (`Model`): The COBRAk model to analyze.
    * `optimization_dict` (`dict[str, float]`): Dictionary of reaction IDs and their optimization values.
    * `tfba_variability_dict` (`dict[str, tuple[float, float]]`): Dictionary of reaction IDs and their TFBA variability (lower and upper bounds).
    * `with_kappa` (`bool`, optional): Whether to include κ saturation terms. Defaults to `True`.
    * `with_gamma` (`bool`, optional): Whether to include γ thermodynamic terms. Defaults to `True`.
    * `with_iota` (`bool`, optional): Whether to include ι inhibition terms. Defaults to `False` and untested!
    * `with_alpha` (`bool`, optional): Whether to include α activation terms. Defaults to `False` and untested!
    * `active_reactions` (`list[str]`, optional): List of active reaction IDs. Defaults to `[]`.
    * `min_active_flux` (`float`, optional): Minimum flux value for active reactions. Defaults to `1e-5`.
    * `calculate_reacs` (`bool`, optional): Whether to calculate reaction flux variability. Defaults to `True`.
    * `calculate_concs` (`bool`, optional): Whether to calculate metabolite concentration variability. Defaults to `True`.
    * `calculate_rest` (`bool`, optional): Whether to calculate variability of other variables (e.g., enzyme delivery, κ, γ). Defaults to `True`.
    * `strict_mode` (`bool`, optional): Whether to use strict mode. Defaults to `False`.
    * `min_mdf` (`float`, optional): Minimum MDF value. Defaults to `STANDARD_MIN_MDF`.
    * `min_flux_cutoff` (`float`, optional): Minimum flux cutoff value. Defaults to `1e-8`.
    * `solver` (Solver, optional): Used NLP solver. Defaults to IPOPT.
    * `do_not_delete_with_z_var_one` (`bool`, optional): Whether to delete reactions with Z variable equal to one. Defaults to `False`.
    * `parallel_verbosity_level` (`int`, optional): Verbosity level for parallel processing. Defaults to `0`.
    * `approximation_value` (`float`, optional): Approximation value for κ, γ, ι, and α terms. Defaults to `0.0001`. This value is the
       minimal value for κ, γ, ι, and α terms, and can lead to an overapproximation in this regard.

    # Returns
    * `dict[str, tuple[float, float]]`: A dictionary of variable IDs and their variability (lower and upper bounds).
    """
    cobrak_model = deepcopy(cobrak_model)
    cobrak_model = delete_unused_reactions_in_optimization_dict(
        cobrak_model=cobrak_model,
        optimization_dict=optimization_dict,
        do_not_delete_with_z_var_one=do_not_delete_with_z_var_one,
    )

    for active_reaction in active_reactions:
        cobrak_model.reactions[active_reaction].min_flux = min_active_flux

    model: ConcreteModel = get_nlp_from_cobrak_model(
        cobrak_model=deepcopy(cobrak_model),
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        approximation_value=approximation_value,
        variability_data=deepcopy(tfba_variability_dict),
        strict_mode=strict_mode,
        irreversible_mode_min_mdf=min_mdf,
    )
    model_var_names = get_model_var_names(model)

    min_values: dict[str, float] = {}
    max_values: dict[str, float] = {}
    objective_targets: list[tuple[int, str]] = []

    """
    min_flux_sum_result = perform_nlp_irreversible_optimization(
        deepcopy(cobrak_model),
        objective_target=FLUX_SUM_VAR_ID,
        objective_sense=-1,
        with_kappa=with_kappa,
        with_gamma=with_gamma,
        with_iota=with_iota,
        with_alpha=with_alpha,
        approximation_value=approximation_value,
        variability_dict=deepcopy(tfba_variability_dict),
        strict_mode=strict_mode,
        min_mdf=min_mdf,
        with_flux_sum_var=True,
        solver=solver,
    )
    """

    if calculate_concs or calculate_rest:
        min_mdf_result = perform_nlp_irreversible_optimization(
            deepcopy(cobrak_model),
            objective_target=MDF_VAR_ID,
            objective_sense=-1,
            with_kappa=with_kappa,
            with_gamma=with_gamma,
            with_iota=with_iota,
            with_alpha=with_alpha,
            approximation_value=approximation_value,
            variability_dict=deepcopy(tfba_variability_dict),
            strict_mode=strict_mode,
            solver=solver,
        )

    if calculate_concs:
        for met_id, metabolite in cobrak_model.metabolites.items():
            met_var_name = f"{LNCONC_VAR_PREFIX}{met_id}"
            if met_var_name in model_var_names:
                min_mdf_conc = min_mdf_result[met_var_name]
                max_mdf_conc = min_mdf_result[met_var_name]
                if metabolite.log_min_conc in (min_mdf_conc, max_mdf_conc):
                    min_values[met_var_name] = metabolite.log_min_conc
                else:
                    objective_targets.append((-1, met_var_name))
                if metabolite.log_max_conc in (min_mdf_conc, max_mdf_conc):
                    max_values[met_var_name] = metabolite.log_max_conc
                else:
                    objective_targets.append((+1, met_var_name))

    for reac_id, reaction in cobrak_model.reactions.items():
        # min_flux_sum_flux = min_flux_sum_result[reac_id]
        if calculate_reacs:
            # if reaction.min_flux in (min_flux_sum_flux,):
            #    min_values[reac_id] = (
            #        reaction.min_flux if reaction.min_flux >= min_flux_cutoff else 0.0
            #    )
            # else:
            # if reaction.max_flux in (min_flux_sum_flux,):
            #    max_values[reac_id] = reaction.max_flux
            # else:
            objective_targets.extend(((-1, reac_id), (+1, reac_id)))

        if not calculate_rest:
            continue

        kappa_var_name = f"{KAPPA_VAR_PREFIX}{reac_id}"
        gamma_var_name = f"{GAMMA_VAR_PREFIX}{reac_id}"
        if kappa_var_name in model_var_names:
            objective_targets.extend(((-1, kappa_var_name), (+1, kappa_var_name)))
        if gamma_var_name in model_var_names:
            objective_targets.extend(((-1, gamma_var_name), (+1, gamma_var_name)))
        if reaction.enzyme_reaction_data is not None:
            full_enzyme_id = get_full_enzyme_id(
                reaction.enzyme_reaction_data.identifiers
            )
            if full_enzyme_id:
                enzyme_delivery_var_name = get_reaction_enzyme_var_id(reac_id, reaction)
                # if 0.0 in (min_flux_sum_flux,):
                #    min_values[enzyme_delivery_var_name] = 0.0
                # else:
                objective_targets.extend(
                    ((-1, enzyme_delivery_var_name), (+1, enzyme_delivery_var_name))
                )

    if len(extra_tested_vars_min) > 0:
        for extra_tested_var in extra_tested_vars_max:
            if extra_tested_var in model_var_names:
                objective_targets.append((-1, extra_tested_var))

    if len(extra_tested_vars_max) > 0:
        for extra_tested_var in extra_tested_vars_max:
            if extra_tested_var in model_var_names:
                objective_targets.append((+1, extra_tested_var))

    objectives_data: list[tuple[str, str]] = []
    for obj_sense, target_id in objective_targets:
        print(target_id)
        if obj_sense == -1:
            objective_name = f"MIN_OBJ_{target_id}"
            pyomo_sense = minimize
        else:
            objective_name = f"MAX_OBJ_{target_id}"
            pyomo_sense = maximize
        setattr(
            model,
            objective_name,
            Objective(expr=getattr(model, target_id), sense=pyomo_sense),
        )
        getattr(model, objective_name).deactivate()
        objectives_data.append((objective_name, target_id))

    objectives_data_batches = split_list(
        objectives_data, len(objectives_data)
    )  # cpu_count())

    results_list = Parallel(n_jobs=-1, verbose=parallel_verbosity_level)(
        delayed(_batch_nlp_variability_optimization)(
            batch,
            cobrak_model,
            with_kappa,
            with_gamma,
            with_iota,
            with_alpha,
            approximation_value,
            tfba_variability_dict,
            strict_mode,
            min_mdf,
            solver,
        )
        for batch in objectives_data_batches
    )
    for result in chain(*results_list):
        is_minimization = result[0]
        target_id = result[1]
        result_value = result[2]
        if is_minimization:
            min_values[target_id] = result_value
        else:
            max_values[target_id] = result_value

    for key, min_value in min_values.items():
        if (key in cobrak_model.reactions) or (
            key.startswith(ENZYME_VAR_PREFIX) and (min_value is not None)
        ):
            min_values[key] = min_value if min_value >= min_flux_cutoff else 0.0

    all_target_ids = sorted(
        set(
            list(min_values.keys())
            + list(max_values.keys())
            + [obj_target[1] for obj_target in objective_targets]
        )
    )
    all_target_ids = [x[1] for x in objectives_data]
    variability_dict: dict[str, tuple[float, float]] = {
        target_id: (min_values[target_id], max_values[target_id])
        for target_id in all_target_ids
    }

    return variability_dict
