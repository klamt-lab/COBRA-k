"""COBRAk LPs and MILPs.

This file contains all linear program (LP) and mixed-integer linear program (MILP)
functions that can be used with COBRAk models.
With LP, one can integrate stoichiometric and enzymatic constraints.
With MILP, one can additionally integrate thermodynamic constraints.
For non-linear-programs (NLP), see nlps.py in the same folder.
"""

# IMPORT SECTION #
from copy import deepcopy
from itertools import chain
from math import ceil, floor
from typing import Any

from joblib import Parallel, cpu_count, delayed
from numpy import percentile
from pyomo.common.errors import ApplicationError
from pyomo.environ import (
    Binary,
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Reals,
    TerminationCondition,
    Var,
    exp,
    log,
    maximize,
    minimize,
    value,
)
from pyomo.opt.base.solvers import SolverFactoryClass

from cobrak.pyomo_functionality import add_linear_approximation_to_pyomo_model

from .constants import (
    ALL_OK_KEY,
    BIG_M,
    DF_VAR_PREFIX,
    DG0_VAR_PREFIX,
    ENZYME_VAR_INFIX,
    ENZYME_VAR_PREFIX,
    ERROR_BOUND_LOWER_CHANGE_PREFIX,
    ERROR_BOUND_UPPER_CHANGE_PREFIX,
    ERROR_CONSTRAINT_PREFIX,
    ERROR_SUM_VAR_ID,
    ERROR_VAR_PREFIX,
    FLUX_SUM_VAR_ID,
    KAPPA_PRODUCTS_VAR_PREFIX,
    KAPPA_SUBSTRATES_VAR_PREFIX,
    LNCONC_VAR_PREFIX,
    MDF_VAR_ID,
    PROT_POOL_MET_NAME,
    PROT_POOL_REAC_NAME,
    QUASI_INF,
    STANDARD_MIN_MDF,
    Z_VAR_PREFIX,
    ZB_VAR_PREFIX,
)
from .dataclasses import CorrectionConfig, Model, Reaction, Solver
from .pyomo_functionality import get_model_var_names, get_objective, get_solver
from .standard_solvers import SCIP
from .utilities import (
    add_statuses_to_optimziation_dict,
    apply_variability_dict,
    delete_unused_reactions_in_variability_dict,
    get_base_id,
    get_full_enzyme_id,
    get_full_enzyme_mw,
    get_model_dG0s,
    get_model_kms,
    get_model_max_kcat_times_e_values,
    get_potentially_active_reactions_in_variability_dict,
    get_pyomo_solution_as_dict,
    get_reaction_enzyme_var_id,
    get_reaction_string,
    have_all_unignored_km,
    is_any_error_term_active,
    split_list,
)


# "PRIVATE" FUNCTIONS SECTION #
def _add_concentration_vars_and_constraints(
    model: ConcreteModel,
    cobrak_model: Model,
) -> ConcreteModel:
    """Adds the logarithmized concentration variables and associated constraints.

    For every metabolite in the COBRAk model that is associated with at least one reaction with a ΔG'°,
    a logarithmized concentration variable is added with the name f'{LNCONC_VAR_PREFIX}{met_id}'.
    Additonally, as constraints, the mimimal and maximal logarithmized concentration are set.
    Optionally, if concentration ratio constraints are set in the COBRAk model, these constraints
    are also added.

    Args:
        model (ConcreteModel): The (N/MI)LP to which the variables and constraints shall be added.
        cobrak_model (Model): The COBRAk model from which the new variables and constraints are deduced.

    Returns:
        ConcreteModel: The (N/MI)LP with added concentration variables and constraints.
    """
    # Add metabolite logarithmized concentration variables
    for met_id, metabolite in cobrak_model.metabolites.items():
        # Any enzymatic pseudo-metabolites are not regarded thermodynamically
        if met_id.startswith((ENZYME_VAR_PREFIX, PROT_POOL_MET_NAME)):
            continue

        # Check if metabolite occurs in at least one reaction with ΔG'°
        is_met_in_reac_with_concentration_constraints = False
        for reaction in cobrak_model.reactions.values():
            has_kappa = (reaction.enzyme_reaction_data is not None) and (
                have_all_unignored_km(
                    reaction, cobrak_model.kinetic_ignored_metabolites
                )
            )
            has_gamma = reaction.dG0 is not None
            if (not has_kappa) and (not has_gamma):
                continue
            if met_id in reaction.stoichiometries:
                is_met_in_reac_with_concentration_constraints = True
                break
        if not is_met_in_reac_with_concentration_constraints:
            continue

        # Finally, add the metabolite with given concentration bounds
        met_var_id = f"{LNCONC_VAR_PREFIX}{met_id}"
        setattr(
            model,
            met_var_id,
            Var(
                within=Reals, bounds=(metabolite.log_min_conc, metabolite.log_max_conc)
            ),
        )

    return model


def _add_conc_sum_constraints(
    cobrak_model: Model,
    model: ConcreteModel,
) -> ConcreteModel:
    met_sum_ids: list[str] = []
    for var_id in get_model_var_names(model):
        if not var_id.startswith(LNCONC_VAR_PREFIX):
            continue
        if not any(
            var_id.endswith(suffix) for suffix in cobrak_model.conc_sum_include_suffixes
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
        add_linear_approximation_to_pyomo_model(
            model=model,
            y_function=exp,
            y_function_derivative=exp,
            x_reference_var_id=met_sum_id,
            new_y_var_name=f"exp_{met_sum_id}",
            min_x=getattr(model, met_sum_id).lb,
            max_x=getattr(model, met_sum_id).ub,
            max_rel_difference=cobrak_model.conc_sum_max_rel_error,
            min_abs_error=cobrak_model.conc_sum_min_abs_error,
        )
        conc_sum_expr += getattr(model, f"exp_{met_sum_id}")

    setattr(
        model,
        "met_sum_var",
        Var(within=Reals, bounds=(1e-12, cobrak_model.max_conc_sum)),
    )
    setattr(
        model,
        "met_sum_constraint",
        Constraint(rule=conc_sum_expr <= getattr(model, "met_sum_var")),
    )

    return model


def _add_df_and_dG0_var_for_reaction(
    model: ConcreteModel,
    reac_id: str,
    reaction: Reaction,
    cobrak_model: Model,
    strict_df_equality: bool,
    add_error_term: bool = False,
    max_abs_dG0_correction: float = QUASI_INF,
) -> tuple[ConcreteModel, str]:
    """Add driving force- and ΔG'°-related constraints and variables to pyomo model.

    The following variables are added for the reaction:
    * f'{DG0_VAR_PREFIX}{reac_id}' - Variables representing the reaction ΔG'°, from
        ΔG'°-uncertainty to ΔG'°+uncertainty
    * f'{DF_VAR_PREFIX}{reac_id}' - Variables representing the current driving force -
        either correctly (strict_df_equality==True) or this variable is in (-inf;actual driving force)

    Args:
        model (ConcreteModel): The pyomo model which shall be expanded
        reac_id (str): The reaction's ID
        reaction (Reaction): The reaction from which the thermodyamic properties are read
        cobrak_model (Model): The COBRAk model from which R and T are taken
        strict_df_equality (bool): Whether or not the driving force is exactly modeled (see above)

    Returns:
        tuple[ConcreteModel, str]: The expanded model, and the new driving force variable name
    """
    dG0_value = reaction.dG0
    dG0_uncertainty = (
        reaction.dG0_uncertainty if reaction.dG0_uncertainty is not None else 0.0
    )

    dG0_var_name = f"{DG0_VAR_PREFIX}{reac_id}"
    setattr(
        model,
        dG0_var_name,
        Var(
            within=Reals,
            bounds=(dG0_value - dG0_uncertainty, dG0_value + dG0_uncertainty),
        ),
    )

    f_var_name = f"{DF_VAR_PREFIX}{reac_id}"
    setattr(model, f_var_name, Var(within=Reals, bounds=(-QUASI_INF, QUASI_INF)))

    # e.g.,
    #     RT * ln(([S]*[T])/([A]²*[B]))
    # <=> RT*ln([S]) + RT*ln([T]) - 2*RT*ln([A]) - RT*ln([B])
    f_expression_lhs = -getattr(model, dG0_var_name)
    if add_error_term:
        error_var_id = f"{ERROR_VAR_PREFIX}_dG0_{reac_id}"
        setattr(
            model,
            error_var_id,
            Var(within=Reals, bounds=(0.0, max_abs_dG0_correction)),
        )
        f_expression_lhs += getattr(model, error_var_id)
        f_bound_correction_var_id = f"{ERROR_BOUND_UPPER_CHANGE_PREFIX}{f_var_name}"
        setattr(model, f_bound_correction_var_id, Var())
        getattr(model, f_bound_correction_var_id).fix(max_abs_dG0_correction)
    for met_id in reaction.stoichiometries:
        stoichiometry = reaction.stoichiometries[met_id]
        if met_id.startswith(ENZYME_VAR_PREFIX) or (met_id == PROT_POOL_MET_NAME):
            continue
        f_expression_lhs += (
            (-1)
            * cobrak_model.R
            * cobrak_model.T
            * stoichiometry
            * getattr(model, f"{LNCONC_VAR_PREFIX}{met_id}")
        )
    if strict_df_equality:
        setattr(
            model,
            f"{DF_VAR_PREFIX}_constraint_{reac_id}",
            Constraint(rule=f_expression_lhs == getattr(model, f_var_name)),
        )
    else:
        setattr(
            model,
            f"{DF_VAR_PREFIX}_constraint_{reac_id}",
            Constraint(rule=f_expression_lhs >= getattr(model, f_var_name)),
        )
    return model, f_var_name


def _add_error_sum_to_model(
    model: ConcreteModel,
    cobrak_model: Model,
    correction_config: CorrectionConfig,
) -> Model:
    """Adds an error sum to the model, which can be either a quadratic or linear sum of error variables,
       optionally weighted based on certain parameters from the COBRA model.

    Args:
        model (ConcreteModel): The Pyomo model to which the error sum will be added.
        cobrak_model (Model): The COBRA model containing the necessary parameters for weighting.
        correction_config (CorrectionConfig): The correction configuration determining how the sum is built.

    Returns:
        Model: The modified model with the added error sum variable and constraint.
    """
    error_model_var_ids = [
        var_name
        for var_name in get_model_var_names(model)
        if var_name.startswith(ERROR_VAR_PREFIX)
    ]

    if correction_config.use_weights:
        kcat_times_e_weight: float = percentile(
            sorted(get_model_max_kcat_times_e_values(cobrak_model)),
            correction_config.weight_percentile,
        ) / len(cobrak_model.enzymes)
        dG0_weight: float = percentile(
            sorted(get_model_dG0s(cobrak_model, abs_values=True)),
            correction_config.weight_percentile,
        )
        km_weight: float = percentile(
            sorted([abs(log(k_m)) for k_m in get_model_kms(cobrak_model)]),
            correction_config.weight_percentile,
        )

    error_expr: Expression = 0.0
    for error_model_var_id in error_model_var_ids:
        if not correction_config.use_weights:
            weight = 1.0
        else:
            if "kcat_" in error_model_var_id:
                weight = 1 / kcat_times_e_weight
            elif error_model_var_id.endswith(("substrate", "product")):
                weight = 1 / km_weight
            elif "dG0" in error_model_var_id:
                weight = 1 / dG0_weight
            elif (
                error_model_var_id[len(ERROR_VAR_PREFIX) + 1 :].split("_origstart")[0]
                in correction_config.extra_weights
            ):
                weight = abs(
                    correction_config.extra_weights[
                        error_model_var_id[len(ERROR_VAR_PREFIX) + 1 :].split(
                            "_origstart"
                        )[0]
                    ]
                )
            else:
                weight = 1.0

        if correction_config.error_sum_as_qp:
            error_expr += weight * getattr(model, error_model_var_id) ** 2
        else:
            error_expr += weight * getattr(model, error_model_var_id)
    setattr(model, ERROR_SUM_VAR_ID, Var(within=Reals, bounds=(1e-6, QUASI_INF)))
    setattr(
        model,
        f"constraint_{ERROR_SUM_VAR_ID}",
        Constraint(
            rule=getattr(model, ERROR_SUM_VAR_ID) >= error_expr,
        ),
    )
    return model


def _add_extra_linear_constraints_to_lp(
    model: ConcreteModel,
    cobrak_model: Model,
) -> ConcreteModel:
    """Adds extra linear constraints from the COBRAk Model to the Pyomo model.

    This function iterates through each extra linear constraint defined in the COBRAk Model.
    For each constraint, it checks if all required variables exist in the current Pyomo model.
    If any variable is missing, the constraint is skipped. Otherwise, it adds the constraint
    to the model, setting either a lower bound, an upper bound, or both, based on the values
    specified in the COBRAk model.

    Args:
        model (ConcreteModel): The Pyomo model to which the constraints will be added.
        cobrak_model (Model): The COBRAk Model containing the extra linear constraints.

    Returns:
        ConcreteModel: The updated Pyomo model with the added extra linear constraints.
    """
    for extra_constraint_counter, extra_linear_constraint in enumerate(
        cobrak_model.extra_linear_constraints
    ):
        missing_var = False
        extra_constraint_lhs = 0.0
        for var_id in extra_linear_constraint.stoichiometries:
            if var_id not in get_model_var_names(model):
                missing_var = True
                continue
            extra_constraint_lhs += extra_linear_constraint.stoichiometries[
                var_id
            ] * getattr(model, var_id)

        if missing_var:
            continue

        base_extra_constraint_name = f"Extra_constraint_{extra_constraint_counter}_"

        if extra_linear_constraint.lower_value is not None:
            setattr(
                model,
                f"{base_extra_constraint_name}LB",
                Constraint(
                    expr=extra_constraint_lhs >= extra_linear_constraint.lower_value
                ),
            )
        if extra_linear_constraint.upper_value is not None:
            setattr(
                model,
                f"{base_extra_constraint_name}UB",
                Constraint(
                    expr=extra_constraint_lhs <= extra_linear_constraint.upper_value
                ),
            )
    return model


def _add_enzyme_constraints_to_lp(
    model: ConcreteModel,
    cobrak_model: Model,
    ignored_reacs: list[str] = [],
    add_error_term: bool = False,
    error_cutoff: float = 0.1,
    max_rel_correction: float = QUASI_INF,
) -> ConcreteModel:
    """Add MOMENT [1]-like enzyme constraint to (N/MI)LP.

    This function adds enzyme constraints to a given (N/MI)LP Pyomo model based on the
    COBRAk model's enzyme data. It follows the MOMENT approach described in [1]
    to incorporate enzyme usage into metabolic models.

    Args:
        model (ConcreteModel): The Pyomo instance of the (N/MI)LP model.
        cobrak_model (Model): The associated metabolic model containing enzyme data.
        ignored_reacs (list[str], optional): A list of reaction IDs for which enzyme constraints
                                             should not be included. Defaults to an empty list.

    Returns:
        ConcreteModel: The modified Pyomo model with added enzyme constraints.

    [1] Adadi et al. PLoS computational biology 8.7 (2012): e1002575. https://doi.org/10.1371/journal.pcbi.1002575
    """
    # Add the protein pool delivery pseudo-reaction (flux in [g/gDW])
    setattr(
        model,
        PROT_POOL_REAC_NAME,
        Var(within=Reals, bounds=(0.0, cobrak_model.max_prot_pool)),
    )

    # Collect all kcats and get error-eligible reactions
    if add_error_term:
        all_max_kcat_times_e_values = get_model_max_kcat_times_e_values(cobrak_model)
        max_kcat_times_e_lowbound = all_max_kcat_times_e_values[
            : ceil(error_cutoff * len(all_max_kcat_times_e_values))
        ][-1]

    # Expression for the sum which describes thee used protein pool (in [g/gDW])
    prot_pool_sum = 0.0
    # Go through each reaction, checking for the inclusion of enzyme constraints
    for reac_id, reaction in cobrak_model.reactions.items():
        # Ignore reactions where the user specifically said that no enzyme constraints
        # shall be included.
        if reac_id in ignored_reacs:
            continue

        # Check enzyme reaction data...
        enzyme_reaction_data = reaction.enzyme_reaction_data
        # ...with no such data, we can't add enzyme constraints
        # ...and with no associated enzyme, we can't add enzyme constraints
        # ...and with (erronously given, often happens for ATPM) empty enzyme
        # identifiers, we also cannot add enzyme constraints.
        if (
            (enzyme_reaction_data is None)
            or (len(enzyme_reaction_data.identifiers) == 0)
            or ("" in enzyme_reaction_data.identifiers)
        ):
            continue

        # Calculate the sum of the molecular weights of the reaction
        # enzyme's quartery structure (i.e., sum of all MWs of all
        # subenzymes).
        full_enzyme_mw = get_full_enzyme_mw(cobrak_model, reaction)
        min_enzyme_conc = None
        max_enzyme_conc = None
        for enzyme_id in enzyme_reaction_data.identifiers:
            enzyme = cobrak_model.enzymes[enzyme_id]

            # If given, add concentration range constraints for the specific
            # enzyme quartery structure.
            if enzyme.min_conc is not None:
                if min_enzyme_conc is None:
                    min_enzyme_conc = 0.0
                min_enzyme_conc = max(min_enzyme_conc, enzyme.min_conc)
            if enzyme.max_conc is not None:
                if max_enzyme_conc is None:
                    max_enzyme_conc = float("inf")
                max_enzyme_conc = min(max_enzyme_conc, enzyme.max_conc)

        # Collect all further relevant data for the enzyme constraint...
        # ...k_cat [in 1/h]
        k_cat = enzyme_reaction_data.k_cat
        # ...full enzyme name (for constraint naming)
        full_enzyme_id = get_reaction_enzyme_var_id(reac_id, reaction)
        # ...minimal and maximal possible enzyme concentration
        if min_enzyme_conc is None:
            min_enzyme_conc = 0.0
        if max_enzyme_conc is None:
            # More than MW*max_prot_pool cannot be used anyway
            max_enzyme_conc = full_enzyme_mw * cobrak_model.max_prot_pool

        # Add enzyme concentration variable
        setattr(
            model,
            full_enzyme_id,
            Var(within=Reals, bounds=(min_enzyme_conc, max_enzyme_conc)),
        )
        # Add enzyme concentration constraint
        if add_error_term:
            max_k_cat_times_e: float = (
                k_cat
                * cobrak_model.max_prot_pool
                / get_full_enzyme_mw(cobrak_model, reaction)
            )
            if max_k_cat_times_e <= max_kcat_times_e_lowbound:
                kcat_times_e_error_var_id = f"{ERROR_VAR_PREFIX}_kcat_times_e_{reac_id}"
                setattr(
                    model,
                    kcat_times_e_error_var_id,
                    Var(within=Reals, bounds=(0.0, QUASI_INF)),
                )
                enzyme_constraint_expr: Expression = getattr(model, reac_id) <= getattr(
                    model, full_enzyme_id
                ) * k_cat + getattr(model, kcat_times_e_error_var_id)

                setattr(
                    model,
                    f"enzyme_error_bound_constraint_{reac_id}",
                    Constraint(
                        expr=getattr(model, kcat_times_e_error_var_id)
                        <= max_rel_correction * getattr(model, full_enzyme_id) * k_cat
                    ),
                )
            else:
                enzyme_constraint_expr: Expression = (
                    getattr(model, reac_id) <= getattr(model, full_enzyme_id) * k_cat
                )
        else:
            enzyme_constraint_expr: Expression = (
                getattr(model, reac_id) <= getattr(model, full_enzyme_id) * k_cat
            )
        setattr(
            model,
            f"enzyme_constraint_{reac_id}",
            Constraint(expr=enzyme_constraint_expr),
        )

        # Add current enzyme usage to total enzyme pool
        prot_pool_sum += full_enzyme_mw * getattr(model, full_enzyme_id)

    # Finally, set the protein pool
    setattr(
        model,
        f"{PROT_POOL_REAC_NAME}_constraint",
        Constraint(expr=getattr(model, PROT_POOL_REAC_NAME) >= prot_pool_sum),
    )

    return model


def _add_kappa_substrates_and_products_vars(
    model: ConcreteModel,
    reac_id: str,
    reaction: Reaction,
    cobrak_model: Model,
    strict_kappa_products_equality: bool,
    add_error_term: bool = False,
    max_rel_km_correction: float = QUASI_INF,
    kms_lowbound: float = QUASI_INF,
    kms_highbound: float = 0.0,
) -> tuple[ConcreteModel, str, str]:
    """Add kappa variables and constraints to the given model for the specified reaction.

    This function adds kappa substrate and product variables to the model and sets up the corresponding constraints
    based on the stoichiometries of the reaction's metabolites and the kinetic parameters.

    Kappa substrates is the sum of logarithmized substrate concentrations minus the sum of logrithmized associated kms.
    Kappa products is the sum of logarithmized products concentrations minus the sum of logrithmized associated kms.
    Kinetically ignored metabolites (as set in the model) are not evaluated.

    Parameters:
    - model (ConcreteModel): The Pyomo model to which the variables and constraints will be added.
    - reac_id (str): The identifier for the reaction.
    - reaction (Reaction): The COBRAk Reaction object containing stoichiometry and enzyme reaction data.
    - cobrak_model (Model): The COBRAk Model object containing kinetic ignored metabolites.
    - strict_kappa_products_equality (bool): Whether the kappa products constraint should be an equality or >= inequality.
                                             The latter can make kappa only lower, not higher than with an equality.

    Returns:
    - tuple[ConcreteModel, str, str]: The updated model, kappa substrates variable ID, and kappa products variable ID.

    Raises:
        ValueError: If the Reaction has no EnzymeReactionData instance
    """
    if reaction.enzyme_reaction_data is None:
        raise ValueError

    # Add kappa variables
    kappa_substrates_var_id = f"{KAPPA_SUBSTRATES_VAR_PREFIX}{reac_id}"
    kappa_products_var_id = f"{KAPPA_PRODUCTS_VAR_PREFIX}{reac_id}"

    setattr(
        model,
        kappa_substrates_var_id,
        Var(within=Reals, bounds=(-QUASI_INF, QUASI_INF)),
    )
    setattr(
        model,
        kappa_products_var_id,
        Var(within=Reals, bounds=(-QUASI_INF, QUASI_INF)),
    )

    if add_error_term:
        max_reac_product_changes = 0.0
        max_reac_substrate_changes = 0.0
    kappa_substrates_lhs: Expression = -1.0 * getattr(model, kappa_substrates_var_id)
    kappa_substrates_sum = 0.0
    kappa_products_lhs: Expression = -1.0 * getattr(model, kappa_products_var_id)
    kappa_products_sum = 0.0
    for reac_met_id, stoichiometry in reaction.stoichiometries.items():
        if reac_met_id in cobrak_model.kinetic_ignored_metabolites:
            continue
        if reac_met_id.startswith(ENZYME_VAR_PREFIX):
            continue
        k_m = reaction.enzyme_reaction_data.k_ms[reac_met_id]
        if stoichiometry > 0.0:  # Product
            kappa_products_lhs += stoichiometry * getattr(
                model, f"{LNCONC_VAR_PREFIX}{reac_met_id}"
            )
            kappa_products_sum -= stoichiometry * log(k_m)

            if add_error_term and k_m <= kms_lowbound:
                km_product_error_var_id = (
                    f"{ERROR_VAR_PREFIX}_{reac_id}____{reac_met_id}_product"
                )
                max_product_change = abs(
                    log(k_m) - log((1 + max_rel_km_correction) * k_m)
                )
                setattr(
                    model,
                    km_product_error_var_id,
                    Var(
                        within=Reals,
                        bounds=(
                            0.0,
                            max_product_change,
                        ),
                    ),
                )
                kappa_products_sum -= abs(stoichiometry) * getattr(
                    model, km_product_error_var_id
                )
                max_reac_product_changes += max_product_change
        else:  # Educt
            kappa_substrates_lhs += abs(stoichiometry) * getattr(
                model, f"{LNCONC_VAR_PREFIX}{reac_met_id}"
            )
            kappa_substrates_sum -= abs(stoichiometry) * log(k_m)

            if add_error_term and k_m >= kms_highbound:
                km_substrate_error_var_id = (
                    f"{ERROR_VAR_PREFIX}_{reac_id}____{reac_met_id}_substrate"
                )
                max_substrate_change = abs(
                    log(k_m) - log((1 - max_rel_km_correction) * k_m)
                )
                setattr(
                    model,
                    km_substrate_error_var_id,
                    Var(
                        within=Reals,
                        bounds=(
                            0.0,
                            max_substrate_change,
                        ),
                    ),
                )
                kappa_substrates_sum += abs(stoichiometry) * getattr(
                    model, km_substrate_error_var_id
                )
                max_reac_substrate_changes += max_substrate_change

    if add_error_term:
        products_bound_change_var_id = (
            f"{ERROR_BOUND_UPPER_CHANGE_PREFIX}{kappa_products_var_id}"
        )
        setattr(model, products_bound_change_var_id, Var())
        getattr(model, products_bound_change_var_id).fix(max_reac_product_changes)

        substrates_bound_change_var_id = (
            f"{ERROR_BOUND_LOWER_CHANGE_PREFIX}{kappa_products_var_id}"
        )
        setattr(model, substrates_bound_change_var_id, Var())
        getattr(model, substrates_bound_change_var_id).fix(max_reac_substrate_changes)

    setattr(
        model,
        f"{KAPPA_SUBSTRATES_VAR_PREFIX}_constraint_{reac_id}",
        Constraint(rule=kappa_substrates_lhs == -kappa_substrates_sum),
    )
    kappa_products_constraint_id = f"{KAPPA_PRODUCTS_VAR_PREFIX}_constraint_{reac_id}"
    if strict_kappa_products_equality:
        setattr(
            model,
            kappa_products_constraint_id,
            Constraint(rule=kappa_products_lhs == -kappa_products_sum),
        )
    else:
        setattr(
            model,
            kappa_products_constraint_id,
            Constraint(rule=kappa_products_lhs <= -kappa_products_sum),
        )

    return model, kappa_substrates_var_id, kappa_products_var_id


def _add_thermodynamic_constraints_to_lp(
    model: ConcreteModel,
    cobrak_model: Model,
    add_thermobottleneck_analysis_vars: bool,
    min_mdf: float,
    strict_kappa_products_equality: bool = False,
    add_dG0_error_term: bool = False,
    dG0_error_cutoff: float = 1.0,
    max_abs_dG0_correction: float = QUASI_INF,
    add_km_error_term: bool = False,
    km_error_cutoff: float = 1.0,
    max_rel_km_correction: float = QUASI_INF,
    ignored_reacs: list[str] = [],
) -> ConcreteModel:
    """Add thermodynamic constraints to a (N/MI)LP model.

    This function incorporates thermodynamic constraints into a given (N/MI)LP Pyomo model
    based on the COBRAk model's thermodynamic data. It follows the TFBA and OptMDFpathway methods described in
    [1] and [2] to ensure thermodynamic feasibility of metabolic pathways.

    [1] Soh et al. Metabolic flux analysis: methods and protocols (2014): 49-63. https://doi.org/10.1007/978-1-4939-1170-7_3
    [2] Hädicke et al. PLoS computational biology 14.9 (2018): e1006492. https://doi.org/10.1371/journal.pcbi.1006492

    Args:
        model (ConcreteModel): The Pyomo instance of the (N/MI)LP model.
        cobrak_model (Model): The associated metabolic model containing thermodynamic data.
        add_thermobottleneck_analysis_vars (bool): Whether to add variables for thermodynamic bottleneck analysis.
        min_mdf (float): The minimum metabolic driving force (MDF) to be enforced.
        strict_kappa_products_equality (bool, optional): Whether to enforce strict equality for kappa products.
                                                         Defaults to False.

    Returns:
        ConcreteModel: The modified Pyomo model with added thermodynamic constraints.
    """
    cobrak_model = deepcopy(cobrak_model)

    model = _add_concentration_vars_and_constraints(model, cobrak_model)

    # Set OptMDFpathway constraints and kappa vars
    if add_thermobottleneck_analysis_vars:
        zb_sum_expression = 0.0

    if add_dG0_error_term:
        dG0_highbound = _get_dG0_highbound(cobrak_model, dG0_error_cutoff)
    else:
        dG0_highbound = 0.0

    if add_km_error_term:
        kms_lowbound, kms_highbound = _get_km_bounds(cobrak_model, km_error_cutoff)
    else:
        kms_lowbound, kms_highbound = 0.0, 0.0

    setattr(model, MDF_VAR_ID, Var(within=Reals, bounds=(min_mdf, QUASI_INF)))
    for reac_id, reaction in cobrak_model.reactions.items():
        if reac_id in ignored_reacs:
            continue
        has_kappa = True
        if (reaction.enzyme_reaction_data is None) or (
            not have_all_unignored_km(
                reaction, cobrak_model.kinetic_ignored_metabolites
            )
        ):
            has_kappa = False

        if reaction.dG0 is not None:
            model, f_var_name = _add_df_and_dG0_var_for_reaction(
                model,
                reac_id,
                reaction,
                cobrak_model,
                strict_df_equality=True,
                add_error_term=add_dG0_error_term and (reaction.dG0 >= dG0_highbound),
                max_abs_dG0_correction=max_abs_dG0_correction,
            )

            if add_thermobottleneck_analysis_vars:
                zb_varname = f"zb_var_{reac_id}"
                setattr(model, zb_varname, Var(within=Binary))
                zb_sum_expression += getattr(model, zb_varname)

            z_varname = f"{Z_VAR_PREFIX}{reac_id}"
            setattr(model, z_varname, Var(within=Binary))

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
            if not add_thermobottleneck_analysis_vars:
                bigm_optmdfpathway_1_constraint = getattr(model, f_var_name) + (
                    1 - getattr(model, z_varname)
                ) * BIG_M >= getattr(model, MDF_VAR_ID)
            else:
                bigm_optmdfpathway_1_constraint = getattr(model, f_var_name) + (
                    1 - getattr(model, z_varname)
                ) * BIG_M + 2 * BIG_M * getattr(model, zb_varname) >= getattr(
                    model, MDF_VAR_ID
                )
            setattr(
                model,
                f"bigm_optmdfpathway_1_{reac_id}",
                Constraint(rule=bigm_optmdfpathway_1_constraint),
            )

        if not has_kappa:
            continue

        # Add kappa variables
        model, _, _ = _add_kappa_substrates_and_products_vars(
            model,
            reac_id,
            reaction,
            cobrak_model,
            strict_kappa_products_equality,
            add_km_error_term,
            max_rel_km_correction=max_rel_km_correction,
            kms_lowbound=kms_lowbound,
            kms_highbound=kms_highbound,
        )

    if add_thermobottleneck_analysis_vars:
        zb_sum_var_id = "zb_sum"
        setattr(model, zb_sum_var_id, Var(within=Reals))
        setattr(
            model,
            "zb_sum_constraint",
            Constraint(rule=getattr(model, zb_sum_var_id) == zb_sum_expression),
        )

    return model


def _apply_error_scenario(
    model: ConcreteModel,
    cobrak_model: Model,
    correction_config: CorrectionConfig,
) -> None:
    """Applies an error scenario to the model by introducing error terms for reactions and metabolites.

    This function processes an error scenario dictionary to apply errors to either reactions
    (flux error terms) or metabolites (concentration error terms) based on the provided flags.
    For each variable in the scenario, it creates positive and negative error variables and
    adds constraints to the model to account for the error margins. The function then calls
    `_apply_scenario` to apply any remaining scenario items that do not require error terms.

    Args:
        model (ConcreteModel): The Pyomo model to which the error scenario will be applied.
        cobrak_model (Model): The COBRA model containing information about reactions and metabolites.
        correction_config (CorrectionConfig): The corrrection configuration determining which errors are appliued.

    Raises:
        ValueError: If the CorrectionConfig's ```var_lb_ub_application```member variable has an invalid value.
    """
    model_var_ids = get_model_var_names(model)
    error_scenario = deepcopy(correction_config.error_scenario)
    for var_id in list(error_scenario.keys()):
        if (
            (
                correction_config.add_flux_error_term
                and (var_id in cobrak_model.reactions)
            )
            or (
                correction_config.add_met_logconc_error_term
                and (var_id[len(LNCONC_VAR_PREFIX) :] in cobrak_model.metabolites)
            )
            or (
                correction_config.add_enzyme_conc_error_term
                and (
                    var_id[len(ENZYME_VAR_PREFIX) :].split(ENZYME_VAR_INFIX)[0]
                    in cobrak_model.enzymes
                )
            )
        ):
            pass
        else:
            continue

        if var_id not in model_var_ids:
            continue

        lb, ub = (
            error_scenario[var_id][0],
            error_scenario[var_id][1],
        )
        base_var_id = f"{ERROR_VAR_PREFIX}_{var_id}_origstart_{str(error_scenario[var_id][0]).replace('.', '-')}__{str(error_scenario[var_id][1]).replace('.', '-')}_origend"
        plus_var_id = f"{base_var_id}_plus"
        minus_var_id = f"{base_var_id}_minus"
        setattr(
            model,
            plus_var_id,
            Var(within=Reals, bounds=(0.0, QUASI_INF)),
        )
        setattr(
            model,
            minus_var_id,
            Var(within=Reals, bounds=(0.0, QUASI_INF)),
        )
        error_var_plus: Var = getattr(model, plus_var_id)
        error_var_minus: Var = getattr(model, minus_var_id)
        model_var: Var = getattr(model, var_id)
        match correction_config.var_lb_ub_application:
            case "":
                lb_expr = model_var >= lb + error_var_plus - error_var_minus
                ub_expr = model_var <= ub + error_var_plus - error_var_minus
            case "exp":
                lb_expr = exp(model_var) >= exp(lb) + error_var_plus - error_var_minus
                ub_expr = exp(model_var) <= exp(ub) + error_var_plus - error_var_minus
            case "log":
                lb_expr = (
                    log(model_var + 1e-8)
                    >= log(lb + 1e-8) + error_var_plus - error_var_minus
                )
                ub_expr = (
                    log(model_var + 1e-8)
                    <= log(ub + 1e-8) + error_var_plus - error_var_minus
                )
            case _:
                raise ValueError

        setattr(
            model,
            f"{ERROR_CONSTRAINT_PREFIX}{var_id}_lb",
            Constraint(expr=lb_expr),
        )
        setattr(
            model,
            f"{ERROR_CONSTRAINT_PREFIX}{var_id}_ub",
            Constraint(expr=ub_expr),
        )

        del error_scenario[var_id]

    _apply_scenario(model, error_scenario)


def _apply_scenario(
    model: ConcreteModel, scenario: dict[str, tuple[float, float]]
) -> None:
    """Applies a scenario to the Pyomo model by setting variable bounds.

    This function takes a dictionary of variable IDs and their corresponding (lower bound,
    upper bound) tuples, and applies these bounds to the variables in the model.

    Args:
        model (ConcreteModel): The Pyomo model to which the scenario will be applied.
        scenario (dict[str, tuple[float, float]]): A dictionary where keys are variable IDs
            and values are tuples specifying (lower bound, upper bound) for the variables.
    """
    model_var_names = get_model_var_names(model)
    for var_id, (lb, ub) in scenario.items():
        if var_id not in model_var_names:
            continue
        var: Var = getattr(model, var_id)
        var.setlb(lb)
        var.setub(ub)


def _batch_variability_optimization(
    pyomo_solver: SolverFactoryClass,
    model: ConcreteModel,
    batch: list[tuple[str, str]],
    solve_extra_options: dict[str, Any] = {},
) -> list[tuple[bool, str, float | None]]:
    """Perform batch flux variability optimization on a given model. Used in (parallelized) Flux Variability analysis function.

    This function iterates over a batch of objective-target pairs, activates the corresponding objective in the model for each batch member,
    performs the optimization for each member, and collects the results.

    When possible, warmstart is used.

    Parameters:
    - pyomo_solver: The pyomo solver instance used to solve the model.
    - model (Model): The pyomo model on which the optimization will be performed.
    - batch (list[tuple[str, str]]): A list of tuples where each tuple contains an objective name and a target variable ID.

    Returns:
    - list[tuple[bool, str, float | None]]: A list of tuples containing:
        - A boolean indicating if the objective name starts with "MIN_OBJ_".
        - The target variable ID.
        - The optimization result (or None if the solver reached the maximum time limit).

    Example:
    ```
    results = _batch_variability_optimization(solver, model, [('OBJ_1', 'var1'), ('OBJ_2', 'var2')])
    ```
    """
    resultslist: list[tuple[bool, str, float | None]] = []
    for objective_name, target_id in batch:
        getattr(model, objective_name).activate()
        try:
            results = pyomo_solver.solve(
                model, tee=False, warmstart=True, **solve_extra_options
            )
        except Exception:
            results = pyomo_solver.solve(model, tee=False, **solve_extra_options)
        if results.solver.termination_condition == TerminationCondition.maxTimeLimit:
            result = None
        else:
            result = value(getattr(model, target_id))
        getattr(model, objective_name).deactivate()
        resultslist.append((objective_name.startswith("MIN_OBJ_"), target_id, result))
    return resultslist


def _get_dG0_highbound(cobrak_model: Model, dG0_error_cutoff: float) -> float:
    """Calculate the high bound for dG0 values based on a specified error cutoff.

    This function retrieves all dG0 values from the COBRA model, sorts them, and determines
    the high bound by selecting the value at the percentile defined by $1 - dG0_error_cutoff$.

    Args:
        cobrak_model (Model): The COBRA model containing the dG0 values.
        dG0_error_cutoff (float): The fraction of dG0 values to consider for the high bound.

    Returns:
        float: The high bound value for dG0 based on the specified cutoff.
    """
    all_dG0s: list[float] = get_model_dG0s(cobrak_model)
    all_dG0s.sort()
    return all_dG0s[floor((1 - dG0_error_cutoff) * len(all_dG0s)) :][0]


def _get_km_bounds(cobrak_model: Model, km_error_cutoff: float) -> tuple[float, float]:
    """Determine the low and high bounds for km values based on a specified error cutoff.

    This function collects all km values from the reactions in the COBRA model, sorts them,
    and calculates both the low and high bounds by selecting the values at the percentiles
    defined by km_error_cutoff and $1 - km_error_cutoff$, respectively.

    Args:
        cobrak_model (Model): The COBRA model containing the enzyme reaction data with km values.
        km_error_cutoff (float): The fraction of km values to consider for the bounds.

    Returns:
        tuple[float, float]: A tuple containing the low bound and high bound for km values.
    """
    all_kms: list[float] = []
    for reaction in cobrak_model.reactions.values():
        if reaction.enzyme_reaction_data is None:
            continue
        for km_value in reaction.enzyme_reaction_data.k_ms.values():
            all_kms.append(km_value)
    all_kms.sort()
    kms_lowbound = all_kms[: ceil(km_error_cutoff * len(all_kms))][-1]
    kms_highbound = all_kms[floor((1 - km_error_cutoff) * len(all_kms)) :][0]
    return kms_lowbound, kms_highbound


def _get_steady_state_lp_from_cobrak_model(
    cobrak_model: Model,
    ignored_reacs: list[str] = [],
) -> ConcreteModel:
    """Returns the basic linear steady-state constraints as a pyomo model.

    These constraints were first introduced in [1] and constitute:
    * N*r = 0
    * Reaction min and max fluxes

    The extra linear constraints (A*d <= b) are not yet included.

    [1] Watson, M. R. (1984). Metabolic maps for the Apple II. https://doi.org/10.1042/bst0121093

    Returns:
        _description_
    """
    model = ConcreteModel()

    for reac_id, reaction in cobrak_model.reactions.items():
        if reac_id in ignored_reacs:
            continue
        setattr(
            model,
            reac_id,
            Var(within=Reals, bounds=(reaction.min_flux, reaction.max_flux)),
        )

    for met_id in cobrak_model.metabolites:
        constraint_name = f"Steady-state_of_{met_id}"

        constraint_lhs = 0.0
        for reac_id, reaction in cobrak_model.reactions.items():
            if reac_id in ignored_reacs:
                continue
            if met_id in reaction.stoichiometries:
                constraint_lhs += reaction.stoichiometries[met_id] * getattr(
                    model, reac_id
                )

        setattr(model, constraint_name, Constraint(expr=constraint_lhs == 0))

    return model


# "PUBLIC" FUNCTIONS SECTION #
def add_flux_sum_var(model: ConcreteModel, cobrak_model: Model) -> ConcreteModel:
    """Add a flux sum variable to a (N/MI)LP model.

    This function introduces a flux sum variable to a given (N/MI)LP Pyomo model, which represents
    the total sum of absolute fluxes across all reactions in the COBRAk model. The methodology is based on
    the pFBA (Parsimonious Flux Balance Analysis) approach [1].

    [1] Lewis et al. Molecular systems biology 6.1 (2010): 390. https://doi.org/10.1038/msb.2010.47

    Args:
        model (ConcreteModel): The Pyomo instance of the (N/MI)LP model.
        cobrak_model (Model): The associated metabolic model containing reaction data.

    Returns:
        ConcreteModel: The modified Pyomo model with the added flux sum variable and constraint.
    """
    flux_sum_expr = 0.0
    for reac_id in cobrak_model.reactions:
        try:
            flux_sum_expr += getattr(model, reac_id)
        except AttributeError:
            continue

    setattr(model, FLUX_SUM_VAR_ID, Var(within=Reals, bounds=(0.0, 1e9)))
    setattr(
        model,
        "FLUX_SUM_CONSTRAINT",
        Constraint(rule=getattr(model, FLUX_SUM_VAR_ID) == flux_sum_expr),
    )

    return model


def add_loop_constraints_to_lp(
    model: ConcreteModel,
    cobrak_model: Model,
    only_nonthermodynamic: bool,
    ignored_reacs: list[str] = [],
) -> ConcreteModel:
    """Add mixed-integer loop constraints to a (N/MI)LP model to prevent thermodynamically infeasible cycles.

    This function incorporates loop constraints into a given (N/MI)LP Pyomo model based on the COBRAk model's
    reaction data. It follows the ll-COBRA methodology described in [1] to prevent the formation
    of thermodynamically infeasible cycles in metabolic networks.

    [1] Schellenberger et al. (2011). Biophysical journal, 100(3), 544-553. https://doi.org/10.1016/j.bpj.2010.12.3707

    Args:
        model (ConcreteModel): The Pyomo instance of the (N/MI)LP model.
        cobrak_model (Model): The associated metabolic model containing reaction data.
        only_nonthermodynamic (bool): If True, only add constraints to reactions without thermodynamic data.

    Returns:
        ConcreteModel: The modified Pyomo model with added loop constraints.
    """
    base_id_constraints: dict[str, Expression] = {}
    num_elements_per_constraint = {}
    for reac_id, reaction in cobrak_model.reactions.items():
        if reac_id in ignored_reacs:
            continue
        if (only_nonthermodynamic) and (reaction.dG0 is not None):
            continue

        base_id = get_base_id(reac_id, cobrak_model.fwd_suffix, cobrak_model.rev_suffix)
        if base_id not in base_id_constraints:
            base_id_constraints[base_id] = 0.0
            num_elements_per_constraint[base_id] = 0

        zv_var_id = "zV_var_" + reac_id
        setattr(model, zv_var_id, Var(within=Binary))
        setattr(
            model,
            reac_id + "_base",
            Constraint(
                rule=getattr(model, reac_id) <= BIG_M * getattr(model, zv_var_id)
            ),
        )

        base_id_constraints[base_id] += getattr(model, zv_var_id)
        num_elements_per_constraint[base_id] += 1

    for base_id, constraint_lhs in base_id_constraints.items():
        if num_elements_per_constraint[base_id] > 1:
            setattr(
                model,
                base_id + "_base_constraint",
                Constraint(rule=constraint_lhs <= 1.0),
            )

    return model


def get_lp_from_cobrak_model(
    cobrak_model: Model,
    with_enzyme_constraints: bool,
    with_thermodynamic_constraints: bool,
    with_loop_constraints: bool,
    with_flux_sum_var: bool = False,
    ignored_reacs: list[str] = [],
    min_mdf: float = STANDARD_MIN_MDF,
    add_thermobottleneck_analysis_vars: bool = False,
    strict_kappa_products_equality: bool = False,
    add_extra_linear_constraints: bool = True,
    correction_config: CorrectionConfig = CorrectionConfig(),
) -> ConcreteModel:
    """Construct a linear programming (LP) model from a COBRAk model with various constraints and configurations.

    This function creates a steady-state LP model from the provided COBRAk Model and enhances it with
    different types of constraints and variables based on the specified parameters. It allows for the
    inclusion of enzyme constraints, thermodynamic constraints, loop constraints, and additional
    linear constraints. Furthermore, it supports the addition of flux sum variables and error handling
    configurations.

    See the following chapters of COBRAk's documentation for more on these constraints:

    * Steady-state and extra linear constraints ⇒ Chapter "Linear Programs"
    * Enzyme constraints ⇒ Chapter "Linear Programs"
    * Thermodynamic constraints ⇒ Chapter "Mixed-Integer Linear Programs"

    Parameters
    ----------
    cobrak_model : Model
        The COBRAk Model from which to construct the LP model.
    with_enzyme_constraints : bool
        If True, adds enzyme-pool constraints to the model.
    with_thermodynamic_constraints : bool
        If True, adds thermodynamic MILP constraints to the model, ensuring that reaction fluxes are
        thermodynamically feasible by considering Gibbs free energy changes.
    with_loop_constraints : bool
        If True, adds loop constraints to prevent or control flux loops in the metabolic network.
        This constraint makes the LP a MILP as a binary variable controls whether either the
        forward or the reverse reaction is running.
    with_flux_sum_var : bool, optional
        If True, adds a flux sum variable to the model, which aggregates the total flux through
        all reactions for optimization or analysis purposes. Defaults to False.
    ignored_reacs : list[str], optional
        List of reaction IDs to ignore in the model, which will be excluded. Defaults to [].
    min_mdf : float, optional
        Minimum value for Max-Min Driving Force (MDF). Only relevant with thermodynamic
        constraints. Defaults to STANDARD_MIN_MDF.
    add_thermobottleneck_analysis_vars : bool, optional
        If True, adds variables for thermodynamic bottleneck analysis, helping to identify
        potential bottlenecks in the metabolic network where thermodynamic constraints might limit
        flux. Defaults to False.
    strict_kappa_products_equality : bool, optional
        If True, enforces strict equality for kappa products, ensuring consistency in
        thermodynamic parameters related to reaction products. Defaults to False.
    add_extra_linear_constraints : bool, optional
        If True, adds extra linear constraints from the COBRAk Mmodel, allowing for additional
        linear constraints. Defaults to True.
    correction_config : CorrectionConfig, optional
        Configuration for parameter correction handling in the model, allowing for the inclusion of error terms
        in constraints related to enzyme activity, thermodynamics, etc. Defaults to CorrectionConfig().

    Returns
    -------
    ConcreteModel
        The constructed LP model with the specified constraints and configurations.
    """
    # Initialize the steady-state LP model from the COBRA model, ignoring specified reactions
    model: ConcreteModel = _get_steady_state_lp_from_cobrak_model(
        cobrak_model=cobrak_model,
        ignored_reacs=ignored_reacs,
    )

    # Add enzyme constraints if enabled
    if with_enzyme_constraints:
        model = _add_enzyme_constraints_to_lp(
            model=model,
            cobrak_model=cobrak_model,
            ignored_reacs=ignored_reacs,
            add_error_term=correction_config.add_kcat_times_e_error_term,
            error_cutoff=correction_config.kcat_times_e_error_cutoff,
            max_rel_correction=correction_config.max_rel_kcat_times_e_correction,
        )

    # Add thermodynamic constraints if enabled
    if with_thermodynamic_constraints:
        model = _add_thermodynamic_constraints_to_lp(
            model=model,
            cobrak_model=cobrak_model,
            add_thermobottleneck_analysis_vars=add_thermobottleneck_analysis_vars,
            min_mdf=min_mdf,
            strict_kappa_products_equality=strict_kappa_products_equality,
            add_dG0_error_term=correction_config.add_dG0_error_term,
            dG0_error_cutoff=correction_config.dG0_error_cutoff,
            max_abs_dG0_correction=correction_config.max_abs_dG0_correction,
            add_km_error_term=correction_config.add_km_error_term,
            km_error_cutoff=correction_config.km_error_cutoff,
            max_rel_km_correction=correction_config.max_rel_km_correction,
            ignored_reacs=ignored_reacs,
        )

        if cobrak_model.max_conc_sum < float("inf"):
            model = _add_conc_sum_constraints(cobrak_model, model)

    # Add loop constraints if enabled
    if with_loop_constraints:
        model = add_loop_constraints_to_lp(
            model,
            cobrak_model,
            only_nonthermodynamic=with_thermodynamic_constraints,
            ignored_reacs=ignored_reacs,
        )

    # Add flux sum variable if enabled
    if with_flux_sum_var:
        model = add_flux_sum_var(
            model,
            cobrak_model,
        )

    # Apply error scenarios and add error sum term if error handling is configured
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

    # Add extra linear constraints if enabled
    if add_extra_linear_constraints:
        model = _add_extra_linear_constraints_to_lp(
            model,
            cobrak_model,
        )

    return model


def perform_lp_min_active_reactions_analysis(
    cobrak_model: Model,
    with_enzyme_constraints: bool,
    variability_dict: dict[str, tuple[float, float]],
    min_mdf: float = 0.0,
    verbose: bool = False,
    solver: Solver = SCIP,
) -> float:
    """Run a mixed-integer linear program to determine the minimum number of active reactions.

    This function constructs and solves a mixed-integer linear programming model to find the minimum number of
    reactions that need to be active to satisfy the given variability constraints. It uses a binary
    variable for each reaction to indicate whether it is active, and the objective is to minimize
    the sum of these binary variables. The model includes constraints based on enzyme data,
    thermodynamic feasibility, and loop prevention, depending on the specified parameters.

    Parameters
    ----------
    cobrak_model : Model
        The COBRA model containing the metabolic network and reaction data.
    with_enzyme_constraints : bool
        If True, includes enzyme-pool constraints in the model.
    variability_dict : dict[str, tuple[float, float]]
        A dictionary where keys are reaction IDs and values are tuples specifying (lower bound,
        upper bound) for reaction fluxes.
    min_mdf : float, optional
        Minimum value for Min-Max Driving Force (MDF), setting a lower bound on fluxes.
        Defaults to 0.0.
    verbose : bool, optional
        If True, enables solver output. Defaults to False.
    solver: Solver
        The MILP solver used for this analysis. Defaults to SCIP.

    Returns
    -------
    float
        The minimum number of active reactions required to satisfy the constraints.
    """
    # Create a deep copy of the COBRAk model to avoid modifying the original model
    cobrak_model = deepcopy(cobrak_model)

    # Remove reactions that are not present in the variability dictionary
    minz_cobrak_model = delete_unused_reactions_in_variability_dict(
        cobrak_model, variability_dict
    )

    # Construct the LP model with the specified constraints
    minz_model, _ = get_lp_from_cobrak_model(
        minz_cobrak_model,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=True,
        with_loop_constraints=False,
        min_mdf=min_mdf,
    )

    # Initialize the sum of binary variables to zero
    extrazsum_expression = 0.0

    # Iterate over all potentially active reactions
    for reac_id in get_potentially_active_reactions_in_variability_dict(
        cobrak_model, variability_dict
    ):
        # Create a binary variable for each reaction to indicate activity
        extraz_varname = f"extraz_var_{reac_id}"
        setattr(minz_model, extraz_varname, Var(within=Binary))

        # Add a constraint to relate reaction flux to the binary variable
        setattr(
            minz_model,
            f"extraz_const_{reac_id}",
            Constraint(
                rule=getattr(minz_model, reac_id)
                <= BIG_M * getattr(minz_model, extraz_varname)
            ),
        )

        # Accumulate the binary variables in the sum expression
        extrazsum_expression += getattr(minz_model, extraz_varname)

    # Add a variable to represent the total sum of active reactions
    setattr(minz_model, "extrazsum", Var(within=Reals))

    # Add a constraint to equate the sum variable to the sum expression
    setattr(
        minz_model,
        "extrazsum_const",
        Constraint(rule=getattr(minz_model, "extrazsum") == extrazsum_expression),
    )

    # Set the objective function to minimize the number of active reactions
    minz_model.obj = get_objective(minz_model, "extrazsum", minimize)

    # Initialize the solver with the specified options and attributes
    solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)

    # Solve the LP model
    solver.solve(minz_model, tee=verbose, **solver.solve_extra_options)

    # Retrieve the solution as a dictionary
    minz_dict = get_pyomo_solution_as_dict(minz_model)

    # Return the minimum number of active reactions
    return minz_dict["extrazsum"]


def perform_lp_optimization(
    cobrak_model: Model,
    objective_target: str | dict[str, float],
    objective_sense: int,
    with_enzyme_constraints: bool = False,
    with_thermodynamic_constraints: bool = False,
    with_loop_constraints: bool = False,
    variability_dict: dict[str, tuple[float, float]] = {},
    ignored_reacs: list[str] = [],
    min_mdf: float = STANDARD_MIN_MDF,
    verbose: bool = False,
    with_flux_sum_var: bool = False,
    solver: Solver = SCIP,
    correction_config: CorrectionConfig = CorrectionConfig(),
) -> dict[str, float]:
    """Perform linear programming optimization on a COBRAk model to determine flux distributions.

    This function constructs and solves an LP problem for the given metabolic model using specified constraints,
    variables, and solver options. It supports various types of constraints such as enzyme constraints, thermodynamic
    constraints, and loop constraints. Additionally, it can handle variability dictionaries and ignored reactions.

    Parameters:
        cobrak_model (Model): A COBRAk Model object representing the metabolic network.
        objective_target (str | dict[str, float]): The target for optimization. Can be a reaction ID if optimizing a single
            reaction or a dictionary specifying flux values for multiple reactions.
        objective_sense (int): The sense of the optimization problem (+1: maximize, -1: minimize).
        with_enzyme_constraints (bool, optional): Whether to include enzyme constraints in the model. Defaults to False.
        with_thermodynamic_constraints (bool, optional): Whether to include thermodynamic constraints in the model.
            Defaults to False.
        with_loop_constraints (bool, optional): Whether to include loop closure constraints in the model. Defaults to False.
        variability_dict (dict[str, tuple[float, float]], optional): A dictionary specifying variable bounds for reactions
            or metabolites. Defaults to an empty dict.
        ignored_reacs (list[str], optional): List of reaction IDs to deactivate during optimization. Defaults to an empty list.
        min_mdf (float, optional): Minimum metabolic distance factor threshold for thermodynamic constraints. Defaults to STANDARD_MIN_MDF.
        verbose (bool, optional): Whether to print solver output information. Defaults to False.
        with_flux_sum_var (bool, optional): Whether to include flux sum variable in the model. Defaults to False.
        solver (Solver, optional): Solver used for LP. Default is SCIP.
        correction_config (CorrectionConfig, optional): Configuration for handling prameter corrections and scenarios during optimization.

    Returns:
        dict[str, float]: A dictionary containing the flux distribution results for each reaction in the model.
    """
    optimization_cobrak_model = deepcopy(cobrak_model)
    if variability_dict != {}:
        optimization_cobrak_model = delete_unused_reactions_in_variability_dict(
            cobrak_model,
            variability_dict,
        )
    optimization_model = get_lp_from_cobrak_model(
        cobrak_model=optimization_cobrak_model,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=with_thermodynamic_constraints,
        with_loop_constraints=with_loop_constraints,
        with_flux_sum_var=with_flux_sum_var,
        min_mdf=min_mdf,
        correction_config=correction_config,
    )

    for deactivated_reaction in set(ignored_reacs):
        try:
            setattr(
                optimization_model,
                f"DEACTIVATE_{deactivated_reaction}",
                Constraint(
                    expr=getattr(optimization_model, deactivated_reaction) == 0.0
                ),
            )
        except AttributeError:
            continue

    optimization_model = apply_variability_dict(
        optimization_model,
        cobrak_model,
        variability_dict,
        correction_config.error_scenario,
    )
    optimization_model.obj = get_objective(
        optimization_model, objective_target, objective_sense
    )

    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)
    results = pyomo_solver.solve(
        optimization_model, tee=verbose, **solver.solve_extra_options
    )

    fba_dict = get_pyomo_solution_as_dict(optimization_model)

    return add_statuses_to_optimziation_dict(fba_dict, results)


def perform_lp_thermodynamic_bottleneck_analysis(
    cobrak_model: Model,
    with_enzyme_constraints: bool = False,
    min_mdf: float = STANDARD_MIN_MDF,
    verbose: bool = False,
    solver: Solver = SCIP,
) -> tuple[list[str], dict[str, float]]:
    """Perform thermodynamic bottleneck analysis on a COBRAk model using mixed-integer linear programming.

    This function identifies a minimal set of thermodynamic bottlenecks in a COBRAk model by minimizing the sum of
    newly introduced binary variables that indicate bottleneck reactions, i.e. reactions that do not allow the
    max-min driving force (MDF) to become at least the set min_mdf.
    This methology was first described in [1]. Keep in mind that results from this function are optimal, but not
    neccessarily unique!

    [1] Bekiaris et al. (2023). Nature Communications, 14(1), 4660.  https://doi.org/10.1038/s41467-023-40297-8

    Args:
        cobrak_model (Model): The COBRAk model to analyze for thermodynamic bottlenecks.
        with_enzyme_constraints (bool): Whether to include enzyme constraints in the analysis.
        min_mdf (float, optional): Minimum max-min driving force (MDF) to be enforced. Defaults to STANDARD_MIN_MDF.
        verbose (bool, optional): If True, print detailed information about identified bottlenecks. Defaults to False.
        solver_name (str, optional): Name of the solver to use for optimization. Defaults to "scip".
        solver_options (dict[str, float | int | str], optional): Options for the solver, such as number of threads
                                                                 and LP method. Defaults to an empty dictionary.

    Returns:
        list[str], dict[str, float]: A list of reaction IDs identified as thermodynamic bottlenecks, and the associated solution dict.
    """
    cobrak_model = deepcopy(cobrak_model)
    thermo_constraint_lp = get_lp_from_cobrak_model(
        cobrak_model,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=True,
        with_loop_constraints=False,
        add_thermobottleneck_analysis_vars=True,
        min_mdf=min_mdf,
    )

    thermo_constraint_lp.obj = get_objective(
        thermo_constraint_lp,
        "zb_sum",
        objective_sense=-1,
    )
    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)
    try:
        pyomo_solver.solve(
            thermo_constraint_lp, tee=verbose, **solver.solve_extra_options
        )
        solution_dict = get_pyomo_solution_as_dict(thermo_constraint_lp)
    except (ApplicationError, AttributeError, ValueError):
        solution_dict = {ALL_OK_KEY: False}

    bottleneck_counter = 1
    bottleneck_reactions = []
    for var_id, var_value in solution_dict.items():
        if not var_id.startswith(ZB_VAR_PREFIX):
            continue
        if var_value <= 0.01:
            continue
        bottleneck_reac_id = var_id.replace(ZB_VAR_PREFIX, "")
        bottleneck_reactions.append(bottleneck_reac_id)
        if verbose:
            bottleneck_dG0 = cobrak_model.reactions[bottleneck_reac_id].dG0
            if bottleneck_dG0 is not None:
                printed_dG0 = round(bottleneck_dG0, 3)
            printed_string = get_reaction_string(cobrak_model, bottleneck_reac_id)
            print(
                f"#{bottleneck_counter}: {bottleneck_reac_id} with ΔG'° of {printed_dG0} kJ/mol, {printed_string}"
            )
        bottleneck_counter += 1

    return bottleneck_reactions, solution_dict


def perform_lp_variability_analysis(
    cobrak_model: Model,
    with_enzyme_constraints: bool = False,
    with_thermodynamic_constraints: bool = False,
    active_reactions: list[str] = [],
    min_active_flux: float = 1e-3,
    calculate_reacs: bool = True,
    calculate_concs: bool = True,
    calculate_rest: bool = True,
    further_tested_vars: list[str] = [],
    min_mdf: float = STANDARD_MIN_MDF,
    min_flux_cutoff: float = 1e-5,
    abs_df_cutoff: float = 1e-5,
    min_enzyme_cutoff: float = 1e-5,
    max_active_enzyme_cutoff: float = 1e-4,
    solver: Solver = SCIP,
    parallel_verbosity_level: int = 0,
) -> dict[str, tuple[float, float]]:
    """Perform linear programming variability analysis on a COBRAk model.

    This function conducts a variability analysis on a COBRAk model using linear programming (LP).
    It evaluates the range of possible flux values for each reaction and all other occuring variables in the model,
    considering optional enzyme and thermodynamic constraints. The methodology is based on the approach
    described by [1] and parallelized as outlined in [2].

    [1] Mahadevan & Schilling. (2003). Metabolic engineering, 5(4), 264-276. https://doi.org/10.1016/j.ymben.2003.09.002
    [2] Gudmundsson & Thiele. BMC Bioinformatics 11, 489 (2010). https://doi.org/10.1186/1471-2105-11-489

    Args:
        cobrak_model (Model): The COBRAk model to analyze.
        with_enzyme_constraints (bool): Whether to include enzyme constraints in the analysis.
        with_thermodynamic_constraints (bool): Whether to include thermodynamic constraints in the analysis.
        active_reactions (list[str], optional): List of reactions to be set as active with a minimum flux.
                                                Defaults to an empty list.
        min_active_flux (float, optional): Minimum flux value for active reactions. Defaults to 1e-5.
        calculate_reacs (bool, optional): If True, analyze reaction fluxes. Default: True.
        calculate_concs (bool, optional): If True, analyze concentrations. Default: True.
        calculate_rest (bool, optional): If True, analyze all other parameters (e.g. kappa values and driving forces). Default: True.
        min_mdf (float, optional): Minimum metabolic driving force (MDF) to be enforced. Defaults to 0.0.
        min_flux_cutoff (float, optional): Minimum flux cutoff for considering a reaction active. Defaults to 1e-8.
        solver (Solver, optional): MILP solver used for variability analysis. Default is SCIP, recommended is CPLEX_FOR_VARIABILITY_ANALYSIS
                                   or GUROBI_FOR_VARIABILITY_ANALYSIS if you have a CPLEX or Gurobi license.
        parallel_verbosity_level (int, optional): Sets the verbosity level for the analysis parallelization. The higher,
                                                  the value, the more is printed. Default: 0.

    Returns:
        dict[str, tuple[float, float]]: A dictionary mapping variable IDs to their minimum and maximum values
                                        determined by the variability analysis.
    """
    cobrak_model = deepcopy(cobrak_model)
    for active_reaction in active_reactions:
        cobrak_model.reactions[active_reaction].min_flux = min_active_flux

    model = get_lp_from_cobrak_model(
        cobrak_model=cobrak_model,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=with_thermodynamic_constraints,
        with_loop_constraints=True,
        min_mdf=min_mdf,
        strict_kappa_products_equality=True,
    )
    model_var_names = get_model_var_names(model)

    min_values: dict[str, float] = {}
    max_values: dict[str, float] = {}
    objective_targets: list[tuple[int, str]] = []

    max_flux_sum_result = perform_lp_optimization(
        cobrak_model,
        objective_target=FLUX_SUM_VAR_ID,
        objective_sense=+1,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=True,
        with_loop_constraints=True,
        with_flux_sum_var=True,
        solver=solver,
    )
    min_flux_sum_result = perform_lp_optimization(
        cobrak_model,
        objective_target=FLUX_SUM_VAR_ID,
        objective_sense=-1,
        with_enzyme_constraints=with_enzyme_constraints,
        with_thermodynamic_constraints=True,
        with_loop_constraints=True,
        with_flux_sum_var=True,
        solver=solver,
    )

    if (calculate_concs or calculate_rest) and with_thermodynamic_constraints:
        min_mdf_result = perform_lp_optimization(
            cobrak_model,
            objective_target=MDF_VAR_ID,
            objective_sense=-1,
            with_enzyme_constraints=with_enzyme_constraints,
            with_thermodynamic_constraints=True,
            with_loop_constraints=True,
            solver=solver,
        )
        max_mdf_result = perform_lp_optimization(
            cobrak_model,
            objective_target=MDF_VAR_ID,
            objective_sense=+1,
            with_enzyme_constraints=with_enzyme_constraints,
            with_thermodynamic_constraints=True,
            with_loop_constraints=True,
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
        min_flux_sum_flux = min_flux_sum_result[reac_id]
        max_flux_sum_flux = max_flux_sum_result[reac_id]

        if calculate_reacs:
            if reaction.min_flux in (min_flux_sum_flux, max_flux_sum_flux):
                min_values[reac_id] = (
                    reaction.min_flux if reaction.min_flux >= min_flux_cutoff else 0.0
                )
            else:
                objective_targets.append((-1, reac_id))
            if reaction.max_flux in (min_flux_sum_flux, max_flux_sum_flux):
                max_values[reac_id] = reaction.max_flux
            else:
                objective_targets.append((+1, reac_id))

        if not calculate_rest:
            continue

        f_var_name = f"{DF_VAR_PREFIX}{reac_id}"
        kappa_substrates_var_name = f"{KAPPA_SUBSTRATES_VAR_PREFIX}{reac_id}"
        kappa_products_var_name = f"{KAPPA_PRODUCTS_VAR_PREFIX}{reac_id}"
        if f_var_name in model_var_names:
            if min_mdf in (min_mdf_result[f_var_name], max_mdf_result[f_var_name]):
                min_values[f_var_name] = min_mdf
            else:
                objective_targets.append((-1, f_var_name))
            objective_targets.append((+1, f_var_name))
        if kappa_substrates_var_name in model_var_names:
            objective_targets.extend(
                ((-1, kappa_substrates_var_name), (+1, kappa_substrates_var_name))
            )
        if kappa_products_var_name in model_var_names:
            objective_targets.extend(
                ((-1, kappa_products_var_name), (+1, kappa_products_var_name))
            )
        if reaction.enzyme_reaction_data is not None and with_enzyme_constraints:
            full_enzyme_id = get_full_enzyme_id(
                reaction.enzyme_reaction_data.identifiers
            )
            if full_enzyme_id:
                enzyme_delivery_var_name = get_reaction_enzyme_var_id(reac_id, reaction)
                if 0.0 in (min_flux_sum_flux, max_flux_sum_flux):
                    min_values[enzyme_delivery_var_name] = 0.0
                else:
                    objective_targets.append((-1, enzyme_delivery_var_name))
                objective_targets.append((+1, enzyme_delivery_var_name))

    for further_tested_var in further_tested_vars:
        objective_targets.append((+1, further_tested_var))
        objective_targets.append((-1, further_tested_var))

    objectives_data: list[tuple[str, str]] = []
    for obj_sense, target_id in objective_targets:
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

    objectives_data_batches = split_list(objectives_data, cpu_count())
    pyomo_solver = get_solver(solver.name, solver.solver_options, solver.solver_attrs)

    results_list = Parallel(n_jobs=-1, verbose=parallel_verbosity_level)(
        delayed(_batch_variability_optimization)(
            pyomo_solver, model, batch, solver.solve_extra_options
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
        if key in cobrak_model.reactions:
            min_values[key] = min_value if min_value >= min_flux_cutoff else 0.0
        if key.startswith(ENZYME_VAR_PREFIX):
            min_values[key] = min_value if min_value >= min_enzyme_cutoff else 0.0
        if key.startswith(DF_VAR_PREFIX):
            min_values[key] = min_value if abs(min_value) >= abs_df_cutoff else 0.0

    enzyme_var_to_reac_id = {
        get_reaction_enzyme_var_id(reac_id, reaction): reac_id
        for reac_id, reaction in cobrak_model.reactions.items()
    }
    for key, max_value in max_values.items():
        if key.startswith(ENZYME_VAR_PREFIX) and (
            (max_values[key] != 0.0) or (max_values[enzyme_var_to_reac_id[key]] > 0.0)
        ):
            max_values[key] = max(max_value, max_active_enzyme_cutoff)
        if key.startswith(DF_VAR_PREFIX):
            max_values[key] = max_value if abs(max_value) >= abs_df_cutoff else 0.0

    all_target_ids = sorted(
        set(
            list(min_values.keys())
            + list(max_values.keys())
            + [obj_target[1] for obj_target in objective_targets]
        )
    )
    variability_dict: dict[str, tuple[float, float]] = {
        target_id: (min_values[target_id], max_values[target_id])
        for target_id in all_target_ids
    }

    return variability_dict
