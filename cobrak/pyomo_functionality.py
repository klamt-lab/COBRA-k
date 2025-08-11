"""Utilities to work with pyomo ConcreteModel instances directly."""

from collections.abc import Callable

from joblib import cpu_count
from numpy import linspace
from pydantic.dataclasses import dataclass
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Reals,
    SolverFactory,
    Var,
    maximize,
    minimize,
)
from pyomo.opt.base.solvers import SolverFactoryClass

from .constants import OBJECTIVE_VAR_NAME, QUASI_INF


@dataclass
class ApproximationPoint:
    """Represents a point in a linear approximation.

    This dataclass is used to store the slope, intercept, and x-coordinate of a point in a linear approximation.

    Attributes:
    - slope (float): The slope of the line passing through this point.
    - intercept (float): The y-intercept of the line passing through this point.
    - x_point (float): The x-coordinate of this point.
    """
    slope: float
    intercept: float
    x_point: float


def add_linear_approximation_to_pyomo_model(
    model: ConcreteModel,
    y_function: Callable[[float], float],
    y_function_derivative: Callable[[float], float],
    x_reference_var_id: str,
    new_y_var_name: str,
    min_x: float,
    max_x: float,
    max_rel_difference: float,
    max_num_segments: int = float("inf"),
    min_abs_error: float = 1e-6,
) -> ConcreteModel:
    """Add a linear approximation of a given function to a Pyomo model.

    This function approximates the provided function `y_function` with a piecewise linear function
    and adds the approximation to the given Pyomo model. The approximation is based on the derivative
    of the function `y_function_derivative`. The approximation is added as a new variable and a set
    of constraints to the model.

    Parameters:
    - model (ConcreteModel): The Pyomo model to which the approximation will be added.
    - y_function (Callable[[float], float]): The function to be approximated.
    - y_function_derivative (Callable[[float], float]): The derivative of the function to be approximated.
    - x_reference_var_id (str): The name of the variable in the model that will be used as the independent variable for the approximation.
    - new_y_var_name (str): The name of the new variable that will be added to the model to represent the approximation.
    - min_x (float): The minimum value of the independent variable for the approximation.
    - max_x (float): The maximum value of the independent variable for the approximation.
    - max_rel_difference (float): The maximum allowed relative difference between the approximation and the original function.
    - max_num_segments (int, optional): The maximum number of segments to use for the piecewise linear approximation. Defaults to infinity.
    - min_abs_error (float, optional): The minimum absolute error allowed between the approximation and the original function. Defaults to 1e-6.

    Returns:
    - ConcreteModel: The Pyomo model with the added approximation.
    """
    # Find fitting approximation
    num_segments = 2
    approximation_points: list[ApproximationPoint] = []
    while True:
        ignored_is = []
        x_points = linspace(min_x, max_x, num_segments)
        approximation_points = [
            ApproximationPoint(
                slope=y_function_derivative(x_point),
                intercept=y_function(x_point)
                - y_function_derivative(x_point) * x_point,
                x_point=x_point,
            )
            for x_point in x_points
        ]

        max_found_min_rel_difference = -float("inf")
        x_midpoints_data: list[tuple[int, int, float]] = []
        for i in range(len(x_points) - 1):
            first_index, second_index = i, i + 1
            if (
                approximation_points[first_index].slope
                - approximation_points[second_index].slope
                == 0
            ):
                continue
            x_midpoint = (
                approximation_points[second_index].intercept
                - approximation_points[first_index].intercept
            ) / (
                approximation_points[first_index].slope
                - approximation_points[second_index].slope
            )
            x_midpoints_data.append((first_index, second_index, x_midpoint))

        for first_index, second_index, x_value in x_midpoints_data:
            real_y = y_function(x_value)
            y_approx_one = (
                approximation_points[first_index].slope * x_value
                + approximation_points[first_index].intercept
            )
            y_approx_two = (
                approximation_points[second_index].slope * x_value
                + approximation_points[second_index].intercept
            )
            errors_absolute = (real_y - y_approx_one, real_y - y_approx_two)
            if max(errors_absolute) < min_abs_error:
                ignored_is.append(first_index)
            errors_relative = (
                abs(errors_absolute[0] / real_y),
                abs(errors_absolute[1] / real_y),
            )
            max_found_min_rel_difference = max(
                max_found_min_rel_difference, min(errors_relative)
            )

        if (max_found_min_rel_difference <= max_rel_difference) or (
            num_segments == max_num_segments
        ):
            break

        num_segments += 1
    # Add approximation to model
    min_approx_y = (
        approximation_points[0].slope * x_points[0] + approximation_points[0].intercept
    )
    max_approx_y = (
        approximation_points[-1].slope * x_points[-1]
        + approximation_points[-1].intercept
    )
    setattr(
        model, new_y_var_name, Var(within=Reals, bounds=(min_approx_y, max_approx_y))
    )
    for approx_i, approximation_point in enumerate(approximation_points):
        if approx_i in ignored_is:
            continue
        setattr(
            model,
            f"{new_y_var_name}_constraint_{approx_i}",
            Constraint(
                rule=getattr(model, new_y_var_name)
                >= approximation_point.slope * getattr(model, x_reference_var_id)
                + approximation_point.intercept
            ),
        )
    return model


def set_target_as_var_and_value(
    model: ConcreteModel,
    target: str | dict[str, float],
    var_name: str,
    constraint_name: str,
) -> tuple[ConcreteModel, Expression]:
    """Set a target as a variable and its value in a Pyomo model.

    This function adds a new variable to the given Pyomo model and sets its value to the provided target.
    The target can be either a single variable name or a dictionary of variable names with their corresponding multipliers.

    Parameters:
    - model (ConcreteModel): The Pyomo model to which the variable and constraint will be added.
    - target (str | dict[str, float]): The target for the new variable. It can be a single variable name or a dictionary of variable names with their corresponding multipliers.
    - var_name (str): The name of the new variable that will be added to the model.
    - constraint_name (str): The name of the new constraint that will be added to the model to set the value of the new variable.

    Returns:
    - tuple[ConcreteModel, Expression]: The Pyomo model with the added variable and constraint, and the expression representing the target.
    """
    if isinstance(target, str):
        expr = getattr(model, target)
    else:
        expr = 0.0
        for target_id, multiplier in target.items():  # type: ignore
            expr += multiplier * getattr(model, target_id)
    setattr(model, var_name, Var(within=Reals, bounds=(-QUASI_INF, QUASI_INF)))
    setattr(
        model,
        constraint_name,
        Constraint(expr=getattr(model, var_name) == expr),
    )
    return model, expr


def add_objective_to_model(
    model: ConcreteModel,
    objective_target: str | dict[str, float],
    objective_sense: int,
    objective_name: str,
    objective_var_name: str = OBJECTIVE_VAR_NAME,
) -> ConcreteModel:
    """Add an objective function to a Pyomo model.

    This function adds an objective function to the given Pyomo model based on the provided target and sense.
    The target can be a single variable name or a dictionary of variable names with their corresponding multipliers.
    The sense can be either maximization (as int, value > 0) or minimization (as int, value < 0).

    Parameters:
    - model (ConcreteModel): The Pyomo model to which the objective function will be added.
    - objective_target (str | dict[str, float]): The target for the objective function. It can be a single variable name or a dictionary of variable names with their corresponding multipliers.
    - objective_sense (int): The sense of the objective function. It can be an integer (positive for maximization, negative for minimization, zero for no objective).
    - objective_name (str): The name of the new objective function that will be added to the model.
    - objective_var_name (str, optional): The name of the new variable that will be added to the model to represent the objective function. Defaults to OBJECTIVE_VAR_NAME.

    Returns:
    - ConcreteModel: The Pyomo model with the added objective function.
    """
    setattr(
        model,
        objective_name,
        get_objective(
            model,
            objective_target,
            objective_sense,
            objective_var_name,
        ),
    )
    return model


def get_objective(
    model: ConcreteModel,
    objective_target: str | dict[str, float],
    objective_sense: int,
    objective_var_name: str = OBJECTIVE_VAR_NAME,
) -> Objective:
    """Create and return a pyomo objective function for the given model.

    Sets up an objective function based on the provided target and sense.
    The target can be a single variable or a weighted sum of multiple variables.
    The sense can be either maximization (as int, value > 0) or minimization (as int, value < 0).

    Parameters:
    - model (ConcreteModel): The Pyomo model to which the objective function will be added.
    - objective_target (str | dict[str, float]): The target for the objective function. It can be a single variable name or a dictionary of
                                                 variable names with their corresponding multipliers.
    - objective_sense (int): The sense of the objective function. It can be an integer
                                        (positive for maximization, negative for minimization, zero for no objective).

    Returns:
    - Objective: The Pyomo Objective object representing the objective function.
    """
    model, expr = set_target_as_var_and_value(
        model,
        objective_target,
        objective_var_name,
        "constraint_of_" + objective_var_name,
    )

    if isinstance(objective_sense, int):
        if objective_sense > 0:
            expr *= objective_sense
            pyomo_sense = maximize
        elif objective_sense < 0:
            expr *= abs(objective_sense)
            pyomo_sense = minimize
        else:  # objective_sense == 0
            expr = 0.0
            pyomo_sense = minimize
    else:
        print(f"ERROR: Objective sense is {objective_sense}, but must be an integer.")
        raise ValueError
    return Objective(expr=expr, sense=pyomo_sense)


def get_model_var_names(model: ConcreteModel) -> list[str]:
    """Extracts and returns a list of names of all variable components from a Pyomo model.

    This function iterates over all variable objects (`Var`) defined in the given Pyomo concrete model instance.
    It collects the name attribute of each variable object and returns these names as a list of strings.

    Parameters:
        model (ConcreteModel): A Pyomo concrete model instance containing various components, including variables.

    Returns:
        list[str]: A list of string names representing all variable objects in the provided Pyomo model.

    Examples:

        >>> from pyomo.environ import ConcreteModel, Var
        >>> m = ConcreteModel()
        >>> m.x = Var(initialize=1.0)
        >>> m.y = Var([1, 2], initialize=lambda m,i: i)  # Creates two variables y[1] and y[2]
        >>> var_names = get_model_var_names(m)
        >>> print(var_names)
        ['x', 'y[1]', 'y[2]']
    """
    return [v.name for v in model.component_objects(Var)]


def get_solver(
    solver_name: str,
    solver_options: dict[str, float | int | str],
    solver_attrs: dict[str, float | int | str],
) -> SolverFactoryClass:
    """Create and configure a solver for the given solver name and options.

    This function returns a Pyomo solver using the specified solver name and applies the provided options to it.

    Parameters:
    - solver_name (str): The name of the solver to be used (e.g., 'glpk', 'cbc').
    - solver_options (dict[str, float | int | str]): A dictionary of solver options where keys are option names and values are the corresponding option values.
    - solver_attrs (dict[str, float | int | str]): A dictionary of solver attributes where keys are attribute names and values are the corresponding attribute values.

    Returns:
    - SolverFactoryClass: The configured solver instance.
    """
    if solver_name == "ipopt" and cpu_count() > 16:
        solver = SolverFactory(solver_name)
    else:
        solver = SolverFactory(solver_name)
    for attr_name, attr_value in solver_attrs.items():
        setattr(solver, attr_name, attr_value)
    if not solver_name.startswith("multistart"):
        for option_name, option_value in solver_options.items():
            solver.options[option_name] = option_value
    return solver
