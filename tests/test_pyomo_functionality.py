"""pytest tests for COBRA-k's module pyomo_functionality"""

import unittest

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Reals,
    Var,
)
from pyomo.solvers.plugins.solvers.GLPK import GLPKSHELL

from cobrak.dataclasses import Solver
from cobrak.pyomo_functionality import get_objective, get_solver


def test_get_objective_single_variable() -> None:  # noqa: D103
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", 1)
    assert isinstance(objective, Objective)


def test_get_objective_weighted_sum() -> None:  # noqa: D103
    model = ConcreteModel()
    model.test_var1 = Var(within=Reals, bounds=(-1, 1))
    model.test_var2 = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, {"test_var1": 1.0, "test_var2": 2.0}, 1)
    assert isinstance(objective, Objective)


def test_get_objective_minimization() -> None:  # noqa: D103
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", -1)
    assert isinstance(objective, Objective)


def test_get_objective_zero_sense() -> None:  # noqa: D103
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", 0)
    assert isinstance(objective, Objective)


def test_get_solver() -> None:  # noqa: D103
    glpk = Solver(
        name="glpk",
        solver_options={"timelimit": 600, "mipgap": 0.01},
    )
    solver = get_solver(glpk)
    assert isinstance(solver, GLPKSHELL)


if __name__ == "__main__":
    unittest.main()
