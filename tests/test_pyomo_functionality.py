import unittest

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Reals,
    Var,
)
from pyomo.solvers.plugins.solvers.GLPK import GLPKSHELL

from cobrak.pyomo_functionality import get_objective, get_solver


def test_get_objective_single_variable():
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", 1)
    assert isinstance(objective, Objective)


def test_get_objective_weighted_sum():
    model = ConcreteModel()
    model.test_var1 = Var(within=Reals, bounds=(-1, 1))
    model.test_var2 = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, {"test_var1": 1.0, "test_var2": 2.0}, 1)
    assert isinstance(objective, Objective)


def test_get_objective_minimization():
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", -1)
    assert isinstance(objective, Objective)


def test_get_objective_zero_sense():
    model = ConcreteModel()
    model.test_var = Var(within=Reals, bounds=(-1, 1))
    objective = get_objective(model, "test_var", 0)
    assert isinstance(objective, Objective)


def test_get_solver():
    solver = get_solver("glpk", {"timelimit": 600, "mipgap": 0.01}, {})
    assert isinstance(solver, GLPKSHELL)


if __name__ == "__main__":
    unittest.main()
