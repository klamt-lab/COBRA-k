"""Tests for cobrak.cobrapy_model_functionality."""

import cobra
import pytest

from cobrak.cobrapy_model_functionality import get_fullsplit_cobra_model


def _get_model(bounds: list[tuple[float, float]]) -> cobra.Model:
    model = cobra.Model("test_model")
    metabolite = cobra.Metabolite("A")
    model.add_metabolites([metabolite])
    for i, (lb, ub) in enumerate(bounds):
        reaction = cobra.Reaction(f"R{i}", lower_bound=lb, upper_bound=ub)
        model.add_reactions([reaction])
        reaction.add_metabolites({metabolite: -1})
    return model


@pytest.mark.parametrize(
    ("lb", "ub", "expected_fwd", "expected_rev"),
    [
        (-10, 5, (0, 5), (0, 10)),
        (-10, 0, (0, 0), (0, 10)),
        (-3, 3, (0, 3), (0, 3)),
        (-10, -2, (0, 0), (2, 10)),
        (-1000, -2, (0, 0), (2, 1000)),
    ],
)
def test_split_reversible_bounds(lb, ub, expected_fwd, expected_rev) -> None:
    split_model = get_fullsplit_cobra_model(_get_model([(lb, ub)]))
    assert split_model.reactions.get_by_id("R0_FWD").bounds == expected_fwd
    assert split_model.reactions.get_by_id("R0_REV").bounds == expected_rev


def test_irreversible_reaction_not_split() -> None:
    split_model = get_fullsplit_cobra_model(_get_model([(0, 5)]))
    assert split_model.reactions.get_by_id("R0").bounds == (0, 5)
    assert "R0_REV" not in split_model.reactions


@pytest.mark.parametrize(
    ("lb", "ub", "expected_fwd", "expected_rev"),
    [
        (-1000, 2000, (0, 100), (0, 100)),
        (-1000, -2, (0, 0), (2, 100)),
        (-1000, -200, (0, 0), (100, 100)),
    ],
)
def test_split_bounds_with_cap(lb, ub, expected_fwd, expected_rev) -> None:
    split_model = get_fullsplit_cobra_model(_get_model([(lb, ub)]), reac_lb_ub_cap=100)
    assert split_model.reactions.get_by_id("R0_FWD").bounds == expected_fwd
    assert split_model.reactions.get_by_id("R0_REV").bounds == expected_rev
