import pytest

from cobrak.constants import (
    ALL_OK_KEY,
    ENZYME_VAR_INFIX,
    OBJECTIVE_VAR_NAME,
    REAC_FWD_SUFFIX,
    REAC_REV_SUFFIX,
    SOLVER_STATUS_KEY,
    TERMINATION_CONDITION_KEY,
)
from cobrak.dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)
from cobrak.utilities import (
    add_statuses_to_optimziation_dict,
    compare_optimization_result_reaction_uses,
    delete_orphaned_metabolites_and_enzymes,
    delete_unused_reactions_in_optimization_dict,
    delete_unused_reactions_in_variability_dict,
    get_active_reacs_from_optimization_dict,
    get_base_id,
    get_base_id_optimzation_result,
    get_cobrak_enzyme_reactions_string,
    get_extra_linear_constraint_string,
    get_full_enzyme_id,
    get_full_enzyme_mw,
    get_fwd_rev_corrected_flux,
    get_metabolites_in_elementary_conservation_relations,
    get_potentially_active_reactions_in_variability_dict,
    get_reaction_enzyme_var_id,
    get_reaction_string,
    get_solver_status_from_pyomo_results,
    get_stoichiometric_matrix,
    get_stoichiometrically_coupled_reactions,
    get_substrate_and_product_exchanges,
    get_termination_condition_from_pyomo_results,
    have_all_unignored_km,
    is_objsense_maximization,
    last_n_elements_equal,
    sort_dict_keys,
)


# Example fixtures and mock objects
@pytest.fixture
def mock_optimization_dict():
    return {
        "R1": 1.0,
        "R2": 2.0,
        "R3": 3.0,
        OBJECTIVE_VAR_NAME: 10.0,
        ALL_OK_KEY: True,
        TERMINATION_CONDITION_KEY: 0.2,
        SOLVER_STATUS_KEY: 0,
    }


@pytest.fixture
def mock_pyomo_results():
    class MockSolverResults:
        solver = type(
            "Solver", (object,), {"status": "ok", "termination_condition": "optimal"}
        )()

    return MockSolverResults()


@pytest.fixture
def mock_model():
    return Model(
        reactions={
            "R1": Reaction(
                stoichiometries={"A": -1, "B": 1},
                min_flux=0,
                max_flux=1000,
                dG0=0.0,
                dG0_uncertainty=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E1", "E1B"],
                    k_cat=1.0,
                    k_ms={"A": 0.1},
                    k_is={},
                    k_as={},
                    special_stoichiometries={},
                ),
            ),
            "R2": Reaction(
                stoichiometries={"B": -1, "C": 1},
                min_flux=0,
                max_flux=1000,
                dG0=0.0,
                dG0_uncertainty=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=[
                        "E2",
                    ],
                    k_cat=1.0,
                    k_ms={"B": 0.2},
                    k_is={},
                    k_as={},
                    special_stoichiometries={},
                ),
            ),
        },
        metabolites={
            "A": Metabolite(name="A"),
            "B": Metabolite(name="B"),
            "C": Metabolite(name="C"),
        },
        enzymes={
            "E1": Enzyme(
                molecular_weight=10.0,
            ),
            "E1B": Enzyme(
                molecular_weight=10.0,
            ),
            "E2": Enzyme(
                molecular_weight=1.0,
            ),
        },
        max_prot_pool=1000.0,
        R=8.314,
        T=298.15,
        kinetic_ignored_metabolites=[],
    )


@pytest.fixture
def mock_variability_dict():
    return {
        "R1": (0.0, 10.0),
        "R2": (-5.0, 5.0),
        "R3": (0.0, 0.0),
    }


@pytest.fixture
def mock_datasets():
    return [
        {
            "R1": EnzymeReactionData(
                identifiers=["E1"],
                k_cat=1.0,
                k_ms={"A": 0.1},
                k_is={},
                k_as={},
            ),
            "R2": None,
        },
        {
            "R1": None,
            "R2": EnzymeReactionData(
                identifiers=["E2"],
                k_cat=1.0,
                k_ms={},
                k_is={},
                k_as={},
            ),
        },
    ]


# Test functions
def test_add_statuses_to_optimziation_dict(mock_optimization_dict, mock_pyomo_results):
    result = add_statuses_to_optimziation_dict(
        mock_optimization_dict, mock_pyomo_results
    )
    assert result[SOLVER_STATUS_KEY] == 0
    assert result[TERMINATION_CONDITION_KEY] == 0.2
    assert result[ALL_OK_KEY]


def test_get_base_id():
    assert get_base_id(f"R1{REAC_FWD_SUFFIX}") == "R1"
    assert get_base_id(f"R2{REAC_REV_SUFFIX}") == "R2"
    assert get_base_id("R3") == "R3"


def test_get_base_id_optimzation_result(mock_model, mock_optimization_dict):
    result = get_base_id_optimzation_result(mock_model, mock_optimization_dict)
    assert result["R1"] == 1.0
    assert result["R2"] == 2.0


def test_get_fwd_rev_corrected_flux():
    # Test case where reverse reaction exists and has greater flux, and custom fwd and rev
    result = {"R1_fwd": 10.0, "R1_rev": 15.0}
    assert (
        get_fwd_rev_corrected_flux(
            "R1_fwd", {"R1_fwd", "R1_rev"}, result, fwd_suffix="_fwd", rev_suffix="_rev"
        )
        == 0.0
    )

    # Test case where reverse reaction exists and has lesser flux
    result = {f"R1_{REAC_FWD_SUFFIX}": 20.0, f"R1_{REAC_REV_SUFFIX}": 5.0}
    assert (
        get_fwd_rev_corrected_flux(
            f"R1_{REAC_FWD_SUFFIX}",
            {f"R1_{REAC_FWD_SUFFIX}", f"R1_{REAC_REV_SUFFIX}"},
            result,
        )
        == 15.0
    )

    # Test case where reverse reaction does not exist
    result = {f"R1_{REAC_FWD_SUFFIX}": 30.0}
    assert (
        get_fwd_rev_corrected_flux(
            f"R1_{REAC_FWD_SUFFIX}", {f"R1_{REAC_FWD_SUFFIX}"}, result
        )
        == 30.0
    )

    # Test case where reverse reaction is not usable
    result = {f"R1_{REAC_FWD_SUFFIX}": 25.0, f"R1_{REAC_REV_SUFFIX}": 10.0}
    assert (
        get_fwd_rev_corrected_flux(
            f"R1_{REAC_FWD_SUFFIX}", {f"R1_{REAC_FWD_SUFFIX}"}, result
        )
        == 25.0
    )


def test_get_solver_status_from_pyomo_results(mock_pyomo_results):
    assert get_solver_status_from_pyomo_results(mock_pyomo_results) == 0


def test_get_termination_condition_from_pyomo_results(mock_pyomo_results):
    assert get_termination_condition_from_pyomo_results(mock_pyomo_results) == 0.2


def test_sort_dict_keys():
    input_dict = {"b": 2, "a": 1, "c": 3}
    result = sort_dict_keys(input_dict)
    assert result == {"a": 1, "b": 2, "c": 3}


def test_compare_optimization_result_reaction_uses(mock_model, mock_optimization_dict):
    results = [
        mock_optimization_dict,
        {"R1": 0.5, "R2": 1.5, "R3": 2.5, "objective": 15.0},
    ]
    compare_optimization_result_reaction_uses(mock_model, results)
    # This function prints results, so we can't directly assert the output.
    # You can manually check the printed output or redirect stdout to capture it.


def test_delete_orphaned_metabolites_and_enzymes(mock_model):
    model = delete_orphaned_metabolites_and_enzymes(mock_model)
    assert "A" in model.metabolites
    assert "B" in model.metabolites
    assert "C" in model.metabolites
    assert "E1" in model.enzymes


def test_delete_unused_reactions_in_optimization_dict(
    mock_model, mock_optimization_dict
):
    model = delete_unused_reactions_in_optimization_dict(
        mock_model, mock_optimization_dict
    )
    assert "R1" in model.reactions
    assert "R2" in model.reactions
    assert "R3" not in model.reactions


def test_delete_unused_reactions_in_variability_dict(mock_model, mock_variability_dict):
    model = delete_unused_reactions_in_variability_dict(
        mock_model, mock_variability_dict
    )
    assert "R1" in model.reactions
    assert "R2" in model.reactions
    assert "R3" not in model.reactions


def test_get_active_reacs_from_optimization_dict(mock_model, mock_optimization_dict):
    result = get_active_reacs_from_optimization_dict(mock_model, mock_optimization_dict)
    assert result == ["R1", "R2"]


def test_get_cobrak_enzyme_reactions_string(mock_model):
    result = get_cobrak_enzyme_reactions_string(mock_model, "E1")
    assert result == "R1"


def test_get_reaction_string(mock_model):
    result = get_reaction_string(mock_model, "R1")
    assert result == "-1.0 A \u21d2 1.0 B"


def test_get_extra_linear_constraint_string():
    constraint = ExtraLinearConstraint(
        lower_value=0.0,
        upper_value=10.0,
        stoichiometries={"var1": 1.0, "var2": -1.0},
    )
    result = get_extra_linear_constraint_string(constraint)
    assert result == "0.0 \u2264  + 1.0 var1 - 1.0 var2\u2264 10.0"


def test_get_full_enzyme_id():
    result = get_full_enzyme_id(["E1", "E2"])
    assert result == "E1_AND_E2"


def test_get_full_enzyme_mw(mock_model):
    result = get_full_enzyme_mw(mock_model, mock_model.reactions["R1"])
    assert result == 20.0


def test_get_metabolites_in_elementary_conservation_relations(mock_model):
    result = get_metabolites_in_elementary_conservation_relations(mock_model)
    assert sorted(result) == ["A", "B", "C"]


def test_get_potentially_active_reactions_in_variability_dict(
    mock_model, mock_variability_dict
):
    result = get_potentially_active_reactions_in_variability_dict(
        mock_model, mock_variability_dict
    )
    assert result == ["R1", "R2"]


def test_get_reaction_enzyme_var_id(mock_model):
    result = get_reaction_enzyme_var_id("R1", mock_model.reactions["R1"])
    assert result == f"enzyme_E1_AND_E1B{ENZYME_VAR_INFIX}R1"


def test_get_stoichiometric_matrix(mock_model):
    result = get_stoichiometric_matrix(mock_model)
    assert result == [
        [-1, 0.0],
        [1, -1],
        [0.0, 1],
    ]


def test_get_stoichiometrically_coupled_reactions(mock_model):
    result = get_stoichiometrically_coupled_reactions(mock_model)
    assert result == [["R1", "R2"]]


def test_get_substrate_and_product_exchanges(mock_model):
    result = get_substrate_and_product_exchanges(mock_model)
    assert result == (tuple(), tuple())


def test_have_all_unignored_km(mock_model):
    result = have_all_unignored_km(
        mock_model.reactions["R1"], mock_model.kinetic_ignored_metabolites
    )
    assert not result


def test_is_objsense_maximization():
    assert is_objsense_maximization(1)
    assert not is_objsense_maximization(-1)


def test_last_n_elements_equal():
    # Test cases where the last n elements are equal
    assert last_n_elements_equal([1, 2, 3, 4, 4, 4], 3)
    assert last_n_elements_equal(["a", "b", "b", "b"], 3)

    # Test cases where the last n elements are not equal
    assert not last_n_elements_equal([1, 2, 3, 4, 5, 6], 3)
    assert not last_n_elements_equal(["a", "b", "c", "d"], 2)

    # Test cases with fewer elements than n
    assert not last_n_elements_equal([1, 2], 3)
    assert not last_n_elements_equal([], 1)

    # Edge case: n = 0 (should always return True)
    assert last_n_elements_equal([1, 2, 3], 0)
    assert last_n_elements_equal([], 0)
