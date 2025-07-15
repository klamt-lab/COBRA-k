from cobrak.cobrapy_model_functionality import (
    create_irreversible_cobrapy_model_from_stoichiometries,
)


def test_create_irreversible_cobrapy_model_from_stoichiometries():
    model = create_irreversible_cobrapy_model_from_stoichiometries(
        {
            "EX_A": {"A": +1},
            "R1": {"A": -1, "B": +1},
            "EX_B": {"B": -1},
        }
    )
    assert len(model.reactions) == 3
    assert len(model.metabolites) == 2
