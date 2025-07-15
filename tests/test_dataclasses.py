from cobrak.dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)


def test_model_creation():
    Model(
        metabolites={
            "X": Metabolite(0.0, 0.0, {"X": "X"}, "X", "X", 0),
        },
        enzymes={
            "X": Enzyme(0.0),
        },
        reactions={
            "X": Reaction(
                {"X": 0.0},
                0.0,
                0.0,
                0.0,
                0.0,
                EnzymeReactionData(
                    ["X"],
                    0.1,
                ),
                {"X": "X"},
                "X",
            ),
        },
        extra_linear_constraints=[
            ExtraLinearConstraint(
                {"X": 0.0},
                0.0,
                0.0,
            )
        ],
    )
