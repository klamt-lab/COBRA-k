"""Currently just contains the toy model form COBRAk's documentation as example model"""

# IMPORT SECTION
from math import log

from .constants import STANDARD_R, STANDARD_T
from .dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)

# EXAMPLE MODEL DEFINITION SECTION
toy_model = Model(
    reactions={
        # Metabolic reactions
        "Glycolysis": Reaction(
            # Stoichiometrically relevant member variables
            stoichiometries={
                "S": -1,  # Negative stoichiometry → S is consumed by Glycolysis
                "M": +1,  # Positive stoichiometry → M is produced by Glycolysis
                "ATP": +2,  # Two ATP molecules are produced by Glycolysis
            },
            min_flux=0.0,  # Minimal flux in mmol⋅gDW⁻¹⋅h⁻¹; should be ≥0 for most analyses
            max_flux=1_000.0,  # Maximal flux in mmol⋅gDW⁻¹⋅h⁻¹
            # Thermodynamically relevant member variables
            # (only neccessary if thermodynamic constraints are used)
            dG0=-10.0,  # Standard Gibb's free energy ΔG'° in kJ⋅mol⁻¹; Default is None (no ΔG'°)
            dG0_uncertainty=None,  # ΔG'° uncertainty in kJ⋅mol⁻¹; Default is None (no uncertainty)
            # Let's set the variable for enzyme-kinetic parameters
            # of the dataclass EnzymeReactionData
            # (Default is None, i.e. no enzyme parameters given)
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=[
                    "E_glyc"
                ],  # Subunit(s) which constitute the reaction's catalyst
                k_cat=70_000.0,  # Turnover number in h⁻¹
                k_ms={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                    "S": 0.0001,  # e.g., K_m of reaction Glycolysis regarding metabolite A
                    "M": 0.0001,
                    "ATP": 0.0001,
                },
                special_stoichiometries={},  # No special stoichiometry, all subunits occur once
            ),
            # Extra information member variables
            annotation={"description": "This is reaction Glycolysis"},  # Default is {}
            name="Reaction Glycolysis",  # Default is ""
        ),
        "Respiration": Reaction(
            stoichiometries={
                "M": -1,
                "C": +1,
                "ATP": +6,
            },
            min_flux=0.0,
            max_flux=1_000.0,
            dG0=-10.0,
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=["E_resp"],
                k_cat=70_000.0,
                k_ms={
                    "M": 0.03,
                    "C": 0.0001,
                    "ATP": 0.0001,
                },
            ),
        ),
        "Overflow": Reaction(
            stoichiometries={
                "M": -1,
                "P": +1,
            },
            min_flux=0.0,
            max_flux=1_000.0,
            dG0=-10.0,
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=["E_over"],
                k_cat=70_000.0,
                k_ms={
                    "M": 0.005,
                    "P": 0.0001,
                },
            ),
        ),
        # Exchange reactions
        "EX_S": Reaction(
            stoichiometries={
                "S": +1,
            },
            min_flux=0.0,
            max_flux=1_000.0,
        ),
        "EX_C": Reaction(
            stoichiometries={
                "C": -1.0,
            },
            min_flux=0.0,
            max_flux=1_000.0,
        ),
        "EX_P": Reaction(
            stoichiometries={
                "P": -1,
            },
            min_flux=0.0,
            max_flux=1_000.0,
        ),
        "EX_ATP": Reaction(
            stoichiometries={
                "ATP": -1,
            },
            min_flux=0.0,
            max_flux=1_000.0,
        ),
    },
    metabolites={
        "S": Metabolite(
            log_min_conc=log(
                1e-6
            ),  # optional, minimal ln(concentration); Default is ln(1e-6 M)
            log_max_conc=log(
                0.02
            ),  # optional, maximal ln(concentration); Default is ln(0.02 M)
            annotation={
                "description": "This is metabolite S"
            },  # optional, default is ""
            name="Metabolite S",  # optional, default is ""
            formula="X",  # optional, default is ""
            charge=0,  # optional, default is 0
        ),
        "M": Metabolite(),
        "C": Metabolite(),
        "P": Metabolite(),
        "ATP": Metabolite(),
    },
    enzymes={
        "E_glyc": Enzyme(
            molecular_weight=1_000.0,  # Molecular weight in kDa
            min_conc=None,  # Optional concentration in mmol⋅gDW⁻¹; Default is None (minimum is 0)
            max_conc=None,  # Optional maximal concentration in mmol⋅gDW⁻¹; Default is None (only protein pool restricts)
            annotation={"description": "Enzyme of Glycolysis"},  # Default is {}
            name="Glycolysis enzyme",  # Default is ""
        ),
        "E_resp": Enzyme(molecular_weight=2_500.0),
        "E_over": Enzyme(molecular_weight=500.0),
    },
    max_prot_pool=0.25,  # In g⋅gDW⁻¹; This value is used for our analyses with enzyme constraints
    # We set the following two constraints:
    # 1.0 * EX_A - 1.0 * Glycolysis ≤ 0.0
    # and
    # 1.0 * EX_A + 1.0 * Glycolysis ≥ 0.0
    # in other words, effectively,
    # 1.0 * EX_A = 1.0 * Glycolysis
    extra_linear_constraints=[
        ExtraLinearConstraint(
            stoichiometries={
                "EX_S": -1.0,
                "Glycolysis": 1.0,
            },
            lower_value=0.0,
            upper_value=0.0,
        )
    ],  # Keep in mind that this is a list as multiple extra flux constraints are possible
    kinetic_ignored_metabolites=[],
    R=STANDARD_R,
    T=STANDARD_T,
    max_conc_sum=float("inf"),
    annotation={"description": "COBRAk toy model"},
)
