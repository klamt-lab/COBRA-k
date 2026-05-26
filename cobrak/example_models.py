"""Contains the toy model from COBRAk's documentation and publication as example model as well as iCH360_cobrak."""

# IMPORT SECTION
import importlib.resources as r
from math import log
from typing import Any

from .constants import STANDARD_R, STANDARD_T
from .dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)
from .io import json_load


# LAZY (ONLY LOADED ON CALL) MODEL DEFINITIONS
def __getattr__(name: str) -> Any:  # called only for missing attrs  # noqa: ANN401
    """Called by Python when an attribute that does not yet exist in the
    module namespace is requested.

    We use it to load the JSON the first time ``iCH360_cobrak`` is asked
    for, then store the result in ``globals()`` so the load happens only
    once.
    """
    if name == "iCH360_cobrak":
        data = json_load(
            str(r.files("cobrak").joinpath("data/iCH360_cobrak.json"))
        )  # <-load the file now
        globals()[name] = data  # cache for future accesses
        return data
    raise AttributeError(name)


# DIRECT EXAMPLE MODEL DEFINITION SECTION
toy_model = Model(
    reactions={
        # Metabolic reactions
        "Glycolysis": Reaction(
            # Stoichiometrically relevant member variables
            stoichiometries={
                "S": -1,  # Negative stoichiometry → S is consumed by Glycolysis
                "ADP": -2,
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
                k_cat=140_000.0,  # Turnover number in h⁻¹
                k_ms={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                    "S": 0.0001,  # e.g., K_m of reaction Glycolysis regarding metabolite A
                    "ADP": 0.0001,
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
                "ADP": -4,
                "C": +1,
                "ATP": +4,
            },
            min_flux=0.0,
            max_flux=1_000.0,
            dG0=-10.0,
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=["E_resp"],
                k_cat=140_000.0,
                k_ms={
                    "ADP": 0.00027,
                    "M": 0.00027,
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
                k_cat=140_000.0,
                k_ms={
                    "M": 0.001,
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
        # ATP "maintenance" reaction
        "ATP_Consumption": Reaction(
            stoichiometries={
                "ATP": -1,
                "ADP": +1,
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
        "ADP": Metabolite(),
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
    max_prot_pool=0.4,  # In g⋅gDW⁻¹; This value is used for our analyses with enzyme constraints
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
    kinetic_ignored_metabolite_exceptions=[],
    R=STANDARD_R,
    T=STANDARD_T,
    max_conc_sum=float("inf"),
    annotation={"description": "COBRA-k toy model"},
)


data_toy_model = Model(
    metabolites={
        # Most parameter searches need the metabolites to use BiGG
        # database IDs, which are already used in many existing
        # metabolite networks. If not, you can look them up (and
        # download the database itself) from https://bigg.ucsd.edu/
        "g6p_c": Metabolite(
            # ...additionally, for the search of ΔG'° values, you
            # have to use eQuilibrator-API-compatible
            # annotations, i.e. metabolite IDs from a multitude of databases
            # Lets define three examples (for a full list of supported
            # identifiers, check out the USED_IDENTIFIERS_FOR_EQUILIBRATOR
            # constant in cobrak.constants; note: Sometimes, INCHI strings
            # and keys cannot be read out correctly)
            annotation={
                "bigg.metabolite": "g6p",  # From the BiGG database (https://bigg.ucsd.edu/)
                "kegg.compound": "C00092",  # From the KEGG database (https://www.genome.jp/kegg/)
                "metanetx.chemical": "MNXM160",  # From MetaNetX (https://www.metanetx.org/)
                # ...again, this extra annotation is fully optional as long
                # as the eQuilibrator-API can read out your metabolite ID
            },
        ),
        # Also for f6p, we have to define an eQuilibrator-API-compatible annotation explicitly
        "f6p_c": Metabolite(
            annotation={
                "bigg.metabolite": "f6p",  # From the BiGG database (https://bigg.ucsd.edu/)
            }
        ),
    },
    reactions={
        # While this reaction uses a BiGG ID, the reaction ID can actually
        # be of any format (unlike many other IDs, as explained here)
        "PGI_fw": Reaction(
            stoichiometries={
                # Note again that the metabolites use BiGG IDs
                "g6p_c": -1.0,
                "f6p_c": 1.0,
            },
            annotation={
                # For k_cat, k_M and other kinetic reaction parameters,
                # it is *neccessary* to give the reaction a valid Enyme Commission
                # (EC) number through such an 'ec-code' annotation, which is already
                # included in many published metabolic models.
                # If you do not know the EC code of your reaction, you can try to
                # look it up through databases such as, amongst others,
                # EXPASY ENZYME (https://enzyme.expasy.org/) or also BiGG (https://bigg.ucsd.edu/)
                "ec-code": "5.3.1.9",
            },
            # As always with COBRA-k, reactions have to be *irreversible*,
            # reversible reactions have to be split up beforehand,
            # which you can e.g. automatically do for SBML models with the COBRA-k function
            # ```load_annotated_sbml_model_as_cobrak_model``` in ```cobrak.io````
            # while keeping the ```do_model_fullsplit``` argument at ```True```.
            min_flux=0.0,
            max_flux=1000.0,
            # And, again, as always with COBRA-k, not only reversible reactions have to be split
            # into seperate ones, but also reactions that are catalyzed by multiple enzymes (isozymes).
            # In this case, the reaction is split into as many variants as there are enzymes that catalyze it.
            # Again, you can do this automatically for SBML models with the COBRA-k function
            # ```load_annotated_sbml_model_as_cobrak_model``` in ```cobrak.io````
            # while keeping the ```do_model_fullsplit``` argument at ```True```.
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=["b4025"],
            ),
        ),
        # To make this model work, we'll also add pseudo-functions that deliver
        # g6p_c and take up f6p_c into the environment. As these reactions are not
        # mass-balanced (mass of substrates ≠ mass of products), thermokinetic parameters
        # and constraints do not make any sense here, so that we can omit any annotations.
        "EX_g6p_c": Reaction(
            stoichiometries={"g6p_c": +1.0}
        ),  # Produce substrate glucose-6-phosphate
        "EX_f6p_c": Reaction(
            stoichiometries={"f6p_c": -1.0}
        ),  # Take up product fructose-6-phosphate
    },
    enzymes={
        # For the automated collection of molecular weights, enzymes need an ID and/or
        # name that can be found in Uniprot (https://www.uniprot.org/). Make sure
        # that the ID and/or name is not ambiguous in your modeled organism (the
        # automatic routine only searches for enzymes of the modeled organism), so that
        # the right enzyme can be chosen.
        "b4025": Enzyme(
            name="pgi",
        ),
    },
)
