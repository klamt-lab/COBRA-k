from cobrak.constants import LNCONC_VAR_PREFIX
from cobrak.dataclasses import Metabolite, Model, Reaction
from cobrak.tellurium_functionality import (
    _get_numbersafe_id,
    _get_reaction_string_of_cobrak_reaction,
    get_tellurium_string_from_cobrak_model,
)


def test_get_reaction_string():
    cobrak_model = Model(
        reactions={
            "reaction_id": Reaction(
                min_flux=0.0,
                max_flux=1.0,
                stoichiometries={"metabolite_id": 1.0},
                dG0=1.0,
                dG0_uncertainty=0.1,
                enzyme_reaction_data=None,
            ),
        },
        metabolites={"metabolite_id": Metabolite()},
        enzymes={},
        max_prot_pool=1.0,
        extra_linear_constraints=[],
        kinetic_ignored_metabolites=[],
        R=1.0,
        T=1.0,
    )
    reac_id = "reaction_id"
    cobrak_reaction = cobrak_model.reactions[reac_id]
    e_conc = 1.0
    met_concs = {"metabolite_id": 1.0}
    reac_flux = 1.0
    nlp_results = {LNCONC_VAR_PREFIX + "metabolite_id": 0.0}
    kinetic_ignored_metabolites = []
    unoptimized_reactions = {}
    reac_string = _get_reaction_string_of_cobrak_reaction(
        cobrak_model,
        reac_id,
        cobrak_reaction,
        e_conc,
        met_concs,
        reac_flux,
        nlp_results,
        kinetic_ignored_metabolites,
        unoptimized_reactions,
    )
    assert isinstance(reac_string, str)


def test_get_tellurium_string():
    cobrak_model = Model(
        reactions={
            "reaction_id": Reaction(
                min_flux=0.0,
                max_flux=1.0,
                stoichiometries={"metabolite_id": 1.0},
                dG0=1.0,
                dG0_uncertainty=0.1,
                enzyme_reaction_data=None,
            ),
        },
        metabolites={"metabolite_id": Metabolite()},
        enzymes={},
        max_prot_pool=1.0,
        extra_linear_constraints=[],
        kinetic_ignored_metabolites=[],
        R=1.0,
        T=1.0,
    )
    cell_density = 1.0
    e_concs = {"reaction_id": 1.0}
    met_concs = {"metabolite_id": 1.0}
    nlp_results = {"reaction_id": 1.0, LNCONC_VAR_PREFIX + "metabolite_id": 1e-6}
    tellurium_string = get_tellurium_string_from_cobrak_model(
        cobrak_model,
        cell_density,
        e_concs,
        met_concs,
        nlp_results,
    )
    assert isinstance(tellurium_string, str)


def test_get_numbersafe_id():
    met_id = "metabolite_id"
    numbersafe_id = _get_numbersafe_id(met_id)
    assert numbersafe_id == met_id

    met_id = "1metabolite_id"
    numbersafe_id = _get_numbersafe_id(met_id)
    assert numbersafe_id == "x1metabolite_id"


if __name__ == "__main__":
    test_get_reaction_string()
    test_get_tellurium_string()
    test_get_numbersafe_id()
    print("All tests passed.")
