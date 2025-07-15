"""This script is a wrapper for the ΔG'° determination with the eQuilibrator API.

This wrapper intends to work with BiGG-styled cobrapy metabolic models.
"""

# IMPORTS SECTION #
from typing import Any

import cobra
from equilibrator_api import Q_, ComponentContribution, Reaction

from .constants import USED_IDENTIFIERS_FOR_EQUILIBRATOR


## PUBLIC FUNCTIONS ##
def equilibrator_get_model_dG0_and_uncertainty_values(
    sbml_path: str,
    inner_to_outer_compartments: list[str],
    phs: dict[str, float],
    pmgs: dict[str, float],
    ionic_strengths: dict[str, float],
    potential_differences: dict[tuple[str, str], float],
    exclusion_prefixes: list[str] = [],
    exclusion_inner_parts: list[str] = [],
    ignore_uncertainty: bool = False,
    max_uncertainty: float = 1_000.0,
) -> tuple[dict[str, float], dict[str, float]]:
    """Cobrapy model wrapper for the ΔG'° determination of reactions using the eQuilibrator-API.

    Reactions are identified according to all annotation (in the cobrapy reaction's annotation member variable)
    given in this modules global USED_IDENTIFIERS list.

    Args:
        sbml_path (str): The path to the SBML-encoded constraint-based metabolic model for which ΔG'° values are determined.
        inner_to_outer_compartments (List[str]): A list with compartment IDs going from inner (e.g., in E. coli,
            the cytosol or 'c' in iML1515) to outer (e.g., the extracellular component or 'e' in iML1515). Used
            for the ΔG'° calculation in multi-compartmental reactions.
        phs (Dict[str, float]): A dictionary with compartment IDs as keys and the compartment pHs as values.
        pmgs (Dict[str, float]): A dictionary with compartment IDs as keys and the compartment pMgs as values.
        ionic_strengths (Dict[str, float]): A dictionary with compartment IDs as keys and the ionic strengths as values.
        potential_differences (Dict[Tuple[str, str], float]): A dictionary containing tuples with 2 elements describing
            the ID of an innter and outer compartment, and the potential difference between them.
        max_uncertainty (float): The maximal accepted uncertainty value (defaults to 1000 kJ⋅mol⁻¹). If a calculated uncertainty
            is higher than this value, the associated ΔG'° is *not* used (i.e., the specific reaction gets no ΔG'°).

    Returns:
        Dict[str, Dict[str, float]]: A dictionary with the reaction IDs as keys, and dictionaries as values which,
            in turn, contain the ΔG'° of a reaction under the key 'dG0' and the calculated uncertainty as 'uncertainty'.
    """
    cobra_model = cobra.io.read_sbml_model(sbml_path)

    reaction_dG0s: dict[str, float] = {}
    reaction_dG0_uncertainties: dict[str, float] = {}
    cc = ComponentContribution()
    for reaction_x in cobra_model.reactions:
        reaction: cobra.Reaction = reaction_x

        stop = False
        for exclusion_prefix in exclusion_prefixes:
            if reaction.id.startswith(exclusion_prefix):
                stop = True
        for exclusion_inner_part in exclusion_inner_parts:
            if exclusion_inner_part in reaction.id:
                stop = True
        if stop:
            continue

        stoichiometries: list[float] = []
        compartments: list[str] = []
        identifiers: list[str] = []
        identifier_keys: list[str] = []
        for metabolite_x in reaction.metabolites.keys():
            metabolite: cobra.Metabolite = metabolite_x
            stoichiometries.append(reaction.metabolites[metabolite])
            compartments.append(metabolite.compartment)
            identifier = ""
            for used_identifier in USED_IDENTIFIERS_FOR_EQUILIBRATOR:
                if used_identifier not in metabolite.annotation.keys():
                    continue
                metabolite_identifiers = metabolite.annotation[used_identifier]
                identifier_temp = ""
                if isinstance(metabolite_identifiers, list):
                    identifier_temp = metabolite_identifiers[0]
                elif isinstance(metabolite_identifiers, str):
                    identifier_temp = metabolite_identifiers
                if used_identifier == "inchi":
                    compound = cc.get_compound_by_inchi(identifier_temp)
                elif used_identifier == "inchi_key":
                    compound_list = cc.search_compound_by_inchi_key(identifier_temp)
                    if len(compound_list) > 0:
                        compound = compound_list[0]
                    else:
                        compound = None
                else:
                    identifier_temp = used_identifier + ":" + identifier_temp
                    compound = cc.get_compound(identifier_temp)
                if compound is not None:
                    identifier_key = used_identifier
                    identifier = identifier_temp
                    break
            if identifier == "":
                break
            identifier_keys.append(identifier_key)
            identifiers.append(identifier)

        if identifier == "":
            print(
                f"ERROR: Metabolite {metabolite_x.id} has no identifier of the given types!"
            )
            print(metabolite_x.annotation)
            continue

        # Check for three cases:
        # 1: Single-compartment reaction
        # 2: Double-compartment reaction
        # 3: Multi-compartment reaction (not possible)
        unique_reaction_compartments = list(set(compartments))
        num_compartments = len(unique_reaction_compartments)
        if num_compartments == 1:
            # Set compartment conditions
            compartment = unique_reaction_compartments[0]
            cc.p_h = Q_(phs[compartment])
            cc.p_mg = Q_(pmgs[compartment])
            cc.ionic_strength = Q_(str(ionic_strengths[compartment]) + "mM")

            # Build together reaction
            reaction_dict: dict[Any, float] = {}
            for i in range(len(stoichiometries)):
                identifier_string = identifiers[i]
                identifier_key = identifier_keys[i]
                stoichiometry = stoichiometries[i]
                if identifier_key == "inchi":
                    compound = cc.get_compound_by_inchi(identifier_string)
                elif identifier_key == "inchi_key":
                    compound = cc.search_compound_by_inchi_key(identifier_string)[0]
                else:
                    compound = cc.get_compound(identifier_string)
                reaction_dict[compound] = stoichiometry
            cc_reaction = Reaction(reaction_dict)

            # Check whether or not the reaction is balanced and...
            if not cc_reaction.is_balanced():
                print(f"INFO: Reaction {reaction.id} is not balanced")
                continue

            standard_dg_prime = cc.standard_dg_prime(cc_reaction)
            uncertainty = standard_dg_prime.error.m_as("kJ⋅mol⁻¹")
            if uncertainty < max_uncertainty:
                dG0 = standard_dg_prime.value.m_as("kJ⋅mol⁻¹")
                reaction_dG0s[reaction.id] = dG0
                if ignore_uncertainty:
                    reaction_dG0_uncertainties[reaction.id] = 0.0
                else:
                    reaction_dG0_uncertainties[reaction.id] = abs(uncertainty)

                print(
                    f"No error with reaction {reaction.id}, ΔG'° succesfully calculated!"
                )
            else:
                print(
                    f"INFO: Reaction {reaction.id} uncertainty is too high with {uncertainty} kJ⋅mol⁻¹; ΔG'° not assigned for this reaction"
                )
        elif num_compartments == 2:
            index_zero = inner_to_outer_compartments.index(
                unique_reaction_compartments[0]
            )
            index_one = inner_to_outer_compartments.index(
                unique_reaction_compartments[1]
            )

            if index_one > index_zero:
                outer_compartment = unique_reaction_compartments[1]
                inner_compartment = unique_reaction_compartments[0]
            else:
                outer_compartment = unique_reaction_compartments[0]
                inner_compartment = unique_reaction_compartments[1]

            ph_inner = Q_(phs[inner_compartment])
            ph_outer = Q_(phs[outer_compartment])
            ionic_strength_inner = Q_(str(ionic_strengths[inner_compartment]) + " mM")
            ionic_strength_outer = Q_(str(ionic_strengths[outer_compartment]) + " mM")
            pmg_inner = Q_(pmgs[inner_compartment])
            pmg_outer = Q_(pmgs[outer_compartment])

            if (inner_compartment, outer_compartment) in potential_differences:
                potential_difference = Q_(
                    str(potential_differences[(inner_compartment, outer_compartment)])
                    + " V"
                )
            elif (outer_compartment, inner_compartment) in potential_differences:
                potential_difference = Q_(
                    str(potential_differences[(outer_compartment, inner_compartment)])
                    + " V"
                )
            else:
                print("ERROR")
                continue

            inner_reaction_dict: dict[str, float] = {}
            outer_reaction_dict: dict[str, float] = {}
            for i in range(len(stoichiometries)):
                key = identifiers[i]
                stoichiometry = stoichiometries[i]
                try:
                    compound_key = cc.get_compound(key)
                except Exception:  # sqlalchemy.orm.exc.MultipleResultsFound
                    print("ERROR")
                    continue

                if compound_key is None:
                    print("NONE in compound")
                    continue

                if compartments[i] == inner_compartment:
                    inner_reaction_dict[compound_key] = stoichiometry
                else:
                    outer_reaction_dict[compound_key] = stoichiometry

            cc_inner_reaction = Reaction(inner_reaction_dict)
            cc_outer_reaction = Reaction(outer_reaction_dict)

            cc.p_h = ph_inner
            cc.ionic_strength = ionic_strength_inner
            cc.p_mg = pmg_inner
            try:
                standard_dg_prime = cc.multicompartmental_standard_dg_prime(
                    cc_inner_reaction,
                    cc_outer_reaction,
                    e_potential_difference=potential_difference,
                    p_h_outer=ph_outer,
                    p_mg_outer=pmg_outer,
                    ionic_strength_outer=ionic_strength_outer,
                )
                uncertainty = standard_dg_prime.error.m_as("kJ⋅mol⁻¹")
                if uncertainty < max_uncertainty:
                    dG0 = standard_dg_prime.value.m_as("kJ⋅mol⁻¹")
                    reaction_dG0s[reaction.id] = dG0
                    if ignore_uncertainty:
                        reaction_dG0_uncertainties[reaction.id] = 0.0
                    else:
                        reaction_dG0_uncertainties[reaction.id] = abs(uncertainty)
            except ValueError:
                print("ERROR: Multi-compartmental reaction is not balanced")
                continue
        else:
            print("ERROR: More than two compartments are not possible")
            continue

    return reaction_dG0s, reaction_dG0_uncertainties
