"""Contains methods that directly apply on COBRApy models."""

# IMPORTS SECTION #
from copy import deepcopy

import cobra

from .constants import REAC_ENZ_SEPARATOR, REAC_FWD_SUFFIX, REAC_REV_SUFFIX


# FUNCTIONS SECTION #
def get_fullsplit_cobra_model(
    cobra_model: cobra.Model,
) -> cobra.Model:
    """Return a COBRApy model where reactions are split according to reversibility and enzymes.

    "Reversibility" means that, if a reaction i can run in both directions (α_i<0), then it is split as follows:
    Ri: A<->B [-50;100]=> Ri_FWD: A->B [0;100]; Ri_REV: B->A [0;50]
    where the ending "FWD" and "REV" are set in COBRAk's constants REAC_FWD_SUFFIX and REAC_REV_SUFFIX.

    "enzymes" means that, if a reaction i can be catalyzed by multiple enzymes (i.e., at least one OR block in the
    reaction's gene-protein rule), then it is split for each reaction. Say, for example,
    Rj: A->B [0;100]
    has the following gene-protein rule:
    (E1 OR E2)
    ...then, Rj is split into:
    Rj_ENZ_E1: A->B [0;100]
    Rj_ENZ_E2: A->B [0;100]
    where the infix "_ENZ_" is set in COBRAk's constants REAC_ENZ_SEPARATOR.

    Args:
        cobra_model (cobra.Model): The COBRApy model that shall be 'fullsplit'.

    Returns:
        cobra.Model: The 'fullsplit' COBRApy model.
    """
    fullsplit_cobra_model = cobra.Model(cobra_model.id)

    fullsplit_cobra_model.add_metabolites(cobra_model.metabolites)

    for gene in cobra_model.genes:
        fullsplit_cobra_model.genes.add(deepcopy(gene))

    for reaction_x in cobra_model.reactions:
        reaction: cobra.Reaction = reaction_x

        is_reversible = False
        if reaction.lower_bound < 0.0:  # and (not reaction.id.startswith("EX_")):
            is_reversible = True

        single_enzyme_blocks = (
            reaction.gene_reaction_rule.replace("(", "").replace(")", "").split(" or ")
        )

        for single_enzyme_block in single_enzyme_blocks:
            if single_enzyme_block:
                new_reac_base_id = (
                    reaction.id
                    + REAC_ENZ_SEPARATOR
                    + single_enzyme_block.replace(" ", "_")
                )
            else:
                new_reac_base_id = reaction.id
            new_reaction_1 = cobra.Reaction(
                id=new_reac_base_id,
                lower_bound=reaction.lower_bound,
                upper_bound=reaction.upper_bound,
            )
            new_reaction_1.annotation = deepcopy(reaction.annotation)
            if single_enzyme_block:
                new_reaction_1.gene_reaction_rule = single_enzyme_block
            new_reaction_1_met_addition = {}
            for met, stoichiometry in reaction.metabolites.items():
                new_reaction_1_met_addition[met] = stoichiometry
            new_reaction_1.add_metabolites(new_reaction_1_met_addition)

            if is_reversible:
                original_lb = new_reaction_1.lower_bound
                new_reaction_2 = cobra.Reaction(
                    id=new_reac_base_id,
                )
                new_reaction_2.annotation = deepcopy(reaction.annotation)
                if single_enzyme_block:
                    new_reaction_2.gene_reaction_rule = single_enzyme_block
                new_reaction_1.id += REAC_FWD_SUFFIX
                new_reaction_1.lower_bound = 0
                new_reaction_2.id += REAC_REV_SUFFIX
                new_reaction_2.lower_bound = 0
                new_reaction_2.upper_bound = abs(original_lb)

                new_reaction_2_met_addition = {}
                for met, stoichiometry in new_reaction_1.metabolites.items():
                    new_reaction_2_met_addition[met] = -stoichiometry
                new_reaction_2.add_metabolites(new_reaction_2_met_addition)
                fullsplit_cobra_model.add_reactions([new_reaction_2])
            fullsplit_cobra_model.add_reactions([new_reaction_1])

    return fullsplit_cobra_model


def create_irreversible_cobrapy_model_from_stoichiometries(
    stoichiometries: dict[str, dict[str, float]],
) -> cobra.Model:
    """Create an irreversible COBRApy model out of the given dictionary.

    E.g., if the following dict is the argument:
    {
        "EX_A": { "A": +1 },
        "R1": { "A": -1, "B": +1 },
        "EX_B": { "B": -1 },
    }
    ...then, the following three irreversible (i.e, flux from 0 to 1_000) reactions
    are created and returned as a single COBRApy model:
    EX_A: -> A
    R1: A -> B
    EX_B: B ->

    Args:
        stoichiometries (dict[str, dict[str, float]]): The model-describing dictionary

    Returns:
        cobra.Model: The resulting COBRApy model with the given reactions and metabolites
    """
    cobra_model: cobra.Model = cobra.Model()
    reac_ids = stoichiometries.keys()
    metabolite_ids_list = []
    for stoichiometry_entry in stoichiometries.values():
        metabolite_ids_list.extend(list(stoichiometry_entry.keys()))
    metabolite_ids = set(metabolite_ids_list)
    cobra_model.add_metabolites(
        [cobra.Metabolite(id=met_id, compartment="c") for met_id in metabolite_ids]
    )
    cobra_model.add_reactions(
        [
            cobra.Reaction(
                id=reac_id,
                name=reac_id,
                lower_bound=0.0,
                upper_bound=1000.0,
            )
            for reac_id in reac_ids
        ]
    )
    for reac_id in reac_ids:
        reaction: cobra.Reaction = cobra_model.reactions.get_by_id(reac_id)
        reaction.add_metabolites(
            {
                cobra_model.metabolites.get_by_id(met_id): stoichiometry
                for met_id, stoichiometry in stoichiometries[reac_id].items()
            }
        )

    return cobra_model
