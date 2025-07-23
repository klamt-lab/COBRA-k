"""Contains methods that directly apply on COBRApy models."""

# IMPORTS SECTION #
from copy import deepcopy
from dataclasses import asdict

import cobra

from .constants import (
    REAC_ENZ_SEPARATOR,
    REAC_FWD_SUFFIX,
    REAC_REV_SUFFIX,
    STANDARD_R,
    STANDARD_T,
)
from .dataclasses import ExtraLinearConstraint


# FUNCTIONS SECTION #
def get_fullsplit_cobra_model(
    cobra_model: cobra.Model,
    fwd_suffix: str = REAC_FWD_SUFFIX,
    rev_suffix: str = REAC_REV_SUFFIX,
    add_cobrak_sbml_annotation: bool = False,
    cobrak_default_min_conc: float = 1e-6,
    cobrak_default_max_conc: float = 0.2,
    cobrak_extra_linear_constraints: list[ExtraLinearConstraint] = [],
    cobrak_kinetic_ignored_metabolites: list[str] = [],
    cobrak_no_extra_versions: bool = False,
    reac_lb_ub_cap: float = float("inf"),
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

    if add_cobrak_sbml_annotation:
        settings_reac = cobra.Reaction(
            id="cobrak_global_settings",
            name="Global COBRA-k settings",
            lower_bound=0.0,
            upper_bound=0.0,
        )
        settings_reac.annotation["cobrak_max_prot_pool"] = 1000.0
        settings_reac.annotation["cobrak_R"] = STANDARD_R
        settings_reac.annotation["cobrak_T"] = STANDARD_T
        settings_reac.annotation["cobrak_kinetic_ignored_metabolites"] = {}
        settings_reac.annotation["cobrak_reac_rev_suffix"] = (
            rev_suffix  # A "special" suffix to show that this is added
        )
        settings_reac.annotation["cobrak_reac_fwd_suffix"] = fwd_suffix
        settings_reac.annotation["cobrak_reac_enz_separator"] = REAC_ENZ_SEPARATOR
        settings_reac.annotation["cobrak_extra_linear_constraints"] = str(
            [asdict(x) for x in cobrak_extra_linear_constraints]
        )
        settings_reac.annotation["cobrak_kinetic_ignored_metabolites"] = str(
            cobrak_kinetic_ignored_metabolites
        )

        fullsplit_cobra_model.add_reactions([settings_reac])

    fullsplit_cobra_model.add_metabolites(cobra_model.metabolites)

    for gene in cobra_model.genes:
        fullsplit_cobra_model.genes.add(deepcopy(gene))

    for reaction_x in cobra_model.reactions:
        reaction: cobra.Reaction = reaction_x

        if add_cobrak_sbml_annotation:
            for old_name, new_name in (
                ("dG0", "cobrak_dG0"),
                ("dG0_uncertainty", "cobrak_dG0_uncertainty"),
            ):
                if old_name in reaction.annotation:
                    reaction.annotation[new_name] = reaction.annotation[old_name]

            fwd_dG0 = (
                float(reaction.annotation["cobrak_dG0"])
                if "cobrak_dG0" in reaction.annotation
                else None
            )
            dG0_uncertainty = (
                abs(float(reaction.annotation["cobrak_dG0_uncertainty"]))
                if "cobrak_dG0_uncertainty" in reaction.annotation
                else None
            )

        is_reversible = False
        if reaction.lower_bound < 0.0:
            is_reversible = True

        single_enzyme_blocks = (
            reaction.gene_reaction_rule.replace("(", "").replace(")", "").split(" or ")
        )
        current_reac_version = 0
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
                upper_bound=min(reac_lb_ub_cap, reaction.upper_bound),
            )
            new_reaction_1.annotation = deepcopy(reaction.annotation)
            if add_cobrak_sbml_annotation:
                if fwd_dG0 is not None:
                    new_reaction_1.annotation[f"cobrak_dG0_V{current_reac_version}"] = (
                        fwd_dG0
                    )
                if dG0_uncertainty is not None:
                    new_reaction_1.annotation[
                        f"cobrak_dG0_uncertainty_V{current_reac_version}"
                    ] = dG0_uncertainty
                new_reaction_1.annotation[f"cobrak_id_V{current_reac_version}"] = (
                    new_reaction_1.id + (fwd_suffix if is_reversible else "")
                )
            if single_enzyme_block:
                new_reaction_1.gene_reaction_rule = single_enzyme_block
            new_reaction_1_met_addition = {}
            for met, stoichiometry in reaction.metabolites.items():
                new_reaction_1_met_addition[met] = stoichiometry
            new_reaction_1.add_metabolites(new_reaction_1_met_addition)

            if is_reversible:
                current_reac_version += 1

                original_lb = new_reaction_1.lower_bound
                new_reaction_2 = cobra.Reaction(
                    id=new_reac_base_id,
                )
                new_reaction_2.annotation = deepcopy(reaction.annotation)
                if add_cobrak_sbml_annotation:
                    if fwd_dG0 is not None:
                        new_reaction_2.annotation[
                            f"cobrak_dG0_V{current_reac_version}"
                        ] = -fwd_dG0
                    if dG0_uncertainty is not None:
                        new_reaction_2.annotation[
                            f"cobrak_dG0_uncertainty_V{current_reac_version}"
                        ] = dG0_uncertainty
                    new_reaction_2.annotation[f"cobrak_id_V{current_reac_version}"] = (
                        new_reaction_2.id + rev_suffix
                    )
                if single_enzyme_block:
                    new_reaction_2.gene_reaction_rule = single_enzyme_block
                new_reaction_1.id += fwd_suffix
                new_reaction_1.lower_bound = 0
                new_reaction_2.id += rev_suffix
                new_reaction_2.lower_bound = 0
                new_reaction_2.upper_bound = min(reac_lb_ub_cap, abs(original_lb))

                new_reaction_2_met_addition = {}
                for met, stoichiometry in new_reaction_1.metabolites.items():
                    new_reaction_2_met_addition[met] = -stoichiometry
                new_reaction_2.add_metabolites(new_reaction_2_met_addition)

                fullsplit_cobra_model.add_reactions([new_reaction_2])
            fullsplit_cobra_model.add_reactions([new_reaction_1])
            current_reac_version += 1
            if cobrak_no_extra_versions and (
                ("cobrak_k_cat_V0" not in reaction.annotation)
                or ("cobrak_k_cat" not in reaction.annotation)
            ):
                break

    for metabolite in fullsplit_cobra_model.metabolites:
        for old_name, new_name in (("Cmin", "cobrak_Cmin"), ("Cmax", "cobrak_Cmax")):
            if old_name in metabolite.annotation:
                metabolite.annotation[new_name] = metabolite.annotation[old_name]
        if "cobrak_Cmin" not in metabolite.annotation:
            metabolite.annotation["cobrak_Cmin"] = cobrak_default_min_conc
        if "cobrak_Cmax" not in metabolite.annotation:
            metabolite.annotation["cobrak_Cmax"] = cobrak_default_max_conc

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
