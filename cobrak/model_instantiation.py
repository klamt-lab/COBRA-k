"""This module contains the most convenient ways to create new Model instances from COBRApy models."""

import tempfile

# IMPORT SECTION #
from copy import deepcopy

import cobra
from cobra.manipulation import remove_genes
from numpy import log

from .brenda_functionality import brenda_select_enzyme_kinetic_data_for_model
from .cobrapy_model_functionality import get_fullsplit_cobra_model
from .constants import (
    REAC_ENZ_SEPARATOR,
    REAC_FWD_SUFFIX,
    REAC_REV_SUFFIX,
    STANDARD_R,
    STANDARD_T,
)
from .dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)
from .equilibrator_functionality import (
    equilibrator_get_model_dG0_and_uncertainty_values,
)
from .io import (
    get_files,
    json_load,
    json_write,
    load_unannotated_sbml_as_cobrapy_model,
    standardize_folder,
)
from .sabio_rk_functionality import sabio_select_enzyme_kinetic_data_for_model
from .uniprot_functionality import uniprot_get_enzyme_molecular_weights
from .utilities import (
    combine_enzyme_reaction_datasets,
    delete_orphaned_metabolites_and_enzymes,
    get_full_enzyme_mw,
    parse_external_resources,
)


# DATACLASS-USING FUNCTIONS #
def delete_enzymatically_suboptimal_reactions_in_fullsplit_cobrapy_model(
    cobra_model: cobra.Model,
    enzyme_reaction_data: dict[str, EnzymeReactionData | None],
    enzyme_molecular_weights: dict[str, float],
    fwd_suffix: str = REAC_FWD_SUFFIX,
    rev_suffix: str = REAC_REV_SUFFIX,
    reac_enz_separator: str = REAC_ENZ_SEPARATOR,
    special_enzyme_stoichiometries: dict[str, dict[str, float]] = {},
) -> cobra.Model:
    """Removes enzymatically suboptimal reactions from a fullsplit COBRApy model.

    This function identifies and deletes reactions in a COBRApy model that are enzymatically suboptimal based on
    enzyme reaction data and molecular weights. I.e., it retains only the reactions with the minimum molecular weight
    to k_cat (MW/k_cat) ratio for each base reaction. "base" reaction stands for any originally identical reaction, e.g.,
    if there are somehow now multiple phosphoglucokinase (PGK) reactions due to an enzyme fullsplit, only one of these
    PGK variants will be retained in the returned model.

    Args:
        cobra_model (cobra.Model): The COBRA model from which suboptimal reactions will be removed.
        enzyme_reaction_data (dict[str, EnzymeReactionData | None]): A dictionary mapping reaction IDs to
            ```EnzymeReactionData``` objects or ```None``` if the data is missing.
        enzyme_molecular_weights (dict[str, float]): A dictionary mapping enzyme identifiers to their molecular weights.

    Returns:
        cobra.Model: The modified COBRA model with suboptimal reactions removed.
    """
    reac_ids: list[str] = [reaction.id for reaction in cobra_model.reactions]
    ignored_reac_ids: list[str] = []
    base_reacs_to_min_mw_by_k_cat: dict[str, tuple[str, float]] = {}
    for reac_id in reac_ids:
        if reac_enz_separator not in reac_id:
            continue
        if reac_id not in enzyme_reaction_data:
            ignored_reac_ids.append(reac_id)
            continue

        current_enzyme_reaction_data = enzyme_reaction_data[reac_id]
        if current_enzyme_reaction_data is None:
            ignored_reac_ids.append(reac_id)
            continue

        mw = 0.0
        for identifier in current_enzyme_reaction_data.identifiers:
            if reac_id in special_enzyme_stoichiometries:
                if identifier in special_enzyme_stoichiometries[reac_id]:
                    stoichiometry = special_enzyme_stoichiometries[reac_id][identifier]
                else:
                    stoichiometry = 1.0
            else:
                stoichiometry = 1.0
            mw += stoichiometry * enzyme_molecular_weights[identifier]
        k_cat = current_enzyme_reaction_data.k_cat
        mw_by_k_cat = mw / k_cat

        if reac_id.endswith(fwd_suffix):
            direction_addition = fwd_suffix
        elif reac_id.endswith(rev_suffix):
            direction_addition = rev_suffix
        else:
            direction_addition = ""
        base_id = reac_id.split(reac_enz_separator)[0] + direction_addition

        if (
            base_id not in base_reacs_to_min_mw_by_k_cat
            or mw_by_k_cat < base_reacs_to_min_mw_by_k_cat[base_id][1]
        ):
            base_reacs_to_min_mw_by_k_cat[base_id] = (reac_id, mw_by_k_cat)
    enz_reacs_to_keep = [entry[0] for entry in base_reacs_to_min_mw_by_k_cat.values()]
    reacs_to_delete = [
        reac_id
        for reac_id in reac_ids
        if (reac_enz_separator in reac_id)
        and (reac_id not in enz_reacs_to_keep)
        and (reac_id not in ignored_reac_ids)
    ]
    cobra_model.remove_reactions(reacs_to_delete)
    return cobra_model


def delete_enzymatically_suboptimal_reactions_in_cobrak_model(
    cobrak_model: Model,
) -> Model:
    """Delete enzymatically suboptimal reactions in a COBRA model, similar to the idea in sMOMENT/AutoPACMEN [1].

    This function processes each reaction in the provided COBRA model to
    determine if it is enzymatically suboptimal based on Molecular Weight by k_cat (MW/kcat).
    Suboptimal reactions are identified by comparing their MW/kcat value with that of other reactions
    sharing the same base identifier, retaining only those with the lowest MW/kcat.
    The function then removes these suboptimal reactions from the model and cleans up orphaned metabolites.

    - The function assumes that the 'enzyme_reaction_data' attribute of each reaction includes
      identifiers and k_cat information for enzyme-catalyzed reactions. If not, those reactions are skipped.
    - Reactions with identical base IDs (but different directional suffixes) are considered as variants of the same reaction.
    - After removing suboptimal reactions, the function calls `delete_orphaned_metabolites_and_enzymes` to clean up any orphaned metabolites and enzymes that may have been left behind.

    [1] https://doi.org/10.1186/s12859-019-3329-9

    Parameters:
        cobrak_model (cobra.Model): A COBRA model containing biochemical reactions.

    Returns:
        cobra.Model: The updated COBRA model after removing enzymatically suboptimal reactions.
    """
    reac_id_to_mw_by_kcat: dict[str, float] = {}
    reac_id_to_base_id: dict[str, str] = {}
    base_id_to_min_mw_by_kcat: dict[str, float] = {}
    for reac_id, reac_data in cobrak_model.reactions.items():
        if reac_data.enzyme_reaction_data is None:
            continue
        if reac_data.enzyme_reaction_data.identifiers == []:
            continue
        if reac_data.enzyme_reaction_data.identifiers == [""]:
            continue

        mw_by_kcat = (
            get_full_enzyme_mw(cobrak_model, reac_data)
            / reac_data.enzyme_reaction_data.k_cat
        )

        reac_id_to_mw_by_kcat[reac_id] = mw_by_kcat

        if reac_id.endswith(cobrak_model.fwd_suffix):
            direction_addition = cobrak_model.fwd_suffix
        elif reac_id.endswith(cobrak_model.rev_suffix):
            direction_addition = cobrak_model.rev_suffix
        else:
            direction_addition = ""
        base_id = reac_id.split(cobrak_model.reac_enz_separator)[0] + direction_addition

        reac_id_to_base_id[reac_id] = base_id
        if base_id not in base_id_to_min_mw_by_kcat:
            base_id_to_min_mw_by_kcat[base_id] = mw_by_kcat
        else:
            base_id_to_min_mw_by_kcat[base_id] = min(
                base_id_to_min_mw_by_kcat[base_id], mw_by_kcat
            )

    reacs_to_delete = [
        reac_id
        for reac_id, base_id in reac_id_to_base_id.items()
        if reac_id_to_mw_by_kcat[reac_id] != base_id_to_min_mw_by_kcat[base_id]
    ]
    for reac_to_delete in reacs_to_delete:
        del cobrak_model.reactions[reac_to_delete]

    return delete_orphaned_metabolites_and_enzymes(cobrak_model)


def get_cobrak_model_from_sbml_and_thermokinetic_data(
    sbml_path: str,
    extra_linear_constraints: list[ExtraLinearConstraint],
    dG0s: dict[str, float],
    dG0_uncertainties: dict[str, float],
    conc_ranges: dict[str, tuple[float, float]],
    enzyme_molecular_weights: dict[str, float],
    enzyme_reaction_data: dict[str, EnzymeReactionData | None],
    max_prot_pool: float,
    kinetic_ignored_metabolites: list[str],
    enzyme_conc_ranges: dict[str, tuple[float, float] | None] = {},
    do_model_fullsplit: bool = False,
    do_delete_enzymatically_suboptimal_reactions: bool = True,
    R: float = STANDARD_R,
    T: float = STANDARD_T,
    fwd_suffix: str = REAC_FWD_SUFFIX,
    rev_suffix: str = REAC_REV_SUFFIX,
    reac_enz_separator: str = REAC_ENZ_SEPARATOR,
    omitted_metabolites: list[str] = [],
    keep_parameter_refs: bool = False,
) -> Model:
    """Creates a COBRAk model from an SBML and given further thermokinetic (thermodynamic and enzymatic) data.

    This function constructs a `Model` by integrating thermokinetic data and additional constraints
    into an existing COBRA model. It allows for the specification of concentration ranges, enzyme molecular weights, and
    reaction data, among other parameters.

    Args:
        sbml_path (str): The SBML model to be converted.
        extra_linear_constraints (list[ExtraLinearConstraint]): Additional linear constraints to be applied to the model.
        dG0s (dict[str, float]): Standard Gibbs free energy changes for reactions.
        dG0_uncertainties (dict[str, float]): Uncertainties in the standard Gibbs free energy changes.
        conc_ranges (dict[str, tuple[float, float]]): Concentration ranges for metabolites.
        enzyme_molecular_weights (dict[str, float]): Molecular weights of enzymes.
        enzyme_reaction_data (dict[str, EnzymeReactionData | None]): Enzyme reaction data for reactions.
        max_prot_pool (float): Maximum protein pool constraint.
        kinetic_ignored_metabolites (list[str]): Metabolites to be ignored in kinetic calculations.
        enzyme_conc_ranges (dict[str, tuple[float, float] | None], optional): Concentration ranges for enzymes. Defaults to {}.
        do_model_fullsplit (bool, optional): Whether to perform a full split of the model. Defaults to True.
        do_delete_enzymatically_suboptimal_reactions (bool, optional): Whether to delete enzymatically suboptimal reactions. Defaults to True.
        R (float, optional): Universal gas constant. Defaults to STANDARD_R.
        T (float, optional): Temperature in Kelvin. Defaults to STANDARD_T.

    Raises:
        ValueError: If a concentration range for a metabolite is not provided and no default is set.

    Returns:
        Model: The constructed `Model` with integrated thermokinetic data and constraints.
    """
    cobra_model = cobra.io.read_sbml_model(sbml_path)

    if do_model_fullsplit:
        cobra_model = get_fullsplit_cobra_model(cobra_model)

    cobrak_model = Model(
        reactions={},
        metabolites={},
        enzymes={},
        max_prot_pool=max_prot_pool,
        extra_linear_constraints=extra_linear_constraints,
        kinetic_ignored_metabolites=kinetic_ignored_metabolites,
        R=R,
        T=T,
        fwd_suffix=fwd_suffix,
        rev_suffix=rev_suffix,
        reac_enz_separator=reac_enz_separator,
    )

    for metabolite in cobra_model.metabolites:
        if metabolite.id in omitted_metabolites:
            continue

        if metabolite.id in conc_ranges:
            min_conc = conc_ranges[metabolite.id][0]
            max_conc = conc_ranges[metabolite.id][1]
        elif "DEFAULT" in conc_ranges:
            min_conc = conc_ranges["DEFAULT"][0]
            max_conc = conc_ranges["DEFAULT"][1]
        else:
            print(f"ERROR: No concentration range for metabolite {metabolite.id}.")
            print("Fixes: 1) Set its specific range; 2) Set a 'DEFAULT' range.")
            raise ValueError

        cobrak_model.metabolites[metabolite.id] = Metabolite(
            log_min_conc=log(min_conc),
            log_max_conc=log(max_conc),
            annotation=metabolite.annotation,
            name=metabolite.name,
            formula="" if not metabolite.formula else metabolite.formula,
            charge=metabolite.charge,
        )

    for reaction in cobra_model.reactions:
        dG0 = dG0s.get(reaction.id)

        dG0_uncertainty = dG0_uncertainties.get(reaction.id)

        if reaction.id in enzyme_reaction_data:
            used_enzyme_reaction_data = enzyme_reaction_data[reaction.id]
            if (used_enzyme_reaction_data is not None) and (not keep_parameter_refs):
                used_enzyme_reaction_data.k_cat_references = []
                used_enzyme_reaction_data.k_m_references = {}
                used_enzyme_reaction_data.k_i_references = {}
                used_enzyme_reaction_data.k_a_references = {}
        else:
            used_enzyme_reaction_data = None

        cobrak_model.reactions[reaction.id] = Reaction(
            min_flux=reaction.lower_bound,
            max_flux=reaction.upper_bound,
            stoichiometries={
                metabolite.id: value
                for (metabolite, value) in reaction.metabolites.items()
                if metabolite.id not in omitted_metabolites
            },
            dG0=dG0,
            dG0_uncertainty=dG0_uncertainty,
            enzyme_reaction_data=used_enzyme_reaction_data,
            annotation=reaction.annotation,
            name=reaction.name,
        )

    cobra_gene_ids = [gene.id for gene in cobra_model.genes]
    for enzyme_id, molecular_weight in enzyme_molecular_weights.items():
        min_enzyme_conc = None
        max_enzyme_conc = None
        if enzyme_id in enzyme_conc_ranges:
            conc_range = enzyme_conc_ranges[enzyme_id]
            if conc_range is not None:
                min_enzyme_conc = conc_range[0]
                max_enzyme_conc = conc_range[1]
        if enzyme_id in cobra_gene_ids:
            name = cobra_model.genes.get_by_id(enzyme_id).id
            annotation = cobra_model.genes.get_by_id(enzyme_id).annotation
        else:
            name = ""
            annotation = {}
        cobrak_model.enzymes[enzyme_id] = Enzyme(
            molecular_weight=molecular_weight,
            min_conc=min_enzyme_conc,
            max_conc=max_enzyme_conc,
            name=name,
            annotation=annotation,
        )

    if do_delete_enzymatically_suboptimal_reactions:
        cobrak_model = delete_enzymatically_suboptimal_reactions_in_cobrak_model(
            cobrak_model
        )

    return cobrak_model


def get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
    sbml_path: str,
    path_to_external_resources: str,
    folder_of_sabio_database: str,
    brenda_version: str,
    prefer_brenda: bool,
    base_species: str,
    max_prot_pool: float,
    conc_ranges: dict[str, tuple[float, float]],
    inner_to_outer_compartments: list[str],
    phs: dict[str, float],
    pmgs: dict[str, float],
    ionic_strenghts: dict[str, float],
    potential_differences: dict[tuple[str, str], float],
    kinetic_ignored_enzymes: list[str],
    custom_kms_and_kcats: dict[str, EnzymeReactionData | None],
    kinetic_ignored_metabolites: list[str],
    do_model_fullsplit: bool = True,
    do_delete_enzymatically_suboptimal_reactions: bool = True,
    ignore_dG0_uncertainty: bool = True,
    enzyme_conc_ranges: dict[str, tuple[float, float] | None] = {},
    dG0_exclusion_prefixes: list[str] = [],
    dG0_exclusion_inner_parts: list[str] = [],
    dG0_corrections: dict[str, float] = {},
    extra_linear_constraints: list[ExtraLinearConstraint] = [],
    data_cache_folder: str = "",
    R: float = STANDARD_R,
    T: float = STANDARD_T,
    keep_parameter_refs: bool = False,
    enzymes_to_delete: list[str] = [],
    max_taxonomy_level: float = float("inf"),
) -> Model:
    """This functions creates a Model out of an SBML model and automatically collects kinetic and thermodynamic data.

    To make this function work, the functions needs in the COBRApy model...
    ...for the reactions: An EC number annotation [for enzymatic constraints]
    ...for the metabolites: BiGG IDs
    ...and all the other non-optional data in its arguments

    Args:
        sbml_path (str): The SBML model file path from which the Model and associated data retrieval shall be made
        path_to_external_resources (str): _description_
        brenda_version (str): _description_
        base_species (str): _description_
        max_prot_pool (float): _description_
        conc_ranges (dict[str, tuple[float, float]]): _description_
        inner_to_outer_compartments (list[str]): _description_
        phs (dict[str, float]): _description_
        pmgs (dict[str, float]): _description_
        ionic_strenghts (dict[str, float]): _description_
        potential_differences (dict[tuple[str, str], float]): _description_
        kinetic_ignored_enzymes (list[str]): _description_
        custom_kms_and_kcats (dict[str, EnzymeReactionData  |  None]): _description_
        kinetic_ignored_metabolites (list[str]): _description_
        do_model_fullsplit (bool, optional): _description_. Defaults to True.
        do_delete_enzymatically_suboptimal_reactions (bool, optional): _description_. Defaults to True.
        ignore_dG0_uncertainty (bool, optional): _description_. Defaults to True.
        enzyme_conc_ranges (dict[str, tuple[float, float]  |  None], optional): _description_. Defaults to {}.
        dG0_exclusion_prefixes (list[str], optional): _description_. Defaults to [].
        dG0_exclusion_inner_parts (list[str], optional): _description_. Defaults to [].
        extra_linear_constraints (list[ExtraLinearConstraint], optional): _description_. Defaults to [].
        data_cache_folder (str, optional): _description_. Defaults to "".
        R (float, optional): _description_. Defaults to STANDARD_R.
        T (float, optional): _description_. Defaults to STANDARD_T.

    Returns:
        Model: _description_
    """
    cobra_model = load_unannotated_sbml_as_cobrapy_model(sbml_path)
    remove_genes(
        model=cobra_model,
        gene_list=enzymes_to_delete,
        remove_reactions=False,
    )

    path_to_external_resources = standardize_folder(path_to_external_resources)
    if data_cache_folder:
        data_cache_folder = standardize_folder(data_cache_folder)
        data_cache_files = get_files(data_cache_folder)

    parse_external_resources(path_to_external_resources, brenda_version)

    if do_model_fullsplit:
        fullsplit_model = get_fullsplit_cobra_model(cobra_model)
    else:
        fullsplit_model = deepcopy(cobra_model)

    enzyme_reaction_data: dict[str, EnzymeReactionData | None] = {}
    if (not data_cache_folder) or (
        (data_cache_folder)
        and (
            ("_cache_dG0.json" not in data_cache_files)
            or ("_cache_dG0_uncertainties.json" not in data_cache_files)
        )
    ):
        brenda_enzyme_reaction_data = brenda_select_enzyme_kinetic_data_for_model(
            cobra_model=fullsplit_model,
            brenda_json_targz_file_path=f"{path_to_external_resources}brenda_2024_1.json.tar.gz",
            bigg_metabolites_json_path=f"{path_to_external_resources}bigg_models_metabolites.json",
            brenda_version=brenda_version,
            base_species=base_species,
            ncbi_parsed_json_path=f"{path_to_external_resources}parsed_taxdmp.json",
            kinetic_ignored_metabolites=kinetic_ignored_metabolites,
            kinetic_ignored_enzyme_ids=kinetic_ignored_enzymes,
            custom_enzyme_kinetic_data=custom_kms_and_kcats,
            max_taxonomy_level=max_taxonomy_level,
        )
        sabio_enzyme_reaction_data = sabio_select_enzyme_kinetic_data_for_model(
            cobra_model=fullsplit_model,
            sabio_target_folder=folder_of_sabio_database,
            base_species=base_species,
            ncbi_parsed_json_path=f"{path_to_external_resources}parsed_taxdmp.json",
            bigg_metabolites_json_path=f"{path_to_external_resources}bigg_models_metabolites.json",
            kinetic_ignored_metabolites=kinetic_ignored_metabolites,
            kinetic_ignored_enzyme_ids=kinetic_ignored_enzymes,
            custom_enzyme_kinetic_data=custom_kms_and_kcats,
            max_taxonomy_level=max_taxonomy_level,
        )

        enzyme_reaction_data = combine_enzyme_reaction_datasets(
            [
                (
                    brenda_enzyme_reaction_data
                    if prefer_brenda
                    else sabio_enzyme_reaction_data
                ),
                (
                    sabio_enzyme_reaction_data
                    if prefer_brenda
                    else brenda_enzyme_reaction_data
                ),
            ]
        )

        if data_cache_folder:
            json_write(
                f"{data_cache_folder}_cache_enzyme_reaction_data.json",
                enzyme_reaction_data,
            )
    else:
        enzyme_reaction_data = json_load(
            f"{data_cache_folder}_cache_enzyme_reaction_data.json",
            dict[str, EnzymeReactionData | None],
        )

    enzyme_molecular_weights = uniprot_get_enzyme_molecular_weights(
        model=fullsplit_model,
        cache_basepath=data_cache_folder,
    )

    if do_delete_enzymatically_suboptimal_reactions:
        fullsplit_model = (
            delete_enzymatically_suboptimal_reactions_in_fullsplit_cobrapy_model(
                fullsplit_model,
                enzyme_reaction_data,
                enzyme_molecular_weights,
            )
        )

    if (not data_cache_folder) or (
        (data_cache_folder)
        and (
            ("_cache_dG0.json" not in data_cache_files)
            or ("_cache_dG0_uncertainties.json" not in data_cache_files)
        )
    ):
        with tempfile.TemporaryDirectory() as tmpdict:
            cobra.io.write_sbml_model(fullsplit_model, tmpdict + "temp.xml")
        dG0s, dG0_uncertainties = equilibrator_get_model_dG0_and_uncertainty_values(
            tmpdict + "temp.xml",
            inner_to_outer_compartments,
            phs,
            pmgs,
            ionic_strenghts,
            potential_differences,
            dG0_exclusion_prefixes,
            dG0_exclusion_inner_parts,
            ignore_dG0_uncertainty,
        )
        if data_cache_folder:
            json_write(f"{data_cache_folder}_cache_dG0.json", dG0s)
            json_write(
                f"{data_cache_folder}_cache_dG0_uncertainties.json", dG0_uncertainties
            )
    else:
        dG0s = json_load(f"{data_cache_folder}_cache_dG0.json", dict[str, float])
        dG0_uncertainties = json_load(
            f"{data_cache_folder}_cache_dG0_uncertainties.json",
            dict[str, float],
        )

        dG0_keys = list(dG0s.keys())
        for dG0_key in dG0_keys:
            if any(
                dG0_key.startswith(dG0_exclusion_prefix)
                for dG0_exclusion_prefix in dG0_exclusion_prefixes
            ) or any(
                dG0_exclusion_inner_part in dG0_key
                for dG0_exclusion_inner_part in dG0_exclusion_inner_parts
            ):
                del dG0s[dG0_key]
                if dG0_key in dG0_uncertainties:
                    del dG0_uncertainties[dG0_key]

    for key, value in dG0_corrections.items():
        dG0s[key] += value

    with tempfile.TemporaryDirectory() as tmpdict:
        cobra.io.write_sbml_model(fullsplit_model, tmpdict + "temp.xml")
        return get_cobrak_model_from_sbml_and_thermokinetic_data(
            sbml_path=tmpdict + "temp.xml",
            extra_linear_constraints=extra_linear_constraints,
            dG0s=dG0s,
            dG0_uncertainties=dG0_uncertainties,
            conc_ranges=conc_ranges,
            enzyme_molecular_weights=enzyme_molecular_weights,
            enzyme_reaction_data=enzyme_reaction_data,
            max_prot_pool=max_prot_pool,
            kinetic_ignored_metabolites=kinetic_ignored_metabolites,
            enzyme_conc_ranges=enzyme_conc_ranges,
            R=R,
            T=T,
            keep_parameter_refs=keep_parameter_refs,
            do_delete_enzymatically_suboptimal_reactions=False,
        )
