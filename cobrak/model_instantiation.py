"""This module contains the most convenient ways to create new Model instances from COBRApy models."""

import logging
import tempfile

# IMPORT SECTION #
from copy import deepcopy
from os.path import exists

import cobra
from cobra.manipulation import remove_genes
from numpy import log

from cobrak.molmass_functionality import add_molar_masses_to_model_metabolites

from .brenda_functionality import brenda_select_enzyme_kinetic_data_for_sbml
from .cobrapy_model_functionality import get_fullsplit_cobra_model
from .constants import (
    EC_INNER_TO_OUTER_COMPARTMENTS,
    EC_IONIC_STRENGTHS,
    EC_PHS,
    EC_PMGS,
    EC_POTENTIAL_DIFFERENCES,
    REAC_ENZ_SEPARATOR,
    REAC_FWD_SUFFIX,
    REAC_REV_SUFFIX,
    STANDARD_CONC_RANGES,
    STANDARD_MAX_PROT_POOL,
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
    equilibrator_get_model_dG0_and_uncertainty_values_for_sbml,
)
from .expasy_functionality import get_ec_number_transfers
from .io import (
    get_files,
    json_load,
    json_write,
    load_unannotated_sbml_as_cobrapy_model,
    standardize_folder,
)
from .sabio_rk_functionality import sabio_select_enzyme_kinetic_data_for_sbml
from .uniprot_functionality import (
    uniprot_get_enzyme_molecular_weights_for_sbml,
    uniprot_get_enzyme_sequences_for_sbml,
)
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
        cobra_model (cobra.Model): The COBRA-k model from which suboptimal reactions will be removed.
        enzyme_reaction_data (dict[str, EnzymeReactionData | None]): A dictionary mapping reaction IDs to
            ```EnzymeReactionData``` objects or ```None``` if the data is missing.
        enzyme_molecular_weights (dict[str, float]): A dictionary mapping enzyme identifiers to their molecular weights.

    Returns:
        cobra.Model: The modified COBRA-k model with suboptimal reactions removed.
    """
    reac_ids: list[str] = [reaction.id for reaction in cobra_model.reactions]
    ignored_reac_ids_with_mws: list[str] = []
    base_reacs_to_min_mw_by_k_cat: dict[str, tuple[str, float]] = {}
    for reac_id in reac_ids:
        if reac_enz_separator not in reac_id:
            continue
        if reac_id not in enzyme_reaction_data:
            enzyme_ids = reac_id.split(reac_enz_separator)[1].split("_and_")
            enzyme_ids[-1] = (
                enzyme_ids[-1].replace(rev_suffix, "").replace(fwd_suffix, "")
            )
            if not all(
                enzyme_id in enzyme_molecular_weights for enzyme_id in enzyme_ids
            ):
                ignored_reac_ids_with_mws.append((reac_id, 0.0))
                continue
            enzyme_reaction_data[reac_id] = EnzymeReactionData(identifiers=enzyme_ids)

        try:
            current_enzyme_reaction_data = enzyme_reaction_data[reac_id]
        except KeyError:
            logging.warning(f"The dict enzyme_reaction_data does not have {reac_id}")  # noqa: G004, LOG015
            continue
        if current_enzyme_reaction_data is None:
            ignored_reac_ids_with_mws.append((reac_id, 0.0))
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
            try:
                mw += stoichiometry * enzyme_molecular_weights[identifier]
            except KeyError:
                logging.warning(f"Cannot find {identifier} in enzyme_molecular_weights")  # noqa: G004, LOG015

        k_cat = current_enzyme_reaction_data.k_cat
        if k_cat > 1e19:
            ignored_reac_ids_with_mws.append((reac_id, mw))
            continue

        if reac_id.endswith(fwd_suffix):
            direction_addition = fwd_suffix
        elif reac_id.endswith(rev_suffix):
            direction_addition = rev_suffix
        else:
            direction_addition = ""
        base_id = reac_id.split(reac_enz_separator)[0] + direction_addition

        mw_by_k_cat = mw / k_cat
        if (
            base_id not in base_reacs_to_min_mw_by_k_cat
            or mw_by_k_cat < base_reacs_to_min_mw_by_k_cat[base_id][1]
        ):
            base_reacs_to_min_mw_by_k_cat[base_id] = (reac_id, mw_by_k_cat)
    enz_reacs_to_keep = [entry[0] for entry in base_reacs_to_min_mw_by_k_cat.values()]

    # Remove superfluous reactions
    extra_reacs_to_delete = _extra_reacs_to_delete(
        ignored_reac_ids=ignored_reac_ids_with_mws,
        enz_reacs_to_keep=enz_reacs_to_keep,
        rev_suffix=rev_suffix,
        fwd_suffix=fwd_suffix,
        reac_enz_separator=reac_enz_separator,
    )

    reacs_to_delete = [
        reac_id
        for reac_id in reac_ids
        if (reac_enz_separator in reac_id)
        and (reac_id not in enz_reacs_to_keep)
        and (reac_id not in [item[0] for item in ignored_reac_ids_with_mws])
    ] + extra_reacs_to_delete
    cobra_model.remove_reactions(reacs_to_delete)
    return cobra_model


def _extra_reacs_to_delete(
    ignored_reac_ids: list[tuple[str, float]],
    enz_reacs_to_keep: list[str],
    rev_suffix: str,
    fwd_suffix: str,
    reac_enz_separator: str,
) -> list[str]:
    ignored_base_id_suffix_to_reac_ids_and_mws: dict[
        tuple[str, bool], list[tuple[str, float]]
    ] = {}
    for reac_id, mw in ignored_reac_ids:
        base_id = reac_id.split(reac_enz_separator)[0]
        is_rev = reac_id.endswith(rev_suffix)
        tuple_ = (base_id, is_rev)
        if tuple_ not in ignored_base_id_suffix_to_reac_ids_and_mws:
            ignored_base_id_suffix_to_reac_ids_and_mws[tuple_] = []
        ignored_base_id_suffix_to_reac_ids_and_mws[tuple_].append((reac_id, mw))
    extra_reacs_to_delete: list[str] = []
    for reacidlist in ignored_base_id_suffix_to_reac_ids_and_mws.values():
        if len(reacidlist) <= 1:
            continue
        allowed_indices = []
        for i, (singleid, _) in enumerate(reacidlist):
            if (singleid.replace(rev_suffix, fwd_suffix) in enz_reacs_to_keep) or (
                reac_id.replace(fwd_suffix, rev_suffix) in enz_reacs_to_keep
            ):
                allowed_indices.append(i)
        if allowed_indices == []:
            mwlist = [item[1] for item in reacidlist]
            allowed_indices = [mwlist.index(min(mwlist))]
        for i, (singleid, _) in enumerate(reacidlist):
            if i not in allowed_indices:
                extra_reacs_to_delete.append(singleid)
    return extra_reacs_to_delete


def delete_enzymatically_suboptimal_reactions_in_cobrak_model(
    cobrak_model: Model,
    ignored_ids: list[str] = ["s0001"],
    enz_reacs_to_keep: list[str] = [],
) -> Model:
    """Delete enzymatically suboptimal reactions in a COBRA-k model, similar to the idea in sMOMENT/AutoPACMEN [1].

    This function processes each reaction in the provided COBRA-k model to
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
        cobrak_model (cobra.Model): A COBRA-k model containing biochemical reactions.

    Returns:
        cobra.Model: The updated COBRA-k model after removing enzymatically suboptimal reactions.
    """
    reac_id_to_mw_by_kcat: dict[str, float] = {}
    reac_id_to_base_id: dict[str, str] = {}
    base_id_to_min_mw_by_kcat: dict[str, float] = {}
    ignored_reac_ids_with_mws: list[str] = []
    for reac_id, reac_data in cobrak_model.reactions.items():
        if (
            reac_data.enzyme_reaction_data is None
            or reac_data.enzyme_reaction_data.identifiers in ([], [""])
        ):
            ignored_reac_ids_with_mws.append((reac_id, 0.0))
            continue
        if reac_data.enzyme_reaction_data.k_cat >= 1e19:
            ignored_reac_ids_with_mws.append(
                (reac_id, get_full_enzyme_mw(cobrak_model, reac_data))
            )
            continue
        if any(
            ignored_id in reac_data.enzyme_reaction_data.identifiers
            for ignored_id in ignored_ids
        ):
            ignored_reac_ids_with_mws.append((reac_id, 0.0))
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
    extra_reacs_to_delete = _extra_reacs_to_delete(
        ignored_reac_ids=ignored_reac_ids_with_mws,
        enz_reacs_to_keep=enz_reacs_to_keep,
        rev_suffix=cobrak_model.rev_suffix,
        fwd_suffix=cobrak_model.fwd_suffix,
        reac_enz_separator=cobrak_model.reac_enz_separator,
    )
    for reac_to_delete in reacs_to_delete + extra_reacs_to_delete:
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
    max_prot_pool: float = STANDARD_MAX_PROT_POOL,
    kinetic_ignored_metabolites: list[str] = [],
    kinetic_ignored_metabolite_exceptions: list[tuple[str, str]] = [],
    enzyme_conc_ranges: dict[str, tuple[float, float] | None] = {},
    do_model_fullsplit: bool = False,
    do_delete_enzymatically_suboptimal_reactions: bool = True,
    R: float = STANDARD_R,
    T: float = STANDARD_T,
    fwd_suffix: str = REAC_FWD_SUFFIX,
    rev_suffix: str = REAC_REV_SUFFIX,
    reac_enz_separator: str = REAC_ENZ_SEPARATOR,
    omitted_metabolites: list[str] = [],
    ignored_enzyme_ids: str = ["s0001"],
    remove_enzyme_reaction_data_if_no_kcat_set: bool = False,
    sequences: dict[str, str] = {},
    add_molar_masses: bool = True,
) -> Model:
    """Creates a COBRAk model from an SBML and given further thermokinetic (thermodynamic and enzymatic) data.

    This function constructs a `Model` by integrating thermokinetic data and additional constraints
    into an existing COBRA-k model. It allows for the specification of concentration ranges, enzyme molecular weights, and
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
        omitted_metabolites (list[str], optional): Metabolites that shall not be included in the model. Their stoichiometries
         will be jsut deleted. Useful to e.g. delete enzyme-constraint pseudo-metabolites. Defauls to [].
        ignored_enzyme_ids (list[str], optional): Enzymes that shall not be included if their ID occurs in any identifiers part. Defaults to ["s0001"],
         i.e. spontaneously occurring reactions.
        remove_enzyme_reaction_data_if_no_kcat_set (bool, optional): If no $k_{cat}$ is set for a reaction, shall its EnzymeReactionData
         be set to None? If False, the default EnzymeReactionData with a very high (effectively non-existing) $k_{cat}$ is used. Defaults to False.
         sequences (dict[str, str], optional): Data for protein sequences
        add_molar_masses: bool, default True
            Whether or not to calculate molar masses for all metabolites through their formula member variable
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
        kinetic_ignored_metabolite_exceptions=kinetic_ignored_metabolite_exceptions,
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
            annotation={
                key: value
                for key, value in metabolite.annotation.items()
                if not key.startswith("cobrak_")
            },
            name=metabolite.name,
            formula="" if not metabolite.formula else metabolite.formula,
            charge=metabolite.charge,
            compartment=metabolite.compartment,
        )

    for reaction in cobra_model.reactions:
        dG0 = dG0s.get(reaction.id)

        dG0_uncertainty = dG0_uncertainties.get(reaction.id)

        used_enzyme_reaction_data = enzyme_reaction_data.get(reaction.id, None)
        if used_enzyme_reaction_data is None:
            identifiers = reaction.gene_reaction_rule.split(" and ")
            used_enzyme_reaction_data = (
                EnzymeReactionData(
                    identifiers=identifiers,
                )
                if identifiers != [""]
                else None
            )

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
            annotation={
                key: value
                for key, value in reaction.annotation.items()
                if not key.startswith("cobrak_")
            },
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
        sequence = sequences.get(enzyme_id, "")
        cobrak_model.enzymes[enzyme_id] = Enzyme(
            molecular_weight=molecular_weight,
            min_conc=min_enzyme_conc,
            max_conc=max_enzyme_conc,
            name=name,
            annotation=annotation,
            sequence=sequence,
        )

    if do_delete_enzymatically_suboptimal_reactions:
        cobrak_model = delete_enzymatically_suboptimal_reactions_in_cobrak_model(
            cobrak_model,
            ignored_ids=ignored_enzyme_ids,
        )

    if remove_enzyme_reaction_data_if_no_kcat_set:
        for reaction in cobrak_model.reactions.values():
            if reaction.enzyme_reaction_data is None:
                continue
            if reaction.enzyme_reaction_data.k_cat > 1e19:
                reaction.enzyme_reaction_data = None

    if add_molar_masses:
        cobrak_model = add_molar_masses_to_model_metabolites(cobrak_model)

    return cobrak_model


def get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
    sbml_path: str,
    database_data_folder: str,
    brenda_version: str,
    base_species: str,
    prefer_brenda: bool = False,
    use_ec_number_transfers: bool = True,
    max_prot_pool: float = STANDARD_MAX_PROT_POOL,
    conc_ranges: dict[str, tuple[float, float]] = STANDARD_CONC_RANGES,
    inner_to_outer_compartments: list[str] = EC_INNER_TO_OUTER_COMPARTMENTS,
    phs: dict[str, float] = EC_PHS,
    pmgs: dict[str, float] = EC_PMGS,
    ionic_strenghts: dict[str, float] = EC_IONIC_STRENGTHS,
    potential_differences: dict[tuple[str, str], float] = EC_POTENTIAL_DIFFERENCES,
    kinetic_ignored_enzymes: list[str] = [],
    custom_kms_and_kcats: dict[str, EnzymeReactionData | None] = {},
    kinetic_ignored_metabolites: list[str] = [],
    kinetic_ignored_metabolite_exceptions: list[tuple[str, str]] = [],
    do_model_fullsplit: bool = True,
    do_delete_enzymatically_suboptimal_reactions: bool = True,
    ignore_dG0_uncertainty: bool = True,
    enzyme_conc_ranges: dict[str, tuple[float, float] | None] = {},
    dG0_exclusion_prefixes: list[str] = [],
    dG0_exclusion_inner_parts: list[str] = [],
    dG0_corrections: dict[str, float] = {},
    extra_linear_constraints: list[ExtraLinearConstraint] = [],
    R: float = STANDARD_R,
    T: float = STANDARD_T,
    enzymes_to_delete: list[str] = [],
    max_taxonomy_level: float = 1_000.0,
    add_hill_coefficients: bool = True,
    add_protein_sequences: bool = False,
    kis_and_kas_only_for_same_compartments: bool = False,
    add_molar_masses: bool = True,
) -> Model:
    """Build a fully-featured :class:`~cobrak.Model` from an SBML file **and** automatically
    retrieve all required kinetic and thermodynamic data from the local
    ``database_data_folder`` (or download it on-the-fly if missing).

    The function orchestrates a multi-step pipeline:

    1. **Load the SBML** as an un-annotated COBRApy model and optionally delete
       user-specified enzymes (genes) from the model.
    2. **Prepare the external data cache** – ensure that the folder structure
       exists, locate cached JSON files, and (re)generate missing caches.
    3. **Parse EC-number transfers** (optional) to allow cross-species mapping of
       enzyme identifiers.
    4. **Create a “full-split” model** where each enzyme-specific reaction variant
       is represented as a separate COBRApy reaction (controlled by
       ``do_model_fullsplit``).
    5. **Collect enzyme kinetic parameters** from BRENDA and SABIO-RK, optionally
       preferring one source over the other, and combine the two datasets.
    6. **Fetch enzyme molecular weights** from UniProt (cached for future runs).
    7. **Optionally prune sub-optimal enzyme reactions** based on the
       8. **Compute standard Gibbs free energies** (ΔG⁰) and their uncertainties
       using eQuilibrator, applying user-defined compartment, pH, ionic-strength,
       and membrane-potential settings, as well as any exclusion rules.
    9. **Apply user-provided ΔG⁰ corrections** (e.g. literature adjustments).
    10. **Assemble the final COBRA-k model** by calling
        :func:`get_cobrak_model_from_sbml_and_thermokinetic_data` with all
        gathered data, then clean up orphaned metabolites/enzymes.

    Parameters
    ----------
    sbml_path : str
        Path to the SBML file that will be converted into a COBRA-k model.
    database_data_folder : str
        Root folder containing cached kinetic, thermodynamic and annotation data.
        The function will create the folder if it does not exist.
    brenda_version : str
        Version identifier of the BRENDA JSON archive (e.g. ``"2023.1"``).
    base_species : str
        NCBI taxonomy identifier (or scientific name) of the organism for which
        kinetic data should be retrieved.
    prefer_brenda : bool, optional
        If ``True`` BRENDA data are used preferentially when both BRENDA and
        SABIO-RK contain information for the same reaction; otherwise SABIO-RK
        is preferred. Default: ``False``.
    use_ec_number_transfers : bool, optional
        Enable mapping of EC numbers between organisms using the
        ``enzyme.rdf`` file from Expasy. Default: ``True``.
    max_prot_pool : float, optional
        Upper bound on the total protein mass (g·gDW⁻¹) that can be allocated to
        enzymes. Default: :data:`STANDARD_MAX_PROT_POOL`.
    conc_ranges : dict[str, tuple[float, float]], optional
        Log-linear concentration bounds for metabolites (in M). Keys are metabolite
        IDs; the special key ``"DEFAULT"`` provides a fallback range. Default:
        :data:`STANDARD_CONC_RANGES`.
    inner_to_outer_compartments : list[str], optional
        Mapping of inner to outer compartments required by eQuilibrator for
        ΔG⁰ calculations. Default: :data:`EC_INNER_TO_OUTER_COMPARTMENTS`.
    phs : dict[str, float], optional
        pH values for each compartment. Default: :data:`EC_PHS`.
    pmgs : dict[str, float], optional
        Magnesium concentrations (M) for each compartment. Default: :data:`EC_PMGS`.
    ionic_strenghts : dict[str, float], optional
        Ionic strength (M) for each compartment. Default: :data:`EC_IONIC_STRENGTHS`.
    potential_differences : dict[tuple[str, str], float], optional
        Membrane potential differences (V) between compartment pairs. Default:
        :data:`EC_POTENTIAL_DIFFERENCES`.
    kinetic_ignored_enzymes : list[str], optional
        Enzyme identifiers that should be ignored when extracting kinetic data.
        Default: ``[]``.
    custom_kms_and_kcats : dict[str, EnzymeReactionData | None], optional
        User-provided kinetic parameters that override any database values.
        Default: ``{}``.
    kinetic_ignored_metabolites : list[str], optional
        Metabolite IDs that shall be excluded from kinetic calculations
        (e.g., pseudo-metabolites). Default: ``[]``.
    do_model_fullsplit : bool, optional
        Whether to split reactions per enzyme before further processing.
        Default: ``True``.
    do_delete_enzymatically_suboptimal_reactions : bool, optional
        If ``True`` remove reactions that are not optimal with respect to the
        ``MW/k_cat`` criterion. Default: ``True``.
    ignore_dG0_uncertainty : bool, optional
        When ``True`` discard ΔG⁰ uncertainty values after they have been computed.
        Default: ``True``.
    enzyme_conc_ranges : dict[str, tuple[float, float] | None], optional
        Optional concentration bounds for enzymes (in M). ``None`` means no bound.
        Default: ``{}``.
    dG0_exclusion_prefixes : list[str], optional
        Reaction IDs starting with any of these prefixes are removed from the
        ΔG⁰ dataset. Default: ``[]``.
    dG0_exclusion_inner_parts : list[str], optional
        Sub-strings that, if present anywhere in a reaction ID, cause its ΔG⁰
        entry to be removed. Default: ``[]``.
    dG0_corrections : dict[str, float], optional
        Additive corrections (in kJ·mol⁻¹) to specific ΔG⁰ values after they have
        been computed. Default: ``{}``.
    extra_linear_constraints : list[ExtraLinearConstraint], optional
        Additional linear constraints (e.g., flux bounds) to be added to the model.
        Default: ``[]``.
    R : float, optional
        Universal gas constant (kJ·mol⁻¹·K⁻¹). Default: :data:`STANDARD_R`.
    T : float, optional
        Temperature in Kelvin for thermodynamic calculations. Default:
        :data:`STANDARD_T`.
    enzymes_to_delete : list[str], optional
        Gene identifiers that should be removed from the initial COBRApy model
        before any further processing. Default: ``[]``.
    max_taxonomy_level : float, optional
        Upper bound on the NCBI taxonomy distance used when selecting kinetic
        data from related organisms. Default: ``1_000.0``.
    add_hill_coefficients : bool, optional
        If ``True`` include Hill coefficients from SABIO-RK where available.
        Default: ``True``.
    add_protein_sequences: bool, optional
        Whether to add protein sequences or not to Enzyme instances. Default: ``False``
    kis_and_kas_only_for_same_compartments: bool, default False
        If True, kis and kas can only be attributed to a reaction if the affected metabolite has
        shares one of the reaction metabolite's compartments
    add_molar_masses: bool, default True
        Whether or not to calculate molar masses for all metabolites through their formula member variable

    Returns
    -------
    Model
        A fully populated :class:`~cobrak.Model` instance containing e.g.:
        * Metabolite objects with concentration bounds,
        * Reaction objects with flux bounds, ΔG⁰ values, and enzyme reaction data,
        * Enzyme objects with molecular weights and concentration bounds,
        * Any extra linear constraints supplied by the user,
        * The global protein pool constraint.

    Raises
    ------
    FileNotFoundError
        If ``sbml_path`` does not exist or required external files (e.g.
        ``enzyme.rdf`` when ``use_ec_number_transfers`` is ``True``) are missing.
    ValueError
        When a required concentration range for a metabolite is not provided and
        no ``"DEFAULT"`` range exists.
    RuntimeError
        If any of the external data retrieval steps (BRENDA, SABIO-RK,
        UniProt, eQuilibrator) fail unexpectedly.

    Notes
    -----
    * The function heavily relies on caching to avoid repeated expensive web
      queries. Cache files are stored alongside ``database_data_folder`` with
      ``_cache_`` prefixes.
    * The returned model is already cleaned of orphaned metabolites and enzymes
      via :func:`delete_orphaned_metabolites_and_enzymes`.
    * Users can bypass the full pipeline by providing pre-computed cache files;
      in that case the function will simply load the cached data.
    """
    cobra_model = load_unannotated_sbml_as_cobrapy_model(sbml_path)
    remove_genes(
        model=cobra_model,
        gene_list=enzymes_to_delete,
        remove_reactions=False,
    )

    database_data_folder = standardize_folder(database_data_folder)
    data_cache_files = get_files(database_data_folder)

    parse_external_resources(database_data_folder, brenda_version)
    if use_ec_number_transfers:
        transfer_json_path = f"{database_data_folder}ec_number_transfers.json"
        if not exists(transfer_json_path):
            if not exists(f"{database_data_folder}enzyme.rdf"):
                print(
                    f"ERROR: Argument use_ec_number_transfers is True, but no necessary enzyme.rdf can be found in {database_data_folder}"
                )
                print(
                    "You may download it from https://ftp.expasy.org/databases/enzyme/"
                )
                print(
                    f"After downloading, put it into the folder {database_data_folder}"
                )
            ec_number_transfers = get_ec_number_transfers(
                f"{database_data_folder}enzyme.rdf"
            )
            json_write(transfer_json_path, ec_number_transfers)
    else:
        transfer_json_path = ""

    fullsplit_model = (
        get_fullsplit_cobra_model(cobra_model)
        if do_model_fullsplit
        else deepcopy(cobra_model)
    )

    enzyme_reaction_data: dict[str, EnzymeReactionData | None] = {}
    if (not database_data_folder) or (
        (database_data_folder)
        and (
            ("_cache_dG0.json" not in data_cache_files)
            or ("_cache_dG0_uncertainties.json" not in data_cache_files)
            or ("_cache_enzyme_reaction_data.json" not in data_cache_files)
        )
    ):
        with tempfile.TemporaryDirectory() as tmpdict:
            temp_sbml_path = tmpdict + "temp.xml"
            cobra.io.write_sbml_model(fullsplit_model, temp_sbml_path)

            brenda_enzyme_reaction_data = brenda_select_enzyme_kinetic_data_for_sbml(
                sbml_path=temp_sbml_path,
                brenda_json_targz_file_path=f"{database_data_folder}brenda_{brenda_version}.json.tar.gz",
                bigg_metabolites_json_path=f"{database_data_folder}bigg_models_metabolites.json",
                brenda_version=brenda_version,
                base_species=base_species,
                ncbi_parsed_json_path=f"{database_data_folder}parsed_taxdmp.json",
                kinetic_ignored_metabolites=kinetic_ignored_metabolites,
                kinetic_ignored_metabolite_exceptions=kinetic_ignored_metabolite_exceptions,
                kinetic_ignored_enzyme_ids=kinetic_ignored_enzymes,
                custom_enzyme_kinetic_data=custom_kms_and_kcats,
                max_taxonomy_level=max_taxonomy_level,
                transfered_ec_number_json=transfer_json_path,
                kis_and_kas_only_for_same_compartments=kis_and_kas_only_for_same_compartments,
            )
            sabio_enzyme_reaction_data = sabio_select_enzyme_kinetic_data_for_sbml(
                sbml_path=temp_sbml_path,
                sabio_target_folder=database_data_folder,
                base_species=base_species,
                ncbi_parsed_json_path=f"{database_data_folder}parsed_taxdmp.json",
                bigg_metabolites_json_path=f"{database_data_folder}bigg_models_metabolites.json",
                kinetic_ignored_metabolites=kinetic_ignored_metabolites,
                kinetic_ignored_metabolite_exceptions=kinetic_ignored_metabolite_exceptions,
                kinetic_ignored_enzyme_ids=kinetic_ignored_enzymes,
                custom_enzyme_kinetic_data=custom_kms_and_kcats,
                max_taxonomy_level=max_taxonomy_level,
                add_hill_coefficients=add_hill_coefficients,
                transfered_ec_number_json=transfer_json_path,
                kis_and_kas_only_for_same_compartments=kis_and_kas_only_for_same_compartments,
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

        if database_data_folder:
            json_write(
                f"{database_data_folder}_cache_enzyme_reaction_data.json",
                enzyme_reaction_data,
            )
    else:
        enzyme_reaction_data = json_load(
            f"{database_data_folder}_cache_enzyme_reaction_data.json",
            dict[str, EnzymeReactionData | None],
        )

    with tempfile.TemporaryDirectory() as tmpdict:
        sbml_path = tmpdict + "temp.xml"
        cobra.io.write_sbml_model(fullsplit_model, sbml_path)
        enzyme_molecular_weights = uniprot_get_enzyme_molecular_weights_for_sbml(
            sbml_path=sbml_path,
            cache_basepath=database_data_folder,
            base_species=base_species,
        )
        if database_data_folder:
            json_write(
                f"{database_data_folder}_cache_uniprot_molecular_weights.json",
                enzyme_molecular_weights,
            )

        if add_protein_sequences:
            sequences = uniprot_get_enzyme_sequences_for_sbml(
                sbml_path=sbml_path,
                cache_basepath=database_data_folder,
                base_species=base_species,
            )
            if database_data_folder:
                json_write(
                    f"{database_data_folder}_cache_uniprot_sequences.json",
                    sequences,
                )
        else:
            sequences = {}

    if do_delete_enzymatically_suboptimal_reactions:
        fullsplit_model = (
            delete_enzymatically_suboptimal_reactions_in_fullsplit_cobrapy_model(
                fullsplit_model,
                enzyme_reaction_data,
                enzyme_molecular_weights,
            )
        )

    if (not database_data_folder) or (
        (database_data_folder)
        and (
            ("_cache_dG0.json" not in data_cache_files)
            or ("_cache_dG0_uncertainties.json" not in data_cache_files)
        )
    ):
        with tempfile.TemporaryDirectory() as tmpdict:
            cobra.io.write_sbml_model(fullsplit_model, tmpdict + "temp.xml")
            dG0s, dG0_uncertainties = (
                equilibrator_get_model_dG0_and_uncertainty_values_for_sbml(
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
            )
        if database_data_folder:
            json_write(f"{database_data_folder}_cache_dG0.json", dG0s)
            json_write(
                f"{database_data_folder}_cache_dG0_uncertainties.json",
                dG0_uncertainties,
            )
    else:
        dG0s = json_load(f"{database_data_folder}_cache_dG0.json", dict[str, float])
        dG0_uncertainties = json_load(
            f"{database_data_folder}_cache_dG0_uncertainties.json",
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
        return delete_orphaned_metabolites_and_enzymes(
            get_cobrak_model_from_sbml_and_thermokinetic_data(
                sbml_path=tmpdict + "temp.xml",
                extra_linear_constraints=extra_linear_constraints,
                dG0s=dG0s,
                dG0_uncertainties=dG0_uncertainties
                if not ignore_dG0_uncertainty
                else {},
                conc_ranges=conc_ranges,
                enzyme_molecular_weights=enzyme_molecular_weights,
                enzyme_reaction_data=enzyme_reaction_data,
                max_prot_pool=max_prot_pool,
                kinetic_ignored_metabolites=kinetic_ignored_metabolites,
                kinetic_ignored_metabolite_exceptions=kinetic_ignored_metabolite_exceptions,
                enzyme_conc_ranges=enzyme_conc_ranges,
                R=R,
                T=T,
                do_delete_enzymatically_suboptimal_reactions=False,
                sequences=sequences,
                add_molar_masses=add_molar_masses,
            )
        )
