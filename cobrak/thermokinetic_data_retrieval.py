"""Functions for directly retrieving thermokinetic data for and into a COBRA-k Model instance."""

from os.path import exists
from tempfile import TemporaryDirectory

from pydantic import NonNegativeInt, validate_call

from .brenda_functionality import brenda_select_enzyme_kinetic_data_for_sbml
from .constants import (
    EC_INNER_TO_OUTER_COMPARTMENTS,
    EC_IONIC_STRENGTHS,
    EC_PHS,
    EC_PMGS,
    EC_POTENTIAL_DIFFERENCES,
)
from .dataclasses import EnzymeReactionData, Model
from .equilibrator_functionality import (
    equilibrator_get_model_dG0_and_uncertainty_values_for_sbml,
)
from .expasy_functionality import get_ec_number_transfers
from .io import (
    json_write,
    save_cobrak_model_as_annotated_sbml_model,
    standardize_folder,
)
from .sabio_rk_functionality import sabio_select_enzyme_kinetic_data_for_sbml
from .uniprot_functionality import uniprot_get_enzyme_molecular_weights_for_sbml
from .utilities import combine_enzyme_reaction_datasets, parse_external_resources


@validate_call(validate_return=True)
def add_thermokinetic_data_to_cobrak_model(
    cobrak_model: Model,
    mws: dict[str, float] = {},
    kcats: dict[str, float] = {},
    kms: dict[str, dict[str, float]] = {},
    kis: dict[str, dict[str, float]] = {},
    kas: dict[str, dict[str, float]] = {},
    dG0s: dict[str, float] = {},
    conc_ranges: dict[str, tuple[float, float]] = {},
    delete_old_dG0s: bool = False,
    overwrite_existing_dG0s: bool = True,
    overwrite_existing_enzyme_reaction_data: bool = True,
) -> Model:
    """Populate a COBRA‑K model with thermodynamic and kinetic parameters.

    Parameters
    ----------
    cobrak_model : Model
        The model to be updated.
    mws : dict[str, float], optional
        Molecular weights for enzymes keyed by enzyme ID.
    kcats : dict[str, float], optional
        kcat values keyed by reaction ID.
    kms : dict[str, dict[str, float]], optional
        Michaelis‑Menten constants keyed by reaction ID then metabolite ID.
    kis : dict[str, dict[str, float]], optional
        Inhibition constants keyed by reaction ID then metabolite ID.
    kas : dict[str, dict[str, float]], optional
        Activation constants keyed by reaction ID then metabolite ID.
    dG0s : dict[str, float], optional
        Standard Gibbs free energies keyed by reaction ID.
    conc_ranges : dict[str, tuple[float, float]], optional
        Log‑concentration bounds for metabolites keyed by metabolite ID.
    delete_old_dG0s : bool, default False
        If True, remove any existing dG0 values before adding new ones.
    overwrite_existing_dG0s : bool, default True
        Overwrite existing dG0 values when a new value is supplied.
    overwrite_existing_enzyme_reaction_data : bool, default True
        Overwrite existing enzyme reaction data when new data is supplied.

    Returns
    -------
    Model
        The updated model instance.
    """
    # Molecular weights
    for enzyme_id, enzyme in cobrak_model.enzymes.items():
        if enzyme_id in mws:
            enzyme.molecular_weight = mws[enzyme_id]

    # dG0s, kcats, kms, kis and kas
    for reac_id, reaction in cobrak_model.reactions.items():
        # dG0s
        if delete_old_dG0s:
            reaction.dG0 = None
        elif (overwrite_existing_dG0s or reaction.dG0 is None) and (reac_id in dG0s):
            reaction.dG0 = dG0s[reac_id]

        # enzyme_reaction_data
        if delete_old_dG0s:
            reaction.enzyme_reaction_data = None
        elif (
            overwrite_existing_enzyme_reaction_data
            or reaction.enzyme_reaction_data is None
        ):
            if reaction.enzyme_reaction_data is None:
                continue
            if reac_id in kcats:
                reaction.enzyme_reaction_data.k_cat = kcats[reac_id]
            if reac_id in kms:
                for met_id, value in kms[reac_id].items():
                    reaction.enzyme_reaction_data.kms[met_id] = value
            if reac_id in kis:
                for met_id, value in kis[reac_id].items():
                    reaction.enzyme_reaction_data.kis[met_id] = value
            if reac_id in kas:
                for met_id, value in kas[reac_id].items():
                    reaction.enzyme_reaction_data.kas[met_id] = value

    # concentration ranges
    for met_id, (min_log_conc, max_log_conc) in conc_ranges.items():
        if met_id not in cobrak_model.metabolites:
            continue
        cobrak_model.metabolites[met_id].min_log_conc = min_log_conc
        cobrak_model.metabolites[met_id].max_log_conc = max_log_conc

    return cobrak_model


@validate_call(validate_return=True)
def add_enzyme_reaction_data_to_cobrak_model(
    cobrak_model: Model,
    enzyme_reaction_data: dict[str, EnzymeReactionData],
    delete_old_enzyme_reaction_data: bool = False,
    overwrite_existing_enzyme_reaction_data: bool = True,
) -> Model:
    """Insert pre‑computed :class:`EnzymeReactionData` objects into a model

    Model
        The model to be updated.
    enzyme_reaction_data : dict[str, EnzymeReactionData]
        Mapping from reaction IDs to enzyme reaction data objects.
    delete_old_enzyme_reaction_data : bool, default False
        If True, remove any existing data before inserting new data.
    overwrite_existing_enzyme_reaction_data : bool, default True
        Overwrite existing data when a matching reaction ID is found.

    Returns
    -------
    Model
        The model with updated enzyme reaction data.
    """
    for reac_id, reaction in cobrak_model.reactions.items():
        if (reaction.enzyme_reaction_data is not None) and (
            not overwrite_existing_enzyme_reaction_data
        ):
            continue
        if (reaction.enzyme_reaction_data is not None) and (
            delete_old_enzyme_reaction_data
        ):
            reaction.enzyme_reaction_data = None
        if reac_id in enzyme_reaction_data:
            reaction.enzyme_reaction_data = enzyme_reaction_data[reac_id]
    return cobrak_model


@validate_call(validate_return=True)
def automatically_add_database_thermokinetic_data_too_cobrak_model(
    cobrak_model: Model,
    database_data_path: str,
    use_brenda: bool = True,
    use_sabio_rk: bool = True,
    prefer_brenda: bool = True,
    use_ec_number_transfers: bool = True,
    max_taxonomy_level: int = 1_000,
    kinetic_ignored_enzyme_ids: list[str] = ["s0001"],
) -> Model:
    """Retrieve kinetic and thermodynamic data from external databases and add them to a model.

    Parameters
    ----------
    cobrak_model : Model
        The model to be enriched.
    database_data_path : str
        Path to the folder containing the required database files.
    use_brenda : bool, default True
        Include BRENDA data if True.
    use_sabio_rk : bool, default True
        Include SABIO‑RK data if True.
    prefer_brenda : bool, default True
        When both databases provide data, give precedence to BRENDA.
    use_ec_number_transfers : bool, default True
        Use EC‑number transfer mappings when searching for data.
    max_taxonomy_level : int, default 1000
        Maximum taxonomic distance allowed for data transfer.
    kinetic_ignored_enzyme_ids : list[str], default ["s0001"]
        Enzyme IDs to be ignored during kinetic data retrieval.

    Returns
    -------
    Model
        The model populated with database‑derived thermokinetic data.
    """
    database_data_path = standardize_folder(database_data_path)

    enzyme_reaction_data = get_database_kcats_kms_kis_and_kas_for_cobrak_model(
        cobrak_model=cobrak_model,
        database_data_path=database_data_path,
        use_brenda=use_brenda,
        use_sabio_rk=use_sabio_rk,
        prefer_brenda=prefer_brenda,
        use_ec_number_transfers=use_ec_number_transfers,
        max_taxonomy_level=max_taxonomy_level,
        kinetic_ignored_enzyme_ids=kinetic_ignored_enzyme_ids,
    )
    cobrak_model = add_enzyme_reaction_data_to_cobrak_model(
        cobrak_model=cobrak_model,
        enzyme_reaction_data=enzyme_reaction_data,
    )
    mws = get_database_mws_for_cobrak_model(
        cobrak_model=cobrak_model,
        database_data_path=database_data_path,
    )
    return add_thermokinetic_data_to_cobrak_model(
        cobrak_model=cobrak_model,
        mws=mws,
    )


@validate_call(validate_return=True)
def get_database_kcats_kms_kis_and_kas_for_cobrak_model(
    cobrak_model: Model,
    database_data_path: str,
    use_brenda: bool = True,
    use_sabio_rk: bool = True,
    prefer_brenda: bool = False,
    use_ec_number_transfers: bool = True,
    max_taxonomy_level: NonNegativeInt = 1_000,
    kinetic_ignored_enzyme_ids: list[str] = ["s0001"],
) -> dict[str, EnzymeReactionData]:
    """Query BRENDA and/or SABIO‑RK for kinetic parameters and return (if given) a unified dataset.

    Parameters
    ----------
    cobrak_model : Model
        The model for which kinetic data are required.
    database_data_path : str
        Directory containing the database files.
    use_brenda : bool, default True
        Retrieve data from BRENDA if True.
    use_sabio_rk : bool, default True
        Retrieve data from SABIO‑RK if True.
    prefer_brenda : bool, default False
        When both sources contain data for a reaction, keep BRENDA's values.
    use_ec_number_transfers : bool, default True
        Apply EC‑number transfer mappings when searching for data.
    max_taxonomy_level : NonNegativeInt, default 1000
        Maximum allowed taxonomic distance for data transfer.
    kinetic_ignored_enzyme_ids : list[str], default ["s0001"]
        Enzyme IDs to be excluded from kinetic data retrieval.

    Returns
    -------
    dict[str, EnzymeReactionData]
        Mapping from reaction IDs to populated :class:`EnzymeReactionData` objects.
    """
    database_data_path = standardize_folder(database_data_path)
    if not use_brenda and not use_sabio_rk:
        print(
            "ERROR: Arguments use_brenda and use_sabio_rk are both False, but at least one of the databases has to be used"
        )
        raise ValueError

    if use_ec_number_transfers:
        transfer_json_path = f"{database_data_path}ec_number_transfers.json"
        if not exists(transfer_json_path):
            if not exists(f"{database_data_path}enzyme.rdf"):
                print(
                    f"ERROR: Argument use_ec_number_transfers is True, but no necessary enzyme.rdf can be found in {database_data_path}"
                )
                print(
                    "You may download it from https://ftp.expasy.org/databases/enzyme/"
                )
                print(f"After downloading, put it into the folder {database_data_path}")
            ec_number_transfers = get_ec_number_transfers(
                f"{database_data_path}enzyme.rdf"
            )
            json_write(transfer_json_path, ec_number_transfers)
    else:
        transfer_json_path = ""

    parse_external_resources(
        path=database_data_path,
        brenda_version="2024_1",
        parse_brenda=use_brenda,
    )

    with TemporaryDirectory() as tmpdict:
        sbml_path = tmpdict + "temp.xml"
        save_cobrak_model_as_annotated_sbml_model(
            cobrak_model=cobrak_model,
            filepath=sbml_path,
        )

        brenda_enzyme_reaction_data = brenda_select_enzyme_kinetic_data_for_sbml(
            sbml_path=sbml_path,
            brenda_json_targz_file_path=f"{database_data_path}brenda_2024_1.json.tar.gz",
            bigg_metabolites_json_path=f"{database_data_path}bigg_models_metabolites.json",
            brenda_version="2024_1",
            base_species="Escherichia coli",
            ncbi_parsed_json_path=f"{database_data_path}parsed_taxdmp.json",
            kinetic_ignored_metabolites=cobrak_model.kinetic_ignored_metabolites,
            kinetic_ignored_enzyme_ids=kinetic_ignored_enzyme_ids,
            transfered_ec_number_json=transfer_json_path,
            max_taxonomy_level=max_taxonomy_level,
        )

        sabio_enzyme_reaction_data = sabio_select_enzyme_kinetic_data_for_sbml(
            sbml_path=sbml_path,
            sabio_target_folder="examples/common_needed_external_resources",
            bigg_metabolites_json_path=f"{database_data_path}bigg_models_metabolites.json",
            base_species="Escherichia coli",
            ncbi_parsed_json_path=f"{database_data_path}parsed_taxdmp.json",
            kinetic_ignored_metabolites=cobrak_model.kinetic_ignored_metabolites,
            kinetic_ignored_enzyme_ids=kinetic_ignored_enzyme_ids,
            transfered_ec_number_json=transfer_json_path,
            max_taxonomy_level=max_taxonomy_level,
            add_hill_coefficients=False,
        )

    if use_brenda and use_sabio_rk:
        return combine_enzyme_reaction_datasets(
            [brenda_enzyme_reaction_data, sabio_enzyme_reaction_data]
            if prefer_brenda
            else [sabio_enzyme_reaction_data, brenda_enzyme_reaction_data],
        )
    if use_brenda:
        return sabio_enzyme_reaction_data
    return brenda_enzyme_reaction_data


@validate_call(validate_return=True)
def get_database_mws_for_cobrak_model(
    cobrak_model: Model,
    database_data_path: str = "",
) -> dict[str, float]:
    """Retrieve enzyme molecular weights from UniProt for a given model.

    Parameters
    ----------
    cobrak_model : Model
        The model whose enzymes require molecular weights.
    database_data_path : str, optional
        Base path for caching UniProt queries (default empty string).

    Returns
    -------
    dict[str, float]
        Mapping from enzyme IDs to molecular weight values (Daltons).
    """
    database_data_path = standardize_folder(database_data_path)
    with TemporaryDirectory() as tmpdict:
        sbml_path = tmpdict + "temp.xml"
        save_cobrak_model_as_annotated_sbml_model(
            cobrak_model=cobrak_model,
            filepath=sbml_path,
        )
        return uniprot_get_enzyme_molecular_weights_for_sbml(
            sbml_path=sbml_path,
            cache_basepath=database_data_path,
        )


@validate_call(validate_return=True)
def get_database_dG0s_for_cobrak_model(
    cobrak_model: Model,
    inner_to_outer_compartments: list[str] = EC_INNER_TO_OUTER_COMPARTMENTS,
    phs: dict[str, float] = EC_PHS,
    pmgs: dict[str, float] = EC_PMGS,
    ionic_strenghts: dict[str, float] = EC_IONIC_STRENGTHS,
    potential_differences: dict[str, float] = EC_POTENTIAL_DIFFERENCES,
    calculate_multicompartmental: bool = True,
    exclusion_prefixes: list[str] = [],
    exclusion_inner_parts: list[str] = [],
    ignore_uncertainty: bool = False,
    max_uncertainty: float = 1_000.0,
) -> tuple[dict[str, float], dict[str, float]]:
    """Compute standard Gibbs free energies (and uncertainties) for all reactions in a model.

    Parameters
    ----------
    cobrak_model : Model
        The model for which dG⁰ values are to be calculated.
    inner_to_outer_compartments : list[str], optional
        Mapping of inner to outer compartments for multi‑compartment calculations.
    phs : dict[str, float], optional
        pH values per compartment.
    pmgs : dict[str, float], optional
        Proton‑motives per compartment.
    ionic_strenghts : dict[str, float], optional
        Ionic strength per compartment.
    potential_differences : dict[str, float], optional
        Electrical potential differences per compartment.
    calculate_multicompartmental : bool, default True
        Whether to compute dG⁰ for reactions spanning multiple compartments.
    exclusion_prefixes : list[str], optional
        Reaction ID prefixes to exclude from calculation.
    exclusion_inner_parts : list[str], optional
        Inner compartment identifiers to exclude.
    ignore_uncertainty : bool, default False
        If True, uncertainties are not calculated.
    max_uncertainty : float, default 1000.0
        Upper bound for acceptable uncertainty; reactions exceeding this are omitted.

    Returns
    -------
    tuple[dict[str, float], dict[str, float]]
        Two dictionaries mapping reaction IDs to dG⁰ values and to uncertainties,
        respectively.
    """
    with TemporaryDirectory() as tmpdict:
        sbml_path = tmpdict + "temp.xml"
        save_cobrak_model_as_annotated_sbml_model(
            cobrak_model=cobrak_model,
            filepath=sbml_path,
        )
        return equilibrator_get_model_dG0_and_uncertainty_values_for_sbml(
            sbml_path=sbml_path,
            inner_to_outer_compartments=inner_to_outer_compartments,
            phs=phs,
            pmgs=pmgs,
            ionic_strengths=ionic_strenghts,
            potential_differences=potential_differences,
            exclusion_prefixes=exclusion_prefixes,
            exclusion_inner_parts=exclusion_inner_parts,
            ignore_uncertainty=ignore_uncertainty,
            max_uncertainty=max_uncertainty,
            calculate_multicompartmental=calculate_multicompartmental,
        )
