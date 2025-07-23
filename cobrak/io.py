"""General (COBRAk-independent) helper functions, primarily for I/O tasks such as pickle and JSON file handlings."""

# IMPORTS SECTION #
import contextlib
import json
import os
import pickle
import tempfile
import zipfile
from ast import literal_eval
from dataclasses import asdict, is_dataclass
from typing import Any, TypeVar
from zipfile import ZipFile

import cobra
from numpy import exp, log
from pydantic import BaseModel, TypeAdapter

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

# "PRIVATE" FUNCTIONS SECTION #
T = TypeVar("T")


def _add_annotation_to_cobra_reaction(
    cobra_reaction: cobra.Reaction,
    reac_id: str,
    reac_data: Reaction,
    version: str,
) -> None:
    """Adds annotations from a COBRAk Reaction object to a COBRApy Reaction object.

    This function updates the annotation dictionary of the COBRA Reaction object
    with data from the Reaction object, including thermodynamic and kinetic data.
    It also updates the gene_reaction_rule attribute of the COBRA Reaction object
    if enzyme reaction data is available.

    Parameters
    ----------
    cobra_reaction : cobra.Reaction
        The COBRA Reaction object to be updated.
    reac_id : str
        The ID of the reaction.
    reac_data : Reaction
        The Reaction object containing the data to be added.
    version : str
        The version of the data being added.

    Returns
    -------
    None (the cobra.Reaction is permanently changed as it given as reference)
    """
    cobra_reaction.annotation |= reac_data.annotation
    cobra_reaction.annotation[f"cobrak_id_{version}"] = reac_id

    if reac_data.dG0 is not None:
        cobra_reaction.annotation[f"cobrak_dG0_{version}"] = str(reac_data.dG0)
    if reac_data.dG0_uncertainty is not None:
        cobra_reaction.annotation[f"cobrak_dG0_uncertainty_{version}"] = str(
            reac_data.dG0_uncertainty
        )
    if reac_data.enzyme_reaction_data is not None:
        cobra_reaction.annotation[f"cobrak_k_cat_{version}"] = str(
            reac_data.enzyme_reaction_data.k_cat
        )
        cobra_reaction.gene_reaction_rule = " and ".join(
            reac_data.enzyme_reaction_data.identifiers
        )
        if reac_data.enzyme_reaction_data.k_ms is not None:
            cobra_reaction.annotation[f"cobrak_k_ms_{version}"] = str(
                reac_data.enzyme_reaction_data.k_ms
            )
        if reac_data.enzyme_reaction_data.k_is is not None:
            cobra_reaction.annotation[f"cobrak_k_is_{version}"] = str(
                reac_data.enzyme_reaction_data.k_is
            )
        if reac_data.enzyme_reaction_data.k_as is not None:
            cobra_reaction.annotation[f"cobrak_k_as_{version}"] = str(
                reac_data.enzyme_reaction_data.k_as
            )
        cobra_reaction.annotation[f"cobrak_special_stoichiometries_{version}"] = str(
            reac_data.enzyme_reaction_data.special_stoichiometries
        )


# "PUBLIC" FUNCTIONS SECTION #
def get_base_id(
    reac_id: str,
    fwd_suffix: str = REAC_FWD_SUFFIX,
    rev_suffix: str = REAC_REV_SUFFIX,
    reac_enz_separator: str = REAC_ENZ_SEPARATOR,
) -> str:
    """Extract the base ID from a reaction ID by removing specified suffixes and separators.

    Processes a reaction ID to remove forward and reverse suffixes
    as well as any enzyme separators, to obtain the base reaction ID.

    Args:
        reac_id (str): The reaction ID to be processed.
        fwd_suffix (str, optional): The suffix indicating forward reactions. Defaults to REAC_FWD_SUFFIX.
        rev_suffix (str, optional): The suffix indicating reverse reactions. Defaults to REAC_REV_SUFFIX.
        reac_enz_separator (str, optional): The separator used between reaction and enzyme identifiers. Defaults to REAC_ENZ_SEPARATOR.

    Returns:
        str: The base reaction ID with specified suffixes and separators removed.
    """
    reac_id_split = reac_id.split(reac_enz_separator)
    return (
        (reac_id_split[0] + "\b")
        .replace(f"{fwd_suffix}\b", "")
        .replace(f"{rev_suffix}\b", "")
        .replace("\b", "")
    )


# FUNCTIONS SECTION #
def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    with contextlib.suppress(FileExistsError):
        os.makedirs(folder)


def ensure_json_existence(path: str) -> None:
    """Ensures that a JSON file exists at the specified path.

    If the file does not exist, it creates an empty JSON file with "{}" as its content.

    Args:
        path (str): The file path where the JSON file should exist.
    """
    if os.path.isfile(path):
        return
    with open(path, "w", encoding="utf-8") as f:
        f.write("{}")


def convert_cobrak_model_to_annotated_cobrapy_model(
    cobrak_model: Model,
    combine_base_reactions: bool = False,
    add_enzyme_constraints: bool = False,
) -> cobra.Model:
    """Converts a COBRAk model to an annotated COBRApy model.

    This function takes a COBRAk model and converts it to a COBRApy model,
    adding annotations and constraints as specified by the input parameters.

    The function adds the following annotation keys to the COBRApy model:

    * `cobrak_Cmin`: The minimum concentration of a metabolite.
    * `cobrak_Cmax`: The maximum concentration of a metabolite.
    * `cobrak_id_<version>`: The ID of the reaction in the COBRAk model.
    * `cobrak_dG0_<version>`: The standard Gibbs free energy change of a reaction.
    * `cobrak_dG0_uncertainty_<version>`: The uncertainty of the standard Gibbs free energy change of a reaction.
    * `cobrak_k_cat_<version>`: The turnover number of an enzyme.
    * `cobrak_k_ms_<version>`: The Michaelis constant of an enzyme.
    * `cobrak_k_is_<version>`: The inhibition constant of an enzyme.
    * `cobrak_k_as_<version>`: The activation constant of an enzyme.
    * `cobrak_special_stoichiometries_<version>`: Special stoichiometries of a reaction.
    * `cobrak_max_prot_pool`: The maximum protein pool size.
    * `cobrak_R`: The gas constant.
    * `cobrak_T`: The temperature.
    * `cobrak_kinetic_ignored_metabolites`: A list of metabolites that are ignored in kinetic simulations.
    * `cobrak_extra_linear_constraints`: A list of extra linear constraints.
    * `cobrak_mw`: The molecular weight of an enzyme.
    * `cobrak_min_conc`: The minimum concentration of an enzyme.
    * `cobrak_max_conc`: The maximum concentration of an enzyme.

    The conversion process also involves the merging of forward and reverse reactions, as well as isomeric alternatives,
    into a single reaction in the COBRApy model. When the combine_base_reactions parameter is set to True,
    the function combines these reactions into a single entity, while still preserving the unique characteristics
    of each original reaction. To achieve this, the function uses a versioning system, denoted by the <version> suffix,
    to differentiate between the annotations of the original reactions. For example, the cobrak_id_<version> annotation
    key will contain the ID of the original reaction, with <version> indicating whether it corresponds to the forward or
    reverse direction, or an isomeric alternative. This versioning system allows the model to retain the distinct properties
    of each reaction, such as their standard Gibbs free energy changes or enzyme kinetics, while still representing them as
    a single, unified reaction. The <version> suffix can take on values such as V0, V1, etc., with each value corresponding
    to a specific original reaction.

    The conversion of a COBRAk model to a COBRApy model also includes the optional direct addition of enzyme constraints
    in the style of GECKO [1] (or expaned sMOMENT [2]),
    which can be enabled through the add_enzyme_constraints parameter. When this parameter is set to True,
    the function introduces new pseudo-metabolites and pseudo-reactions to the model, allowing for the simulation
    of enzyme kinetics and protein expression. Specifically, a protein pool pseudo-metabolite is added, which
    represents the total amount of protein available in the system. Additionally, pseudo-reactions are created
    to deliver enzymes to the protein pool, taking into account the molecular weight and concentration of each enzyme.
    The function also adds pseudo-reactions to form enzyme complexes, which are essential for simulating the k_cat-based kinetics
    of enzymatic reactions.

    [1] https://doi.org/10.15252/msb.20167411
    [2] https://doi.org/10.1186/s12859-019-3329-9

    Parameters
    ----------
    cobrak_model : Model
        The COBRAk model to be converted.
    combine_base_reactions : bool, optional
        Whether to combine base reactions into a single reaction (default: False).
    add_enzyme_constraints : bool, optional
        Whether to add enzyme constraints to the model (default: False).

    Returns
    -------
    cobra.Model
        The converted COBRApy model.

    Raises
    ------
    ValueError
        If combine_base_reactions and add_enzyme_constraints are both True.
    """
    if combine_base_reactions and add_enzyme_constraints:
        print(
            "ERROR: Stoichiometric enzyme constraints do not work with combined base reactions\n"
            "       as for these enzyme constraints, reactions must remain irreversible."
        )
        raise ValueError

    cobra_model = cobra.Model()

    # Add metabolites
    added_metabolites: list[cobra.Metabolite] = []
    for met_id, met_data in cobrak_model.metabolites.items():
        cobra_metabolite: cobra.Metabolite = cobra.Metabolite(
            id=met_id,
            compartment=met_id.split("_")[-1] if "_" in met_id else "c",
            name=met_data.name,
            formula=met_data.formula,
            charge=met_data.charge,
        )
        cobra_metabolite.annotation = met_data.annotation

        # Add full annotation
        cobra_metabolite.annotation["cobrak_Cmin"] = exp(met_data.log_min_conc)
        cobra_metabolite.annotation["cobrak_Cmax"] = exp(met_data.log_max_conc)

        added_metabolites.append(cobra_metabolite)
    cobra_model.add_metabolites(added_metabolites)

    if add_enzyme_constraints:
        enzyme_reacs: list[cobra.Reaction] = []

        # If set: Add protein pool reaction (flux in g⋅gDW⁻¹)
        prot_pool_met = cobra.Metabolite(id="prot_pool", compartment="c")
        cobra_model.add_metabolites([prot_pool_met])

        prot_pool_reac = cobra.Reaction(
            "prot_pool_delivery",
            lower_bound=0.0,
            upper_bound=cobrak_model.max_prot_pool,
        )
        prot_pool_reac.add_metabolites(
            {
                prot_pool_met: 1.0,
            }
        )
        enzyme_reacs.append(prot_pool_reac)

        # If set: Add enzyme concentration delivery reactions (flux in mmol⋅gDW⁻¹)
        for enzyme_id, enzyme_data in cobrak_model.enzymes.items():
            enzyme_met = cobra.Metabolite(id=enzyme_id, compartment="c")

            lower_bound = (
                enzyme_data.min_conc if enzyme_data.min_conc is not None else 0.0
            )
            upper_bound = (
                enzyme_data.max_conc if enzyme_data.max_conc is not None else 100_000.0
            )
            enzyme_reac = cobra.Reaction(
                id="enzyme_delivery_" + enzyme_id,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
            )
            enzyme_reac.add_metabolites(
                {
                    prot_pool_met: -enzyme_data.molecular_weight,
                    enzyme_met: 1.0,
                }
            )
            enzyme_reacs.append(enzyme_reac)

        cobra_model.add_reactions(enzyme_reacs)

    # If set: Add enzyme complex metabolites and delivery reactions
    if add_enzyme_constraints:
        enzyme_complexes: set[tuple[str, ...]] = {
            tuple(
                enzyme_id
                for enzyme_id in reaction_data.enzyme_reaction_data.identifiers
            )
            for reaction_data in cobrak_model.reactions.values()
            if reaction_data.enzyme_reaction_data is not None
        }
        complex_reacs: list[cobra.Reaction] = []
        for enzyme_complex in enzyme_complexes:
            if len(enzyme_complex) <= 1:
                continue
            if enzyme_complex == ("",):
                continue
            complex_met = cobra.Metabolite(id="_".join(enzyme_complex), compartment="c")
            complex_reac = cobra.Reaction(
                id="complex_delivery_" + "_".join(enzyme_complex),
                lower_bound=0.0,
                upper_bound=10_000.0,
            )
            complex_reac.add_metabolites(
                {
                    cobra_model.metabolites.get_by_id(enzyme_id): -1
                    for enzyme_id in enzyme_complex
                    if enzyme_id
                }
            )
            complex_reac.add_metabolites(
                {
                    complex_met: 1.0,
                }
            )
            complex_reacs.append(complex_reac)
        cobra_model.add_reactions(complex_reacs)

    # Add reactions
    added_reactions: list[cobra.Reaction] = []
    if not combine_base_reactions:
        for reac_id, reac_data in cobrak_model.reactions.items():
            cobra_reaction = cobra.Reaction(
                id=reac_id,
                lower_bound=reac_data.min_flux,
                upper_bound=reac_data.max_flux,
                name=reac_data.name,
            )
            cobra_reaction.add_metabolites(
                {
                    cobra_model.metabolites.get_by_id(met_id): stoich
                    for met_id, stoich in reac_data.stoichiometries.items()
                }
            )

            # Add full annotation
            _add_annotation_to_cobra_reaction(cobra_reaction, reac_id, reac_data, "V0")

            if (
                reac_data.enzyme_reaction_data is not None
                and add_enzyme_constraints
                and reac_data.enzyme_reaction_data.identifiers != []
            ):
                complex_met_id = "_".join(reac_data.enzyme_reaction_data.identifiers)
                if complex_met_id:
                    cobra_reaction.add_metabolites(
                        {
                            cobra_model.metabolites.get_by_id(complex_met_id): -1
                            / reac_data.enzyme_reaction_data.k_cat
                        }
                    )
            added_reactions.append(cobra_reaction)
    else:
        base_id_to_reac_ids: dict[str, list[str]] = {}
        for reac_id in cobrak_model.reactions:
            base_id = get_base_id(
                reac_id,
                cobrak_model.fwd_suffix,
                cobrak_model.rev_suffix,
                cobrak_model.reac_enz_separator,
            )
            if base_id not in base_id_to_reac_ids:
                base_id_to_reac_ids[base_id] = []
            base_id_to_reac_ids[base_id].append(reac_id)

        for base_id, reac_ids in base_id_to_reac_ids.items():
            rev_ids = [
                reac_id
                for reac_id in reac_ids
                if reac_id.endswith(cobrak_model.rev_suffix)
            ]
            fwd_ids = [
                reac_id
                for reac_id in reac_ids
                if not reac_id.endswith(cobrak_model.rev_suffix)
            ]

            if len(rev_ids) > 0:
                min_flux = -max(
                    cobrak_model.reactions[rev_id].max_flux for rev_id in rev_ids
                )
                name = cobrak_model.reactions[rev_ids[0]].name
            else:
                min_flux = max(
                    cobrak_model.reactions[fwd_id].min_flux for fwd_id in fwd_ids
                )
            if len(fwd_ids) > 0:
                max_flux = max(
                    cobrak_model.reactions[fwd_id].max_flux for fwd_id in fwd_ids
                )
                met_stoichiometries = {
                    cobra_model.metabolites.get_by_id(met_id): stoich
                    for met_id, stoich in cobrak_model.reactions[
                        fwd_ids[0]
                    ].stoichiometries.items()
                }
                name = cobrak_model.reactions[fwd_ids[0]].name
            else:
                max_flux = min(
                    cobrak_model.reactions[rev_id].max_flux for rev_id in rev_ids
                )
                met_stoichiometries = {
                    cobra_model.metabolites.get_by_id(met_id): -stoich
                    for met_id, stoich in cobrak_model.reactions[
                        rev_ids[0]
                    ].stoichiometries.items()
                }

            cobra_reaction = cobra.Reaction(
                id=base_id,
                lower_bound=min_flux,
                upper_bound=max_flux,
                name=name,
            )
            cobra_reaction.add_metabolites(met_stoichiometries)
            for number, reac_id in enumerate(reac_ids):
                version = f"V{number}"
                reac_data = cobrak_model.reactions[reac_id]
                _add_annotation_to_cobra_reaction(
                    cobra_reaction, reac_id, reac_data, version
                )

            added_reactions.append(cobra_reaction)

    # Add global information reaction
    added_reactions.append(
        cobra.Reaction(
            id="cobrak_global_settings",
            lower_bound=0.0,
            upper_bound=0.0,
        )
    )
    added_reactions[-1].annotation["cobrak_max_prot_pool"] = cobrak_model.max_prot_pool
    added_reactions[-1].annotation["cobrak_R"] = cobrak_model.R
    added_reactions[-1].annotation["cobrak_T"] = cobrak_model.T
    added_reactions[-1].annotation["cobrak_kinetic_ignored_metabolites"] = str(
        cobrak_model.kinetic_ignored_metabolites
    )
    added_reactions[-1].annotation["cobrak_reac_rev_suffix"] = cobrak_model.rev_suffix
    added_reactions[-1].annotation["cobrak_reac_fwd_suffix"] = cobrak_model.fwd_suffix
    added_reactions[-1].annotation["cobrak_reac_enz_separator"] = (
        cobrak_model.reac_enz_separator
    )
    added_reactions[-1].annotation["cobrak_extra_linear_constraints"] = str(
        [asdict(x) for x in cobrak_model.extra_linear_constraints]
    )

    cobra_model.add_reactions(added_reactions)

    gene_ids = [x.id for x in cobra_model.genes]
    for enzyme_id, enzyme_data in cobrak_model.enzymes.items():
        if enzyme_id not in gene_ids:
            cobra_model.genes.append(cobra.Gene(enzyme_id, name=enzyme_data.name))
        gene = cobra_model.genes.get_by_id(enzyme_id)
        gene.annotation["cobrak_mw"] = enzyme_data.molecular_weight
        if enzyme_data.min_conc is not None:
            gene.annotation["cobrak_min_conc"] = enzyme_data.min_conc
        if enzyme_data.max_conc is not None:
            gene.annotation["cobrak_max_conc"] = enzyme_data.max_conc

    return cobra_model


def save_cobrak_model_as_annotated_sbml_model(
    cobrak_model: Model,
    filepath: str,
    combine_base_reactions: bool = False,
    add_enzyme_constraints: bool = False,
) -> None:
    """Exports a COBRAk model to an annotated SBML file.

    This function converts a `Model` to a COBRApy model and writes it to an SBML file at the specified file path.
    Optionally, stoichiometric GECKO [1]-like enzyme constraints can be added during the conversion.

    [1] Sánchez et al. Molecular systems biology, 13(8), 935. https://doi.org/10.15252/msb.20167411

    Args:
        cobrak_model (Model): The `Model` to be exported.
        filepath (str): The file path where the SBML file will be saved.
        add_enzyme_constraints (bool, optional): Whether to add enzyme constraints during the conversion. Defaults to False.
    """
    cobra.io.write_sbml_model(
        convert_cobrak_model_to_annotated_cobrapy_model(
            cobrak_model,
            combine_base_reactions,
            add_enzyme_constraints,
        ),
        filepath,
    )


def get_files(path: str) -> list[str]:
    """Returns the names of the files in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder of which the file names shall be returned
    """
    files: list[str] = []
    for _, _, filenames in os.walk(path):
        files.extend(filenames)
    return files


def get_folders(path: str) -> list[str]:
    """Returns the names of the folders in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder whose folders shall be returned
    """
    return [
        folder
        for folder in os.listdir(path)
        if os.path.isdir(os.path.join(path, folder))
    ]


def json_load(path: str, dataclass_type: T = Any) -> T:
    """Load JSON data from a file and validate it against a specified dataclass type.

    This function reads the content of a JSON file located at the given `path`, parses it,
    and validates the parsed data against the provided `dataclass_type`. If the data is valid
    according to the dataclass schema, it returns an instance of the dataclass populated with
    the data. Otherwise, it raises an exception.

    Parameters:
    ----------
    path : str
        The file path to the JSON file that needs to be loaded.

    dataclass_type : Type[T]
        A dataclass type against which the JSON data should be validated and deserialized.

    Returns:
    -------
    T
        An instance of the specified `dataclass_type` populated with the data from the JSON file.

    Raises:
    ------
    JSONDecodeError
        If the content of the file is not a valid JSON string.

    ValidationError
        If the parsed JSON data does not conform to the schema defined by `dataclass_type`.

    Examples:
    --------
    >>> @dataclass
    ... class Person:
    ...     name: str
    ...     age: int

    >>> person = json_load('person.json', Person)
    >>> print(person.name, person.age)
    John Doe 30
    """
    with open(path, encoding="utf-8") as f:
        data = f.read()

    return TypeAdapter(dataclass_type).validate_json(data)


def json_write(path: str, json_data: Any) -> None:  # noqa: ANN401
    """Writes a JSON file at the given path with the given data as content.

    Can be also used for any of COBRAk's dataclasses as well as any
    dictionary of the form dict[str, dict[str, T] | None] where
    T stands for a COBRAk dataclass or any other JSON-compatible
    object type.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * json_data: Any ~ The dictionary or list which shalll be the content of
      the created JSON file
    """
    if is_dataclass(json_data):
        json_write(path, asdict(json_data))
    elif isinstance(json_data, BaseModel):
        json_output = json_data.model_dump_json(indent=2)
        with open(path, "w+", encoding="utf-8") as f:
            f.write(json_output)
    elif isinstance(json_data, dict) and sum(
        is_dataclass(value) for value in json_data.values()
    ):
        json_dict: dict[str, dict[str, Any] | None] = {}
        for key, data in json_data.items():
            if data is None:
                json_dict[key] = None
            elif is_dataclass(data):
                json_dict[key] = asdict(data)
            else:
                json_dict[key] = data
        json_write(path, json_dict)
    else:
        json_output = json.dumps(json_data, indent=4)
        with open(path, "w+", encoding="utf-8") as f:
            f.write(json_output)


def json_zip_load(path: str) -> dict:
    """Loads the given zipped JSON file and returns it as json_data (a list
    or a dictionary).

    Arguments
    ----------
    * path: str ~ The path of the JSON file without ".zip" at the end

    Returns
    -------
    dict or list ~ The loaded JSON data
    """
    # Create a temporary directory to extract the zip file contents
    with tempfile.TemporaryDirectory() as temp_dir:
        # Construct the full path to the zip file
        zip_path = f"{path}.zip"

        # Open the zip file and extract its contents to the temporary directory
        with zipfile.ZipFile(zip_path, "r") as zip_file:
            zip_file.extractall(temp_dir)

        # Construct the full path to the JSON file in the temporary directory
        json_path = os.path.join(temp_dir, os.path.basename(path))

        # Open and load the JSON file
        with open(json_path, encoding="utf-8") as json_file:
            json_data = json.load(json_file)

    return json_data


def json_zip_write(
    path: str,
    json_data: Any,
    zip_method: int = zipfile.ZIP_LZMA,  # noqa: ANN401
) -> None:
    """Writes a zipped JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written without ".zip" at the end
    * json_data: Any ~ The dictionary or list which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(json_data, indent=4).encode("utf-8")
    with ZipFile(path + ".zip", "w", compression=zip_method) as zip_file:
        zip_file.writestr(os.path.basename(path), json_output)


def load_annotated_cobrapy_model_as_cobrak_model(
    cobra_model: cobra.Model,
    exclude_enzyme_constraints: bool = True,
) -> Model:
    """Converts a COBRApy model with annotations into a COBRAk Model.

    This function takes a COBRApy model, which may contain specific annotations for metabolites,
    reactions, and genes, and converts it into a COBRAk model. The conversion involves extracting
    relevant annotations and constructing COBRAk-specific data structures for metabolites, reactions,
    and enzymes.

    Parameters:
    - cobra_model (cobra.Model): The COBRApy model to be converted. This model should contain
      annotations that are compatible with the COBRAk model structure.
    - exclude_enzyme_constraints (bool): Whether or not to exclude all stoichiometric enzyme constraint additions.
      Defaults to True.

    Returns:
    - Model: A COBRAk model constructed from the annotated COBRApy model, including metabolites,
      reactions, and enzymes with their respective parameters and constraints.

    Notes:
    - The function assumes that certain annotations (e.g., "cobrak_Cmin", "cobrak_dG0") are present
      in the COBRApy model. Missing annotations will result in default values being used.
    - Reactions with IDs like "prot_pool_delivery" and those starting with "enzyme_delivery_" are ignored.
    - Ensure that the COBRApy model is correctly annotated to fully leverage the conversion process.
    """
    if exclude_enzyme_constraints:
        gene_ids = [gene.id for gene in cobra_model.genes]

    if "cobrak_global_settings" in [x.id for x in cobra_model.reactions]:
        global_settings_reac = cobra_model.reactions.get_by_id("cobrak_global_settings")
        max_prot_pool = float(global_settings_reac.annotation["cobrak_max_prot_pool"])
        kinetic_ignored_metabolites = literal_eval(
            global_settings_reac.annotation["cobrak_kinetic_ignored_metabolites"]
        )
        extra_linear_constraints = [
            ExtraLinearConstraint(**x)
            for x in literal_eval(
                global_settings_reac.annotation["cobrak_extra_linear_constraints"]
            )
        ]
        R = float(global_settings_reac.annotation["cobrak_R"])
        T = float(global_settings_reac.annotation["cobrak_T"])
        reac_fwd_suffix = global_settings_reac.annotation["cobrak_reac_fwd_suffix"]
        reac_rev_suffix = global_settings_reac.annotation["cobrak_reac_rev_suffix"]
        reac_enz_separator = global_settings_reac.annotation[
            "cobrak_reac_enz_separator"
        ]
    else:
        max_prot_pool = 1_000.0
        extra_linear_constraints = []
        kinetic_ignored_metabolites = []
        R = STANDARD_R
        T = STANDARD_T
        reac_fwd_suffix = REAC_FWD_SUFFIX
        reac_rev_suffix = REAC_REV_SUFFIX
        reac_enz_separator = REAC_ENZ_SEPARATOR

    cobrak_metabolites: dict[str, Metabolite] = {}
    for metabolite in cobra_model.metabolites:
        if exclude_enzyme_constraints and sum(
            met_split in gene_ids for met_split in metabolite.id.split("_")
        ):
            continue

        if "cobrak_Cmin" in metabolite.annotation:
            log_min_conc = log(float(metabolite.annotation["cobrak_Cmin"]))
        else:
            log_min_conc = log(1e-6)
        if "cobrak_Cmax" in metabolite.annotation:
            log_max_conc = log(float(metabolite.annotation["cobrak_Cmax"]))
        else:
            log_max_conc = log(0.02)

        cobrak_metabolites[metabolite.id] = Metabolite(
            log_min_conc=log_min_conc,
            log_max_conc=log_max_conc,
            annotation={
                key: literal_eval(value) if "[" in value else value
                for key, value in metabolite.annotation.items()
                if not key.startswith("cobrak_")
            },
            formula=str(metabolite.formula),
            charge=metabolite.charge,
        )

    cobrak_reactions: dict[str, Reaction] = {}
    for reaction in cobra_model.reactions:
        if (
            reaction.id == "prot_pool_delivery"
            or reaction.id.startswith("enzyme_delivery_")
            or reaction.id.startswith("complex_delivery_")
            or reaction.id.startswith("cobrak_global_settings")
        ):
            continue

        version_data = [
            (key.replace("cobrak_id_", ""), reaction.annotation[key])
            for key in reaction.annotation
            if key.startswith("cobrak_id_")
        ]
        for version, version_reac_id in version_data:
            if f"cobrak_dG0_{version}" in reaction.annotation:
                dG0 = float(reaction.annotation[f"cobrak_dG0_{version}"])
            else:
                dG0 = None
            if f"cobrak_dG0_uncertainty_{version}" in reaction.annotation:
                dG0_uncertainty = float(
                    reaction.annotation[f"cobrak_dG0_uncertainty_{version}"]
                )
            else:
                dG0_uncertainty = None

            if f"cobrak_k_cat_{version}" in reaction.annotation:
                if reac_enz_separator in version_reac_id:
                    identifiers = (
                        (
                            version_reac_id.replace("_and", "").split(
                                reac_enz_separator
                            )[1]
                            + "\b"
                        )
                        .replace(f"{reac_fwd_suffix}\b", "")
                        .replace(f"{reac_rev_suffix}\b", "")
                        .replace("\b", "")
                        .split("_")
                    )
                else:
                    identifiers = reaction.gene_reaction_rule.split(" and ")

                k_cat = float(reaction.annotation[f"cobrak_k_cat_{version}"])
                if f"cobrak_k_ms_{version}" in reaction.annotation:
                    k_ms = literal_eval(reaction.annotation[f"cobrak_k_ms_{version}"])
                else:
                    k_ms = None
                if f"cobrak_k_is_{version}" in reaction.annotation:
                    k_is = literal_eval(reaction.annotation[f"cobrak_k_is_{version}"])
                else:
                    k_is = None
                if f"cobrak_k_as_{version}" in reaction.annotation:
                    k_as = literal_eval(reaction.annotation[f"cobrak_k_as_{version}"])
                else:
                    k_as = None
                if f"cobrak_special_stoichiometries_{version}" in reaction.annotation:
                    special_stoichiometries = literal_eval(
                        reaction.annotation[f"cobrak_special_stoichiometries_{version}"]
                    )
                else:
                    special_stoichiometries = {}
                enzyme_reaction_data = EnzymeReactionData(
                    identifiers=identifiers,
                    k_cat=k_cat,
                    k_ms=k_ms,
                    k_is=k_is,
                    k_as=k_as,
                    special_stoichiometries=special_stoichiometries,
                )
            else:
                enzyme_reaction_data = None

            if len(version_data) > 1:
                if version_reac_id.endswith(reac_rev_suffix):
                    stoich_multiplier = -1
                    min_flux = 0.0
                    max_flux = -reaction.lower_bound
                else:
                    stoich_multiplier = +1
                    min_flux = 0.0
                    max_flux = reaction.upper_bound
            else:
                min_flux = reaction.lower_bound
                max_flux = reaction.upper_bound
                stoich_multiplier = +1

            cobrak_reactions[version_reac_id] = Reaction(
                min_flux=min_flux,
                max_flux=max_flux,
                stoichiometries={
                    metabolite.id: stoich_multiplier * value
                    for (metabolite, value) in reaction.metabolites.items()
                    if (not exclude_enzyme_constraints)
                    or (
                        not sum(
                            met_split in gene_ids
                            for met_split in metabolite.id.split("_")
                        )
                    )
                },
                dG0=dG0,
                dG0_uncertainty=dG0_uncertainty,
                enzyme_reaction_data=enzyme_reaction_data,
                annotation={
                    key: literal_eval(value) if "[" in value else value
                    for key, value in reaction.annotation.items()
                    if not key.startswith("cobrak_")
                },
            )

    cobrak_enzymes: dict[str, Enzyme] = {}
    for gene in cobra_model.genes:
        if "cobrak_mw" in gene.annotation:
            if "cobrak_min_conc" in gene.annotation:
                min_conc = float(gene.annotation["cobrak_min_conc"])
            else:
                min_conc = None
            if "cobrak_max_conc" in gene.annotation:
                max_conc = float(gene.annotation["cobrak_max_conc"])
            else:
                max_conc = None
            cobrak_enzymes[gene.id] = Enzyme(
                molecular_weight=float(gene.annotation["cobrak_mw"]),
                min_conc=min_conc,
                max_conc=max_conc,
                name=gene.name,
            )

    return Model(
        reactions=cobrak_reactions,
        metabolites=cobrak_metabolites,
        enzymes=cobrak_enzymes,
        max_prot_pool=max_prot_pool,
        extra_linear_constraints=extra_linear_constraints,
        kinetic_ignored_metabolites=kinetic_ignored_metabolites,
        R=R,
        T=T,
        fwd_suffix=reac_fwd_suffix,
        rev_suffix=reac_rev_suffix,
        reac_enz_separator=reac_enz_separator,
    )


def load_annotated_sbml_model_as_cobrak_model(
    filepath: str,
) -> Model:
    """
    Load an annotated SBML model from a file and convert it into a COBRAk Model.

    This function reads an SBML file containing a metabolic model with specific annotations
    and converts it into a COBRAk Model. It uses the COBRApy library to read the SBML
    file and then uses the `load_annotated_cobrapy_model_as_cobrak_model` function to perform
    the conversion.

    Parameters:
    - filepath (str): The path to the SBML file containing the annotated metabolic model.

    Returns:
    - Model: A COBRAk Model constructed from the annotated SBML model, ready for further
      kinetic and thermodynamic analyses.
    """
    return load_annotated_cobrapy_model_as_cobrak_model(
        cobra.io.read_sbml_model(filepath)
    )


def load_unannotated_sbml_as_cobrapy_model(path: str) -> cobra.Model:
    """Loads an unannotated SBML model from a file into a COBRApy model.

    This function reads an SBML file that contains a metabolic model without specific annotations
    and loads it into a COBRApy model object. It utilizes the COBRApy library's `read_sbml_model`
    function to perform the loading.

    Parameters:
    - path (str): The file path to the SBML file containing the metabolic model.

    Returns:
    - cobra.Model: A COBRApy model object representing the metabolic network described in the SBML file.
    """
    return cobra.io.read_sbml_model(path)


def pickle_load(path: str) -> Any:  # noqa: ANN401
    """Returns the value of the given pickle file.

    Arguments
    ----------
    * path: str ~ The path to the pickle file.
    """
    pickle_file = open(path, "rb")
    pickled_object = pickle.load(pickle_file)
    pickle_file.close()
    return pickled_object


def pickle_write(path: str, pickled_object: Any) -> None:  # noqa: ANN401
    """Writes the given object as pickled file with the given path

    Arguments
    ----------
    * path: str ~ The path of the pickled file that shall be created
    * pickled_object: Any ~ The object which shall be saved in the pickle file
    """
    pickle_file = open(path, "wb")
    pickle.dump(pickled_object, pickle_file)
    pickle_file.close()


def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.
    If the given folder path is empty (''), it returns just ''.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Catch empty folders as they don't need to be standardized
    if not folder:
        return ""

    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder
