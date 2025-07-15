"""Functions and associated dataclasses for retrieving kinetic data from SABIO-RK"""

# IMPORT SECTION #
import threading
from dataclasses import dataclass
from io import StringIO
from math import sqrt
from os.path import exists
from shutil import rmtree
from statistics import median
from time import sleep
from zipfile import ZIP_LZMA, ZipFile

import cobra
import requests
from dataclasses_json import dataclass_json

from .dataclasses import EnzymeReactionData, ParameterReference
from .io import (
    ensure_folder_existence,
    get_files,
    json_load,
    json_zip_load,
    standardize_folder,
)
from .ncbi_taxonomy_functionality import (
    get_taxonomy_dict_from_nbci_taxonomy,
    get_taxonomy_scores,
)


# DATACLASSES SECTION #
@dataclass_json
@dataclass
class SabioEntry:
    """Represents the COBRAk-relevant data retrieved from a single SABIO-RK entry.

    Of which type this entry is (k_cat, k_m, k_i) is not determined here. This is done
    in the dataclass SabioDict.
    """

    entry_id: int
    """The entry's ID number"""
    is_recombinant: bool
    """Whether or not the entry is from a recombinant enzyme"""
    kinetics_mechanism_type: str
    """The reaction's kinetic mechanism (e.g., "Michaelis-Menten")"""
    organism: str
    """The organism (latin-greek name) associated with this entry"""
    temperature: float | None
    """[None if not given] The measurement's temperature in °C"""
    ph: float | None
    """[None if not given] The measurement's pH"""
    parameter_value: float
    """The value of the parameter"""
    parameter_unit: str
    """The unit of the value"""
    parameter_associated_species: str
    """The species (metabolite) associated with the parameter"""
    substrates: list[str]
    """The list of substrate names"""
    products: list[str]
    """The list of product names"""
    chebi_ids: list[str]
    """The list of all CHEBI IDs"""


@dataclass_json
@dataclass
class SabioDict:
    """Includes all retrieved SabioEntry instances and shows of which type they are"""

    kcat_entries: dict[str, list[SabioEntry]]
    """Turnover number entries"""
    km_entries: dict[str, list[SabioEntry]]
    """Michaelis-Menten constant entries"""
    ki_entires: dict[str, list[SabioEntry]]
    """Inhibition constant entries"""
    ka_entires: dict[str, list[SabioEntry]]
    """Activation constant entries"""
    hill_entires: dict[str, list[SabioEntry]]
    """Hill number entries"""


# CLASS SECTION #
class SabioThread(threading.Thread):
    """Represents a single Sabio-RK connection, ready for multi-threading (on one CPU core) using the threading module"""

    def __init__(self, temp_folder: str, start_number: int, end_number: int) -> None:
        """Initializes a SabioThread instance.

        Args:
            temp_folder (str): The path to the temporary folder where the results will be saved.
            start_number (int): The starting number for the query range.
            end_number (int): The ending number for the query range.
        """
        super().__init__()

        self.temp_folder = standardize_folder(temp_folder)
        self.start_number = start_number
        self.end_number = end_number

    def run(self) -> None:
        """Executes the thread's SABIO-RK data request

        Constructs a query string, sends a POST request to the SABIO-RK web service,
        and writes the response to a file in the temporary folder.
        """
        query_numbers = " OR ".join(
            [str(i + 1) for i in range(self.start_number, self.end_number + 1)]
        )
        query_dict = {"EntryID": f"({query_numbers})"}
        query_string = " AND ".join([f"{k}:{v}" for k, v in query_dict.items()])
        query_string += ' AND Parametertype:("activation constant" OR "Ki" OR "kcat" OR "km" OR "Hill coefficient") AND EnzymeType:"wildtype"'
        query = {
            "fields[]": [
                "EntryID",
                "Organism",
                "IsRecombinant",
                "ECNumber",
                "KineticMechanismType",
                "SabioCompoundID",
                "ChebiID",
                "Parameter",
                "Substrate",
                "Product",
                "Temperature",
                "pH",
            ],
            "q": query_string,
        }
        request = requests.post(
            "http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv",
            params=query,
            timeout=1e6,
        )
        request.raise_for_status()
        with open(
            f"{self.temp_folder}zzz{self.start_number}.txt", "w", encoding="utf-8"
        ) as f:
            f.write(request.text)


# "PRIVATE" FUNCTIONS SECTION #
def _get_ec_code_entries(
    sabio_entries: dict[str, list[SabioEntry]],
    reac_ec_codes: list[str],
    min_ph: float,
    max_ph: float,
    accept_nan_ph: bool,
    min_temperature: float,
    max_temperature: float,
    accept_nan_temperature: bool,
    substrate_bigg_ids: list[str],
    product_bigg_ids: list[str],
    name_to_bigg_id_dict: dict[str, str],
) -> list[SabioEntry]:
    """Filters SABIO-RK entries based on EC codes, pH, temperature, and substrate/product BIGG IDs.

    Args:
        sabio_entries (dict[str, list[SabioEntry]]): Dictionary of SabioEntry instances keyed by EC code.
        reac_ec_codes (list[str]): List of reaction EC codes to filter by.
        min_ph (float): Minimum pH value for filtering.
        max_ph (float): Maximum pH value for filtering.
        accept_nan_ph (bool): Whether to accept entries with NaN pH values.
        min_temperature (float): Minimum temperature value for filtering (in °C).
        max_temperature (float): Maximum temperature value for filtering  (in °C).
        accept_nan_temperature (bool): Whether to accept entries with NaN temperature values.
        substrate_bigg_ids (list[str]): List of substrate BIGG IDs to filter by.
        product_bigg_ids (list[str]): List of product BIGG IDs to filter by.
        name_to_bigg_id_dict (dict[str, str]): Dictionary mapping compound names to BIGG IDs.

    Returns:
        list[SabioEntry]: List of filtered SabioEntry instances.
    """
    ec_code_entries = []
    for ec_code, sabio_entry_list in sabio_entries.items():
        for sabio_entry in sabio_entry_list:
            if ec_code not in reac_ec_codes:
                continue

            # Check temperature
            if not accept_nan_temperature:  # noqa: SIM102
                if (
                    sabio_entry.temperature is None
                    or sabio_entry.temperature > max_temperature
                    or sabio_entry.temperature < min_temperature
                ):
                    continue

            # Check pH
            if not accept_nan_ph:  # noqa: SIM102
                if (
                    sabio_entry.ph is None
                    or sabio_entry.ph > max_ph
                    or sabio_entry.ph < min_ph
                ):
                    continue

            # Check correction
            correct_direction = False
            for substrate in sabio_entry.substrates:
                substrate_bigg_id = _search_metname_in_bigg_ids(
                    met_id=substrate.lower(),
                    bigg_id="",
                    entry=sabio_entry,
                    name_to_bigg_id_dict=name_to_bigg_id_dict,
                )
                if substrate_bigg_id in substrate_bigg_ids:
                    correct_direction = True
                    break
            if not correct_direction:
                for product in sabio_entry.products:
                    product_bigg_id = _search_metname_in_bigg_ids(
                        met_id=product.lower(),
                        bigg_id="",
                        entry=sabio_entry,
                        name_to_bigg_id_dict=name_to_bigg_id_dict,
                    )
                    if product_bigg_id in product_bigg_ids:
                        correct_direction = True
                        break
            if not correct_direction:
                continue

            # If all checks ok: Add this entry :-)
            ec_code_entries.append(sabio_entry)
    return ec_code_entries


def _get_sabio_tsv_str(target_folder: str) -> str:
    """Retrieves the SABIO-RK TSV string from the target folder.

    If the a zipped TSV cache does not exist, it downloads the data in threaded SabioThread instances,
    processes it, and stores it in a zip file cache.

    Args:
        target_folder (str): The path to the folder where the TSV file or zip file should be stored.

    Returns:
        str: The content of the SABIO-RK TSV file as a string.

    Example:
        target_folder = "/path/to/target/folder"
        tsv_str = _get_sabio_tsv_str(target_folder)
        print(tsv_str)
    """
    target_folder = standardize_folder(target_folder)
    ensure_folder_existence(target_folder)
    zip_filename = f"{target_folder}sabio_single_tsvs.zip"
    tsv_zip_filename = f"{target_folder}sabio_full_tsv.zip"
    sabio_tsv_filename = "sabio.tsv"

    starts = [0 + 250 * i for i in range(80_000 // 250)]

    if not exists(zip_filename):
        print("READING OUT SABIO-RK ONLINE...")

        temp_folder = f"{target_folder}sabiotemp/"
        ensure_folder_existence(temp_folder)

        start_and_end_numbers = [(start, start + 249) for start in starts]

        threads = [
            SabioThread(temp_folder, start_number, end_number)
            for start_number, end_number in start_and_end_numbers
        ]
        for thread in threads:
            thread.start()
            sleep(2.0)
        for thread in threads:
            thread.join()

        tsv_filenames = get_files(temp_folder)
        with ZipFile(zip_filename, "w", ZIP_LZMA) as zipf:
            for filename in tsv_filenames:
                zipf.write(temp_folder + filename, arcname=filename)

        rmtree(temp_folder)

    if not exists(tsv_zip_filename):
        full_tsv_lines = []
        first = True
        with ZipFile(zip_filename, "r") as zipf:
            file_names = zipf.namelist()
            for file_name in file_names:
                with zipf.open(file_name) as file:
                    lines = file.readlines()
                    for line in lines:
                        decoded_line = line.decode("utf-8")
                        if (not first) and ("EntryID	Organism" in decoded_line):
                            continue
                        if len(decoded_line) <= 1:
                            continue
                        full_tsv_lines.append(decoded_line.replace("\n", ""))
                first = False

        with ZipFile(tsv_zip_filename, "w", ZIP_LZMA) as zipf:
            zipf.writestr(
                sabio_tsv_filename, StringIO("\n".join(full_tsv_lines)).getvalue()
            )

    with ZipFile(tsv_zip_filename, "r") as zipf:  # noqa: SIM117
        with zipf.open(sabio_tsv_filename) as file:
            tsv_content = file.read().decode("utf-8")

    return tsv_content.replace("\r", "")


def _search_metname_in_bigg_ids(
    met_id: str,
    bigg_id: str,
    entry: SabioEntry,
    name_to_bigg_id_dict: dict[str, str],
) -> str:
    entry_bigg_id: str = ""
    for missing_addition_end in ["", "1", "__L", "__D"]:
        for missing_addition_front in ["", "alpha-", "beta-", "l-", "d-"]:
            addition_name = f"{missing_addition_front}{met_id}{missing_addition_end}"
            if addition_name in name_to_bigg_id_dict:
                entry_bigg_id = name_to_bigg_id_dict[addition_name]
                break
    if bigg_id:
        for chebi_id in entry.chebi_ids:
            if chebi_id in name_to_bigg_id_dict:
                entry_bigg_id = name_to_bigg_id_dict[chebi_id]
                if entry_bigg_id == bigg_id:
                    break
    return entry_bigg_id


# "PUBLIC" FUNCTIONS SECTION #
def get_full_sabio_dict(sabio_target_folder: str) -> SabioDict:
    """Parses a SABIO-RK web query TSV file from the target folder to create a SabioDict instance containing SABIO-RK entries.

    Args:
        sabio_target_folder (str): The path to the folder containing the TSV file.

    Returns:
        SabioDict: A SabioDict instance whichm in turn, contains SabioEntry instances
    """
    tsv_str = _get_sabio_tsv_str(sabio_target_folder)

    tsv_lines = tsv_str.split("\n")
    titles = tsv_lines[0].split("\t")
    del tsv_lines[0]

    sabio_dict = SabioDict({}, {}, {}, {}, {})
    for tsv_line in tsv_lines:
        line = tsv_line.split("\t")

        parameter_value_str = line[titles.index("parameter.startValue")]
        if not parameter_value_str:
            continue
        parameter_value = float(parameter_value_str)
        if parameter_value <= 0.0:
            continue  # There is no kinetic parameter that is just 0 or below

        parameter_type_str = line[titles.index("parameter.type")]
        match parameter_type_str.lower():
            case "kcat":
                sabio_dict_pointer = sabio_dict.kcat_entries
            case "km":
                sabio_dict_pointer = sabio_dict.km_entries
            case "ki":
                sabio_dict_pointer = sabio_dict.ki_entires
            case "activation constant":
                sabio_dict_pointer = sabio_dict.ka_entires
            case "hill coefficient":
                sabio_dict_pointer = sabio_dict.hill_entires
            case _:
                continue

        ec_number = line[titles.index("ECNumber")]
        entry_id = int(line[titles.index("EntryID")])
        organism = line[titles.index("Organism")]
        is_recombinant = line[titles.index("IsRecombinant")].lower() == "true"
        kinetics_mechanism_type = line[titles.index("KineticMechanismType")]
        parameter_unit = line[titles.index("parameter.unit")]
        parameter_associated_species = line[titles.index("parameter.associatedSpecies")]
        substrates = line[titles.index("Substrate")].split(";")
        products = line[titles.index("Product")].split(";")
        chebi_ids = line[titles.index("ChebiID")].split(";")
        try:
            temperature = float(line[titles.index("Temperature")])
        except (ValueError, IndexError):
            temperature = None
        try:
            ph = float(line[titles.index("pH")])
        except (ValueError, IndexError):
            ph = None

        if ec_number not in sabio_dict_pointer:
            sabio_dict_pointer[ec_number] = []
        sabio_dict_pointer[ec_number].append(
            SabioEntry(
                entry_id=entry_id,
                is_recombinant=is_recombinant,
                kinetics_mechanism_type=kinetics_mechanism_type,
                organism=organism,
                temperature=temperature,
                ph=ph,
                parameter_unit=parameter_unit,
                parameter_value=parameter_value,
                parameter_associated_species=parameter_associated_species,
                substrates=substrates,
                products=products,
                chebi_ids=chebi_ids,
            )
        )
    return sabio_dict


def sabio_select_enzyme_kinetic_data_for_model(
    cobra_model: cobra.Model,
    sabio_target_folder: str,
    base_species: str,
    ncbi_parsed_json_path: str,
    bigg_metabolites_json_path: str,
    kinetic_ignored_metabolites: list[str] = [],
    kinetic_ignored_enzyme_ids: list[str] = [],
    custom_enzyme_kinetic_data: dict[str, EnzymeReactionData | None] = {},
    min_ph: float = -float("inf"),
    max_ph: float = float("inf"),
    accept_nan_ph: bool = True,
    min_temperature: float = -float("inf"),
    max_temperature: float = float("inf"),
    accept_nan_temperature: bool = True,
    kcat_overwrite: dict[str, float] = {},
    transfered_ec_number_json: str = "",
    max_taxonomy_level: int = float("inf"),
) -> dict[str, EnzymeReactionData | None]:
    """Selects enzyme kinetic data for a given COBRA model using SABIO-RK data.

    If this data cannot be found, an internet connection is built to SABIO-RK and the relevant
    data is downloaded, which may take some time in the order of dozens of minutes.
    If you want to download the full SABIO-RK data beforehand, run get_full_sabio_dict() from
    this module beforehand, with the same sabio_target_folder.

    Collected data includes k_cat, k_m, k_i, k_a and Hill coefficients for all EC numbers that
    occur in the model's BiGG-compliant EC number annotation.

    Args:
        cobra_model (cobra.Model): The COBRA model for which enzyme kinetic data is to be selected.
        sabio_target_folder (str): The path to the folder containing SABIO-RK data.
        base_species (str): The base species for taxonomy comparison.
        ncbi_parsed_json_path (str): The path to the NCBI parsed JSON file.
        bigg_metabolites_json_path (str): The path to the BIGG metabolites JSON file.
        kinetic_ignored_metabolites (list[str], optional): List of metabolites to ignore. Defaults to [].
        kinetic_ignored_enzyme_ids (list[str], optional): List of enzyme IDs to ignore. Defaults to [].
        custom_enzyme_kinetic_data (dict[str, EnzymeReactionData | None], optional): Custom enzyme kinetic data. Defaults to {}.
        min_ph (float, optional): Minimum pH value for filtering. Defaults to -float("inf").
        max_ph (float, optional): Maximum pH value for filtering. Defaults to float("inf").
        accept_nan_ph (bool, optional): Whether to accept entries with NaN pH values. Defaults to True.
        min_temperature (float, optional): Minimum temperature value for filtering. Defaults to -float("inf").
        max_temperature (float, optional): Maximum temperature value for filtering. Defaults to float("inf").
        accept_nan_temperature (bool, optional): Whether to accept entries with NaN temperature values. Defaults to True.
        kcat_overwrite (dict[str, float], optional): Dictionary to overwrite kcat values. Defaults to {}.

    Returns:
        dict[str, EnzymeReactionData | None]: A dictionary mapping reaction IDs to enzyme kinetic data.

    Example:
        cobra_model = cobra.io.read_sbml_model("model.xml")
        sabio_target_folder = "/path/to/sabio/data"
        base_species = "Escherichia coli"
        ncbi_parsed_json_path = "/path/to/ncbi.json"
        bigg_metabolites_json_path = "/path/to/bigg_metabolites.json"
        enzyme_kinetic_data = sabio_select_enzyme_kinetic_data_for_model(
            cobra_model,
            sabio_target_folder,
            base_species,
            ncbi_parsed_json_path,
            bigg_metabolites_json_path
        )
        print(enzyme_kinetic_data)
    """
    sabio_dict = get_full_sabio_dict(
        sabio_target_folder,
    )
    ncbi_parsed_json_data = json_zip_load(ncbi_parsed_json_path)
    name_to_bigg_id_dict: dict[str, str] = json_load(
        bigg_metabolites_json_path, dict[str, str]
    )

    # Get reaction<->enzyme reaction data mapping
    enzyme_reaction_data: dict[str, EnzymeReactionData | None] = {}
    transfered_ec_codes: dict[str, str] = (
        json_load(transfered_ec_number_json, dict[str, str])
        if transfered_ec_number_json
        else {}
    )
    for reaction in cobra_model.reactions:
        if "ec-code" not in reaction.annotation:
            continue

        enzyme_identifiers = reaction.gene_reaction_rule.split(" and ")
        has_found_ignored_enzyme = False
        for enzyme_identifier in enzyme_identifiers:
            if enzyme_identifier in kinetic_ignored_enzyme_ids:
                has_found_ignored_enzyme = True
                break
        if has_found_ignored_enzyme:
            continue

        reac_met_ids = [met.id for met in reaction.metabolites]
        substrate_bigg_ids = [
            met_id[: met_id.rfind("_")]
            for met_id in reac_met_ids
            if reaction.metabolites[cobra_model.metabolites.get_by_id(met_id)] < 0
        ]
        product_bigg_ids = [
            met_id[: met_id.rfind("_")]
            for met_id in reac_met_ids
            if reaction.metabolites[cobra_model.metabolites.get_by_id(met_id)] > 0
        ]

        ec_codes = reaction.annotation["ec-code"]
        if isinstance(ec_codes, str):
            ec_codes = [ec_codes]
        reaction_transfered_ec_codes = [
            transfered_ec_codes[ec_code]
            for ec_code in ec_codes
            if ec_code in transfered_ec_codes
        ]
        ec_codes += reaction_transfered_ec_codes

        all_entries = (
            (
                "kcat",
                _get_ec_code_entries(
                    sabio_dict.kcat_entries,
                    ec_codes,
                    min_ph,
                    max_ph,
                    accept_nan_ph,
                    min_temperature,
                    max_temperature,
                    accept_nan_temperature,
                    substrate_bigg_ids,
                    product_bigg_ids,
                    name_to_bigg_id_dict,
                ),
            ),
            (
                "km",
                _get_ec_code_entries(
                    sabio_dict.km_entries,
                    ec_codes,
                    min_ph,
                    max_ph,
                    accept_nan_ph,
                    min_temperature,
                    max_temperature,
                    accept_nan_temperature,
                    substrate_bigg_ids,
                    product_bigg_ids,
                    name_to_bigg_id_dict,
                ),
            ),
            (
                "ki",
                _get_ec_code_entries(
                    sabio_dict.ki_entires,
                    ec_codes,
                    min_ph,
                    max_ph,
                    accept_nan_ph,
                    min_temperature,
                    max_temperature,
                    accept_nan_temperature,
                    substrate_bigg_ids,
                    product_bigg_ids,
                    name_to_bigg_id_dict,
                ),
            ),
            (
                "ka",
                _get_ec_code_entries(
                    sabio_dict.ka_entires,
                    ec_codes,
                    min_ph,
                    max_ph,
                    accept_nan_ph,
                    min_temperature,
                    max_temperature,
                    accept_nan_temperature,
                    substrate_bigg_ids,
                    product_bigg_ids,
                    name_to_bigg_id_dict,
                ),
            ),
        )

        # {'mol', 'katal*g^(-1)', 'M', 'M^2', 'g', 'mol/mol', 'J/mol', '-',
        # 's^(-1)', 's^(-1)*g^(-1)', 'mg/ml', 'mol*s^(-1)*mol^(-1)', 'M^(-1)', 'Pa',
        # 'M^(-1)*s^(-1)', 'mol*s^(-1)*g^(-1)', 'katal'}
        k_cat_per_tax_score: dict[int, list[float]] = {}
        k_cat_refs_per_tax_score: dict[int, list[ParameterReference]] = {}
        k_ms_per_tax_score: dict[str, dict[int, list[float]]] = {}
        k_m_refs_per_tax_score: dict[str, dict[int, list[ParameterReference]]] = {}
        k_is_per_tax_score: dict[str, dict[int, list[float]]] = {}
        k_i_refs_per_tax_score: dict[str, dict[int, list[ParameterReference]]] = {}
        k_as_per_tax_score: dict[str, dict[int, list[float]]] = {}
        k_a_refs_per_tax_score: dict[str, dict[int, list[ParameterReference]]] = {}
        for entries_type, entries in all_entries:
            if entries_type == "kcat":  # Reaction-wide search
                for entry in entries:
                    match entry.parameter_unit:
                        case "s^(-1)":
                            multiplier = 3_600
                        case _:
                            continue

                    taxonomy_dict = get_taxonomy_dict_from_nbci_taxonomy(
                        [base_species, entry.organism], ncbi_parsed_json_data
                    )
                    taxonomy_score = get_taxonomy_scores(base_species, taxonomy_dict)[
                        entry.organism
                    ]
                    if taxonomy_score > max_taxonomy_level:
                        continue
                    if taxonomy_score not in k_cat_per_tax_score:
                        k_cat_per_tax_score[taxonomy_score] = []
                        k_cat_refs_per_tax_score[taxonomy_score] = []
                    k_cat_per_tax_score[taxonomy_score].append(
                        entry.parameter_value * multiplier
                    )
                    k_cat_refs_per_tax_score[taxonomy_score].append(
                        ParameterReference(
                            database="SABIO-RK",
                            comment="SabioEntryID: " + str(entry.entry_id),
                            species=entry.organism,
                            substrate=entry.parameter_associated_species,
                            value=entry.parameter_value * multiplier,
                        )
                    )
            else:  # Metabolite-wide search
                match entries_type:
                    case "ka":
                        values_pointer = k_as_per_tax_score
                        ref_pointer = k_a_refs_per_tax_score
                    case "ki":
                        values_pointer = k_is_per_tax_score
                        ref_pointer = k_i_refs_per_tax_score
                    case "km":
                        values_pointer = k_ms_per_tax_score
                        ref_pointer = k_m_refs_per_tax_score
                    case _:
                        raise ValueError
                for met in cobra_model.metabolites:
                    if met.id in kinetic_ignored_metabolites:
                        continue
                    if (entries_type == "km") and met not in reaction.metabolites:
                        continue
                    bigg_id = met.id[: met.id.rfind("_")]
                    for entry in entries:
                        entry_met_id = (
                            entry.parameter_associated_species.lower().strip()
                        )
                        if entry_met_id in name_to_bigg_id_dict:
                            entry_bigg_id = name_to_bigg_id_dict[entry_met_id]
                        else:
                            entry_bigg_id = _search_metname_in_bigg_ids(
                                met_id=entry_met_id,
                                bigg_id="",
                                entry=entry,
                                name_to_bigg_id_dict=name_to_bigg_id_dict,
                            )
                            if not entry_bigg_id:
                                continue
                        if entry_bigg_id != bigg_id:
                            continue

                        match entry.parameter_unit:
                            case "M^2":
                                applier = sqrt
                            case "M^(-1)":
                                applier = lambda x: 1 / x  # noqa: E731
                            case "M":
                                applier = lambda x: x  # noqa: E731
                            case _:
                                continue
                        taxonomy_dict = get_taxonomy_dict_from_nbci_taxonomy(
                            [base_species, entry.organism], ncbi_parsed_json_data
                        )
                        taxonomy_score = get_taxonomy_scores(
                            base_species, taxonomy_dict
                        )[entry.organism]
                        if taxonomy_score > max_taxonomy_level:
                            continue

                        if met.id not in values_pointer:
                            values_pointer[met.id] = {}
                            ref_pointer[met.id] = {}
                        if taxonomy_score not in values_pointer[met.id]:
                            values_pointer[met.id][taxonomy_score] = []
                            ref_pointer[met.id][taxonomy_score] = []
                        values_pointer[met.id][taxonomy_score].append(
                            applier(entry.parameter_value)
                        )
                        ref_pointer[met.id][taxonomy_score].append(
                            ParameterReference(
                                database="SABIO-RK",
                                comment="SabioEntryID: " + str(entry.entry_id),
                                species=entry.organism,
                                substrate=entry.parameter_associated_species,
                                tax_distance=taxonomy_score,
                                value=applier(entry.parameter_value),
                            )
                        )

        if reaction.id in kcat_overwrite:
            k_cat = kcat_overwrite[reaction.id]
            k_cat_references = [
                ParameterReference(database="OVERWRITE", tax_distance=-1)
            ]
        elif (
            (reaction.id not in kcat_overwrite) and (kcat_overwrite != {})
        ) or not k_cat_per_tax_score:
            continue
        else:
            min_k_cat_tax_score = min(k_cat_per_tax_score.keys())
            k_cat = median(k_cat_per_tax_score[min_k_cat_tax_score])
            k_cat_references = k_cat_refs_per_tax_score[min_k_cat_tax_score]

        k_ms: dict[str, float] = {}
        k_m_references: dict[str, list[ParameterReference]] = {}
        for met_id, k_m_per_tax_score in k_ms_per_tax_score.items():
            k_ms[met_id] = median(k_m_per_tax_score[min(k_m_per_tax_score.keys())])
            k_m_references[met_id] = k_m_refs_per_tax_score[met_id][
                min(k_m_per_tax_score.keys())
            ]

        k_is: dict[str, float] = {}
        k_i_references: dict[str, list[ParameterReference]] = {}
        for met_id, k_i_per_tax_score in k_is_per_tax_score.items():
            k_is[met_id] = median(k_i_per_tax_score[min(k_i_per_tax_score.keys())])
            k_i_references[met_id] = k_i_refs_per_tax_score[met_id][
                min(k_i_per_tax_score.keys())
            ]

        k_as: dict[str, float] = {}
        k_a_references: dict[str, list[ParameterReference]] = {}
        for met_id, k_a_per_tax_score in k_as_per_tax_score.items():
            k_as[met_id] = median(k_a_per_tax_score[min(k_a_per_tax_score.keys())])
            k_a_references[met_id] = k_a_refs_per_tax_score[met_id][
                min(k_a_per_tax_score.keys())
            ]

        enzyme_reaction_data[reaction.id] = EnzymeReactionData(
            identifiers=enzyme_identifiers,
            k_cat=k_cat,
            k_cat_references=k_cat_references,
            k_ms=k_ms,
            k_m_references=k_m_references,
            k_is=k_is,
            k_i_references=k_i_references,
            k_as=k_as,
            k_a_references=k_a_references,
        )

    enzyme_reaction_data = {**enzyme_reaction_data, **custom_enzyme_kinetic_data}

    for reac_id in kcat_overwrite:  # noqa: PLC0206
        if reac_id not in enzyme_reaction_data:
            reaction = cobra_model.reactions.get_by_id(reac_id)
            enzyme_identifiers = reaction.gene_reaction_rule.split(" and ")
            enzyme_reaction_data[reac_id] = EnzymeReactionData(
                identifiers=enzyme_identifiers,
                k_cat=kcat_overwrite[reac_id],
                k_cat_references=[
                    ParameterReference(database="OVERWRITE", tax_distance=-1)
                ],
                k_ms={},
                k_is={},
            )

    return enzyme_reaction_data
