"""Contains functions which allow one to create a model-specific and BRENDA-depending kinetic data database."""

# IMPORTS SECTION #
import contextlib
import copy
import json
import tarfile
from math import isnan
from statistics import median
from typing import Any

import cobra

from .dataclasses import EnzymeReactionData, ParameterReference
from .io import json_load, json_zip_load
from .ncbi_taxonomy_functionality import (
    get_taxonomy_dict_from_nbci_taxonomy,
    most_taxonomic_similar,
)
from .sabio_rk_functionality import _search_metname_in_bigg_ids


# "PRIVATE" FUNCTIONS SECTION #
def _brenda_get_all_enzyme_kinetic_data_for_model(
    model: cobra.Model,
    brenda_json_targz_file_path: str,
    bigg_metabolites_json_path: str,
    brenda_version: str,
    min_ph: float = -float("inf"),
    max_ph: float = float("inf"),
    accept_nan_ph: bool = True,
    min_temperature: float = -float("inf"),
    max_temperature: float = float("inf"),
    accept_nan_temperature: bool = True,
    transfered_ec_codes: dict[str, str] = {},
) -> dict[str, Any]:
    """Reads out a BRENDA JSON file created with parse_brenda_textfile and creates a model-specific JSON.

    Arguments
    ----------
    * sbml_path: str ~ The path of the SBML model of which a specific BRENDA JSON kcat database
      shall be created

    Output
    ----------
    A JSON in the given folder and the name 'kcat_database_brenda.json', and with the following structure:
    <pre>
    {
        '$EC_NUMBER': {
            '$BIGG_ID_METABOLITE': {
                '$ORGANISM': [
                    kcat_list: float
                ],
                (...)
            },
            (...)
        },
        (...)
    }
    </pre>
    """
    brenda_kinetics_database_original = _brenda_parse_full_json(
        brenda_json_targz_file_path=brenda_json_targz_file_path,
        bigg_metabolites_json_path=bigg_metabolites_json_path,
        brenda_version=brenda_version,
        min_ph=min_ph,
        max_ph=max_ph,
        accept_nan_ph=accept_nan_ph,
        min_temperature=min_temperature,
        max_temperature=max_temperature,
        accept_nan_temperature=accept_nan_temperature,
    )

    # Get EC numbers of the model's reactions
    ec_numbers_of_model: list[str] = []
    for reaction in model.reactions:
        if "ec-code" not in reaction.annotation:
            continue

        ec_numbers_of_reaction = reaction.annotation["ec-code"]
        if isinstance(ec_numbers_of_reaction, str):
            ec_numbers_of_reaction = [ec_numbers_of_reaction]

        reaction_transfered_ec_codes = []
        for ec_code in ec_numbers_of_reaction:
            if ec_code in transfered_ec_codes:
                reaction_transfered_ec_codes.append(transfered_ec_codes[ec_code])
        ec_numbers_of_reaction += reaction_transfered_ec_codes

        ec_numbers_of_model += ec_numbers_of_reaction

    ec_numbers_of_model = list(set(ec_numbers_of_model))

    # Get EC number entries for each EC number of the model
    brenda_database_for_model = {}
    for ec_number in ec_numbers_of_model:
        entry_error = False
        if ec_number in brenda_kinetics_database_original:
            ec_number_entry = copy.deepcopy(
                brenda_kinetics_database_original[ec_number]
            )
            if "ERROR" in ec_number_entry:
                entry_error = True
            else:
                ec_number_entry["WILDCARD"] = False
                brenda_database_for_model[ec_number] = copy.deepcopy(ec_number_entry)

        if (ec_number not in brenda_kinetics_database_original) or entry_error:
            eligible_ec_number_entries: list[dict[str, Any]] = []
            for wildcard_level in range(1, 5):
                for database_ec_number in list(
                    brenda_kinetics_database_original.keys()
                ):
                    if _is_fitting_ec_numbers(
                        ec_number, database_ec_number, wildcard_level
                    ):
                        database_ec_number_entry = copy.deepcopy(
                            brenda_kinetics_database_original[database_ec_number]
                        )
                        if "ERROR" not in database_ec_number_entry:
                            eligible_ec_number_entries.append(database_ec_number_entry)
                if len(eligible_ec_number_entries) > 0:
                    break
            ec_number_entry = {}
            for eligible_ec_number_entry in eligible_ec_number_entries:
                for metabolite_key in eligible_ec_number_entry:
                    metabolite_entry = copy.deepcopy(
                        eligible_ec_number_entry[metabolite_key]
                    )
                    if metabolite_key not in ec_number_entry:
                        ec_number_entry[metabolite_key] = metabolite_entry
                    else:
                        ec_number_entry[metabolite_key] = {
                            **ec_number_entry[metabolite_key],
                            **metabolite_entry,
                        }
            ec_number_entry["WILDCARD"] = True
            brenda_database_for_model[ec_number] = ec_number_entry

    return brenda_database_for_model


def _brenda_parse_full_json(
    brenda_json_targz_file_path: str,
    bigg_metabolites_json_path: str,
    brenda_version: str,
    min_ph: float = -float("inf"),
    max_ph: float = float("inf"),
    accept_nan_ph: bool = True,
    min_temperature: float = -float("inf"),
    max_temperature: float = float("inf"),
    accept_nan_temperature: bool = True,
) -> dict[str, dict[str, Any]]:
    """Goes through a BRENDA database JSON and converts it into a machine-readable dictionary.

    The JSON includes kcats for found organisms and substrates.
    As of Sep 24 2024, the BRENDA database can be downloaded as JSON under
    https://www.brenda-enzymes.org/download.php

    The BRENDA database is not in a completely standardized format, that's why this function
    contains many convoluted checks and circumventions of non-standardized data.

    k_cat values from mutated enzymes are excluded (if this information can be successfully
    read out by this function).

    Output
    ----------
    * A dictionary containing the BRENDA JSON kinetic data in the following machine-readable format:
    <pre>
        {
            "$EC_NUMBER": {
                "$SUBSTRATE_WITH_BIGG_ID_1": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                },
                (...),
                "REST": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                }
            }
            (...),
        }
    </pre>
    'REST' stands for a substrate without found BIGG ID.
    """
    # Load BIGG ID <-> metabolite name mapping :D
    name_to_bigg_id_dict: dict[str, str] = json_load(
        bigg_metabolites_json_path, dict[str, str]
    )

    with tarfile.open(brenda_json_targz_file_path, "r:gz") as tar:
        json_filename = f"brenda_{brenda_version}.json"
        json_file = tar.extractfile(json_filename)
        if json_file is None:
            print(f"ERROR: BRENDA JSON in {json_filename} not found :O")
            raise FileNotFoundError
        brenda_json = json.load(json_file)

    result_json: dict[str, dict[str, dict[str, list[Any]]]] = {}
    for ec_number, ec_data in brenda_json["data"].items():
        if "turnover_number" not in ec_data:
            continue

        result_json[ec_number] = {}

        protein_to_organism_dict: dict[str, str] = {}
        protein_to_references_dict: dict[str, list[str]] = {}
        for key, valuedict in ec_data["protein"].items():
            protein_to_organism_dict[key] = valuedict["organism"]
            protein_to_references_dict[key] = valuedict["references"]

        for target_entry in ("turnover_number", "km_value", "ki_value"):
            if (target_entry != "turnover_number") and (target_entry not in ec_data):
                continue

            for kinetics_entry in ec_data[target_entry]:
                if "value" not in kinetics_entry:
                    continue

                try:
                    kinetic_value: float = float(
                        kinetics_entry["value"].split(" ")[0].lstrip().rstrip()
                    )
                except ValueError:
                    continue
                if kinetic_value <= 0.0:
                    # There is no realistic kinetic parameter that is just 0 or below
                    continue

                ph = float("nan")
                temperature = float("nan")
                if "comment" in kinetics_entry:
                    comment = kinetics_entry["comment"]
                    if ("mutant" in comment) or ("mutated" in comment):
                        continue
                    if comment.lower().count("ph ") == 1:
                        pH_split = comment.lower().split("ph ")
                        number = (
                            pH_split[1].split(" ")[0].replace(",", "").replace(")", "")
                        )
                        if len(number) > 0 and number[-1] == ".":
                            number = number[:-1]
                        with contextlib.suppress(ValueError):
                            ph = float(number)
                    if comment.lower().count("°c") == 1:
                        temperature_split = (
                            comment.lower().replace(" °c", "°c").split(" ")
                        )
                        temperature_parts = [x for x in temperature_split if "°c" in x]
                        temperature_string = (
                            temperature_parts[0]
                            .split(",")[0]
                            .replace("â", "")
                            .replace("/room", "")
                            .replace("(", "")
                            .replace(")", "")
                            .replace("/", "")
                        )
                        if temperature_string[0] != "-":
                            temperature_string = temperature_string.split("-")[0]
                        with contextlib.suppress(ValueError):
                            temperature = float(temperature_string.replace("°c", ""))

                if isnan(temperature):
                    if not accept_nan_temperature:
                        continue
                else:
                    if temperature > max_temperature:
                        continue
                    if temperature < min_temperature:
                        continue
                if isnan(ph):
                    if not accept_nan_ph:
                        continue
                else:
                    if ph > max_ph:
                        continue
                    if ph < min_ph:
                        continue

                substrate_raw = (
                    kinetics_entry["value"]
                    .split(" ")[1]
                    .replace("{", "")
                    .replace("}", "")
                    .lower()
                )
                substrate = ""
                if substrate_raw in name_to_bigg_id_dict:
                    substrate = name_to_bigg_id_dict[substrate_raw]
                else:
                    found = False
                    for redundant_addition in [
                        "",
                        "5'-",
                        "3-",
                        "deamido-",
                        "n-",
                        "n ",
                        "+",
                        "/in",
                        "+/in",
                        "[side 2]",
                        "-",
                        "l-",
                        "2'-",
                        "d-",
                        "d-(+)-",
                    ]:
                        replaced_substrate_raw = substrate_raw.replace(
                            redundant_addition, ""
                        )
                        replaced_substrate_raw = replaced_substrate_raw.strip()

                        for missing_addition_end in ["", "1", "__L", "__D"]:
                            for missing_addition_front in [
                                "",
                                "alpha-",
                                "beta-",
                                "l-",
                                "d-",
                            ]:
                                full_replaced_substrate_raw = (
                                    missing_addition_front
                                    + replaced_substrate_raw
                                    + missing_addition_end
                                )

                                if full_replaced_substrate_raw in name_to_bigg_id_dict:
                                    substrate = name_to_bigg_id_dict[
                                        full_replaced_substrate_raw
                                    ]
                                    found = True
                                    break
                            if found:
                                break
                        if found:
                            break

                    manual_replacements = [
                        ("beta-ala", "balabala"),
                        ("deamido-nad+", "nad"),
                        ("itaconate", "itacon"),
                        ("butyrate", "but"),
                        ("propionate", "but"),
                        ("n na+/in", "na1"),
                        ("butanal", "btal"),
                        ("butanol", "btoh"),
                        ("butan-1-ol", "btoh"),
                        ("butan-2-ol", "ibtol"),
                        ("propan-2-ol", "2ppoh"),
                        ("xylitol", "xylt"),
                        ("l-arabinitol", "abt"),
                        ("d-galactosamine", "galam"),
                        ("lactose", "lcts"),
                        ("xylitol", "xylt"),
                        ("l-sorbose", "srb__L"),
                        ("3',5'-cgmp", "35cgmp"),
                    ]
                    for manual_replacement in manual_replacements:
                        if substrate_raw == manual_replacement[0]:
                            substrate = name_to_bigg_id_dict[manual_replacement[1]]
                            found = True

                    if not found:
                        if target_entry == "kcat_km":
                            substrate = "REST"
                        else:
                            continue
                if substrate not in result_json[ec_number]:
                    result_json[ec_number][substrate] = {}

                if not substrate:
                    print(f"INFO: No substrate found for {ec_number}")
                    continue

                ref_num_to_pub: dict[str, dict[str, Any]] = ec_data["reference"]
                read_references = []
                if "references" in kinetics_entry:
                    for ref_number in kinetics_entry["references"]:
                        if ref_number not in ref_num_to_pub:
                            continue
                        read_references.append(
                            "PMID: " + str(ref_num_to_pub[ref_number]["pmid"])
                            if "pmid" in ref_num_to_pub[ref_number]
                            else ",".join(ref_num_to_pub[ref_number]["authors"])[:20]
                            + "..., "
                            + ref_num_to_pub[ref_number]["title"][:30]
                            + "..."
                            + str(ref_num_to_pub[ref_number]["year"])
                        )
                organisms = [
                    ec_data["protein"][protein_num]["organism"]
                    for protein_num in kinetics_entry["proteins"]
                ]
                substrates_list = []
                if "substrates_products" in ec_data:
                    for natural_products_entry in ec_data["substrates_products"]:
                        if set(kinetics_entry["proteins"]).isdisjoint(
                            set(natural_products_entry["proteins"])
                        ):
                            continue
                        reac_string = natural_products_entry["value"]
                        if "?" in reac_string:
                            continue
                        substrates_list = reac_string.split(" = ")[0].split(" + ")
                        for substrate_id in substrates_list.copy():
                            bigg_id = _search_metname_in_bigg_ids(
                                substrate_id.lower(),
                                bigg_id="",
                                entry=None,
                                name_to_bigg_id_dict=name_to_bigg_id_dict,
                            )
                            if bigg_id:
                                substrates_list.append(substrate_id)

                for organism in organisms:
                    if organism not in result_json[ec_number][substrate]:
                        result_json[ec_number][substrate][organism] = []
                    result_json[ec_number][substrate][organism].append(
                        [
                            target_entry,
                            kinetic_value,
                            read_references,
                            kinetics_entry.get("comment", ""),
                            kinetics_entry.get("value", ""),
                            set(substrates_list),
                        ]
                    )
    return result_json


def _is_fitting_ec_numbers(
    ec_number_one: str, ec_number_two: str, wildcard_level: int
) -> bool:
    """Check whether the EC numbers are the same under the used wildcard level.

    Arguments
    ----------
    * ec_number_one: str ~ The first given EC number.
    * ec_number_two: str ~ The second given EC number.
    * wildcard_level: int ~ The wildcard level.
    """
    if wildcard_level == 0:
        ec_number_one_full_numbers = ec_number_one.split(".")
        ec_number_two_full_numbers = ec_number_two.split(".")
    else:
        ec_number_one_full_numbers = ec_number_one.split(".")[:-wildcard_level]
        ec_number_two_full_numbers = ec_number_two.split(".")[:-wildcard_level]

    return ec_number_one_full_numbers == ec_number_two_full_numbers


# "PUBLIC" FUNCTIONS SECTION #
def brenda_select_enzyme_kinetic_data_for_model(
    cobra_model: cobra.Model,
    brenda_json_targz_file_path: str,
    bigg_metabolites_json_path: str,
    brenda_version: str,
    base_species: str,
    ncbi_parsed_json_path: str,
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
    """Select and assign enzyme kinetic data for each reaction in a COBRApy model based on BRENDA
    database entries and taxonomic similarity.

    This function retrieves enzyme kinetic data from a compressed BRENDA JSON file, merges it
    with BiGG metabolite translation data and taxonomy information from NCBI. It then iterates
    over the reactions in the provided COBRApy model to:

      - Filter reactions that have EC code annotations.
      - Identify eligible EC codes (ignoring those with hyphens).
      - Collect kinetic entries (e.g., turnover numbers, KM values, KI values) for each metabolite
        involved in the reaction.
      - Choose the best kinetic parameters (k_cat, k_ms, k_is) based on taxonomic similarity to
        a base species.
      - Apply conversion factors (e.g., s⁻¹ to h⁻¹ for k_cat, mM to M for KM and KI).
      - Respect ignore lists for metabolites and enzymes.
      - Override computed k_cat values if provided in the kcat_overwrite dictionary.
      - Merge with any custom enzyme kinetic data provided.

    Parameters:
        cobra_model (cobra.Model): A COBRA model object representing the metabolic network.
        brenda_json_targz_file_path (str): Path to the compressed JSON file containing
            BRENDA enzyme kinetic data.
        bigg_metabolites_json_path (str): Path to the JSON file mapping metabolite IDs to
            BiGG identifiers.
        brenda_version (str): String identifier for the BRENDA database version.
        base_species (str): Species identifier used as the reference for taxonomic similarity.
        ncbi_parsed_json_path (str): Path to the parsed JSON file containing NCBI taxonomy data.
        kinetic_ignored_metabolites (list[str], optional): List of metabolite IDs to exclude
            from kinetic parameter selection. Defaults to an empty list.
        kinetic_ignored_enzyme_ids (list[str], optional): List of enzyme identifiers to ignore
            when considering a reaction. Defaults to an empty list.
        custom_enzyme_kinetic_data (dict[str, EnzymeReactionData | None], optional):
            Dictionary of custom enzyme kinetic data to override or supplement computed data.
            The keys are reaction IDs and the values are EnzymeReactionData instances or None.
            Defaults to an empty dictionary.
        min_ph (float, optional): The minimum pH value for kinetic data inclusion. Defaults
            to negative infinity.
        max_ph (float, optional): The maximum pH value for kinetic data inclusion. Defaults
            to positive infinity.
        accept_nan_ph (bool, optional): If True, kinetic entries with NaN pH values are accepted.
            Defaults to True.
        min_temperature (float, optional): The minimum temperature value (e.g., in Kelvin) for
            kinetic data inclusion. Defaults to negative infinity.
        max_temperature (float, optional): The maximum temperature value for kinetic data inclusion.
            Defaults to positive infinity.
        accept_nan_temperature (bool, optional): If True, kinetic entries with NaN temperature values
            are accepted. Defaults to True.
        kcat_overwrite (dict[str, float], optional): Dictionary mapping reaction IDs to k_cat values
            that should override computed values. Defaults to an empty dictionary.

    Returns:
        dict[str, EnzymeReactionData | None]:
            A dictionary mapping reaction IDs (str) from the COBRA model to their corresponding
            EnzymeReactionData instances. If no suitable kinetic data are found (or if the enzyme
            is in the ignore list), the value will be None for that reaction.

    Notes:
        - Kinetic values are converted to standardized units:
            - k_cat values are converted from s⁻¹ to h⁻¹.
            - KM and KI values are converted from mM to M.
        - The function leverages taxonomic similarity (using NCBI TAXONOMY data)
          to select the most relevant kinetic values.
        - Custom enzyme kinetic data and k_cat overrides will replace any computed values.
    """
    transfered_ec_codes: dict[str, str] = (
        json_load(transfered_ec_number_json, dict[str, str])
        if transfered_ec_number_json
        else {}
    )
    brenda_database_for_model = _brenda_get_all_enzyme_kinetic_data_for_model(
        cobra_model,
        brenda_json_targz_file_path,
        bigg_metabolites_json_path,
        brenda_version,
        min_ph,
        max_ph,
        accept_nan_ph,
        min_temperature,
        max_temperature,
        accept_nan_temperature,
        transfered_ec_codes=transfered_ec_codes,
    )
    ncbi_parsed_json_data = json_zip_load(ncbi_parsed_json_path)

    bigg_metabolites_data: dict[str, str] = json_load(
        bigg_metabolites_json_path,
        dict[str, str],
    )

    # Get reaction<->enzyme reaction data mapping
    enzyme_reaction_data: dict[str, EnzymeReactionData | None] = {}
    for reaction in cobra_model.reactions:
        if "ec-code" not in reaction.annotation:
            continue

        substrate_names_and_ids = []
        for metabolite, stoichiometry in reaction.metabolites.items():
            if stoichiometry < 0:
                substrate_names_and_ids.extend((metabolite.id, metabolite.name.lower()))
                for checked_string in (metabolite.id, metabolite.name.lower()):
                    bigg_id = _search_metname_in_bigg_ids(
                        checked_string,
                        bigg_id="",
                        entry=None,
                        name_to_bigg_id_dict=bigg_metabolites_data,
                    )
                    if bigg_id:
                        substrate_names_and_ids.append(bigg_id)
        substrate_names_and_ids_set = set(substrate_names_and_ids)

        reaction_ec_codes = reaction.annotation["ec-code"]
        if isinstance(reaction_ec_codes, str):
            reaction_ec_codes = [reaction_ec_codes]
        eligible_reaction_ec_codes = [
            ec_code
            for ec_code in reaction_ec_codes
            if (ec_code in brenda_database_for_model) and ("-" not in ec_code)
        ]

        reaction_transfered_ec_codes = []
        for ec_code in eligible_reaction_ec_codes:
            if ec_code in transfered_ec_codes:
                single_transfered_ec_code = transfered_ec_codes[ec_code]
                if single_transfered_ec_code in brenda_database_for_model:
                    reaction_transfered_ec_codes.append(single_transfered_ec_code)
        eligible_reaction_ec_codes += reaction_transfered_ec_codes

        metabolite_entries: dict[str, dict[str, Any]] = {}
        for ec_code in eligible_reaction_ec_codes:
            ec_code_entry = brenda_database_for_model[ec_code]
            for met_id in ec_code_entry:
                if met_id == "WILDCARD":
                    continue
                if met_id not in metabolite_entries:
                    metabolite_entries[met_id] = {}
                for organism in ec_code_entry[met_id]:
                    if organism not in metabolite_entries[met_id]:
                        metabolite_entries[met_id][organism] = []
                    metabolite_entries[met_id][organism] += ec_code_entry[met_id][
                        organism
                    ]

        # Choose kcats and kms taxonomically
        best_kcat_taxonomy_level = float("inf")
        best_km_taxonomy_levels = {
            metabolite.id: float("inf") for metabolite in cobra_model.metabolites
        }
        best_ki_taxonomy_levels = {
            metabolite.id: float("inf") for metabolite in cobra_model.metabolites
        }
        taxonomically_best_kcats: list[float] = []
        taxonomically_best_kms: dict[str, list[float]] = {}
        taxonomically_best_kis: dict[str, list[float]] = {}
        k_cat_references: list[ParameterReference] = []
        k_m_references: dict[str, list[ParameterReference]] = {}
        k_i_references: dict[str, list[ParameterReference]] = {}
        for metabolite in cobra_model.metabolites:
            idx_last_underscore = metabolite.id.rfind("_")
            met_id = metabolite.id[:idx_last_underscore]
            if metabolite.id in kinetic_ignored_metabolites:
                continue
            if met_id not in metabolite_entries:
                continue
            organisms = list(metabolite_entries[met_id].keys())
            if base_species not in organisms:
                organisms.append(base_species)
            taxonomy_dict = get_taxonomy_dict_from_nbci_taxonomy(
                organisms, ncbi_parsed_json_data
            )
            taxonomy_similarities = most_taxonomic_similar(base_species, taxonomy_dict)
            highest_taxonomy_level = max(taxonomy_similarities.values())
            for taxonomy_level in range(highest_taxonomy_level + 1):
                if taxonomy_level > max_taxonomy_level:
                    continue

                level_organisms = [
                    organism
                    for organism in organisms
                    if taxonomy_similarities[organism] == taxonomy_level
                ]
                for level_organism in level_organisms:
                    if (level_organism not in metabolite_entries[met_id]) and (
                        level_organism == base_species
                    ):  # I.e., if it is the base species
                        continue
                    kinetic_entries = metabolite_entries[met_id][level_organism]
                    if taxonomy_level <= best_kcat_taxonomy_level:
                        kcat_entries = [
                            km_kcat_entry
                            for km_kcat_entry in kinetic_entries
                            if (km_kcat_entry[0] == "turnover_number")
                            and not (
                                substrate_names_and_ids_set.isdisjoint(km_kcat_entry[5])
                            )
                        ]

                        if len(kcat_entries) > 0:
                            if (
                                best_kcat_taxonomy_level > taxonomy_level
                            ):  # "Erase" if we find a better level
                                taxonomically_best_kcats = []
                            best_kcat_taxonomy_level = min(
                                taxonomy_level, best_kcat_taxonomy_level
                            )
                            if taxonomy_level <= best_kcat_taxonomy_level:
                                for kcat_entry in kcat_entries:
                                    taxonomically_best_kcats.append(
                                        kcat_entry[1] * 3_600
                                    )  # convert from s⁻¹ to h⁻¹
                                    k_cat_references.append(
                                        ParameterReference(
                                            database="BRENDA",
                                            comment=kcat_entry[3],
                                            species=level_organism,
                                            pubs=kcat_entry[2],
                                            substrate=kcat_entry[4],
                                            tax_distance=taxonomy_level,
                                            value=kcat_entry[1] * 3_600,
                                        )
                                    )

                    if taxonomy_level <= best_km_taxonomy_levels[metabolite.id]:
                        km_entries = [
                            km_kcat_entry
                            for km_kcat_entry in kinetic_entries
                            if km_kcat_entry[0] == "km_value"
                            and not (
                                substrate_names_and_ids_set.isdisjoint(km_kcat_entry[5])
                            )
                        ]
                        if len(km_entries) > 0:
                            if metabolite.id not in taxonomically_best_kms:
                                taxonomically_best_kms[metabolite.id] = []
                                k_m_references[metabolite.id] = []
                            if (
                                best_km_taxonomy_levels[metabolite.id] > taxonomy_level
                            ):  # "Erase" if we find a better level
                                taxonomically_best_kms[metabolite.id] = []
                            best_km_taxonomy_levels[metabolite.id] = min(
                                taxonomy_level, best_km_taxonomy_levels[metabolite.id]
                            )
                            if taxonomy_level <= best_km_taxonomy_levels[metabolite.id]:
                                for km_entry in km_entries:
                                    taxonomically_best_kms[metabolite.id].append(
                                        km_entry[1] / 1_000
                                    )  # convert from mM to M
                                    k_m_references[metabolite.id].append(
                                        ParameterReference(
                                            database="BRENDA",
                                            comment=km_entry[3],
                                            species=level_organism,
                                            pubs=km_entry[2],
                                            substrate=km_entry[4],
                                            tax_distance=taxonomy_level,
                                            value=km_entry[1] / 1_000,
                                        )
                                    )

                    if taxonomy_level <= best_ki_taxonomy_levels[metabolite.id]:
                        ki_entries = [
                            kinetic_entry
                            for kinetic_entry in kinetic_entries
                            if kinetic_entry[0] == "ki_value"
                            and not (
                                substrate_names_and_ids_set.isdisjoint(kinetic_entry[5])
                            )
                        ]
                        if len(ki_entries) > 0:
                            if metabolite.id not in taxonomically_best_kis:
                                taxonomically_best_kis[metabolite.id] = []
                                k_i_references[metabolite.id] = []
                            if (
                                best_ki_taxonomy_levels[metabolite.id] > taxonomy_level
                            ):  # "Erase" if we find a better level
                                taxonomically_best_kis[metabolite.id] = []
                            best_ki_taxonomy_levels[metabolite.id] = min(
                                taxonomy_level, best_ki_taxonomy_levels[metabolite.id]
                            )
                            if taxonomy_level <= best_ki_taxonomy_levels[metabolite.id]:
                                for ki_entry in ki_entries:
                                    taxonomically_best_kis[metabolite.id].append(
                                        ki_entry[1] / 1_000
                                    )  # convert from mM to M
                                    k_i_references[metabolite.id].append(
                                        ParameterReference(
                                            database="BRENDA",
                                            comment=ki_entry[3],
                                            species=level_organism,
                                            pubs=ki_entry[2],
                                            substrate=ki_entry[4],
                                            value=ki_entry[1] / 1_000,
                                            tax_distance=taxonomy_level,
                                        )
                                    )

        if reaction.id in kcat_overwrite:
            taxonomically_best_kcats = [kcat_overwrite[reaction.id]]
            k_cat_references = [
                ParameterReference(database="OVERWRITE", tax_distance=-1)
            ]
        else:
            taxonomically_best_kcats = []

        reaction_kms = {}
        for met_id, values in taxonomically_best_kms.items():
            if met_id not in [x.id for x in reaction.metabolites]:
                continue
            reaction_kms[met_id] = median(values)

        reaction_kis = {}
        for met_id, values in taxonomically_best_kis.items():
            if met_id not in taxonomically_best_kis:
                continue
            reaction_kis[met_id] = median(taxonomically_best_kis[met_id])

        enzyme_identifiers = reaction.gene_reaction_rule.split(" and ")
        has_found_ignored_enzyme = False
        for enzyme_identifier in enzyme_identifiers:
            if enzyme_identifier in kinetic_ignored_enzyme_ids:
                has_found_ignored_enzyme = True
                break

        if (len(taxonomically_best_kcats) > 0) and (not has_found_ignored_enzyme):
            reaction_kcat = median(taxonomically_best_kcats)  # or max(), min(), ...
            enzyme_reaction_data[reaction.id] = EnzymeReactionData(
                identifiers=enzyme_identifiers,
                k_cat=reaction_kcat,
                k_cat_references=k_cat_references,
                k_ms=reaction_kms,
                k_m_references=k_m_references,
                k_is=reaction_kis,
                k_i_references=k_i_references,
            )
        else:
            enzyme_reaction_data[reaction.id] = None

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
