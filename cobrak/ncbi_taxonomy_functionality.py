"""ncbi_taxonomy.py

This module contains functions which can access NCBI TAXONOMY.
"""

# IMPORTS #
import os
from typing import Any
from zipfile import ZipFile

from .io import json_zip_write, standardize_folder


# PUBLIC FUNCTIONS #
def parse_ncbi_taxonomy(
    ncbi_taxdmp_zipfile_path: str,
    ncbi_parsed_json_path: str,
) -> None:
    """Parses NCBI taxonomy data from a taxdump zip file and saves it as a JSON file.

    This function extracts the necessary files from the NCBI taxdump zip archive, parses the taxonomy data,
    and writes the parsed data to a JSON file. The parsed data includes mappings from taxonomy numbers to names
    and vice versa, as well as the taxonomy tree structure.

    Args:
        ncbi_taxdmp_zipfile_path (str): The file path to the NCBI taxdump zip archive.
        ncbi_parsed_json_path (str): The file path where the parsed JSON data will be saved.
    """
    old_wd = os.getcwd()
    folder = standardize_folder(os.path.dirname(ncbi_taxdmp_zipfile_path))
    filename = os.path.basename(ncbi_taxdmp_zipfile_path)
    os.chdir(folder)

    with ZipFile(filename, "r") as zipfile:
        zipfile.extract("names.dmp")
        zipfile.extract("nodes.dmp")

    with open("names.dmp", encoding="utf-8") as f:
        name_lines = f.readlines()
    with open("nodes.dmp", encoding="utf-8") as f:
        node_lines = f.readlines()

    os.remove("names.dmp")
    os.remove("nodes.dmp")
    os.chdir(old_wd)

    parsed_json_data = {}

    number_to_names_dict: dict[str, Any] = {}
    names_to_number_dict = {}
    for line in name_lines:
        if ("scientific name" not in line) and ("synonym" not in line):
            continue
        number = line.split("|")[0].lstrip().rstrip()
        name = line.split("|")[1].lstrip().rstrip()
        if number not in number_to_names_dict:
            number_to_names_dict[number] = []
        number_to_names_dict[number].append(name)
        names_to_number_dict[name] = number

    parsed_json_data["number_to_names_dict"] = number_to_names_dict
    parsed_json_data["names_to_number_dict"] = names_to_number_dict

    nodes_dict = {}
    for line in node_lines:
        begin = line.split("|")[0].lstrip().rstrip()
        end = line.split("|")[1].lstrip().rstrip()
        if begin == end:
            nodes_dict[begin] = "END"
        else:
            nodes_dict[begin] = end
    parsed_json_data["nodes_dict"] = nodes_dict
    json_zip_write(ncbi_parsed_json_path, parsed_json_data)


def get_taxonomy_dict_from_nbci_taxonomy(
    organisms: list[str],
    parsed_json_data: dict[str, Any],
) -> dict[str, list[str]]:
    """Generates a taxonomy dictionary from NCBI taxonomy data.

    This function constructs a dictionary mapping each organism to its taxonomy path based on the provided NCBI taxonomy data.

    Args:
        organisms (list[str]): A list of organism names for which taxonomy paths are to be retrieved.
        parsed_json_data (dict[str, Any]): Parsed JSON data containing taxonomy information, including:
            - "number_to_names_dict": A dictionary mapping taxonomy numbers to names.
            - "names_to_number_dict": A dictionary mapping organism names to taxonomy numbers.
            - "nodes_dict": A dictionary representing the taxonomy tree structure.

    Returns:
        dict[str, list[str]]: A dictionary where each key is an organism name and the value is a list of taxonomy names
        representing the path from the organism to the root of the taxonomy tree.
    """
    number_to_names_dict = parsed_json_data["number_to_names_dict"]
    names_to_number_dict = parsed_json_data["names_to_number_dict"]
    nodes_dict = parsed_json_data["nodes_dict"]

    organism_to_taxonomy_dicts = {}
    for organism in organisms:
        try:
            node_train = [names_to_number_dict[organism]]
        except KeyError:
            organism_to_taxonomy_dicts[organism] = [organism, "all"]
            continue
        current_number = names_to_number_dict[organism]
        while True:
            next_number = nodes_dict[current_number]
            if next_number == "END":
                break
            node_train.append(next_number)
            current_number = next_number
        node_train_names = [number_to_names_dict[x][0] for x in node_train]
        organism_to_taxonomy_dicts[organism] = node_train_names
    return organism_to_taxonomy_dicts


def get_taxonomy_scores(
    base_species: str,
    taxonomy_dict: dict[str, list[str]],
) -> dict[str, int]:
    """Returns a dictionary with a taxonomic distance from the given organism.

    e.g. if base_species is "Escherichia coli" and taxonomy_dict is
    <pre>
    {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    </pre>
    this function would return
    <pre>
    {
        "Escherichia coli": 0,
        "Pseudomonas": 1,
        "Homo sapiens": 4,
    }
    </pre>

    Arguments
    ----------
    * base_species: str ~ The species to which a relation is made.
    * taxonomy_dict: dict[str, list[str]] ~ A dictionary with organism names as keys and
      their taxonomic levels (sorted from nearest to farthest) as string list.
    """
    base_species_taxonomy = taxonomy_dict[base_species]
    taxonomy_scores: dict[str, int] = {
        base_species: 0,
    }
    for other_species_name, other_species_taxonomy in taxonomy_dict.items():
        score = 0
        for taxonomy_part in base_species_taxonomy:
            if taxonomy_part in other_species_taxonomy:
                break
            score += 1
        taxonomy_scores[other_species_name] = score

    return taxonomy_scores


def most_taxonomic_similar(
    base_species: str, taxonomy_dict: dict[str, list[str]]
) -> dict[str, int]:
    """Returns a dictionary with a score of taxonomic distance from the given organism.

    e.g. if base_species is "Escherichia coli" and taxonomy_dict is
    <pre>
    {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    </pre>
    this function would return
    <pre>
    {
        "Escherichia coli": 0,
        "Pseudomonas": 1,
        "Homo sapiens": 2,
    }
    </pre>

    Arguments
    ----------
    * base_species: str ~ The species to which a relation is made.
    * taxonomy_dict: dict[str, list[str]] ~ A dictionary with organism names as keys and
      their taxonomic levels (sorted from nearest to farthest) as string list.
    """
    base_taxonomy = taxonomy_dict[base_species]
    level_dict: dict[str, int] = {}
    for level, taxonomic_level in enumerate(base_taxonomy):
        level_dict[taxonomic_level] = level

    score_dict: dict[str, int] = {}
    for species, taxonomic_levels in taxonomy_dict.items():
        for taxonomic_level in taxonomic_levels:
            if taxonomic_level in level_dict:
                score_dict[species] = level_dict[taxonomic_level]
                break

    return score_dict
