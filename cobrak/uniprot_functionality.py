"""get_protein_mass_mapping.py

Functions for the generation of a model's mapping of its proteins and their masses.
"""

# IMPORTS SECTION #
import time
from typing import Any

import cobra
import requests

from .io import ensure_folder_existence, json_load, json_write, standardize_folder


# FUNCTIONS SECTION #
def uniprot_get_enzyme_molecular_weights(
    model: cobra.Model,
    cache_basepath: str,
    multiplication_factor: float = 1 / 1000,
) -> dict[str, float]:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are taken  from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format

    Output
    ----------
    A JSON file with the path project_folder+project_name+'_protein_id_mass_mapping.json'
    and the following structure:
    <pre>
    {
        "$PROTEIN_ID": $PROTEIN_MASS_IN_KDA,
        (...),
    }
    </pre>
    """

    # GET UNIPROT ID - PROTEIN MAPPING
    uniprot_id_protein_id_mapping: dict[str, list[str]] = {}
    for gene in model.genes:
        # Without a UniProt ID, no mass mapping can be found
        if "uniprot" not in gene.annotation:
            uniprot_id_protein_id_mapping[gene.id] = [gene.id]
            continue
        uniprot_id = gene.annotation["uniprot"]
        if uniprot_id in uniprot_id_protein_id_mapping:
            uniprot_id_protein_id_mapping[uniprot_id].append(gene.id)
        else:
            uniprot_id_protein_id_mapping[uniprot_id] = [gene.id]

    # GET UNIPROT ID<->PROTEIN MASS MAPPING
    uniprot_id_protein_mass_mapping: dict[str, float] = {}
    # The cache stored UniProt masses for already searched
    # UniProt IDs (each file in the cache folder has the name
    # of the corresponding UniProt ID). This prevents searching
    # UniProt for already found protein masses. :-)
    cache_basepath = standardize_folder(cache_basepath)
    ensure_folder_existence(cache_basepath)
    cache_filepath = f"{cache_basepath}_cache_uniprot_molecular_weights.json"
    try:
        cache_json = json_load(cache_filepath, Any)
    except Exception:
        cache_json = {}
    # Go through each batch of UniProt IDs (multiple UniProt IDs
    # are searched at once in order to save an amount of UniProt API calls)
    # and retrieve the amino acid sequences and using these sequences, their
    # masses.
    print("Starting UniProt ID<->Protein mass search using UniProt API...")
    uniprot_ids = list(uniprot_id_protein_id_mapping.keys())
    batch_size = 15
    batch_start = 0
    while batch_start < len(uniprot_ids):
        # Create the batch with all UniProt IDs
        prebatch = uniprot_ids[batch_start : batch_start + batch_size]
        batch = []
        # Remove all IDs which are present in the cache (i.e.,
        # which were searched for already).
        # The cache consists of pickled protein mass floats, each
        # onein a file with the name of the associated protein.
        for uniprot_id in prebatch:
            if uniprot_id not in cache_json:
                batch.append(uniprot_id)
            else:
                uniprot_id_protein_mass_mapping[uniprot_id] = cache_json[uniprot_id]
                # print(uniprot_id + ":", uniprot_id_protein_mass_mapping[uniprot_id])

        # If all IDs could be found in the cache, continue with the next batch.
        if len(batch) == 0:
            batch_start += batch_size
            continue

        # Create the UniProt query for the batch
        # With 'OR', all given IDs are searched, and subsequently in this script,
        # the right associated masses are being picked.
        query = " OR ".join(batch)
        uniprot_query_url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=tsv&fields=accession,id,mass"
        print(f"UniProt batch search for: {query}")

        # Call UniProt's API :-)
        uniprot_data = requests.get(uniprot_query_url, timeout=1e6).text.split("\n")
        # Wait in order to cool down their server :-)
        time.sleep(2.0)

        # Read out the API-returned lines
        found_ids = []
        for line in uniprot_data[1:]:
            if not line:
                continue
            accession_id = line.split("\t")[0].lstrip().rstrip()
            entry_id = line.split("\t")[1].lstrip().rstrip()
            mass_string = line.split("\t")[2].lstrip().rstrip()
            try:
                # Note that the mass entry from UniProt uses a comma as a thousand separator, so it has to be removed before parsing
                mass = float(mass_string.replace(",", ""))
            except ValueError:  # We may also risk the entry is missing
                # print(f"No protein mass obtainable for protein ID {uniprot_id}")
                continue
            uniprot_id_protein_mass_mapping[accession_id] = float(mass)
            uniprot_id_protein_mass_mapping[entry_id] = float(mass)
            found_ids.extend((accession_id, entry_id))

        # Create the pickled cache files for the searched protein masses
        for uniprot_id in found_ids:
            if (
                uniprot_id in uniprot_id_protein_mass_mapping
            ):  # Takes into account that we may fail to obtain a UniProt ID
                cache_json[uniprot_id] = uniprot_id_protein_mass_mapping[uniprot_id]

        # Continue with the next batch :D
        batch_start += batch_size

    # Create the final protein ID <-> mass mapping
    protein_id_mass_mapping: dict[str, float] = {}
    for uniprot_id in list(uniprot_id_protein_mass_mapping.keys()):
        try:
            protein_ids = uniprot_id_protein_id_mapping[uniprot_id]
        except Exception:
            continue
        for protein_id in protein_ids:
            protein_id_mass_mapping[protein_id] = (
                uniprot_id_protein_mass_mapping[uniprot_id] * multiplication_factor
            )

    # Return protein mass list JSON :D
    json_write(cache_filepath, cache_json)
    return protein_id_mass_mapping
