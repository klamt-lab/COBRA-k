"""bigg_parse_metabolites_file.py

This module contains a function which transforms a BIGG metabolites .txt list
into an machine-readable JSON.
"""

# IMPORTS SECTION #
from .io import json_write


# PUBLIC FUNCTIONS SECTION #
def bigg_parse_metabolites_file(
    bigg_metabolites_txt_path: str,
    bigg_metabolites_json_path: str,
) -> None:
    """Parses a BIGG metabolites text file and returns a dictionary for this file.

    As of Sep 14 2024, a BIGG metabolites list of all BIGG-included metabolites
    is retrievable under http://bigg.ucsd.edu/data_access

    Arguments
    ----------
    * bigg_metabolites_file_path: str ~ The file path to the BIGG metabolites file.
      The usual file name (which has to be included too in this argument) is
      bigg_models_metabolites.txt
    * output_folder: str ~ The folder in which the JSON including the parsed BIGG
      metabolites file data is stored with the name 'bigg_id_name_mapping.json'

    Output
    ----------
    * A JSON file with the name 'bigg_id_name_mapping.json' in the given output folder,
      with the following structure:
    <pre>
     {
         "$BIGG_ID": "$CHEMICAL_OR_USUAL_NAME",
         (...),
         "$BIGG_ID": "$BIGG_ID",
         (...),
     }
    </pre>
    The BIGG ID <-> BIGG ID mapping is done for models which already use the BIGG IDs.
    """
    # Open the BIGG metabolites file as string list, and remove all newlines
    with open(bigg_metabolites_txt_path, encoding="utf-8") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "") for x in lines if len(x) > 0]

    # Mapping variable which will store the BIGG ID<->
    bigg_id_name_mapping = {}
    # Go through each BIGG metabolites file line (which is a tab-separated file)
    # and retrieve the BIGG ID and the name (if there is a name for the given BIGG
    # ID)
    for line in lines:
        bigg_id = line.split("\t")[1]
        bigg_id_name_mapping[bigg_id] = bigg_id

        # Exception to check if there is no name :O
        try:
            name = line.split("\t")[2].lower()
        except Exception:
            continue
        bigg_id_name_mapping[name] = bigg_id

        try:
            database_links = line.split("\t")[4]
        except Exception:
            continue
        for database_link_part in database_links.split(": "):
            if "CHEBI:" not in database_link_part:
                continue
            subpart = database_link_part.split("CHEBI:")[1].strip()
            chebi_id = subpart.split("; ")[0] if "; " in subpart else subpart
            bigg_id_name_mapping[chebi_id] = bigg_id

    # Write the JSON in the given folder :D
    json_write(bigg_metabolites_json_path, bigg_id_name_mapping)
