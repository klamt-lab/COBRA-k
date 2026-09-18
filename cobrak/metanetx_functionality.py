"""Functionalities for reading out MetaNetX files."""

import os

from pydantic import ConfigDict, validate_call

from .dataclasses import Model
from .io import ensure_folder_existence, json_zip_load, json_zip_write


@validate_call(config=ConfigDict(arbitrary_types_allowed=True))
def clean_and_compress_mnx_files(
    chem_prop_filepath: str,
    chem_xref_filepath: str,
    output_dir: str,
) -> None:
    """
    Cleans data from two MetaNetX TSV files (chem_prop and chem_xref) and
    saves the cleaned versions as compressed JSON (.json.zip) files in a specified
    output directory.

    These cleaned versions are small enough to be stored in a GitHub repository :-)
    and can be directly used with COBRA-k's other MetaNetX functions to add SMILES to
    metabolites.

    The two files can be found here (as of Dec 2, 2025):
    https://www.metanetx.org/mnxdoc/mnxref.html

    Args:
        chem_prop_filename: The path to the 'chem_prop.tsv' file.
        chem_xref_filename: The path to the 'chem_xref.tsv' file.
        output_dir: The path to the directory where the cleaned, compressed
                    files will be saved.
    """

    # 1. Create the output directory if it doesn't exist
    ensure_folder_existence(output_dir)

    # --- Processing chem_xref.tsv ---
    chem_xref_output_path = os.path.join(output_dir, "chem_xref.json")
    print(f"Processing '{chem_xref_filepath}'...")

    chem_xref_dict: dict[str, str] = {}
    try:
        with open(chem_xref_filepath, encoding="utf-8") as f:
            for line in f:
                # Skip comments and empty lines
                if line.startswith("#") or len(line.strip()) == 0:
                    continue

                line = line.strip()  # noqa: PLW2901
                linesplit = line.split("\t")

                # Ensure there are enough columns
                if len(linesplit) > 1:
                    external_id = linesplit[0]
                    metanetx_id = linesplit[1]
                    chem_xref_dict[external_id] = metanetx_id

        # Write the cleaned data to a compressed GZ file
        json_zip_write(chem_xref_output_path, chem_xref_dict)
        print(f"Successfully saved compressed file to '{chem_xref_output_path}'")

    except FileNotFoundError:
        print(f"Error: Input file '{chem_xref_filepath}' not found.")
    except Exception as e:
        print(
            f"The following error occurred while processing '{chem_xref_filepath}': {e}"
        )

    # --- Processing chem_prop.tsv ---
    chem_prop_output_path = os.path.join(output_dir, "chem_prop.json")
    print(f"Processing '{chem_prop_filepath}'...")

    chem_prop_dict: dict[str, dict[str, str]] = {}
    try:
        with open(chem_prop_filepath, encoding="utf-8") as f:
            for line in f:
                # Skip comments and empty lines
                if line.startswith("#") or len(line.strip()) == 0:
                    continue

                line = line.strip()  # noqa: PLW2901
                linesplit = line.split("\t")

                # Ensure there are enough columns before accessing them
                if len(linesplit) > 8:
                    metanetx_id = linesplit[0]
                    colloquial_name = linesplit[1]
                    charge = linesplit[4]
                    mass = linesplit[5]
                    smiles = linesplit[8]

                    chem_prop_dict[metanetx_id] = {
                        "colloquial_name": colloquial_name,
                        "charge": charge,
                        "mass": mass,
                        "smiles": smiles,
                    }

        # Write the cleaned data to a compressed GZ file
        json_zip_write(chem_prop_output_path, chem_prop_dict)
        print(f"Successfully saved compressed file to '{chem_prop_output_path}'")

    except FileNotFoundError:
        print(f"Error: Input file '{chem_prop_filepath}' not found.")
    except Exception as e:
        print(
            f"The following error occurred while processing '{chem_prop_filepath}': {e}"
        )

    print("Cleanup and compression of MetaNetX tsv files complete!")


@validate_call(config=ConfigDict(arbitrary_types_allowed=True))
def add_smiles_annotation_to_metabolites(
    cobrak_model: Model,
    chem_prop_json_filepath: str,
    chem_xref_json_filepath: str,
    print_found_smiles: bool = False,
    print_not_found_smiles: bool = False,
    allowed_annotation_keys: list[str] = [],
) -> Model:
    """Annotates metabolites in a COBRA-k model with SMILES strings using preprocessed MetaNetX files.

    The function reads two gzipped JSON files (produced by COBRA-k's
    `clean_and_compress_mnx_files` function - see there):
    1. `chem_xref.json`
    2. `chem_prop.json`
    Note: The JSON's producedby COBRA-k's clean_and_compress_mnx_files are zipped, but
    you must not add the .zip suffix to the given file paths.

    It iterates through the model's metabolites, uses their existing annotations
    to find a matching MetaNetX ID, and then uses the MetaNetX ID to retrieve
    the SMILES string. The SMILES string is then added to the metabolite's
    annotation dictionary under the specified key.

    Args:
        cobrak_model: The COBRA-k `Model` object containing the metabolites to
                      be annotated.
        chem_prop_json_filepath: Path to the zipped JSON file containing
                                 MetaNetX chemical properties (MetaNetX ID to properties),
                                 without the .zip file ending.
        chem_xref_json_filepath: Path to the zipped JSON file containing
                                 MetaNetX cross-references (External ID to MetaNetX ID),
                                 without the .zip file ending.
        print_found_smiles: If True, prints a message for every metabolite
                            where a SMILES string was successfully added.
        print_not_found_smiles: If True, prints a message for every metabolite
                                where no SMILES string could be found.
        allowed_annotation_keys: An optional list of annotation keys (e.g.,
                                 'chebi', 'bigg.metabolite') to be considered only.
                                 If empty, all existing annotation keys are checked.
                                 Note: If a metabolite has multiple eligible annotations,
                                 the very first read out annotation with MetaNetX cross-reference
                                 is used. Thereby, the first annotation key in this list has the
                                 highest precedence. (default: [], i.e. all keys are considered)
        smiles_annotation_key: The key under which the SMILES string should be
                               stored in the metabolite's annotation dictionary
                               (default: 'smiles').

    Returns:
        The updated COBRA-k `Model` object with SMILES annotations added to
        the metabolites.
    """
    chem_xref_dict: dict[tuple[str, str], str] = json_zip_load(chem_xref_json_filepath)
    chem_prop_dict: dict[str, dict[str, str]] = json_zip_load(chem_prop_json_filepath)

    for met_id, met_data in cobrak_model.metabolites.items():
        metanetx_id: str = ""
        eligible_keys: list[str] = allowed_annotation_keys or list(
            met_data.annotation.keys()
        )
        for key in eligible_keys:
            if key not in met_data.annotation:
                continue
            values_unknown_type: list[str] | str = met_data.annotation[key]
            if type(values_unknown_type) is str:
                values: list[str] = [values_unknown_type]
            else:
                values: list[str] = values_unknown_type
            metanetx_id_found = False
            for value in values:
                annotation_id = f"{key}:{value.split(':')[0]}"
                if annotation_id not in chem_xref_dict:
                    continue
                metanetx_id = chem_xref_dict[annotation_id]
                metanetx_id_found = True
                break
            if metanetx_id_found:
                break

        if not metanetx_id or metanetx_id not in chem_prop_dict:
            if print_not_found_smiles:
                print(f"SMILES not found for {met_id}")
            continue
        smiles = chem_prop_dict[metanetx_id]["smiles"]
        met_data.smiles = smiles
        if print_found_smiles:
            print(f"SMILES found for {met_id}: {smiles}")

    return cobrak_model
