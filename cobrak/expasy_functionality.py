"""This module provides functionality to parse Expasy enzyme RDF files and extract EC number transfers."""

import xml.etree.ElementTree as ET


def get_ec_number_transfers(expasy_enzyme_rdf_path: str) -> dict[str, str]:
    """Parses an Expasy enzyme RDF file to extract enzyme EC number transfers.

    Args:
        expasy_enzyme_rdf_path (str): Path to the Expasy enzyme RDF file.

    Returns:
        dict[str, str]: A dictionary where each key is an EC number, and its corresponding value is the EC number it is transferred to.
                        The dictionary includes both directions of the transfer (old to new and new to old).
    """
    tree = ET.parse(expasy_enzyme_rdf_path)
    root = tree.getroot()

    ec_number_transfers: dict[str, str] = {}
    for child in root:
        for subchild in child:
            if "replaces" not in subchild.tag:
                continue
            new_ec_numbers = list(child.attrib.values())
            old_ec_numbers = list(subchild.attrib.values())
            for new_ec_number in new_ec_numbers:
                for old_ec_number in old_ec_numbers:
                    ec_number_transfers[old_ec_number] = new_ec_number
                    ec_number_transfers[new_ec_number] = old_ec_number
    return ec_number_transfers
