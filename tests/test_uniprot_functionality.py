"""pytest tests for COBRA-k's module uniprot_functionality"""

import json
import os
import tempfile
from shutil import rmtree

import cobra

from cobrak.uniprot_functionality import uniprot_get_enzyme_molecular_weights_for_sbml


def test_uniprot_get_enzyme_molecular_weights() -> None:  # noqa: D103
    # Create a test model
    model = cobra.Model("test_model")
    gene = cobra.Gene("gene1")
    gene.annotation = {"uniprot": "P12345"}
    model.genes.append(gene)
    reaction = cobra.Reaction("reaction1")
    reaction.gene_reaction_rule = "gene1"
    model.reactions.append(reaction)

    # Call the function
    cache_basepath = "test_cache"
    with tempfile.TemporaryDirectory() as tmp_dict:
        cobra.io.write_sbml_model(model, tmp_dict + "temp.xml")
        protein_id_mass_mapping = uniprot_get_enzyme_molecular_weights_for_sbml(
            tmp_dict + "temp.xml", cache_basepath
        )

    # Check that the function returns a dictionary
    assert isinstance(protein_id_mass_mapping, dict)

    # Check that the dictionary contains the expected protein ID
    assert "gene1" in protein_id_mass_mapping

    # Check that the protein mass is a float
    assert isinstance(protein_id_mass_mapping["gene1"], float)


def test_uniprot_get_enzyme_molecular_weights_cache() -> None:  # noqa: D103
    # Create a test model
    model = cobra.Model("test_model")
    gene = cobra.Gene("gene1")
    gene.annotation = {"uniprot": "P12345"}
    model.genes.append(gene)
    reaction = cobra.Reaction("reaction1")
    reaction.gene_reaction_rule = "gene1"
    model.reactions.append(reaction)

    # Create a cache file
    cache_basepath = "test_cache"
    cache_filepath = f"{cache_basepath}_cache_uniprot_molecular_weights.json"
    cache_json = {"P12345": 100.0}
    with open(cache_filepath, "w", encoding="utf-8") as f:
        json.dump(cache_json, f)

    # Call the function
    with tempfile.TemporaryDirectory() as tmp_dict:
        cobra.io.write_sbml_model(model, tmp_dict + "temp.xml")
        protein_id_mass_mapping = uniprot_get_enzyme_molecular_weights_for_sbml(
            tmp_dict + "temp.xml", cache_basepath
        )

    # Check that the function returns a dictionary
    assert isinstance(protein_id_mass_mapping, dict)

    # Check that the dictionary contains the expected protein ID
    assert "gene1" in protein_id_mass_mapping

    # Check that the protein mass is a float
    assert isinstance(protein_id_mass_mapping["gene1"], float)


def test_uniprot_get_enzyme_molecular_weights_no_uniprot_id() -> None:  # noqa: D103
    # Create a test model
    model = cobra.Model("test_model")
    gene = cobra.Gene("gene1")
    model.genes.append(gene)
    reaction = cobra.Reaction("reaction1")
    reaction.gene_reaction_rule = "gene1"
    model.reactions.append(reaction)

    # Call the function
    cache_basepath = "test_cache"
    with tempfile.TemporaryDirectory() as tmp_dict:
        cobra.io.write_sbml_model(model, tmp_dict + "temp.xml")
        protein_id_mass_mapping = uniprot_get_enzyme_molecular_weights_for_sbml(
            tmp_dict + "temp.xml", cache_basepath
        )

    # Check that the function returns an empty dictionary
    assert protein_id_mass_mapping == {}


def test_uniprot_get_enzyme_molecular_weights_invalid_cache_basepath() -> None:  # noqa: D103
    # Create a test model
    model = cobra.Model("test_model")
    gene = cobra.Gene("gene1")
    gene.annotation = {"uniprot": "P12345"}
    model.genes.append(gene)
    reaction = cobra.Reaction("reaction1")
    reaction.gene_reaction_rule = "gene1"
    model.reactions.append(reaction)


def test_uniprot_get_enzyme_molecular_weights_cleanup() -> None:  # noqa: D103
    # Create a test model
    model = cobra.Model("test_model")
    gene = cobra.Gene("gene1")
    gene.annotation = {"uniprot": "P12345"}
    model.genes.append(gene)
    reaction = cobra.Reaction("reaction1")
    reaction.gene_reaction_rule = "gene1"
    model.reactions.append(reaction)

    # Call the function
    cache_basepath = "test_cache"
    with tempfile.TemporaryDirectory() as tmp_dict:
        cobra.io.write_sbml_model(model, tmp_dict + "temp.xml")
        uniprot_get_enzyme_molecular_weights_for_sbml(
            tmp_dict + "temp.xml", cache_basepath
        )

    # Check that the cache file was created
    cache_filepath = f"{cache_basepath}_cache_uniprot_molecular_weights.json"
    assert os.path.exists(cache_filepath)

    # Clean up
    rmtree("./test_cache")
    os.remove("test_cache_cache_uniprot_molecular_weights.json")
