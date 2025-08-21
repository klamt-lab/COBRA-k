"""pytest tests for COBRA-k's module ncbi_taxonomy"""

from cobrak.ncbi_taxonomy_functionality import (
    get_taxonomy_scores,
    most_taxonomic_similar,
)


def test_get_taxonomy_scores() -> None:  # noqa: D103
    base_species = "Escherichia coli"
    taxonomy_dict = {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    taxonomy_scores = get_taxonomy_scores(base_species, taxonomy_dict)
    assert taxonomy_scores["Escherichia coli"] == 0
    assert taxonomy_scores["Pseudomonas"] == 1
    assert taxonomy_scores["Homo sapiens"] == 2


def test_most_taxonomic_similar() -> None:  # noqa: D103
    base_species = "Escherichia coli"
    taxonomy_dict = {
        "Escherichia coli": ["Escherichia", "Bacteria", "Organism"],
        "Pseudomonas": ["Pseudomonas", "Bacteria", "Organism"],
        "Homo sapiens": ["Homo", "Mammalia", "Animalia", "Organism"],
    }
    most_similar = most_taxonomic_similar(base_species, taxonomy_dict)
    assert most_similar["Escherichia coli"] == 0
    assert most_similar["Pseudomonas"] == 1
    assert most_similar["Homo sapiens"] == 2
