import json
import os
import pickle

import cobra

from cobrak.dataclasses import (
    Enzyme,
    EnzymeReactionData,
    Metabolite,
    Model,
    Reaction,
)
from cobrak.io import (
    _add_annotation_to_cobra_reaction,
    convert_cobrak_model_to_annotated_cobrapy_model,
    ensure_folder_existence,
    ensure_json_existence,
    get_base_id,
    get_files,
    get_folders,
    json_load,
    json_write,
    json_zip_load,
    json_zip_write,
    pickle_load,
    pickle_write,
    save_cobrak_model_as_annotated_sbml_model,
    standardize_folder,
)


def test_get_base_id():
    reac_id = "reaction_id"
    fwd_suffix = "fwd_suffix"
    rev_suffix = "rev_suffix"
    reac_enz_separator = "reac_enz_separator"

    base_id = get_base_id(reac_id, fwd_suffix, rev_suffix, reac_enz_separator)
    assert base_id == reac_id.replace(f"{fwd_suffix}\b", "").replace(
        f"{rev_suffix}\b", ""
    )


def test_ensure_folder_existence(tmp_path):
    folder = str(tmp_path / "test_folder")
    ensure_folder_existence(folder)
    assert os.path.isdir(folder)


def test_ensure_json_existence(tmp_path):
    path = str(tmp_path / "test_json.json")
    ensure_json_existence(path)
    assert os.path.isfile(path)


def test_add_annotation_to_cobra_reaction():
    cobra_reaction = cobra.Reaction("reaction_id")
    reac_id = "reaction_id"
    reac_data = Reaction(
        min_flux=0.0,
        max_flux=1.0,
        stoichiometries={"metabolite_id": 1.0},
        dG0=1.0,
        dG0_uncertainty=0.1,
        enzyme_reaction_data=EnzymeReactionData(
            identifiers=["enzyme_id"],
            k_cat=1.0,
            k_ms={"A": 1.0},
            k_is={"A": 1.0},
            k_as={"A": 1.0},
            special_stoichiometries={},
        ),
    )
    version = "V0"
    _add_annotation_to_cobra_reaction(cobra_reaction, reac_id, reac_data, version)
    assert cobra_reaction.annotation["cobrak_id_V0"] == reac_id
    assert cobra_reaction.annotation["cobrak_dG0_V0"] == "1.0"
    assert cobra_reaction.annotation["cobrak_dG0_uncertainty_V0"] == "0.1"
    assert cobra_reaction.annotation["cobrak_k_cat_V0"] == "1.0"
    assert cobra_reaction.annotation["cobrak_k_ms_V0"] == "{'A': 1.0}"
    assert cobra_reaction.annotation["cobrak_k_is_V0"] == "{'A': 1.0}"
    assert cobra_reaction.annotation["cobrak_k_as_V0"] == "{'A': 1.0}"
    assert cobra_reaction.annotation["cobrak_special_stoichiometries_V0"] == "{}"


def test_convert_cobrak_model_to_annotated_cobrapy_model(tmp_path):
    cobrak_model = Model(
        reactions={
            "reaction_id": Reaction(
                min_flux=0.0,
                max_flux=1.0,
                stoichiometries={"metabolite_id": 1.0},
                dG0=1.0,
                dG0_uncertainty=0.1,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["enzyme_id"],
                    k_cat=1.0,
                    k_ms={"A": 1.0},
                    k_is={"A": 1.0},
                    k_as={"A": 1.0},
                    special_stoichiometries={},
                ),
            )
        },
        metabolites={"metabolite_id": Metabolite(log_min_conc=0.0, log_max_conc=1.0)},
        enzymes={"enzyme_id": Enzyme(molecular_weight=1.0)},
        max_prot_pool=1.0,
        extra_linear_constraints=[],
        kinetic_ignored_metabolites=[],
        R=1.0,
        T=1.0,
    )
    cobra_model = convert_cobrak_model_to_annotated_cobrapy_model(
        cobrak_model, combine_base_reactions=False, add_enzyme_constraints=False
    )
    assert len(cobra_model.metabolites) == 1


def test_save_cobrak_model_as_annotated_sbml_model(tmp_path):
    cobrak_model = Model(
        reactions={
            "reaction_id": Reaction(
                min_flux=0.0,
                max_flux=1.0,
                stoichiometries={"metabolite_id": 1.0},
                dG0=1.0,
                dG0_uncertainty=0.1,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["enzyme_id"],
                    k_cat=1.0,
                    k_ms={"A": 1.0},
                    k_is={"A": 1.0},
                    k_as={"A": 1.0},
                    special_stoichiometries={},
                ),
            )
        },
        metabolites={"metabolite_id": Metabolite(log_min_conc=0.0, log_max_conc=1.0)},
        enzymes={"enzyme_id": Enzyme(molecular_weight=1.0)},
        max_prot_pool=1.0,
        extra_linear_constraints=[],
        kinetic_ignored_metabolites=[],
        R=1.0,
        T=1.0,
    )
    filepath = str(tmp_path / "test_sbml.xml")
    save_cobrak_model_as_annotated_sbml_model(cobrak_model, filepath)
    assert os.path.isfile(filepath)


def test_get_files(tmp_path):
    folder = str(tmp_path / "test_folder")
    os.makedirs(folder)
    with open(os.path.join(folder, "file1.txt"), "w", encoding="utf-8") as f:
        f.write("Hello, world!")
    with open(os.path.join(folder, "file2.txt"), "w", encoding="utf-8") as f:
        f.write("Goodbye, world!")
    files = get_files(folder)
    assert set(files) == {"file1.txt", "file2.txt"}


def test_get_folders(tmp_path):
    folder = str(tmp_path / "test_folder")
    os.makedirs(folder)
    os.makedirs(os.path.join(folder, "subfolder"))
    folders = get_folders(folder)
    assert folders == ["subfolder"]


def test_json_load(tmp_path):
    path = str(tmp_path / "test_json.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"key": "value"}, f)
    json_data = json_load(path, dict[str, str])
    assert json_data == {"key": "value"}


def test_json_write(tmp_path):
    path = str(tmp_path / "test_json.json")
    json_data = {"key": "value"}
    json_write(path, json_data)
    with open(path) as f:
        json_data_read = json.load(f)
    assert json_data_read == {"key": "value"}


def test_json_zip_write_and_load(tmp_path):
    path = str(tmp_path / "test_json.json")
    json_data = {"key": "value"}
    json_zip_write(path, json_data)
    json_data_read = json_zip_load(path)
    assert json_data_read == {"key": "value"}


def test_pickle_load(tmp_path):
    path = str(tmp_path / "test_pickle.pkl")
    pickle_write(path, {"key": "value"})
    pickled_object = pickle_load(path)
    assert pickled_object == {"key": "value"}


def test_pickle_write(tmp_path):
    path = str(tmp_path / "test_pickle.pkl")
    pickled_object = {"key": "value"}
    pickle_write(path, pickled_object)
    with open(path, "rb") as f:
        pickled_object_read = pickle.load(f)
    assert pickled_object_read == {"key": "value"}


def test_standardize_folder():
    folder = "C:\\test"
    standardized_folder = standardize_folder(folder)
    assert standardized_folder == "C:/test/"
