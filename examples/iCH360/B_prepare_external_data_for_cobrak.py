# IMPORTS SECTION #  # noqa: D100
import cobra
import z_add_path  # noqa: F401

from cobrak.brenda_functionality import brenda_select_enzyme_kinetic_data_for_model
from cobrak.io import ensure_folder_existence, json_write
from cobrak.sabio_rk_functionality import sabio_select_enzyme_kinetic_data_for_model
from cobrak.utilities import combine_enzyme_reaction_datasets, parse_external_resources
from examples.iCH360.A_general_model_constants import (
    kinetic_ignored_enzyme_ids,
    kinetic_ignored_metabolites,
)

# RUNNING SCRIPT SECTION #
common_input_folder = "examples/common_needed_external_resources/"
model_input_folder = "examples/iCH360/external_resources/"
output_folder = "examples/iCH360/prepared_external_resources/"
ensure_folder_existence(output_folder)

cobra_model = cobra.io.read_sbml_model(
    "examples/iCH360/external_resources/EC_iCH360_unadjusted_kcats.xml"
)

reac_ids = [reac.id for reac in cobra_model.reactions]
with open(f"{model_input_folder}drg0_prime_mean.csv", encoding="utf-8") as f:
    dG0_lines = f.readlines()

dG0s: dict[str, float] = {}
for dG0_line in dG0_lines:
    linesplit = dG0_line.split(",")
    reac_id = linesplit[0]
    if reac_id == "":
        continue
    dG0 = float(linesplit[1])
    if reac_id in reac_ids:
        dG0s[reac_id] = dG0
    if reac_id + "_fw" in reac_ids:
        dG0s[reac_id + "_fw"] = +dG0
    if reac_id + "_bw" in reac_ids:
        dG0s[reac_id + "_bw"] = -dG0
json_write(f"{output_folder}dG0s.json", dG0s)

molecular_weights = {}
kcat_per_h = {}
for reac_x in cobra_model.reactions:
    reac: cobra.Reaction = reac_x
    annotation_keys = reac.annotation.keys()
    if "smoment_enzyme" not in annotation_keys:
        continue
    reac.gene_reaction_rule = reac.annotation["smoment_enzyme"]
    molecular_weights[reac.annotation["smoment_enzyme"]] = float(
        reac.annotation["smoment_mw"]
    )
    kcat_per_h[reac.id] = float(reac.annotation["smoment_kcat_per_s"]) * 3_600
json_write(f"{output_folder}mws.json", molecular_weights)
json_write(f"{output_folder}kcats.json", kcat_per_h)

cobra_model.remove_metabolites([cobra_model.metabolites.enzyme_pool])
cobra_model.remove_reactions([cobra_model.reactions.enzyme_pool_supply])
cobra.io.write_sbml_model(
    cobra_model,
    f"{output_folder}EC_iCH360_unadjusted_kcats_cleaned.xml",
)

parse_external_resources(
    path=common_input_folder,
    brenda_version="2024_1",
)

brenda_enzyme_reaction_data = brenda_select_enzyme_kinetic_data_for_model(
    cobra_model=cobra_model,
    brenda_json_targz_file_path=f"{common_input_folder}brenda_2024_1.json.tar.gz",
    bigg_metabolites_json_path=f"{common_input_folder}bigg_models_metabolites.json",
    brenda_version="2024_1",
    base_species="Escherichia coli",
    ncbi_parsed_json_path=f"{common_input_folder}parsed_taxdmp.json",
    kinetic_ignored_metabolites=kinetic_ignored_metabolites,
    kinetic_ignored_enzyme_ids=kinetic_ignored_enzyme_ids,
    kcat_overwrite=kcat_per_h,
    transfered_ec_number_json=f"{common_input_folder}ec_number_transfers.json",
    max_taxonomy_level=6,
)

sabio_enzyme_reaction_data = sabio_select_enzyme_kinetic_data_for_model(
    cobra_model=cobra_model,
    sabio_target_folder="examples/common_needed_external_resources",
    bigg_metabolites_json_path=f"{common_input_folder}bigg_models_metabolites.json",
    base_species="Escherichia coli",
    ncbi_parsed_json_path=f"{common_input_folder}parsed_taxdmp.json",
    kinetic_ignored_metabolites=kinetic_ignored_metabolites,
    kinetic_ignored_enzyme_ids=kinetic_ignored_enzyme_ids,
    kcat_overwrite=kcat_per_h,
    transfered_ec_number_json=f"{common_input_folder}ec_number_transfers.json",
    max_taxonomy_level=6,
    add_hill_coefficients=False,
)

full_enzyme_reaction_data = combine_enzyme_reaction_datasets(
    [sabio_enzyme_reaction_data, brenda_enzyme_reaction_data],
)

json_write(
    f"{output_folder}enzyme_reaction_data_COMBINED.json",
    full_enzyme_reaction_data,
)
