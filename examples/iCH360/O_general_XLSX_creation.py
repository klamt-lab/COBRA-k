from itertools import starmap
from math import exp
from os.path import exists

import z_add_path  # noqa: F401

from cobrak.constants import LNCONC_VAR_PREFIX, OBJECTIVE_VAR_NAME
from cobrak.dataclasses import Model
from cobrak.io import ensure_folder_existence, get_files, json_load, standardize_folder
from cobrak.spreadsheet_functionality import (
    OptimizationDataset,
    VariabilityDataset,
    create_cobrak_spreadsheet,
)
from cobrak.utilities import (
    get_fwd_rev_corrected_flux,
    get_reaction_enzyme_var_id,
    sort_dict_keys,
)

datafolder_settings = (
    # (
    #    "examples/iCH360/RESULTS_GLCUPTAKE",
    #    "maxglc",
    #    "1000",
    #    [2, 3, 4, 5],
    # ),
    (
        "examples/iCH360/RESULTS_MAXAC",
        "_maxac",
        "",
        [2, 3, 4, 5],
    ),
    (
        "examples/iCH360/RESULTS_MAXAC_MINMU",
        "_maxac_minmu0.1",
        "",
        [2, 3, 4, 5],
    ),
)

accombined_tfvas = {}
accombined_tfbas = {}
accombined_extrastatistics = {}
for raw_datafolder, suffix, extrasuffix, other_indices in datafolder_settings:
    datafolder = standardize_folder(raw_datafolder)
    cobrak_model: Model = json_load(
        f"{datafolder}used_cobrak_model__1_{suffix}{extrasuffix}.json",
        Model,
    )
    variability_dict = json_load(
        f"{datafolder}variability_dict__1_{suffix}{extrasuffix}.json",
        dict[str, tuple[float, float]],
    )

    optimization_datasets = {
        final_best_result_path.split("final_best_result__")[1].replace(
            ".json", ""
        ): OptimizationDataset(
            json_load(final_best_result_path),
            with_df=True,
            with_vplus=True,
            with_kappa=True,
            with_gamma=True,
            with_kinetic_differences=True,
        )
        for final_best_result_path in [
            datafolder + filename
            for filename in get_files(datafolder)
            if filename.startswith("final_best_result_")
            and filename.endswith(".json")
            and ("time" not in filename)
            and ("_6_" not in filename)
            and ("_7_" not in filename)
            and ("_8_" not in filename)
        ]
    }
    optimization_datasets = sort_dict_keys(optimization_datasets)

    glcvalues = (
        1.0,
        1.33,
        1.66,
        2.0,
        2.33,
        2.66,
        3.0,
        3.33,
        3.66,
        4.0,
        4.33,
        4.66,
        5.0,
        5.33,
        5.66,
        6.0,
        6.33,
        6.66,
        7.0,
        7.33,
        7.66,
        8.0,
        8.33,
        8.66,
        9.0,
        9.33,
        9.65,
        10.0,
        1000,
        "ALL",
    )

    # Generate conc and flux comparisons
    first_index = 1
    reac_ids = list(cobrak_model.reactions.keys())
    met_ids = [met_id for met_id in cobrak_model.metabolites if met_id.endswith("_c")]
    enzyme_var_ids = list(
        starmap(get_reaction_enzyme_var_id, cobrak_model.reactions.items())
    )

    max_max_flux_difference = -float("inf")
    max_max_met_difference = -float("inf")
    max_max_enzyme_difference = -float("inf")
    max_flux_difference_dict: dict[float, float] = {}
    max_met_difference_dict: dict[float, float] = {}
    max_enzyme_difference_dict: dict[float, float] = {}
    for glcuptake in glcvalues:
        if suffix == "maxglc" and glcuptake == "ALL":
            continue
        if suffix != "maxglc" and glcuptake != "ALL":
            continue
        first_json_path = (
            f"{datafolder}final_best_result__{first_index}_{suffix}{glcuptake}.json"
            if suffix == "maxglc"
            else f"{datafolder}final_best_result__{first_index}_{suffix}.json"
        )
        if not exists(first_json_path):
            continue
        first_dataset = json_load(first_json_path)
        used_first_dataset = {}
        for met_id in met_ids:
            if LNCONC_VAR_PREFIX + met_id in first_dataset:
                used_first_dataset[met_id] = exp(
                    first_dataset[LNCONC_VAR_PREFIX + met_id]
                )
            else:
                used_first_dataset[met_id] = 0.0
        for reac_id in reac_ids:
            if reac_id not in first_dataset:
                used_first_dataset[reac_id] = 0.0
            else:
                used_first_dataset[reac_id] = first_dataset[reac_id]
        for enzyme_var_id in enzyme_var_ids:
            if enzyme_var_id not in first_dataset:
                used_first_dataset[enzyme_var_id] = 0.0
            else:
                used_first_dataset[enzyme_var_id] = first_dataset[enzyme_var_id]

        for other_index in other_indices:
            second_dataset = json_load(
                f"{datafolder}final_best_result__{other_index}_{suffix}{glcuptake}.json"
                if suffix == "maxglc"
                else f"{datafolder}final_best_result__{other_index}_{suffix}.json"
            )
            used_second_dataset = {}
            for met_id in met_ids:
                if LNCONC_VAR_PREFIX + met_id in second_dataset:
                    used_second_dataset[met_id] = exp(
                        second_dataset[LNCONC_VAR_PREFIX + met_id]
                    )
                else:
                    used_second_dataset[met_id] = 0.0
            for reac_id in reac_ids:
                if reac_id not in second_dataset:
                    used_second_dataset[reac_id] = 0.0
                else:
                    used_second_dataset[reac_id] = second_dataset[reac_id]
            for enzyme_var_id in enzyme_var_ids:
                if enzyme_var_id not in second_dataset:
                    used_second_dataset[enzyme_var_id] = 0.0
                else:
                    used_second_dataset[enzyme_var_id] = second_dataset[enzyme_var_id]

            max_met_difference = -float("inf")
            max_flux_difference = -float("inf")
            max_enzyme_difference = -float("inf")
            for key, value in used_first_dataset.items():
                if key in met_ids:
                    old_d = max_met_difference
                    max_met_difference = max(
                        max_met_difference, abs(value - used_second_dataset[key])
                    )
                elif key in reac_ids:
                    corrected_flux_1 = get_fwd_rev_corrected_flux(
                        key,
                        list(used_first_dataset.keys()),
                        used_first_dataset,
                        cobrak_model.fwd_suffix,
                        cobrak_model.rev_suffix,
                    )
                    corrected_flux_2 = get_fwd_rev_corrected_flux(
                        key,
                        list(used_second_dataset.keys()),
                        used_second_dataset,
                        cobrak_model.fwd_suffix,
                        cobrak_model.rev_suffix,
                    )
                    max_flux_difference = max(
                        max_flux_difference, corrected_flux_1 - corrected_flux_2
                    )
                elif key in enzyme_var_ids:
                    max_enzyme_difference = max(
                        max_enzyme_difference, abs(value - used_second_dataset[key])
                    )
            max_flux_difference_dict[glcuptake] = max_flux_difference
            max_met_difference_dict[glcuptake] = max_met_difference
            max_enzyme_difference_dict[glcuptake] = max_enzyme_difference
            max_max_met_difference = max(max_max_met_difference, max_met_difference)
            max_max_flux_difference = max(max_max_flux_difference, max_flux_difference)
            max_max_enzyme_difference = max(
                max_max_enzyme_difference, max_enzyme_difference
            )

    print("===")
    print("Overall maximally occurring differences:")
    print(
        "Maximal overall flux difference:",
        max_max_flux_difference,
    )
    print(
        "Maximal overall metabolite concentration difference:",
        max_max_met_difference,
    )
    print(
        "Maximal overall enzyme concentration difference:",
        max_max_enzyme_difference,
    )
    print("===")

    # Actual XLSX generation
    for glcvalue in glcvalues:
        if glcvalue != "ALL":
            continue
        print(
            "Creating spreadsheets for maximal glucose uptake [g/(mmol*h)] of:",
            glcvalue,
        )

        # Collect data
        glcvalue_datasets = {
            dataset_name: opt_dataset
            for dataset_name, opt_dataset in optimization_datasets.items()
            if (f"{suffix}{glcvalue}" in dataset_name)
            or glcvalue == "ALL"
            and ("time" not in dataset_name)
        }

        if len(glcvalue_datasets) == 0:
            print(" INFO: No data yet for this value")
            continue

        extra_optstatistics_data = {}
        if (glcvalue == "ALL") and any(f"{suffix}" in key for key in glcvalue_datasets):
            try:  # noqa: PLC1901
                sorted_dataset_key = sorted(
                    glcvalue_datasets.keys(),
                    key=lambda x: float(x.split(f"{suffix}")[1]),
                )
            except ValueError:
                sorted_dataset_key = sorted(
                    glcvalue_datasets.keys(),
                )
            glcvalue_datasets = (
                {
                    sorted_dataset_key.replace("1_", ""): glcvalue_datasets[
                        sorted_dataset_key
                    ]
                    for sorted_dataset_key in sorted_dataset_key
                    if sorted_dataset_key.startswith("1_")
                }
                if extrasuffix
                else {
                    sorted_dataset_key.replace("1_", ""): glcvalue_datasets[
                        sorted_dataset_key
                    ]
                    for sorted_dataset_key in sorted_dataset_key
                }
            )
            extra_optstatistics_data = (
                {
                    "Max flux differences (n=5) [mmol⋅gDW⁻¹]": [
                        max_flux_difference_dict[float(key.replace(f"{suffix}", ""))]
                        for key in glcvalue_datasets
                    ],
                    "Max metabolite concentration differences (n=5) [M]": [
                        max_met_difference_dict[float(key.replace(f"{suffix}", ""))]
                        for key in glcvalue_datasets
                    ],
                    "Max enzyme concentration differences (n=5) [mmol⋅gDW⁻¹]": [
                        max_enzyme_difference_dict[float(key.replace(f"{suffix}", ""))]
                        for key in glcvalue_datasets
                    ],
                }
                if extrasuffix
                else {
                    "Max flux differences (n=5) [mmol⋅gDW⁻¹]": [
                        max_flux_difference_dict["ALL"]
                    ],
                    "Max metabolite concentration differences (n=5) [M]": [
                        max_met_difference_dict["ALL"]
                    ],
                    "Max enzyme concentration differences (n=5) [mmol⋅gDW⁻¹]": [
                        max_enzyme_difference_dict["ALL"]
                    ],
                }
            )

        # Create XLSX
        create_cobrak_spreadsheet(
            f"{datafolder}spreadsheet_{glcvalue}{'_' + suffix if glcvalue == 'ALL' else ''}.xlsx",
            cobrak_model,
            {"ecTFVA": VariabilityDataset(variability_dict, with_df=True)},
            glcvalue_datasets,
            is_maximization=True,
            objective_overwrite="Biomass_fw"
            if suffix == "maxglc"
            else OBJECTIVE_VAR_NAME,
            extra_optstatistics_data=extra_optstatistics_data,
        )

        if "maxac" in suffix:
            dataname = "maxac" if suffix == "_maxac" else "maxac_minmu0.1"
            accombined_tfvas["ecTFVA - " + dataname] = VariabilityDataset(
                variability_dict, with_df=True
            )
            accombined_tfbas["NLP - " + dataname] = glcvalue_datasets[
                list(glcvalue_datasets.keys())[0]
            ]
            for key in extra_optstatistics_data:
                if key not in accombined_extrastatistics:
                    accombined_extrastatistics[key] = []
                accombined_extrastatistics[key].extend(extra_optstatistics_data[key])

    print("All done!")

# Create XLSX
accombined_folder = "./examples/iCH360/acetate_spreadsheet/"
ensure_folder_existence(accombined_folder)
create_cobrak_spreadsheet(
    f"{accombined_folder}acetate_spreadsheet.xlsx",
    cobrak_model,
    accombined_tfvas,
    accombined_tfbas,
    is_maximization=True,
    extra_optstatistics_data=accombined_extrastatistics,
)
