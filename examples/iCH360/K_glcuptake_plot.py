import z_add_path  # noqa: D100, F401

from cobrak.io import get_files, json_load, standardize_folder
from cobrak.plotting import dual_axis_plot
from cobrak.spreadsheet_functionality import (
    OptimizationDataset,
)
from cobrak.utilities import sort_dict_keys

datafolder = standardize_folder("examples/iCH360/RESULTS_GLCUPTAKE")

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
        and ("_time" not in filename)
    ]
}
optimization_datasets = sort_dict_keys(optimization_datasets)

xaxis_caption = "Glucose uptake [mmol⋅gDW⁻¹⋅h⁻¹]"
optvalue_id = "Biomass_fw"
left_ids = ["Biomass_fw", "prot_pool_delivery"]
left_colors = ["black", "grey"]
left_captions = ["Growth", "Used enzyme pool"]
right_ids = ["EX_ac_e_fw", "EX_o2_e_bw", "CS_fw"]
right_colors = ["green", "red", "blue"]
right_captions = ["Acetate excretion", "Oxygen uptake", "TCA flux"]
leftaxis_title = "Growth [h⁻¹] / Enzyme pool [g⋅gDW⁻¹]"
rightaxis_title = "Flux [mmol⋅gDW⁻¹⋅h⁻¹]"
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
    # 7.05,
    # 7.1,
    # 7.15,
    # 7.2,
    # 7.25,
    # 7.3,
    # 7.33,
    # 7.35,
    # 7.4,
    # 7.45,
    # 7.5,
    # 7.55,
    # 7.6,
    # 7.66,
    # 7.7,
    # 7.75,
    # 7.8,
    # 7.85,
    # 7.9,
    # 7.95,
    8.0,
    8.33,
    8.66,
    9.0,
    9.33,
    9.65,
    10.0,
    1000,
)
x_id = "EX_glc__D_e_bw"
x_values = []
left_values = [[] for _ in range(len(left_ids))]
right_values = [[] for _ in range(len(right_ids))]
for glcvalue in glcvalues:
    print(f"Collecting {glcvalue}...")
    # Collect data
    glcvalue_datasets = {
        dataset_name: opt_dataset
        for dataset_name, opt_dataset in optimization_datasets.items()
        if (f"maxglc{glcvalue}" in dataset_name) or (glcvalue == "ALL")
    }

    opt_datasets = list(glcvalue_datasets.values())
    if len(opt_datasets) == 0:
        print(" INFO: Not yet calculated")
        continue
    all_objvalues = [opt_dataset.data[optvalue_id] for opt_dataset in opt_datasets]
    opt_opt_dataset = opt_datasets[all_objvalues.index(max(all_objvalues))].data
    for left_i, left_id in enumerate(left_ids):
        left_values[left_i].append(opt_opt_dataset.get(left_id, 0.0))
    for right_i, right_id in enumerate(right_ids):
        right_values[right_i].append(opt_opt_dataset.get(right_id, 0.0))
    x_values.append(opt_opt_dataset[x_id])

print("Create glucose uptake analysis plot...")

extrapoints = [
    (9.65, 6.82, False, "green", "x", r"Measured acetate excretion", 1.03),
    (9.654, 2.9780, False, "blue", "x", r"Measured TCA flux", 0.0906),
    (9.65, 0.65, True, "black", "x", r"Measured growth rate", 0.01),
]

dual_axis_plot(
    xpoints=x_values,
    leftaxis_ypoints_list=left_values,
    rightaxis_ypoints_list=right_values,
    xaxis_caption=xaxis_caption,
    leftaxis_caption=leftaxis_title,
    rightaxis_caption=rightaxis_title,
    leftaxis_colors=left_colors,
    rightaxis_colors=right_colors,
    leftaxis_titles=left_captions,
    rightaxis_titles=right_captions,
    extrapoints=extrapoints,
    has_legend=True,
    legend_direction="best",
    legend_position=(),
    is_leftaxis_logarithmic=False,
    is_rightaxis_logarithmic=False,
    point_style="o",
    line_style="-",
    max_digits_after_comma=1,
    savepath=f"{datafolder}zcomplete_figure_with_protpool.png",
    left_ylim=(0.0, 0.75),
    right_ylim=(0.0, 18.0),
    left_legend_position=[0, 1],
    right_legend_position=[0, 4, 2, 1, 3],
    figure_size_inches=(4, 3),
    special_figure_mode=True,
)

print("All done!")
