from statistics import mean  # noqa: D100

import z_add_path  # noqa: F401

from cobrak.io import get_files, json_load, standardize_folder

roundnumsdict: dict[str, list[int]] = {
    "Raven": [1, 2, 3, 4, 5],
    "MacBook": [6],
}
datafolder = standardize_folder("examples/iCH360/RESULTS_GLCUPTAKE")


for computer_id, roundnums in roundnumsdict.items():
    print(f"==={computer_id} runtime statistics===")

    ectfva_time_filenames = [
        datafolder + filename
        for filename in get_files(datafolder)
        if filename.startswith("variability_dict__")
        and filename.endswith("_time.json")
        and any(f"_{roundnum}_" in filename for roundnum in roundnums)
    ]
    evolution_time_filenames = [
        datafolder + filename
        for filename in get_files(datafolder)
        if filename.startswith("best_evolution_result_")
        and filename.endswith("_time.json")
        and any(f"_{roundnum}_" in filename for roundnum in roundnums)
    ]
    postprocessing_time_filenames = [
        datafolder + filename
        for filename in get_files(datafolder)
        if filename.startswith("final_best_result__")
        and filename.endswith("_time.json")
        and any(f"_{roundnum}_" in filename for roundnum in roundnums)
    ]

    ectfva_times = [
        json_load(ectfva_time_filename)[0] / 60
        for ectfva_time_filename in ectfva_time_filenames
    ]
    evolution_times = [
        json_load(evolution_time_filename)[0] / 60
        for evolution_time_filename in evolution_time_filenames
    ]
    postprocessing_times = [
        json_load(postprocessing_time_filename)[0] / 60
        for postprocessing_time_filename in postprocessing_time_filenames
    ]

    print(
        f"Min/Mean/max ecTFVA time [min]: {round(min(ectfva_times), 2)} / {round(mean(ectfva_times), 2)} / {round(max(ectfva_times), 2)}"
    )
    print(
        f"Min/Mean/max evolution time [min]: {round(min(evolution_times), 2)} / {round(mean(evolution_times), 2)} / {round(max(evolution_times), 2)}"
    )
    print(
        f"Min/Mean/max postprocessing time [min]: {round(min(postprocessing_times), 2)} / {round(mean(postprocessing_times), 2)} / {round(max(postprocessing_times), 2)}"
    )
