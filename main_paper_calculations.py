import contextlib
import json
import os
import subprocess
import tempfile
from dataclasses import asdict, dataclass, is_dataclass
from sys import argv
from typing import Any

from pydantic import Field, TypeAdapter

if argv[-1] == "local":
    print("LOCAL MODE")
    print("In the 'main_paper_calculations_jsons' subfolder,")
    print("the .jsons for all main paper")
    print("are created. From the COBRA-k main folder,")
    print("run examples/iCH360/H_run_calculations.py")
    print("with a JSON path to recacalculate,")
    print(
        "e.g. uv run examples/iCH360/H_run_calculations.py ./main_paper_calculations_jsons/XXX.json"
    )


def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    try:
        os.makedirs(folder)
    except FileExistsError:
        pass


def json_load(path: str, dataclass_type: Any = Any) -> Any:
    """Load JSON data from a file and validate it against a specified dataclass type.

    This function reads the content of a JSON file located at the given `path`, parses it,
    and validates the parsed data against the provided `dataclass_type`. If the data is valid
    according to the dataclass schema, it returns an instance of the dataclass populated with
    the data. Otherwise, it raises an exception.

    Parameters:
    ----------
    path : str
        The file path to the JSON file that needs to be loaded.

    dataclass_type : Type[T]
        A dataclass type against which the JSON data should be validated and deserialized.

    Returns:
    -------
    T
        An instance of the specified `dataclass_type` populated with the data from the JSON file.

    Raises:
    ------
    JSONDecodeError
        If the content of the file is not a valid JSON string.

    ValidationError
        If the parsed JSON data does not conform to the schema defined by `dataclass_type`.

    Examples:
    --------
    >>> @dataclass
    ... class Person:
    ...     name: str
    ...     age: int

    >>> person = json_load('person.json', Person)
    >>> print(person.name, person.age)
    John Doe 30
    """
    with open(path, encoding="utf-8") as f:
        data = f.read()

    return TypeAdapter(dataclass_type).validate_json(data)


def json_write(path: str, json_data: Any) -> None:
    """Writes a JSON file at the given path with the given data as content.

    Can be also used for any of COBRAk's dataclasses as well as any
    dictionary of the form dict[str, dict[str, T] | None] where
    T stands for a COBRAk dataclass or any other JSON-compatible
    object type.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * json_data: Any ~ The dictionary or list which shalll be the content of
      the created JSON file
    """
    if is_dataclass(json_data):
        json_write(path, asdict(json_data))
    elif isinstance(json_data, dict) and sum(
        is_dataclass(value) for value in json_data.values()
    ):
        json_dict: dict[str, dict[str, Any] | None] = {}
        for key, data in json_data.items():
            if data is None:
                json_dict[key] = None
            elif is_dataclass(data):
                json_dict[key] = asdict(data)
            else:
                json_dict[key] = data
        json_write(path, json_dict)
    else:
        json_output = json.dumps(json_data, indent=4)
        with open(path, "w+", encoding="utf-8") as f:
            f.write(json_output)


@dataclass
class RunConfig:
    # Model changes
    manually_changed_kms: dict[str, dict[str, float]]
    manually_changed_kcats: dict[str, dict[str, float]]
    # Folder settings
    results_folder: str
    # ecTFVA settings
    ectfva_active_reacs: list[str]
    # Evolutionary algorithm settings
    round_num: int
    objective_target: str | dict[str, float]
    objective_sense: int
    deactivated_reacs: list[str]
    set_bounds: dict[str, tuple[float, float]]
    working_results: list[dict[str, float]]
    sampling_rounds_per_metaround: int = 2
    sampling_wished_num_feasible_starts: int = 5
    sampling_max_metarounds: int = 3
    evolution_num_gens: int = 150
    pop_size: int = 32
    protein_pool: float | None = None
    uses_bennett_concs: bool = False
    max_conc_sum: float | None = None
    nameaddition: str | None = None
    kicked_reacs: list[str] = Field(default_factory=list)


def create_and_submit_slurm_job(json_path: str) -> None:
    jobname = "PSBR"
    # Define the SLURM script content
    slurm_script_content = f"""#!/bin/bash

#SBATCH -J {jobname} # Your job name
#SBATCH -e {jobname}%j.err
#SBATCH -o {jobname}%j.out
#SBATCH --time=0-10:45:00 # Maximum expected runtime.
#SBATCH --nodes=1                 # Request 1 full node
#SBATCH --ntasks=1                 # Allocate 1 task
#SBATCH --cpus-per-task=72         # Number CPUs per task
#SBATCH --partition general           # Choose Partition (Queue)
#SBATCH --mail-type=FAIL,END       # An email is sent on begin, end, and failure of the job
#SBATCH --mail-user=bekiaris@mpi-magdeburg.mpg.de # E-Mail for notification

export OMP_NUM_THREADS=1

bash -i run_cobrak_python.sh ./examples/iCH360/H_run_calculations.py {json_path}
"""

    # Create a temporary file to store the SLURM script
    with tempfile.NamedTemporaryFile(
        encoding="utf-8", mode="w", delete=False, suffix=".slurm", dir=os.getcwd()
    ) as temp_file:
        temp_file.write(slurm_script_content)
        temp_file_path = temp_file.name

    # Make the SLURM script executable
    os.chmod(temp_file_path, 0o755)

    # Submit the SLURM job using sbatch without blocking the Python script
    subprocess.Popen(["sbatch", temp_file_path])

    print(f"SLURM job script created and submitted: {temp_file_path}")


if __name__ == "__main__":
    if argv[-1] != "local":
        os.chdir("/u/pbekiaris/")

    MANUALLY_CHANGED_KMS = {}
    MANUALLY_CHANGED_KCATS = {
        "NADTRHD_fw": 706056.6952036729 / 10,
        "ME1_fw": 1181225.0119759631 / 100,
    }

    run_configs = []

    for round_num in (1, 2, 3, 4, 5):
        """
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS_8,
                    results_folder="/examples/iCH360/RESULTS_PROTPOOLTESTS/",
                    ectfva_active_reacs=["Biomass_fw"],
                    round_num=round_num,
                    objective_target="Biomass_fw",
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "EX_glc__D_e_bw": (0.0, maxglcvalue[0]),
                        # "Biomass_fw": (0.7, 1.0),
                    },
                    protein_pool=maxglcvalue[1],
                    max_conc_sum=0.4,
                    nameaddition=f"maxglc{maxglcvalue[0]}_protpool{maxglcvalue[1]}_CALIBRATION",
                )
                for maxglcvalue in (
                    (1000.0, 0.221),
                    (1000.0, 0.222),
                    (1000.0, 0.223),
                    (1000.0, 0.224),
                    (1000.0, 0.225),
                    (1000.0, 0.226),
                    (1000.0, 0.227),
                    (1000.0, 0.228),
                    (1000.0, 0.229),
                    (1000.0, 0.230),
                    (1000.0, 0.231),
                    (1000.0, 0.232),
                    (1000.0, 0.233),
                    (1000.0, 0.234),
                    (1000.0, 0.235),
                )
            ]
        )
        """
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS,
                    results_folder="/examples/iCH360/RESULTS_GLCUPTAKE/",
                    ectfva_active_reacs=["Biomass_fw"],
                    round_num=round_num,
                    objective_target={
                        "Biomass_fw": 1.0,
                        "prot_pool_delivery": -0.01,
                    },
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "EX_glc__D_e_bw": (0.0, maxglcvalue),
                    },
                    protein_pool=0.224,
                    max_conc_sum=0.4,
                    nameaddition=f"maxglc{maxglcvalue}",
                    uses_bennett_concs=False,
                    kicked_reacs=[],
                )
                for maxglcvalue in (
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
                )
            ]
        )
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS,
                    results_folder="/examples/iCH360/RESULTS_BENNETT/",
                    ectfva_active_reacs=["Biomass_fw"],
                    round_num=round_num,
                    objective_target={
                        "Biomass_fw": 1.0,
                        "prot_pool_delivery": -0.01,
                    },
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "EX_glc__D_e_bw": (0.0, maxglcvalue),
                    },
                    protein_pool=0.224,
                    max_conc_sum=0.4,
                    nameaddition=f"bennett_maxglc{maxglcvalue}",
                    uses_bennett_concs=True,
                    kicked_reacs=[],
                )
                for maxglcvalue in (9.65,)
            ]
        )
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS,
                    results_folder="/examples/iCH360/RESULTS_MAXAC_MINMU/",
                    ectfva_active_reacs=["Biomass_fw", "EX_ac_e_fw"],
                    round_num=round_num,
                    objective_target={
                        "EX_ac_e_fw": 1.0,
                        "prot_pool_delivery": -0.01,
                    },
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "Biomass_fw": (0.1, 2.0),
                        "EX_glc__D_e_bw": (0.0, 1_000),
                    },
                    protein_pool=0.224,
                    max_conc_sum=0.4,
                    nameaddition="_maxac_minmu0.1",
                    uses_bennett_concs=False,
                    kicked_reacs=[],
                ),
            ]
        )
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS,
                    results_folder="/examples/iCH360/RESULTS_MAXAC/",
                    ectfva_active_reacs=["EX_ac_e_fw"],
                    round_num=round_num,
                    objective_target={
                        "EX_ac_e_fw": 1.0,
                        "prot_pool_delivery": -0.01,
                    },
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "EX_glc__D_e_bw": (0.0, 1_000),
                    },
                    protein_pool=0.224,
                    max_conc_sum=0.4,
                    nameaddition="_maxac",
                    uses_bennett_concs=False,
                    kicked_reacs=[],
                )
            ]
        )
        minmu_ac_ectfva_data = json_load(
            "AutoCOBRAK/examples/iCH360/prepared_external_resources/variability_dict_enforced_growth_and_ac.json"
            if argv[-1] != "local"
            else "./examples/iCH360/prepared_external_resources/variability_dict_enforced_growth_and_ac.json"
        )
        cobrak_model = json_load(
            "AutoCOBRAK/examples/iCH360/prepared_external_resources/iCH360_cobrak.json"
            if argv[-1] != "local"
            else "./examples/iCH360/prepared_external_resources/iCH360_cobrak.json"
        )
        raw_ko_targets = [
            reac_id
            for reac_id, (min_flux, max_flux) in minmu_ac_ectfva_data.items()
            if (min_flux <= 0.0)
            and (max_flux != 0.0)
            and (reac_id in cobrak_model["reactions"])
            and (
                "EX_" not in reac_id
            )  # ("PPS_fw" in reac_id)  # ("pp_" not in reac_id) and ("tex" not in reac_id) and ("EX_" not in reac_id)
        ]
        found_ko_targets = []
        ko_targets = []
        for raw_ko_target in raw_ko_targets:
            if raw_ko_target in found_ko_targets:
                continue

            if raw_ko_target.endswith("_fw"):
                tried_other_id = (raw_ko_target + "\b").replace("_fw\b", "_bw")
            else:
                tried_other_id = (raw_ko_target + "\b").replace("_bw\b", "_fw")
            ko_target_1 = raw_ko_target

            if (tried_other_id in cobrak_model["reactions"]) and (
                tried_other_id not in found_ko_targets
            ):
                ko_target_2 = tried_other_id
            else:
                ko_target_2 = ""
            ko_targets.append((ko_target_1, ko_target_2))
            found_ko_targets.extend((ko_target_1, ko_target_2))
        run_configs.extend(
            [
                RunConfig(
                    manually_changed_kms=MANUALLY_CHANGED_KMS,
                    manually_changed_kcats=MANUALLY_CHANGED_KCATS,
                    results_folder="/examples/iCH360/RESULTS_SINGLES_KOS/",
                    ectfva_active_reacs=["Biomass_fw", "EX_ac_e_fw"],
                    round_num=round_num,
                    objective_target={
                        "Biomass_fw": 1.0,
                        "prot_pool_delivery": -0.01,
                    },
                    objective_sense=+1,
                    deactivated_reacs=[],
                    evolution_num_gens=100,
                    pop_size=64,
                    working_results=[],
                    set_bounds={
                        "EX_glc__D_e_bw": (0.0, 1_000),
                    },
                    kicked_reacs=[
                        ko_target[0],
                        ko_target[1],
                    ],
                    protein_pool=0.224,
                    max_conc_sum=0.4,
                    nameaddition=f"_ko{ko_target[0]}_{ko_target[1]}",
                    uses_bennett_concs=False,
                )
                for ko_target in ko_targets
            ]
        )

    if argv[-1] == "local":
        ensure_folder_existence("./main_paper_calculations_jsons")
        for run_config in run_configs:
            json_path = (
                "./main_paper_calculations_jsons/"
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}.json"
            )
            json_write(json_path, run_config)
    else:
        for run_config in run_configs:
            with contextlib.suppress(FileExistsError):
                os.makedirs("./AutoCOBRAK" + run_config.results_folder)
            json_path = (
                "./AutoCOBRAK"
                + run_config.results_folder
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}.json"
            )
            json_write(json_path, run_config)
            create_and_submit_slurm_job(
                json_path="."
                + run_config.results_folder
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}.json",
            )
