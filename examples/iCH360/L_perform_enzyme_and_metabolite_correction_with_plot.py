# IMPORTS SECTION #
from math import exp, log

import matplotlib.pyplot as plt
import z_add_path  # noqa: F401
from scipy.stats import linregress

from cobrak.constants import (
    ALL_OK_KEY,
    ENZYME_VAR_INFIX,
    ENZYME_VAR_PREFIX,
    ERROR_SUM_VAR_ID,
    LNCONC_VAR_PREFIX,
    OBJECTIVE_VAR_NAME,
)
from cobrak.dataclasses import (
    CorrectionConfig,
    ErrorScenario,
    Model,
    OptResult,
    VarResult,
)
from cobrak.io import ensure_folder_existence, json_load, json_write, standardize_folder
from cobrak.nlps import perform_nlp_irreversible_optimization_with_active_reacs_only
from cobrak.plotting import scatterplot_with_labels
from cobrak.standard_solvers import IPOPT_LONGRUN
from cobrak.utilities import (
    get_full_enzyme_id,
    get_full_enzyme_mw,
    get_reaction_enzyme_var_id,
)

# User variables #
opttarget_id = "Biomass_fw"
datafolder = standardize_folder("examples/iCH360/RESULTS_GLCUPTAKE/")
runname = "_1_maxglc9.65"
correction_folder = standardize_folder("examples/iCH360/correction_results/")

# Commonly used data #
cobrak_model: Model = json_load(
    f"{datafolder}used_cobrak_model_{runname}.json",
    Model,
)
best_result = json_load(
    f"{datafolder}final_best_result_{runname}.json",
    OptResult,
)
variability_dict = json_load(
    f"{datafolder}variability_dict_{runname}.json",
    VarResult,
)
variability_dict[opttarget_id] = (
    best_result[opttarget_id] * 0.995,
    best_result[opttarget_id] * 1.0,
)

# Load Bennett metabolite data #
bennett_data: dict[str, dict[str, float]] = json_load(
    "examples/common_needed_external_resources/Bennett_2009_full_data.json",
)

# Load Schmidt enzyme data #
schmidt_data: dict[str, dict[str, float]] = json_load(
    "examples/common_needed_external_resources/Schmidt_2016_full_data.json",
)
enzyme_id_to_reac_id: dict[str, str] = {}
enzyme_correction_error_scenario: ErrorScenario = {}
enzyme_extra_weights: dict[str, float] = {}
for reac_id, reac_data in cobrak_model.reactions.items():
    if reac_data.enzyme_reaction_data is None:
        continue
    full_enz_id = get_full_enzyme_id(reac_data.enzyme_reaction_data.identifiers)
    if full_enz_id in schmidt_data:
        cobrak_var_id = get_reaction_enzyme_var_id(reac_id, reac_data)
        enzyme_correction_error_scenario[cobrak_var_id] = (
            schmidt_data[full_enz_id],
            schmidt_data[full_enz_id],
        )
        enzyme_extra_weights[cobrak_var_id] = get_full_enzyme_mw(
            cobrak_model, reac_data
        )
        enzyme_id_to_reac_id[full_enz_id] = reac_id

# Create folder #
ensure_folder_existence(correction_folder)


# Correction and plot functions #
def run_metabolite_correction() -> dict[str, float]:
    met_correction_error_scenario: ErrorScenario = {
        f"{LNCONC_VAR_PREFIX}{met_id}": (
            log(bennett_data[met_id]["lb"]),
            log(bennett_data[met_id]["ub"]),
        )
        for met_id in bennett_data
    }

    metcorrection_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
        cobrak_model=cobrak_model,
        objective_target=ERROR_SUM_VAR_ID,
        objective_sense=-1,
        optimization_dict=best_result,
        variability_dict=variability_dict,
        correction_config=CorrectionConfig(
            error_scenario=met_correction_error_scenario,
            add_met_logconc_error_term=True,
            add_flux_error_term=False,
            add_dG0_error_term=False,
            add_km_error_term=False,
            add_kcat_times_e_error_term=False,
            error_sum_as_qp=False,
            use_weights=False,
        ),
        verbose=True,
        solver=IPOPT_LONGRUN,
    )
    print(
        metcorrection_result[ALL_OK_KEY],
        "at",
        metcorrection_result[opttarget_id],
        "1/h :",
        metcorrection_result[OBJECTIVE_VAR_NAME],
        "correction",
    )
    json_write(f"{correction_folder}met_correction_result.json", metcorrection_result)
    return metcorrection_result


def run_enzyme_correction() -> dict[str, float]:
    enzcorrection_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
        cobrak_model=cobrak_model,
        objective_target=ERROR_SUM_VAR_ID,
        objective_sense=-1,
        optimization_dict=best_result,
        variability_dict=variability_dict,
        correction_config=CorrectionConfig(
            error_scenario=enzyme_correction_error_scenario,
            add_met_logconc_error_term=False,
            add_enzyme_conc_error_term=True,
            add_flux_error_term=False,
            add_dG0_error_term=False,
            add_km_error_term=False,
            add_kcat_times_e_error_term=False,
            error_sum_as_qp=False,
            use_weights=False,
            extra_weights={},
            var_lb_ub_application="",
        ),
        verbose=True,
        solver=IPOPT_LONGRUN,
    )
    print(
        enzcorrection_result[ALL_OK_KEY],
        "at",
        enzcorrection_result[opttarget_id],
        "1/h :",
        enzcorrection_result[OBJECTIVE_VAR_NAME],
        "correction",
    )
    json_write(
        f"{correction_folder}enzcorrection_result.json",
        enzcorrection_result,
    )
    return enzcorrection_result


def plot_metabolite_result(
    nlp_result: dict[str, float], ax: plt.Axes, title: str
) -> None:
    variability_nlp: list[tuple[float, float]] = []
    variability_bennett: list[tuple[float, float]] = []
    labels: list[str] = []
    data_nlp: list[float] = []
    data_bennett: list[float] = []
    for met_id in cobrak_model.metabolites:
        nlp_varname = f"{LNCONC_VAR_PREFIX}{met_id}"
        if (nlp_varname not in nlp_result) or (met_id not in bennett_data):
            continue
        if nlp_result[nlp_varname] is None:
            continue
        labels.append(met_id)
        variability_nlp.append(
            (
                exp(nlp_result[nlp_varname]),
                exp(nlp_result[nlp_varname]),
                exp(nlp_result[nlp_varname]),
            )
        )
        variability_bennett.append(
            (
                bennett_data[met_id]["lb"],
                bennett_data[met_id]["ub"],
                bennett_data[met_id]["ub"],
                bennett_data[met_id]["mean"],
            )
        )

        nlp_value = exp(nlp_result[nlp_varname])
        bennett_value = bennett_data[met_id]["mean"]

        data_nlp.append(log(nlp_value))
        data_bennett.append(log(bennett_value))

    _, _, r_value, _, _ = linregress(data_nlp, data_bennett)
    r_squared = r_value**2

    scatterplot_with_labels(
        variability_bennett,
        variability_nlp,
        labels,
        x_label="Measured metabolite concentrations [M]",
        y_label="Predicted metabolite concentrations [M]",
        y_log=True,
        x_log=True,
        add_labels=True,
        ax=ax,
        title=title,
        extratext=f"R²={round(r_squared, 2)}",
        xlim_overwrite=(1.1e-7, 1.5e-1),
        ylim_overwrite=(1.1e-7, 1.5e-1),
    )


def plot_enzyme_result(nlp_result: dict[str, float], ax: plt.Axes, title: str) -> None:
    variability_nlp: list[tuple[float, float]] = []
    variability_schmidt: list[tuple[float, float]] = []
    labels: list[str] = []
    data_nlp: list[float] = []
    data_schmidt: list[float] = []
    for nlp_varname, nlp_value in nlp_result.items():
        if not nlp_varname.startswith(ENZYME_VAR_PREFIX):
            continue

        enzyme_name = nlp_varname[len(ENZYME_VAR_PREFIX) :].split(ENZYME_VAR_INFIX)[0]
        if enzyme_name not in schmidt_data:
            continue
        if enzyme_name not in cobrak_model.enzymes:
            raise ValueError

        predicted_enzyme_conc = (
            nlp_result[nlp_varname] * cobrak_model.enzymes[enzyme_name].molecular_weight
        )
        if predicted_enzyme_conc < 1e-8:  # Quasi deactivated
            continue
        if schmidt_data[enzyme_name] == 0.0:  # Not measured
            continue

        labels.append(enzyme_name)
        variability_nlp.append(
            (
                predicted_enzyme_conc,
                predicted_enzyme_conc,
                predicted_enzyme_conc,
            )
        )
        variability_schmidt.append(
            (
                schmidt_data[enzyme_name],
                schmidt_data[enzyme_name],
                schmidt_data[enzyme_name],
            )
        )
        data_nlp.append(log(predicted_enzyme_conc))
        data_schmidt.append(log(schmidt_data[enzyme_name]))

    _, _, r_value, _, _ = linregress(data_nlp, data_schmidt)
    r_squared = r_value**2

    scatterplot_with_labels(
        variability_schmidt,
        variability_nlp,
        labels,
        x_label="Measured enzyme abundances [g⋅gDW⁻¹]",
        y_label="Predicted enzyme abundances [g⋅gDW⁻¹]",
        y_log=True,
        x_log=True,
        add_labels=False,
        identical_axis_lims=False,
        xlim_overwrite=(1.1e-6, 1e-1),
        ylim_overwrite=(1.1e-6, 1e-1),
        ax=ax,
        title=title,
        extratext=f"R²={round(r_squared, 2)}",
    )


metcorrection_result = run_metabolite_correction()
enzcorrection_result = run_enzyme_correction()
bennett_result = json_load(
    "examples/iCH360/RESULTS_BENNETT/final_best_result__1_bennett_maxglc9.65.json"
)

assert best_result[ALL_OK_KEY]
assert metcorrection_result[ALL_OK_KEY]
assert bennett_result[ALL_OK_KEY]

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
for counter, ax in enumerate(axes.flatten()):
    match counter:
        case 0:
            plot_enzyme_result(best_result, ax, "A")
        case 1:
            plot_metabolite_result(best_result, ax, "B")
        case 2:
            plot_metabolite_result(metcorrection_result, ax, "C")
        case _:
            plot_metabolite_result(bennett_result, ax, "D")

plt.tight_layout()
plt.savefig(f"{correction_folder}correction_plot.png", dpi=500)
plt.show()
