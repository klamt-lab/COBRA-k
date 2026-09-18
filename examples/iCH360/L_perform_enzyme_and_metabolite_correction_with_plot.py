# IMPORTS SECTION #  # noqa: D100
from math import exp, isnan, log
from statistics import mean

import matplotlib.pyplot as plt
import z_add_path  # noqa: F401
from scipy.stats import linregress

from cobrak.constants import (
    ALL_OK_KEY,
    ENZYME_VAR_INFIX,
    ENZYME_VAR_PREFIX,
    ERROR_SUM_VAR_ID,
    KAPPA_VAR_PREFIX,
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
from cobrak.io import ensure_folder_existence, json_load, standardize_folder
from cobrak.nlps import perform_nlp_irreversible_optimization_with_active_reacs_only
from cobrak.plotting import plot_range_bars, scatterplot_with_labels
from cobrak.standard_solvers import IPOPT_LONGRUN
from cobrak.utilities import (
    get_full_enzyme_id,
    get_full_enzyme_mw,
    get_reaction_enzyme_var_id,
    sort_dict_keys,
)


def analyze_metabolite_correction(  # noqa: D103
    model: Model,
    original_result: dict[str, float],
    corrected_result: dict[str, float],
    verbose: bool = False,
    txt_path: str = "",
) -> None:
    met_to_reacs: dict[str, list[tuple[str, float, float, float, float]]] = {}
    met_changes: dict[float, str] = {}
    for reac_id, reac in model.reactions.items():
        if reac_id not in original_result:
            continue
        if original_result[reac_id] < 1e-9:
            continue
        original_kappa = original_result.get(KAPPA_VAR_PREFIX + reac_id)
        corrected_kappa = corrected_result.get(KAPPA_VAR_PREFIX + reac_id)
        original_gamma = original_result.get(KAPPA_VAR_PREFIX + reac_id)
        corrected_gamma = corrected_result.get(KAPPA_VAR_PREFIX + reac_id)
        if None in (original_kappa, original_gamma, corrected_kappa, corrected_gamma):
            continue

        enzyme_id = get_reaction_enzyme_var_id(reac_id, reac)
        mw = get_full_enzyme_mw(cobrak_model=model, reaction=reac)
        for met_id in reac.stoichiometries:
            if met_id not in met_to_reacs:
                met_to_reacs[met_id] = []
            met_to_reacs[met_id].append(
                # (reac_id, original_kappa, corrected_kappa, original_gamma, corrected_gamma)
                (
                    reac_id,
                    original_kappa * original_gamma,
                    corrected_kappa * corrected_gamma,
                    mw * original_result[enzyme_id],
                    mw * corrected_result[enzyme_id],
                )
            )

            if met_id in met_changes:
                continue
            met_var_id = LNCONC_VAR_PREFIX + met_id
            if False in (met_var_id in original_result, met_var_id in corrected_result):
                continue
            met_change = abs(original_result[met_var_id] - corrected_result[met_var_id])
            if met_change in met_changes:
                while met_change in met_changes:
                    met_change += 1e-13
            met_changes[met_change] = (
                met_id,
                original_result[met_var_id],
                corrected_result[met_var_id],
            )
        met_changes = sort_dict_keys(met_changes, reverse=True)

    tested_ids: list[str] = []
    full_text = "[S]=Substrate; [P]=Product\n\n"
    for met_change, (met_id, original_conc, corrected_conc) in met_changes.items():
        if met_id in tested_ids:
            continue
        tested_ids.append(met_id)
        met_text = f"=={met_id}==\n"
        met_text += f"ln(Concentration change): {met_change:.2f}, from {exp(original_conc):.2g} M to {exp(corrected_conc):.2g} M\n"
        met_text += "Affected reactions with kinetic constraints:\n"
        for (
            reac_id,
            original_gk,
            corrected_gk,
            original_mw,
            corrected_mw,
        ) in met_to_reacs[met_id]:
            pos_string = (
                "[S]" if model.reactions[reac_id].stoichiometries[met_id] < 0 else "[P]"
            )
            met_text += f">{pos_string} {reac_id} (ΔG'°={model.reactions[reac_id].dG0:.2f} kJ⋅mol⁻¹) | κ⋅γ from {original_gk:.4f} to {corrected_gk:.4f} // pool usage from {original_mw:.6f} g⋅gDW⁻¹ ({(100 * original_mw / model.max_prot_pool):.4f} %) to {corrected_mw:.6f} g⋅gDW⁻¹ ({(100 * corrected_mw / model.max_prot_pool):.4f} %)\n"
        met_text += "\n"
        if verbose:
            print(met_text)
        full_text += met_text + "\n"
    if txt_path:
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(full_text)


def run_all_analyses_and_plotting(  # noqa: D103
    scenario: str,
) -> None:
    correction_folder: str = standardize_folder("examples/iCH360/correction_results/")
    opttarget_id = "Biomass_fw"
    match scenario:
        case "glucose":
            datafolder = standardize_folder("examples/iCH360/RESULTS_GLCUPTAKE/")
            suffix = "_glucose"
            runname = "_1_maxglc9.65"
            gerosa_name = "Glucose"
        case "glucose_no_manual_kcat_changes":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_GLCUPTAKE_NO_MANUAL_KCAT_CHANGES_DONEPROTADJ/"
            )
            suffix = "_glucose_no_manual_kcat_changes"
            runname = "_1_maxglc9.65"
            gerosa_name = "Glucose"
        case "glucose_glu5k_dG0_lowered":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_GLUANALYSIS_LOWER_GLU5K_DG0_DONEPROTADJ/"
            )
            suffix = "_glucose_glu5k_dG0_lowered"
            runname = "_1_maxglc9.65"
            gerosa_name = "Glucose"
        case "glucose_proline_intake":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_GLUANALYSIS_FREE_PROLINE_INTAKE_DONEPROTADJ/"
            )
            suffix = "_glucose_proline_intake"
            runname = "_1_maxglc9.65"
            gerosa_name = "Glucose"
        case "acetate_protadj":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_ACETATE_PROTADJ/"
            )
            suffix = "_acetate_protadj"
            runname = "_1_maxglc1000"
            gerosa_name = "Acetate"
        case "acetate":
            datafolder = standardize_folder("examples/iCH360/RESULTS_UPTAKE_ACETATE/")
            suffix = "_acetate"
            runname = "_1_maxglc1000"
            gerosa_name = "Acetate"
        case "fructose":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_FRUCTOSE_PROTADJ/"
            )
            suffix = "_fructose"
            runname = "_1_maxglc1000"
            gerosa_name = "Fructose"
        case "glycerol":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_GLYCEROL_PROTADJ/"
            )
            suffix = "_glycerol"
            runname = "_1_maxglc1000"
            gerosa_name = "Glycerol"
        case "gluconate":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_GLUCONATE_PROTADJ/"
            )
            suffix = "_gluconate"
            runname = "_1_maxglc1000"
            gerosa_name = "Gluconate"
        case "pyruvate":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_PYRUVATE_PROTADJ/"
            )
            suffix = "_pyruvate"
            runname = "_1_maxglc1000"
            gerosa_name = "Pyruvate"
        case "succinate":
            datafolder = standardize_folder(
                "examples/iCH360/RESULTS_UPTAKE_SUCCINATE_PROTADJ/"
            )
            suffix = "_succinate"
            runname = "_1_maxglc1000"
            gerosa_name = "Succinate"
        case "ectfva":
            datafolder = standardize_folder("examples/iCH360/ecTFBA_and_ecTFVA_results")
            suffix = "_ectfva"
            gerosa_name = "Glucose"

    # Commonly used data #
    if scenario != "ectfva":
        cobrak_model: Model = json_load(
            f"{datafolder}used_cobrak_model_{runname}.json",
            Model,
        )
        best_result = json_load(
            f"{datafolder}final_best_result_{runname}.json",
            OptResult,
        )
        variability_dict = json_load(
            f"{datafolder.replace('_PROTADJ', '')}variability_dict_{runname}.json",
            VarResult,
        )
        variability_dict[opttarget_id] = (
            best_result[opttarget_id] * (0.995 if scenario != "glycerol" else 0.99525),
            best_result[opttarget_id] * 1.0,
        )
    else:
        cobrak_model: Model = json_load(
            "examples/iCH360/prepared_external_resources/iCH360_cobrak.json",
            Model,
        )
        best_result = json_load(
            "examples/iCH360/ecTFBA_and_ecTFVA_results/ectfba_glc9.65.json",
            OptResult,
        )
        fva_result = json_load(
            "examples/iCH360/ecTFBA_and_ecTFVA_results/ectfva_glc1000.0_mu0.6934087526248266.json",
            VarResult,
        )

    if scenario.startswith(("glucose", "ectfva")):
        # Load Bennett metabolite data #
        in_vivo_met_data: dict[str, dict[str, float]] = json_load(
            "examples/common_needed_external_resources/Bennett_2009_full_data.json",
        )
    else:
        # Load Gerosa metabolite data #
        in_vivo_met_data: dict[str, dict[str, float]] = {}
        with open(
            "examples/iCH360/external_resources/Gerosa_et_al_concentrations.csv",
            encoding="utf-8",
        ) as f:
            csvlines = [x.split("\t") for x in f]
        conc_index = csvlines[0].index(
            gerosa_name
        )  # Index for concentrations of our scenario
        sd_index = conc_index + 1  # Index for standard deviation of our scenario
        for csvline in csvlines[1:]:
            mean_conc = float(csvline[conc_index]) * (
                330 / 1e6
            )  # from µmol⋅gDW⁻¹ to mol⋅l⁻¹ (through 330 g⁻¹)
            sd = abs(float(csvline[sd_index])) * (330 / 1e6)
            in_vivo_met_data[csvline[0]] = {
                "lb": mean_conc - sd,
                "ub": mean_conc + sd,
                "mean": mean_conc,
            }

    # Load Schmidt enzyme data #
    schmidt_data: dict[str, dict[str, float]] = json_load(
        "examples/common_needed_external_resources/Schmidt_2016_full_data.json",
    )
    enzyme_id_to_reac_id: dict[str, str] = {}
    for reac_id, reac_data in cobrak_model.reactions.items():
        if reac_data.enzyme_reaction_data is None:
            continue
        full_enz_id = get_full_enzyme_id(reac_data.enzyme_reaction_data.identifiers)
        if full_enz_id in schmidt_data:
            enzyme_id_to_reac_id[full_enz_id] = reac_id

    # Create folder #
    ensure_folder_existence(correction_folder)

    # Correction and plot functions #
    def run_metabolite_correction() -> dict[str, float]:  # noqa: D103
        met_correction_error_scenario: ErrorScenario = {
            f"{LNCONC_VAR_PREFIX}{met_id}": (
                log(in_vivo_met_data[met_id]["lb"]),
                log(in_vivo_met_data[met_id]["ub"]),
            )
            for met_id in in_vivo_met_data
        }

        if "_PROTADJ" in datafolder:
            cobrak_model.max_prot_pool = -round(best_result[OBJECTIVE_VAR_NAME], 4)

        metcorrection_result = (
            perform_nlp_irreversible_optimization_with_active_reacs_only(
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
        )
        print(
            metcorrection_result[ALL_OK_KEY],
            "at",
            metcorrection_result[opttarget_id],
            "1/h :",
            metcorrection_result[OBJECTIVE_VAR_NAME],
            "correction",
        )
        return metcorrection_result

    def plot_metabolite_result(  # noqa: D103
        nlp_result: dict[str, float],
        ax: plt.Axes,
        title: str,
        variabilities: dict[str, tuple[float, float]] = {},
        extracoords: bool = False,
    ) -> None:
        variability_nlp: list[tuple[float, float]] = []
        variability_bennett: list[tuple[float, float]] = []
        labels: list[str] = []
        data_nlp: list[float] = []
        data_bennett: list[float] = []
        for met_id in cobrak_model.metabolites:
            nlp_varname = f"{LNCONC_VAR_PREFIX}{met_id}"
            if (nlp_varname not in nlp_result) or (met_id not in in_vivo_met_data):
                continue
            if nlp_result[nlp_varname] is None:
                continue
            if isnan(in_vivo_met_data[met_id]["mean"]):
                continue
            labels.append((met_id + "\b").replace("_c\b", "").replace("__", "_"))
            if nlp_varname not in variabilities:
                variability_nlp.append(
                    (
                        exp(nlp_result[nlp_varname]),
                        exp(nlp_result[nlp_varname]),
                        exp(nlp_result[nlp_varname]),
                    )
                )
            else:
                variability_nlp.append(
                    (
                        exp(variabilities[nlp_varname][0]),
                        exp(variabilities[nlp_varname][1]),
                        exp(
                            mean(
                                (
                                    variabilities[nlp_varname][0],
                                    variabilities[nlp_varname][1],
                                )
                            )
                        ),
                    )
                )
            variability_bennett.append(
                (
                    in_vivo_met_data[met_id]["lb"],
                    in_vivo_met_data[met_id]["ub"],
                    in_vivo_met_data[met_id]["mean"],
                )
            )

            nlp_value = exp(nlp_result[nlp_varname])
            bennett_value = in_vivo_met_data[met_id]["mean"]

            data_nlp.append(log(nlp_value))
            data_bennett.append(log(bennett_value))

        _, _, r_value, _, _ = linregress(data_nlp, data_bennett)
        r_squared = r_value**2

        if scenario != "ectfva":
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
                extratext=f"R²={r_squared:.2f}",
                xlim_overwrite=(1.1e-7, 1.75e-1),
                ylim_overwrite=(1.1e-7, 1.75e-1),
                x_labelsize=8,
                y_labelsize=8,
                major_tick_labelsize=8,
                minor_tick_labelsize=8,
                title_labelsize=9,
                extratext_labelsize=8,
                label_fontsize=8,
                labelcoords=(-2.5, 10) if extracoords else (0, 10),
            )
        else:
            plot_range_bars(
                data_captions=[
                    "Feasible metabolite ranges under maximal growth (ecTFBA simulation)",
                    "Measurements with standard deviation",
                ],
                data_labels=labels,
                data_ranges=[
                    [(x[0], x[1]) for x in variability_nlp],
                    [(x[0], x[1]) for x in variability_bennett],
                ],
                data_colors=["red", "black"],
                ax=ax,
                cap_len=0,
                highlight_means=[False, True],
                log_y=True,
                title=title,
                ylabel="Concentration [M]",
                xlabel="Metabolite ID",
                legend_pos="upper left",
                marker_size=40,
                title_labelsize=9,
                axes_labelsize=8,
                ticks_labelsize=6.1,
                legend_labelsize=7.5,
                legend_bbox_to_anchor=(0.0, 0.995),
                ylim=(0.9e-6, 1e-0),
            )

    def plot_enzyme_result(
        nlp_result: dict[str, float],
        ax: plt.Axes,
        title: str,
        variabilities: dict[str, tuple[float, float]] = {},  # noqa: ARG001
    ) -> None:  # noqa: D103
        variability_nlp: list[tuple[float, float]] = []
        variability_schmidt: list[tuple[float, float]] = []
        labels: list[str] = []
        data_nlp: list[float] = []
        data_schmidt: list[float] = []
        for nlp_varname, nlp_value in nlp_result.items():
            if not nlp_varname.startswith(ENZYME_VAR_PREFIX):
                continue

            enzyme_name = nlp_varname[len(ENZYME_VAR_PREFIX) :].split(ENZYME_VAR_INFIX)[
                0
            ]
            if enzyme_name not in schmidt_data:
                continue
            if enzyme_name not in cobrak_model.enzymes:
                raise ValueError

            predicted_enzyme_conc = (
                nlp_result[nlp_varname]
                * cobrak_model.enzymes[enzyme_name].molecular_weight
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
            extratext=f"R²={r_squared:.2f}",
            x_labelsize=8,
            y_labelsize=8,
            major_tick_labelsize=8,
            minor_tick_labelsize=8,
            title_labelsize=9,
            extratext_labelsize=8,
            label_fontsize=8,
        )

    def plot_flux_result(  # noqa: D103
        predicted_data: dict[str, float],
        ax: plt.Axes,
        title: str,
        variabilities: dict[str, tuple[float, float]] = {},  # noqa: ARG001
    ) -> None:
        # Load Gerosa flux data #
        measured_data: dict[str, tuple[float, float, float]] = {}
        with open(
            "examples/iCH360/external_resources/Gerosa_et_al_fluxes.csv",
            encoding="utf-8",
        ) as f:
            csvlines = [x.split("\t") for x in f]
        flux_index = csvlines[0].index(
            gerosa_name
        )  # Index for concentrations of our scenario
        sd_index = flux_index + 1  # Index for standard deviation of our scenario
        for csvline in csvlines[1:]:
            csv_reac_id = csvline[1]
            flux_value = float(csvline[flux_index])
            if flux_value < 0.0:
                if csv_reac_id.endswith("_bw"):
                    csv_reac_id = csv_reac_id.replace("_bw", "_fw")
                else:
                    csv_reac_id = csv_reac_id.replace("_fw", "_bw")
            measured_data[csvline[0]] = (
                abs(float(csvline[sd_index])),
                abs(flux_value),
                csv_reac_id,
            )

        measured_datalist = []
        predicted_datalist = []
        standard_deviation_datalist = []
        labels = []
        for measured_id, (
            standard_deviation,
            measured_flux,
            predicted_id,
        ) in measured_data.items():
            predicted_flux = predicted_data.get(predicted_id, 0.0)

            if not (
                ((measured_flux == 0.0) and (predicted_flux <= 1e-7))
                or ((predicted_flux == 0.0) and (predicted_id == "N/A"))
            ):
                measured_datalist.append(abs(measured_flux))
                predicted_datalist.append(abs(predicted_flux))
                standard_deviation_datalist.append(abs(standard_deviation))
                labels.append(predicted_id)

        _, _, r_value, _, _ = linregress(measured_datalist, predicted_datalist)
        r_squared = r_value**2
        scatterplot_with_labels(
            [
                (
                    measured_datalist[i] - standard_deviation_datalist[i],
                    measured_datalist[i] + standard_deviation_datalist[i],
                    measured_datalist[i],
                )
                for i in range(len(measured_datalist))
            ],
            [(x, x, x) for x in predicted_datalist],
            labels,
            x_label="Measured fluxes [mmol⋅gDW⁻¹⋅h⁻¹]",
            y_label="Predicted fluxes [mmol⋅gDW⁻¹⋅h⁻¹]",
            y_log=False,
            x_log=False,
            add_labels=True,
            identical_axis_lims=False,
            ax=ax,
            title=title,
            extratext=f"R²={r_squared:.2f}",
            xlim_overwrite=(0.0, 17.0),
            ylim_overwrite=(0.0, 17.0),
            x_labelsize=8,
            y_labelsize=8,
            major_tick_labelsize=8,
            minor_tick_labelsize=8,
            title_labelsize=9,
            extratext_labelsize=8,
            label_fontsize=8,
        )

    if scenario != "ectfva":
        metcorrection_result = run_metabolite_correction()
        assert metcorrection_result[ALL_OK_KEY]
        analyze_metabolite_correction(
            cobrak_model,
            best_result,
            metcorrection_result,
            txt_path="examples/iCH360/correction_results/metabolite_correction_report_glucose.txt",
        )

    assert best_result[ALL_OK_KEY]

    if scenario.startswith("glucose"):
        bennett_result = json_load(
            "examples/iCH360/RESULTS_BENNETT/final_best_result__1_bennett_maxglc9.65.json"
        )
        assert bennett_result[ALL_OK_KEY]
        fig, axes = plt.subplots(3, 2, figsize=(7.25, 9.0))
        for counter, ax in enumerate(axes.flatten()):
            ax.set_box_aspect(1)
            match counter:
                case 0:
                    plot_flux_result(best_result, ax, "A")
                case 2:
                    plot_enzyme_result(best_result, ax, "B")
                case 1:
                    plot_metabolite_result(best_result, ax, "C", extracoords=True)
                case 3:
                    plot_metabolite_result(
                        metcorrection_result, ax, "D", extracoords=True
                    )
                case 5:
                    plot_metabolite_result(bennett_result, ax, "E", extracoords=True)
                case _:
                    pass
        fig.delaxes(axes[2, 0])
    elif scenario == "ectfva":
        fig, ax = plt.subplot_mosaic(
            [
                ["A", "B"],
                ["C", "C"],
            ],
            figsize=(7.25, 7.0),
            constrained_layout=True,
        )
        for key in ("A", "B"):
            ax[key].set_box_aspect(1)  # 1:1 width‑to‑height ratio
        plot_flux_result(best_result, ax["A"], "A")
        plot_enzyme_result(best_result, ax["B"], "B")
        plot_metabolite_result(best_result, ax["C"], "C", fva_result)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.0))
        for counter, ax in enumerate(axes.flatten()):
            ax.set_box_aspect(1)
            match counter:
                case 0:
                    plot_flux_result(best_result, ax, "A")
                case 1:
                    plot_metabolite_result(best_result, ax, "B")
                case 2:
                    plot_metabolite_result(metcorrection_result, ax, "C")
                case _:
                    pass
        fig.delaxes(axes[1, 1])

    plt.tight_layout()
    match scenario:
        case "glucose":
            plt.savefig(f"{correction_folder}Figure_3.pdf")
        case "ectfva":
            plt.savefig(f"{correction_folder}Figure_4.pdf")
        case "acetate":
            plt.savefig(f"{correction_folder}Figure_S1.pdf")
        case _:
            pass
    plt.savefig(f"{correction_folder}correction_plot{suffix}.png", dpi=300)
    plt.clf()


run_all_analyses_and_plotting(
    scenario="acetate",
)
run_all_analyses_and_plotting(
    scenario="ectfva",
)
run_all_analyses_and_plotting(
    scenario="glucose",
)
run_all_analyses_and_plotting(
    scenario="glucose_no_manual_kcat_changes",
)
exit(0)
run_all_analyses_and_plotting(
    scenario="glucose_glu5k_dG0_lowered",
)
run_all_analyses_and_plotting(
    scenario="glucose_proline_intake",
)

run_all_analyses_and_plotting(
    scenario="gluconate",
)
run_all_analyses_and_plotting(
    scenario="fructose",
)
run_all_analyses_and_plotting(
    scenario="pyruvate",
)
run_all_analyses_and_plotting(
    scenario="succinate",
)
run_all_analyses_and_plotting(
    scenario="gluconate",
)
run_all_analyses_and_plotting(
    scenario="glycerol",
)
