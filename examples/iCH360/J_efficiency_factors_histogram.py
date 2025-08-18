import matplotlib.pyplot as plt  # noqa: D100
import z_add_path  # noqa: F401

from cobrak.dataclasses import Model
from cobrak.io import ensure_folder_existence, json_load
from cobrak.utilities import get_df_and_efficiency_factors_sorted_lists

resultfolder = "RESULTS_GLCUPTAKE"
runname = "1_maxglc9.65"

cobrak_model = json_load(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak.json", Model
)
result = json_load(f"examples/iCH360/{resultfolder}/final_best_result__{runname}.json")

_, kappas, gammas, _, _, kappa_times_gamma = get_df_and_efficiency_factors_sorted_lists(
    cobrak_model,
    result,
    min_flux=1e-6,
)

# Create a figure and a set of subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # sharex=True; sharey=True

# Iterate over the values and create subplots
for i, value in enumerate(("γ", "κ", "κ⋅γ")):
    match value:
        case "γ":
            data = gammas.values()
        case "κ":
            data = kappas.values()
        case _:
            data = [x[0] for x in kappa_times_gamma.values()]

    # Create the histogram in the corresponding subplot
    axs[i].bar(list(range(len(data))), data, width=2.0)

    match i:
        case 0:
            title = "A"
        case 1:
            title = "B"
        case _:
            title = "C"
    # Add title and labels to the subplot
    axs[i].set_title(title, loc="left", fontweight="bold", fontsize=20)
    axs[i].set_xlabel("Reaction number", fontsize=17)
    axs[i].set_ylabel(f"{value} value", fontsize=17)
    axs[i].set_ylim(0.0, 1.0)
    axs[i].set_xlim(0, len(data))

    axs[i].tick_params(axis="both", which="major", labelsize=16)
    axs[i].tick_params(axis="both", which="minor", labelsize=13)

    # Calculate and add the mean line
    # mean_value = sum(data) / len(data)
    # axs[i].axhline(y=mean_value, color='red', linestyle='--', label=f'Mean = {mean_value:.2f}')
    # axs[i].axhline(y=median(data), color='orange', linestyle='--', label=f'Median = {median(data):.2f}')
    # axs[i].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure
ensure_folder_existence("examples/iCH360/efficiency_factors_data")
plt.savefig(
    "examples/iCH360/efficiency_factors_data/efficiency_factors_histogram.png", dpi=500
)
