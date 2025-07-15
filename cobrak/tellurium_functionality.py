# IMPORTS SECION #
from numpy import exp

from cobrak.constants import (
    DF_VAR_PREFIX,
    GAMMA_VAR_PREFIX,
    KAPPA_VAR_PREFIX,
    LNCONC_VAR_PREFIX,
    STANDARD_R,
    STANDARD_T,
)
from cobrak.dataclasses import Model, Reaction
from cobrak.utilities import (
    delete_orphaned_metabolites_and_enzymes,
    get_unoptimized_reactions_in_nlp_solution,
    have_all_unignored_km,
)


# "PRIVATE" FUNCTIONS SECTION #
def _get_numbersafe_id(met_id: str) -> str:
    match met_id.startswith(("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")):
        case True:
            return f"x{met_id}"
        case False:
            return met_id


def _get_reaction_string_of_cobrak_reaction(
    cobrak_model: Model,
    reac_id: str,
    cobrak_reaction: Reaction,
    e_conc: float,
    met_concs: dict[str, float],
    reac_flux: float,
    nlp_results: dict[str, float],
    kinetic_ignored_metabolites: list[str],
    unoptimized_reactions: dict[str, tuple[float, float]],
) -> str:
    reac_id = _get_numbersafe_id(reac_id)

    has_vplus = (
        (cobrak_reaction.enzyme_reaction_data is not None)
        and (cobrak_reaction.enzyme_reaction_data.k_cat is not None)
        and (cobrak_reaction.enzyme_reaction_data.identifiers != [""])
    )
    has_kappa = has_vplus and have_all_unignored_km(
        cobrak_reaction, kinetic_ignored_metabolites
    )
    has_gamma = cobrak_reaction.dG0 is not None

    # Reaction stoichiometry
    pow_stoich_substrate_strings = []
    pow_stoich_d_km_substrate_strings = []
    pow_stoich_product_strings = []
    pow_stoich_d_km_product_strings = []
    stoich_times_substrate_strings = []
    stoich_times_product_strings = []
    substrates = []
    products = []
    for unsafe_met_id, stoichiometry in cobrak_reaction.stoichiometries.items():
        met_id = _get_numbersafe_id(unsafe_met_id)
        stoich_times_met_string = str(abs(stoichiometry)) + " " + met_id
        pow_stoich_met_string = (
            "(" + met_id + "_molar" + "^" + str(abs(stoichiometry)) + ")"
        )
        pow_stoich_d_km_met_string = (
            "(("
            + met_id
            + "_molar"
            + f" / km_{met_id}_{reac_id}"
            + ")"
            + "^"
            + str(abs(stoichiometry))
            + ")"
        )
        if met_id in kinetic_ignored_metabolites:
            pow_stoich_d_km_met_string = "1.0"
        if stoichiometry < 0:
            stoich_times_substrate_strings.append(stoich_times_met_string)
            pow_stoich_substrate_strings.append(pow_stoich_met_string)
            pow_stoich_d_km_substrate_strings.append(pow_stoich_d_km_met_string)
            substrates.append(met_id)
        else:
            stoich_times_product_strings.append(stoich_times_met_string)
            pow_stoich_product_strings.append(pow_stoich_met_string)
            pow_stoich_d_km_product_strings.append(pow_stoich_d_km_met_string)
            products.append(met_id)
    reac_string = f"{reac_id}: {' + '.join(stoich_times_substrate_strings)} -> {' + '.join(stoich_times_product_strings)};\n"

    kappa_substrates = "(" + " * ".join(pow_stoich_d_km_substrate_strings) + ")"
    gamma_substrates = "(" + " * ".join(pow_stoich_substrate_strings) + ")"
    kappa_products = "(" + " * ".join(pow_stoich_d_km_product_strings) + ")"
    gamma_products = "(" + " * ".join(pow_stoich_product_strings) + ")"

    if reac_id in unoptimized_reactions:
        nlp_flux = unoptimized_reactions[reac_id][0]
        real_flux = unoptimized_reactions[reac_id][1]
        e_conc *= nlp_flux / real_flux

    reac_variables = f"const E_conc_{reac_id} = {e_conc}\n"

    if (not has_vplus) or (not has_kappa and not has_gamma):
        substrate_mult_strlist: list[str] = []
        product_mult_strlist: list[str] = []

        total_mult = 1.0
        for met_id, stoichiometry in cobrak_reaction.stoichiometries.items():
            met_var_id = LNCONC_VAR_PREFIX + met_id
            met_var_id_2 = LNCONC_VAR_PREFIX + met_id[1:]
            if met_id in met_concs:
                used_conc = met_concs[met_id]
            elif met_var_id in nlp_results:
                used_conc = exp(nlp_results[met_var_id])
            elif met_var_id_2 in nlp_results:
                used_conc = exp(nlp_results[met_var_id_2])
            else:
                used_conc = exp(cobrak_model.metabolites[met_id].log_min_conc)

            total_mult *= used_conc**stoichiometry
            if stoichiometry < 0.0:
                substrate_mult_strlist.append(
                    f"({_get_numbersafe_id(met_id)}_molar ^ {abs(stoichiometry)})"
                )
            else:
                product_mult_strlist.append(
                    f"({_get_numbersafe_id(met_id)}_molar ^ {stoichiometry})"
                )

        k_eq = 1.1 * total_mult

        real_flux = 1 - (total_mult / k_eq)
        kappa_mult_str = (
            " * ".join(product_mult_strlist)
            + " / ("
            + " * ".join(substrate_mult_strlist)
            + ") / "
            + str(k_eq)
        )

        reac_assignments = f"{KAPPA_VAR_PREFIX}{reac_id} := {kappa_mult_str}\n"
        menten_multiplier = reac_flux / real_flux
        reac_kinetic = f"{reac_id} = {menten_multiplier} * (1 - {kappa_mult_str})\n"
    else:
        if cobrak_reaction.enzyme_reaction_data is None:
            raise ValueError

        flux_value = 1.0

        # V plus
        reac_variables += (
            f"const k_cat_{reac_id} = {cobrak_reaction.enzyme_reaction_data.k_cat}\n"
        )
        v_plus = f"E_conc_{reac_id} * k_cat_{reac_id}"
        reac_assignments = f"v_plus_{reac_id} := {v_plus}\n"
        flux_value *= cobrak_reaction.enzyme_reaction_data.k_cat * e_conc

        if has_kappa:
            kappa_substrates_value = 1.0
            kappa_products_value = 1.0
            for unsafe_met_id, km in cobrak_reaction.enzyme_reaction_data.k_ms.items():
                original_met_id = unsafe_met_id
                met_id = _get_numbersafe_id(unsafe_met_id)
                reac_variables += f"const km_{met_id}_{reac_id} = {km}\n"

                stoichiometry = cobrak_reaction.stoichiometries[original_met_id]
                kappa_mult = (
                    exp(nlp_results[LNCONC_VAR_PREFIX + original_met_id]) / km
                ) ** abs(stoichiometry)
                if stoichiometry < 0:
                    kappa_substrates_value *= kappa_mult
                else:
                    kappa_products_value *= kappa_mult
            flux_value *= kappa_substrates_value / (
                1 + kappa_substrates_value + kappa_products_value
            )

            kappa = f"{kappa_substrates} / (1 + {kappa_substrates} + {kappa_products})"
            reac_assignments += f"{KAPPA_VAR_PREFIX}{reac_id} := {kappa}\n"
        else:
            kappa = "1"

        if has_gamma:
            used_dG0 = cobrak_reaction.dG0
            if used_dG0 is None:
                raise ValueError

            k_eq = exp(-used_dG0 / (cobrak_model.R * cobrak_model.T))

            reac_variables += f"const dG0_{reac_id} = {used_dG0}\n"
            reac_variables += f"const k_eq_{reac_id} = {k_eq}\n"

            dg = f"(dG0_{reac_id} + R * T * log( ({gamma_products}) ) - R * T * log( ({gamma_substrates}) ))"
            gamma = f"(1 - ({gamma_products} / {gamma_substrates} / k_eq_{reac_id}))"
            reac_assignments += f"{DF_VAR_PREFIX}{reac_id} := -{dg}\n"
            reac_assignments += f"{GAMMA_VAR_PREFIX}{reac_id} := {gamma}\n"

            gamma_products_mult = 1.0
            gamma_substrates_mult = 1.0
            for met_id, stoichiometry in cobrak_reaction.stoichiometries.items():
                if stoichiometry < 0.0:
                    gamma_substrates_mult *= exp(
                        nlp_results[LNCONC_VAR_PREFIX + met_id]
                    ) ** abs(stoichiometry)
                else:
                    gamma_products_mult *= (
                        exp(nlp_results[LNCONC_VAR_PREFIX + met_id]) ** stoichiometry
                    )
            flux_value *= 1 - (gamma_products_mult / gamma_substrates_mult / k_eq)
        else:
            gamma = "1"

        reac_kinetic = f"{reac_id} = ({reac_flux / flux_value} * {v_plus})  *  ({kappa})  *  ({gamma});\n"

    reac_string += reac_kinetic + reac_variables + reac_assignments

    return reac_string


# "PUBLIC" FUNCTIONS SECTION #
def get_tellurium_string_from_cobrak_model(
    cobrak_model: Model,
    cell_density: float,
    e_concs: dict[str, float],
    met_concs: dict[str, float],
    nlp_results: dict[str, float],
) -> str:
    unoptimized_reactions = get_unoptimized_reactions_in_nlp_solution(
        cobrak_model,
        nlp_results,
    )
    tellurium_string = (
        "# General constants\n" + f"R = {STANDARD_R}\n" + f"T = {STANDARD_T}\n"
    )
    for reac_id, cobrak_reaction in cobrak_model.reactions.items():
        if reac_id.endswith(cobrak_model.fwd_suffix):
            reverse_id = reac_id.replace(
                cobrak_model.fwd_suffix, cobrak_model.rev_suffix
            )
        elif reac_id.endswith(cobrak_model.rev_suffix):
            reverse_id = reac_id.replace(
                cobrak_model.rev_suffix, cobrak_model.fwd_suffix
            )
        else:
            reverse_id = ""

        if reverse_id in nlp_results:
            reac_flux = nlp_results[reac_id] - nlp_results[reverse_id]
        else:
            reac_flux = nlp_results[reac_id]

        if reac_flux <= abs(1e-12):
            continue

        e_conc = e_concs.get(reac_id, 1.0)

        tellurium_string += _get_reaction_string_of_cobrak_reaction(
            cobrak_model=cobrak_model,
            reac_id=reac_id,
            cobrak_reaction=cobrak_reaction,
            e_conc=e_conc,
            met_concs=met_concs,
            reac_flux=reac_flux,
            nlp_results=nlp_results,
            kinetic_ignored_metabolites=cobrak_model.kinetic_ignored_metabolites,
            unoptimized_reactions=unoptimized_reactions,
        )

    cobrak_model = delete_orphaned_metabolites_and_enzymes(cobrak_model)

    for unsafe_met_id, metabolite in cobrak_model.metabolites.items():
        original_met_id = unsafe_met_id
        met_id = _get_numbersafe_id(unsafe_met_id)
        if met_id in met_concs:
            tellurium_string += (
                f"\nconst substanceOnly species {met_id} = {met_concs[met_id] * 1_000 / cell_density}"
                f"\n{met_id}_molar := {met_id} * {cell_density / 1_000}\n"
            )
        else:
            if ("x_" + original_met_id in nlp_results) and (
                metabolite.log_min_conc != metabolite.log_max_conc
            ):
                exp_conc = (
                    exp(nlp_results[LNCONC_VAR_PREFIX + original_met_id])
                    * 1_000
                    / cell_density
                )
                tellurium_string += (
                    f"\nsubstanceOnly species {met_id} = {exp_conc}"
                    f"\n{met_id}_molar := {met_id} * {cell_density / 1_000}"
                )
            else:
                prefix = "const " if original_met_id.endswith("_e") else ""
                tellurium_string += (
                    f"\n{prefix}substanceOnly species {met_id} = {exp(metabolite.log_min_conc) * 1_000 / cell_density}"
                    f"\n{met_id}_molar := {met_id} * {cell_density / 1_000}"
                )

    return tellurium_string
