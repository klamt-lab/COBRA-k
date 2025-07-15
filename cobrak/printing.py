"""Pretty-print summaries of optimization and variability results as well as COBRAk Model instances.

For results, its methods generate ```rich``` tables that display flux values and variability information for each category.
For models, its methods generate ```rich```tables that display the model's structure and parameters.
"""

from collections.abc import Callable
from json import dumps
from math import exp
from typing import Any

from rich.columns import Columns
from rich.table import Table

from . import console
from .constants import (
    ALL_OK_KEY,
    DF_VAR_PREFIX,
    ERROR_SUM_VAR_ID,
    ERROR_VAR_PREFIX,
    GAMMA_VAR_PREFIX,
    KAPPA_VAR_PREFIX,
    LNCONC_VAR_PREFIX,
    OBJECTIVE_VAR_NAME,
)
from .dataclasses import Model
from .utilities import (
    get_enzyme_usage_by_protein_pool_fraction,
    get_extra_linear_constraint_string,
    get_metabolite_consumption_and_production,
    get_reaction_enzyme_var_id,
    get_reaction_string,
    get_substrate_and_product_exchanges,
    sort_dict_keys,
)


def _mapcolored(
    value: int | float,
    min_value: int | float,
    max_value: int | float,
    special_value: float = 0.0,
    prefix: str = "",
    suffix: str = "",
) -> str:
    """Print the value of x with a color that depends on its proximity to the minimum and maximum values.

    Args:
        x (float): The value to print.
        min_value (float): The minimum value.
        max_value (float): The maximum value.
    """
    # Map the value to a color between blue and red
    if value != special_value:
        ratio = (value - min_value) / (max_value - min_value)
        r = int(255 * ratio)
        g = 0
        b = int(255 * (1 - ratio))
    else:
        r = g = b = 255

    # Print the value with the corresponding color
    return f"[no highlight][rgb({r},{g},{b})]{prefix}{value}{suffix}[/rgb({r},{g},{b})][/no highlight]"


def _varcolor(key: str, variability: dict[str, tuple[float, float]]) -> tuple[str, str]:
    """Determine the color prefix and suffix for a variable based on its variability.

    Args:
    key (str): The variable key.
    variability (dict[str, tuple[float, float]]): A dictionary containing variability data.

    Returns:
    tuple[str, str]: A tuple containing the color prefix and suffix.
    """
    value_min = variability[key][0]
    value_max = variability[key][1]
    prefix, suffix = _zero_prefix(value_max), _zero_suffix(value_max)
    if value_max == 0.0:
        return f"[no highlight][white]{prefix}", f"{suffix}[/white][/no highlight]"
    if value_min > 0.0:
        return f"[no highlight][red]{prefix}", f"{suffix}[/red][/no highlight]"
    return f"[no highlight][blue]{prefix}", f"{suffix}[/blue][/no highlight]"


def _none_as_na(
    string: str | float | None, rounding: int = 3, prefix: str = "", suffix: str = ""
) -> str:
    """Return a string representation of the input, or 'N/A' if the input is None.

    Args:
    string (str | float | None): The input value to format.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    prefix (str, optional): Prefix to add before the value. Defaults to "".
    suffix (str, optional): Suffix to add after the value. Defaults to "".

    Returns:
    str: The formatted string or 'N/A' if the input is None.
    """
    if string is None:
        return "[grey]N/A[/grey]"
    if isinstance(string, float):
        string = round(string, rounding)
    return prefix + str(string) + suffix


def _get_mapcolored_value_or_na(
    key: str,
    dictionary: dict[str, float],
    min_value: int | float,
    max_value: int | float,
    apply: Callable[[float], float] | None = None,
    special_value: float = 0.0,
    rounding: int = 3,
    prefix: str = "",
    suffix: str = "",
) -> str:
    """Get a value from the dictionary, apply a function if specified, and format it with color coding.

    Args:
    key (str): The key to look up in the dictionary.
    dictionary (dict[str, float]): The dictionary containing the values.
    min_value (int | float): The minimum value for color mapping.
    max_value (int | float): The maximum value for color mapping.
    apply (Callable[[float], float] | None, optional): A function to apply to the value. Defaults to None.
    special_value (float, optional): A value that will be displayed in white. Defaults to 0.0.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    prefix (str, optional): Prefix to add before the value. Defaults to "".
    suffix (str, optional): Suffix to add after the value. Defaults to "".

    Returns:
    str: The formatted string with color coding or 'N/A' if the key is not found.
    """
    if key in dictionary:
        if apply is not None:
            value = _mapcolored(
                round(apply(dictionary[key]), rounding),
                apply(min_value),
                apply(max_value),
                special_value,
                prefix,
                suffix,
            )
        else:
            value = _mapcolored(
                round(dictionary[key], rounding),
                min_value,
                max_value,
                special_value,
                prefix,
                suffix,
            )
    else:
        value = None
    return _none_as_na(value)


def _get_value_or_na(
    key: str,
    dictionary: dict[str, float],
    rounding: int = 3,
    prefix: str = "",
    suffix: str = "",
) -> str:
    """Get a value from the dictionary and format it, or return 'N/A' if the key is not found.

    Args:
    key (str): The key to look up in the dictionary.
    dictionary (dict[str, float]): The dictionary containing the values.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    prefix (str, optional): Prefix to add before the value. Defaults to "".
    suffix (str, optional): Suffix to add after the value. Defaults to "".

    Returns:
    str: The formatted string or 'N/A' if the key is not found.
    """

    value = str(round(dictionary[key], rounding)) if key in dictionary else None
    return _none_as_na(value, rounding, prefix, suffix)


def _get_var_or_na(
    key: str,
    dictionary: dict[str, tuple[float, float]],
    rounding: int = 3,
    prefix: str = "",
    suffix: str = "",
) -> tuple[str, str]:
    """Get a variable range from the dictionary and format it, or return 'N/A' if the key is not found.

    Args:
    key (str): The key to look up in the dictionary.
    dictionary (dict[str, tuple[float, float]]): The dictionary containing the variable ranges.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    prefix (str, optional): Prefix to add before the value. Defaults to "".
    suffix (str, optional): Suffix to add after the value. Defaults to "".

    Returns:
    tuple[str, str]: A tuple containing the formatted minimum and maximum values, or 'N/A' for each if the key is not found.
    """

    if key in dictionary:
        values: tuple[float | None, float | None] = (
            dictionary[key][0],
            dictionary[key][1],
        )
    else:
        values = None, None
    return _none_as_na(values[0], rounding, prefix, suffix), _none_as_na(
        values[1], rounding, prefix, suffix
    )


def _zero_prefix(value: float | int) -> str:
    """Return an opening parenthesis if the value is zero, otherwise return an empty string.

    Args:
    value (float | int): The value to check.

    Returns:
    str: "(" if the value is zero, otherwise "".
    """
    return "(" if value == 0.0 else ""


def _zero_suffix(value: float | int) -> str:
    """Return a closing parenthesis if the value is zero, otherwise return an empty string.

    Args:
    value (float | int): The value to check.

    Returns:
    str: ")" if the value is zero, otherwise "".
    """
    return ")" if value == 0.0 else ""


def print_dict(dictionary: dict[Any, Any], indent: int = 4) -> None:
    """Pretty-print a dictionary in a JSON formatted string with the specified indentation.

    Args:
    dictionary (dict[Any, Any]): The dictionary to print.
    indent (int, optional): The number of spaces for indentation. Defaults to 4.
    """
    console.print(dumps(dictionary, indent=indent))


def print_model(
    cobrak_model: Model,
    print_reacs: bool = True,
    print_enzymes: bool = True,
    print_mets: bool = True,
    print_extra_linear_constraints: bool = True,
    print_settings: bool = True,
    conc_rounding: int = 6,
) -> None:
    """Pretty-print a detailed summary of the model, including reactions, enzymes, metabolites, and settings.

    Args:
    cobrak_model (Model): The model to print.
    print_reacs (bool, optional): Whether to print reactions. Defaults to True.
    print_enzymes (bool, optional): Whether to print enzymes. Defaults to True.
    print_mets (bool, optional): Whether to print metabolites. Defaults to True.
    print_extra_linear_constraints (bool, optional): Whether to print extra linear constraints. Defaults to True.
    print_settings (bool, optional): Whether to print general settings. Defaults to True.
    conc_rounding (int, optional): Number of decimal places to round concentrations to. Defaults to 6.
    """

    console.print("\n[b u]Model[b u]")

    if print_reacs:
        reac_table = Table(title="Reactions", title_justify="left")
        reac_table.add_column("ID")
        reac_table.add_column("String")
        reac_table.add_column("ΔG'°ᵢ")
        reac_table.add_column("kcatᵢ")
        reac_table.add_column("kMᵢ")
        # reac_table.add_column("kIᵢ")
        # reac_table.add_column("kAᵢ")
        reac_table.add_column("Name")
        reac_table.add_column("Annotation")

        for reac_id, reaction in sort_dict_keys(cobrak_model.reactions).items():
            arguments = [
                reac_id,
                get_reaction_string(cobrak_model, reac_id),
                _none_as_na(reaction.dG0),
                (
                    "N/A"
                    if reaction.enzyme_reaction_data is None
                    else str(reaction.enzyme_reaction_data.k_cat)
                ),
                (
                    "N/A"
                    if reaction.enzyme_reaction_data is None
                    else str(reaction.enzyme_reaction_data.k_ms)
                ),
                reaction.name,
                str(reaction.annotation),
            ]
            reac_table.add_row(*arguments)
        console.print(reac_table)

    if print_enzymes and cobrak_model.enzymes != {}:
        enzyme_table = Table(title="Enzymes", title_justify="left")
        enzyme_table.add_column("ID")
        enzyme_table.add_column("MWᵢ")
        enzyme_table.add_column("min([Eᵢ])")
        enzyme_table.add_column("max([Eᵢ])")
        enzyme_table.add_column("Name")
        enzyme_table.add_column("Annotation")

        for enzyme_id, enzyme in sort_dict_keys(cobrak_model.enzymes).items():
            arguments = [
                enzyme_id,
                str(enzyme.molecular_weight),
                _none_as_na(enzyme.min_conc),
                _none_as_na(enzyme.max_conc),
                enzyme.name,
                str(enzyme.annotation),
            ]
            enzyme_table.add_row(*arguments)
        console.print(enzyme_table)

    if print_mets:
        met_table = Table(title="Metabolites", title_justify="left")
        met_table.add_column("ID")
        met_table.add_column("min(cᵢ)")
        met_table.add_column("max(cᵢ)")
        met_table.add_column("Name")
        met_table.add_column("Annotation")

        for met_id, metabolite in sort_dict_keys(cobrak_model.metabolites).items():
            arguments = [
                met_id,
                str(round(exp(metabolite.log_min_conc), conc_rounding)),
                str(round(exp(metabolite.log_max_conc), conc_rounding)),
                metabolite.name,
                str(metabolite.annotation),
            ]
            met_table.add_row(*arguments)
        console.print(met_table)

    if print_extra_linear_constraints and cobrak_model.extra_linear_constraints != []:
        console.print("\n[b u]Extra linear constraints[b u]")
        for extra_linear_constraint in cobrak_model.extra_linear_constraints:
            console.print(get_extra_linear_constraint_string(extra_linear_constraint))

    if print_settings:
        console.print("\n[i]General settings[i]")
        print_strkey_dict_as_table(
            {
                "Protein pool": cobrak_model.T,
                "R [kJ⋅K⁻¹⋅mol⁻¹]": cobrak_model.R,
                "T [K]": cobrak_model.T,
                "Kinetic-ignored mets": ", ".join(
                    cobrak_model.kinetic_ignored_metabolites
                ),
            }
        )


def print_optimization_result(
    cobrak_model: Model,
    optimization_dict: dict[str, float],
    print_exchanges: bool = True,
    print_reactions: bool = True,
    print_enzymes: bool = True,
    print_mets: bool = True,
    print_error_values_if_existing: bool = True,
    add_stoichiometries: bool = False,
    rounding: int = 3,
    conc_rounding: int = 6,
    ignore_unused: bool = False,
    multiple_tables_per_line: bool = True,
    unused_limit: float = 1e-4,
) -> None:
    """Pretty-Print the results of an optimization, including exchanges, reactions, enzymes, and metabolites.

    Args:
    cobrak_model (Model): The model used for optimization.
    optimization_dict (dict[str, float]): A dictionary containing the optimization results.
    print_exchanges (bool, optional): Whether to print exchange reactions. Defaults to True.
    print_reactions (bool, optional): Whether to print non-exchange reactions. Defaults to True.
    print_enzymes (bool, optional): Whether to print enzyme usage. Defaults to True.
    print_mets (bool, optional): Whether to print metabolite concentrations. Defaults to True.
    add_stoichiometries (bool, optional): Whether to include reaction stoichiometries. Defaults to False.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    conc_rounding (int, optional): Number of decimal places to round concentrations to. Defaults to 6.
    ignore_unused (bool, optional): Whether to ignore reactions with zero flux. Defaults to False.
    multiple_tables_per_line (bool, optional): Whether to display multiple tables side by side. Defaults to True.
    """

    table_columns: list[Table] = []

    all_fluxes = [
        optimization_dict[reac_id]
        for reac_id in cobrak_model.reactions
        if reac_id in optimization_dict
    ]
    min_flux = min(all_fluxes)
    max_flux = max(all_fluxes)
    all_dfs = [
        optimization_dict[key]
        for key in optimization_dict
        if key.startswith(DF_VAR_PREFIX)
    ]

    substrate_reac_ids, product_reac_ids = (
        get_substrate_and_product_exchanges(cobrak_model, optimization_dict)
        if print_exchanges
        else ([""], [""])
    )

    if print_exchanges:
        for title, exchange_ids in (
            ("Substrates", substrate_reac_ids),
            ("Products", product_reac_ids),
        ):
            exchange_table = Table(title=title, title_justify="left")
            exchange_table.add_column("ID")
            exchange_table.add_column("Flux")
            for exchange_id in exchange_ids:
                exchange_flux = optimization_dict[exchange_id]
                if ignore_unused and exchange_flux <= unused_limit:
                    continue

                exchange_table.add_row(
                    exchange_id,
                    _mapcolored(
                        round(optimization_dict[exchange_id], rounding),
                        min_flux,
                        max_flux,
                        prefix=_zero_prefix(exchange_flux),
                        suffix=_zero_suffix(exchange_flux),
                    ),
                )
            table_columns.append(exchange_table)

    if print_reactions:
        reac_table = Table(
            title="Non-exchange reactions" if print_exchanges else "Reactions",
            title_justify="left",
        )
        reac_table.add_column("ID")
        reac_table.add_column("vᵢ")
        if add_stoichiometries:
            reac_table.add_column("Stoichiometries")
        reac_table.add_column("dfᵢ")
        reac_table.add_column("κᵢ")
        reac_table.add_column("γᵢ")
        for reac_id in sort_dict_keys(cobrak_model.reactions):
            if ignore_unused and (
                reac_id not in optimization_dict
                or optimization_dict[reac_id] <= unused_limit
            ):
                continue

            if (
                (reac_id not in optimization_dict)
                or (reac_id in product_reac_ids)
                or (reac_id in substrate_reac_ids)
            ):
                continue
            arguments: list[str] = [reac_id]
            if add_stoichiometries:
                arguments.append(get_reaction_string(cobrak_model, reac_id))

            reac_flux = optimization_dict[reac_id]
            prefix, suffix = _zero_prefix(reac_flux), _zero_suffix(reac_flux)

            arguments.extend(
                (
                    _mapcolored(
                        round(reac_flux, rounding),
                        min_flux,
                        max_flux,
                        prefix=prefix,
                        suffix=suffix,
                    ),
                    _get_mapcolored_value_or_na(
                        f"{DF_VAR_PREFIX}{reac_id}",
                        optimization_dict,
                        min(all_dfs) if len(all_dfs) > 0 else 0.0,
                        max(all_dfs) if len(all_dfs) > 0 else 0.0,
                        rounding=rounding,
                        prefix=prefix,
                        suffix=suffix,
                    ),
                    _get_mapcolored_value_or_na(
                        f"{KAPPA_VAR_PREFIX}{reac_id}",
                        optimization_dict,
                        0.0,
                        1.0,
                        rounding=rounding,
                        prefix=prefix,
                        suffix=suffix,
                    ),
                    _get_mapcolored_value_or_na(
                        f"{GAMMA_VAR_PREFIX}{reac_id}",
                        optimization_dict,
                        0.0,
                        1.0,
                        rounding=rounding,
                        prefix=prefix,
                        suffix=suffix,
                    ),
                )
            )
            reac_table.add_row(*arguments)
        table_columns.append(reac_table)

    if print_enzymes:
        enzyme_table = Table(title="Enzyme usage", title_justify="left")
        enzyme_table.add_column("Pool %")
        enzyme_table.add_column("Enzyme IDs")

        enzyme_usage = get_enzyme_usage_by_protein_pool_fraction(
            cobrak_model, optimization_dict
        )
        for pool_fraction, enzyme_ids in enzyme_usage.items():
            if ignore_unused and pool_fraction <= unused_limit:
                continue

            enzyme_table.add_row(
                _mapcolored(
                    round(pool_fraction * 100, rounding),
                    0.0,
                    100.0,
                    prefix=_zero_prefix(pool_fraction),
                    suffix=_zero_suffix(pool_fraction),
                ),
                "; ".join(enzyme_ids),
            )
        table_columns.append(enzyme_table)

    if print_mets:
        met_table = Table(title="Metabolites", title_justify="left")
        met_table.add_column("ID")
        met_table.add_column("Concentration")
        met_table.add_column("Consumption")
        met_table.add_column("Production")
        for met_id, metabolite in sort_dict_keys(cobrak_model.metabolites).items():
            met_var_id = f"{LNCONC_VAR_PREFIX}{met_id}"

            consumption, production = get_metabolite_consumption_and_production(
                cobrak_model, met_id, optimization_dict
            )

            if ignore_unused and production <= unused_limit:
                continue

            prefix, suffix = _zero_prefix(consumption), _zero_suffix(consumption)
            arguments = [met_id]
            arguments.append(
                _get_mapcolored_value_or_na(
                    met_var_id,
                    optimization_dict,
                    metabolite.log_min_conc,
                    metabolite.log_max_conc,
                    apply=exp,
                    special_value=1.0,
                    rounding=conc_rounding,
                    prefix=prefix,
                    suffix=suffix,
                )
            )

            arguments.append(_none_as_na(consumption, prefix=prefix, suffix=suffix))
            arguments.append(_none_as_na(production, prefix=prefix, suffix=suffix))

            met_table.add_row(*arguments)
        table_columns.append(met_table)

    if (
        print_error_values_if_existing
        and sum(
            key.startswith(ERROR_VAR_PREFIX) for key in list(optimization_dict.keys())
        )
        > 0
    ):
        error_table = Table(title="Errors", title_justify="left")
        error_table.add_column("ID")
        sorted_error_values = sort_dict_keys(
            {
                key[len(ERROR_VAR_PREFIX) + 1 :]: value
                for key, value in optimization_dict.items()
                if key.startswith(ERROR_VAR_PREFIX) and key != ERROR_SUM_VAR_ID
            }
        )
        min_error_value = min(list(sorted_error_values.values()))
        max_error_value = max(list(sorted_error_values.values()))
        for error_name, error_value in sorted_error_values.items():
            if ignore_unused and (error_value <= unused_limit):
                continue

            prefix, suffix = _zero_prefix(error_value), _zero_suffix(error_value)
            arguments = []
            arguments.append(error_name)
            arguments.append(
                _get_mapcolored_value_or_na(
                    error_name,
                    sorted_error_values,
                    min_value=min_error_value,
                    max_value=max_error_value,
                    prefix=prefix,
                    suffix=suffix,
                )
            )
            error_table.add_row(*arguments)
        error_table.add_row(*["SUM", str(optimization_dict[ERROR_SUM_VAR_ID])])
        table_columns.append(error_table)

    if multiple_tables_per_line:
        console.print(Columns(table_columns))
    else:
        for table in table_columns:
            console.print(table)

    console.print(
        "OBJECTIVE VALUE:",
        str(optimization_dict[OBJECTIVE_VAR_NAME]),
        "| SOLVE STATUS OK?",
        str(optimization_dict[ALL_OK_KEY]),
    )


def print_strkey_dict_as_table(
    dictionary: dict[str, Any],
    table_title: str = "",
    key_title: str = "",
    value_title: str = "",
) -> None:
    """Print a dictionary as a formatted table.

    Args:
    dictionary (dict[str, Any]): The dictionary to print.
    table_title (str, optional): The title of the table. Defaults to "".
    key_title (str, optional): The title for the key column. Defaults to "".
    value_title (str, optional): The title for the value column. Defaults to "".
    """
    table = Table(title=table_title, title_justify="left", show_header=False)
    table.add_column(key_title, style="cyan", no_wrap=True)
    table.add_column(value_title, style="magenta")
    for key, value in sort_dict_keys(dictionary).items():
        table.add_row(key, str(value))
    console.print(table)


def print_variability_result(
    cobrak_model: Model,
    variability_dict: dict[str, tuple[float, float]],
    print_exchanges: bool = True,
    print_reacs: bool = True,
    print_enzymes: bool = False,
    print_mets: bool = True,
    ignore_unused: bool = False,
    add_stoichiometries: bool = False,
    rounding: int = 3,
    multiple_tables_per_line: bool = True,
) -> None:
    """Print the variability analysis results, including exchanges, reactions, enzymes, and metabolites.

    Args:
    cobrak_model (Model): The model used for variability analysis.
    variability_dict (dict[str, tuple[float, float]]): A dictionary containing the variability results.
    print_exchanges (bool, optional): Whether to print exchange reactions. Defaults to True.
    print_reacs (bool, optional): Whether to print non-exchange reactions. Defaults to True.
    print_enzymes (bool, optional): Whether to print enzyme usage. Defaults to False.
    print_mets (bool, optional): Whether to print metabolite concentrations. Defaults to True.
    ignore_unused (bool, optional): Whether to ignore reactions with zero flux. Defaults to False.
    add_stoichiometries (bool, optional): Whether to include reaction stoichiometries. Defaults to False.
    rounding (int, optional): Number of decimal places to round to. Defaults to 3.
    multiple_tables_per_line (bool, optional): Whether to display multiple tables side by side. Defaults to True.
    """

    table_columns: list[Table] = []

    substrate_reac_ids, product_reac_ids = (
        get_substrate_and_product_exchanges(cobrak_model, variability_dict)
        if print_exchanges
        else ([""], [""])
    )

    reac_columns = [
        "ID",
        "min(vᵢ)",
        "max(vᵢ)",
        "min(dfᵢ)",
        "max(dfᵢ)",
    ]
    if add_stoichiometries:
        reac_columns.insert(1, "Reac string")

    if print_exchanges:
        for title, exchange_ids in (
            ("Substrates", substrate_reac_ids),
            ("Products", product_reac_ids),
        ):
            exchange_table = Table(title=title, title_justify="left")
            if add_stoichiometries:
                exchange_table.add_column("Reac string")
            for reac_column in reac_columns:
                exchange_table.add_column(reac_column)
            for exchange_reac_id in exchange_ids:
                prefix, suffix = _varcolor(exchange_reac_id, variability_dict)
                arguments: list[str] = [
                    exchange_reac_id,
                    *_get_var_or_na(
                        exchange_reac_id, variability_dict, rounding, prefix, suffix
                    ),
                    *_get_var_or_na(
                        f"{DF_VAR_PREFIX}{exchange_reac_id}",
                        variability_dict,
                        rounding,
                        prefix,
                        suffix,
                    ),
                ]
                exchange_table.add_row(*arguments)
            table_columns.append(exchange_table)

    if print_reacs:
        reacs_table = Table(
            title="Non-exchange reactions" if print_exchanges else "Reactions",
            title_justify="left",
        )
        for reac_column in reac_columns:
            reacs_table.add_column(reac_column)
        for reac_id in sort_dict_keys(cobrak_model.reactions):
            if reac_id in [*substrate_reac_ids, *product_reac_ids]:
                continue
            prefix, suffix = _varcolor(reac_id, variability_dict)

            arguments = [
                reac_id,
                *_get_var_or_na(reac_id, variability_dict, rounding, prefix, suffix),
                *_get_var_or_na(
                    f"{DF_VAR_PREFIX}{reac_id}",
                    variability_dict,
                    rounding,
                    prefix,
                    suffix,
                ),
            ]
            reacs_table.add_row(*arguments)
        table_columns.append(reacs_table)

    if print_enzymes:
        enzymes_table = Table(title="Enzymes", title_justify="left")
        enzymes_table.add_column("ID")
        enzymes_table.add_column("min(Eᵢ)")
        enzymes_table.add_column("max(Eᵢ)")

        for reac_id, reaction in sort_dict_keys(cobrak_model.reactions).items():
            if reaction.enzyme_reaction_data is None:
                continue
            enzyme_var_id = get_reaction_enzyme_var_id(reac_id, reaction)
            prefix, suffix = _varcolor(enzyme_var_id, variability_dict)
            reacs_table.add_row(
                *[
                    enzyme_var_id,
                    *_get_var_or_na(
                        enzyme_var_id, variability_dict, rounding, prefix, suffix
                    ),
                ]
            )
        table_columns.append(enzymes_table)

    if print_mets:
        mets_table = Table(title="Metabolites", title_justify="left")
        mets_table.add_column("ID")
        mets_table.add_column("min(cᵢ)")
        mets_table.add_column("max(cᵢ)")
        for met_id in sort_dict_keys(cobrak_model.metabolites):
            min_conc_str, max_conc_str = _get_var_or_na(
                f"{LNCONC_VAR_PREFIX}{met_id}", variability_dict, rounding=1_000
            )
            try:
                min_conc = str(round(exp(float(min_conc_str)), rounding))
                max_conc = str(round(exp(float(max_conc_str)), rounding))
            except ValueError:
                min_conc = min_conc_str
                max_conc = max_conc_str
            color = "[blue]" if min_conc != max_conc else "[red]"
            mets_table.add_row(
                *[
                    met_id,
                    f"{color} {min_conc}",
                    f"{color} {max_conc}",
                ]
            )
        table_columns.append(mets_table)

    if multiple_tables_per_line:
        console.print(Columns(table_columns))
    else:
        for table in table_columns:
            console.print(table)
