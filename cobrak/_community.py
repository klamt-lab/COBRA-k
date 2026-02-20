from copy import deepcopy
from cobrak.io import ExtraLinearWatch, ExtraNonlinearConstraint
from .dataclasses import Model, Reaction, Metabolite, Enzyme, ExtraNonlinearWatch, ExtraLinearConstraint, CommunitySpeciesSetting


def create_multiplied_model(
    models: dict[str, Model],
) -> Model:
    metabolites: dict[str, Metabolite] = {}
    reactions: dict[str, Reaction] = {}
    enzymes: dict[str, Enzyme] = {}
    extra_linear_constraints: list[ExtraLinearConstraint] = []
    extra_nonlinear_constraints: list[ExtraNonlinearConstraint] = []
    extra_linear_watches: dict[str, ExtraLinearWatch] = {}
    extra_nonlinear_watches: dict[str, ExtraNonlinearWatch] = {}
    kinetic_ignored_metabolites: list[str] = []
    for model_id, model in models.items():
        for met_id, metabolite in model.metabolites.items():
            metabolites[f"{met_id}_{model_id}"] = deepcopy(metabolite)
        for reac_id, reaction in model.reactions.items():
            reactions[f"{reac_id}_{model_id}"] = deepcopy(reaction)
            reactions[f"{reac_id}_{model_id}"].stoichiometries = {
                f"{key}_{model_id}": value
                for key, value in reaction.stoichiometries.items()
            }
            if reaction.enzyme_reaction_data:
                reactions[f"{reac_id}_{model_id}"].enzyme_reaction_data.identifiers = [
                    f"{identifier}_{model_id}" for identifier in reactions[f"{reac_id}_{model_id}"].enzyme_reaction_data.identifiers
                ]
                reactions[f"{reac_id}_{model_id}"].enzyme_reaction_data.k_ms = {
                    f"{key}_{model_id}": value
                    for key, value in reaction.enzyme_reaction_data.k_ms.items()
                }
                reactions[f"{reac_id}_{model_id}"].enzyme_reaction_data.k_is = {
                    f"{key}_{model_id}": value
                    for key, value in reaction.enzyme_reaction_data.k_is.items()
                }
                reactions[f"{reac_id}_{model_id}"].enzyme_reaction_data.k_as = {
                    f"{key}_{model_id}": value
                    for key, value in reaction.enzyme_reaction_data.k_as.items()
                }
        for enzyme_id, enzyme in model.enzymes.items():
            enzymes[f"{enzyme_id}_{model_id}"] = deepcopy(enzyme)
        for extra_linear_constraint in model.extra_linear_constraints:
            extra_linear_constraints.append(
                ExtraLinearConstraint(
                    stoichiometries={
                        f"{key}_{model_id}": stoich
                        for key, stoich in extra_linear_constraint.stoichiometries.items()
                    },
                    lower_value=extra_linear_constraint.lower_value,
                    upper_value=extra_linear_constraint.upper_value,
                )
            )
        for extra_nonlinear_constraint in model.extra_nonlinear_constraints:
            extra_nonlinear_constraints.append(
                ExtraNonlinearConstraint(
                    stoichiometries={
                        f"{key}_{model_id}": stoich_and_function
                        for key, stoich_and_function in extra_nonlinear_constraint.stoichiometries.items()
                    },
                    full_application=extra_nonlinear_constraint.full_application,
                    lower_value=extra_nonlinear_constraint.lower_value,
                    upper_value=extra_nonlinear_constraint.upper_value,
                )
            )
        for watch_name, extra_linear_watch in model.extra_linear_watches.items():
            extra_linear_watches[watch_name] = ExtraLinearWatch(
                    stoichiometries={
                        f"{key}_{model_id}": stoich
                        for key, stoich in extra_linear_watch.stoichiometries.items()
                    },
                )
        for watch_name, extra_nonlinear_watch in model.extra_nonlinear_watches.items():
            extra_nonlinear_watches[watch_name] = ExtraNonlinearWatch(
                    stoichiometries={
                        f"{key}_{model_id}": stoich
                        for key, stoich in extra_nonlinear_watch.stoichiometries.items()
                    },
                )

    first_model: Model = list(models.values())[0]
    return Model(
        metabolites=metabolites,
        reactions=reactions,
        enzymes=enzymes,
        extra_linear_watches=extra_linear_watches,
        extra_nonlinear_watches=extra_nonlinear_watches,
        extra_linear_constraints=extra_linear_constraints,
        extra_nonlinear_constraints=extra_nonlinear_constraints,
        kinetic_ignored_metabolites=kinetic_ignored_metabolites,
        R=first_model.R,
        T=first_model.T,
        annotation={},
        reac_enz_separator=first_model.reac_enz_separator,
        fwd_suffix=first_model.fwd_suffix,
        rev_suffix=first_model.rev_suffix,
        max_conc_sum=first_model.max_conc_sum,
        conc_sum_ignore_prefixes=first_model.conc_sum_ignore_prefixes,
        conc_sum_include_suffixes=first_model.conc_sum_include_suffixes,
        conc_sum_max_rel_error=first_model.conc_sum_max_rel_error,
        conc_sum_min_abs_error=first_model.conc_sum_min_abs_error,
        community_species_settings={
            species_id: CommunitySpeciesSetting(
                max_prot_pool=species_model.max_prot_pool,
                max_conc_sum=species_model.max_conc_sum,
                include_mets_in_prot_pool=species_model.include_mets_in_prot_pool,
            )
            for species_id, species_model in models.items()
        }
    )


def _add_exchange_compartment_to_community_model(
    community_model: Model,
    single_models: dict[str, Model],
    exchange_id: str = "xchg",
) -> Model:
    exchange_reactions: dict[str, Reaction] = {}
    exchange_metabolites: dict[str, Metabolite] = {}
    for model_id, model in single_models.items():
        for reaction in model.reactions.values():
            if len(reaction.stoichiometries) != 1:
                continue
            affected_met: str = list(reaction.stoichiometries)[0]
            if f"{affected_met}_{exchange_id}" not in exchange_metabolites:
                exchange_metabolites[f"{affected_met}_{exchange_id}"] = Metabolite(
                    name=f"Community exchange metabolite of {affected_met}",
                )
            met_stoich: int | float = reaction.stoichiometries[affected_met]
            if met_stoich < 0.0:
                exchange_reactions[f"OUT_{model_id}_{affected_met}"] = Reaction(
                    stoichiometries={
                        affected_met: -1.0,
                        f"{affected_met}_{exchange_id}": 1.0,
                    },
                )
                if f"OUT_{model_id}_{affected_met}_{exchange_id}" not in exchange_reactions:
                    exchange_reactions[f"OUT_{model_id}_{affected_met}_{exchange_id}"] = Reaction(
                        stoichiometries={
                            affected_met: -1.0,
                        },
                    )
            else:
                exchange_reactions[f"IN_{model_id}_{affected_met}"] = Reaction(
                    stoichiometries={
                        f"{affected_met}_{exchange_id}": -1.0,
                        affected_met: 1.0,
                    },
                )
                if f"IN_{model_id}_{affected_met}_{exchange_id}" not in exchange_reactions:
                    exchange_reactions[f"IN_{model_id}_{affected_met}_{exchange_id}"] = Reaction(
                        stoichiometries={
                            affected_met: 1.0,
                        },
                    )
    for reac_id, reac in exchange_reactions.items():
        community_model.reactions[reac_id] = reac
    for met_id, met in exchange_metabolites.items():
        community_model.metabolites[met_id] = met

    return community_model


def create_community_model_with_fixed_growth(
    models_and_biomass_reacs: dict[str, tuple[Model, str]],
    growth_rate: float,
    community_growth_reac_id: str = "GROWTH",
    community_growth_met_id: str = "BIOMASS",
    max_considered_flux: float = 1000.0,
) -> Model:
    models: dict[str, Model] = {
        key: value[0]
        for key, value in models_and_biomass_reacs.items()
    }
    community_model = create_multiplied_model(
        models=models,
    )
    community_model = _add_exchange_compartment_to_community_model(
        single_models=models,
        community_model=community_model,
    )

    # Biomass reaction handling
    biomass_reac_ids: list[str] = [
        f"{biomass_reac_id}_{model_id}" for model_id, (_, biomass_reac_id) in models_and_biomass_reacs.items()
    ]
    community_model.reactions[community_growth_reac_id] = Reaction(
        stoichiometries={community_growth_met_id: -1.0},
        min_flux=growth_rate,
        max_flux=growth_rate,
        name="Community fixed growth reaction",
    )
    for biomass_reac_id in biomass_reac_ids:
        community_model.reactions[biomass_reac_id].stoichiometries[community_growth_met_id] = 1.0
        community_model.kinetic_ignored_metabolites.append(community_growth_met_id)

    # Inhomogenous constraint handling
    for model_id, (model, biomass_reac_id) in models_and_biomass_reacs.items():
        for reac_id, reaction in model.reactions.items():
            if reaction.min_flux > 0.0:
                pseudo_met_id = f"Rsnake_LOWER_{reac_id}_{model_id}"
                community_model.metabolites[pseudo_met_id] = Metabolite()
                community_model.reactions[f"{reac_id}_{model_id}"].stoichiometries[pseudo_met_id] = 1.0
                community_model.reactions[f"{biomass_reac_id}_{model_id}"].stoichiometries[pseudo_met_id] = -reaction.min_flux / growth_rate
                community_model.reactions[f"EX_{pseudo_met_id}"] = Reaction(
                    stoichiometries={pseudo_met_id: -1.0},
                    min_flux=0.0,
                    max_flux=max_considered_flux,
                )
                community_model.kinetic_ignored_metabolites.append(pseudo_met_id)
            if reaction.max_flux < max_considered_flux:
                pseudo_met_id = f"Rsnake_UPPER_{reac_id}_{model_id}"
                community_model.metabolites[pseudo_met_id] = Metabolite()
                community_model.reactions[f"{reac_id}_{model_id}"].stoichiometries[pseudo_met_id] = 1.0
                community_model.reactions[f"{biomass_reac_id}_{model_id}"].stoichiometries[pseudo_met_id] = -reaction.max_flux / growth_rate
                community_model.reactions[f"IN_{pseudo_met_id}"] = Reaction(
                    stoichiometries={pseudo_met_id: +1.0},
                    min_flux=0.0,
                    max_flux=max_considered_flux,
                )
                community_model.kinetic_ignored_metabolites.append(pseudo_met_id)

    return community_model


def create_community_model_with_fixed_species_fractions(
    models_and_fractions_and_biomass_reac_ids: dict[str, tuple[Model, float, str]],
    community_growth_reac_id: str = "GROWTH",
    community_growth_met_id: str = "BIOMASS",
    max_considered_flux: float = 1_000.0,
) -> Model:
    models: dict[str, Model] = {
        key: value[0]
        for key, value in models_and_fractions_and_biomass_reac_ids.items()
    }
    community_model = create_multiplied_model(
        models=models,
    )
    community_model = _add_exchange_compartment_to_community_model(
        single_models=models,
        community_model=community_model,
    )

    biomass_reac_ids: list[str] = [
        f"{biomass_reac_id}_{model_id}" for model_id, (_, _, biomass_reac_id) in models_and_fractions_and_biomass_reac_ids.items()
    ]
    community_model.reactions[community_growth_reac_id] = Reaction(
        stoichiometries={community_growth_met_id: -1.0},
        min_flux=0.0,
        max_flux=max_considered_flux,
        name="Community growth reaction",
    )
    for biomass_reac_id in biomass_reac_ids:
        community_model.reactions[biomass_reac_id].stoichiometries[community_growth_met_id] = 1.0

    for model_id, (model, fraction, _) in models_and_fractions_and_biomass_reac_ids.items():
        for reac_id, reaction in model.reactions.items():
            if reaction.min_flux > 0.0:
                community_model.reactions[f"{reac_id}_{model_id}"].min_flux *= fraction
            if reaction.max_flux < max_considered_flux:
                community_model.reactions[f"{reac_id}_{model_id}"].max_flux *= fraction

    return community_model


def create_community_model_with_no_growth(
    models_and_fractions_and_biomass_reac_ids: dict[str, tuple[Model, float, str]],
) -> Model:
    models: dict[str, Model] = {
        key: value[0]
        for key, value in models_and_fractions_and_biomass_reac_ids.items()
    }
    community_model = create_multiplied_model(
        models=models,
    )
    return _add_exchange_compartment_to_community_model(
        single_models=models,
        community_model=community_model,
    )


def remove_community_suffix(
    community_suffixes: list[str],
    var_id: str,
) -> str:
    for community_suffix in community_suffixes:
        if var_id.endswith(f"_{community_suffix}"):
            return var_id[:-len(f"_{community_suffix}")]
    return var_id
