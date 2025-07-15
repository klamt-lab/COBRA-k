"""Contains all dataclasses (and enums) used by COBRAk to define a metabolic model and its extra constraints and optimization objective.

Dataclasses are similar to structs in C: They are not intended to have member functions, only other types of member variables.
The main dataclass used by COBRAk is Model, which contains the full information about the metabolic model. As member variables,
a Model contains further dataclasses (such as Reaction, Metabolite, ...).
As dataclass_json is also invoked, it is possible to store and load the COBRAk dataclasses as JSON.
"""

# IMPORT SECTION #
from copy import deepcopy
from math import log
from typing import Any, Literal, TypeAlias

from pydantic import Field, FiniteFloat, NonNegativeInt, PositiveFloat
from pydantic.dataclasses import dataclass

from .constants import (
    QUASI_INF,
    REAC_ENZ_SEPARATOR,
    REAC_FWD_SUFFIX,
    REAC_REV_SUFFIX,
    STANDARD_R,
    STANDARD_T,
)


# DATACLASSES SECTION #
@dataclass
class Enzyme:
    """Represents an enzyme in a metabolic model.

    Members:
        molecular_weight (float):
            The enzyme's molecular weight in kDa.
        min_conc (float | None):
            [Optional] If wanted, one can set a special minimal concentration
            for the enzyme.
            Defaults to None, i.e., no given concentration value (i.e., only the total
            enzyme pool is the limit).
        max_conc (float | None):
            [Optional] If wanted, one can set a special maximal concentration
            for the enzyme.
            Defaults to None, i.e., no given concentration value (i.e., only the total
            enzyme pool is the limit).
        annotation (dict[str, str | list[str]]):
            [Optional] Dictionary containing additional enzyme annotation,
            e.g., {"UNIPROT_ID": "b12345"}.
            Defaults to '{}'.
        name: str:
            [Optional] Colloquial name of enzyme
    """

    molecular_weight: FiniteFloat
    """The enzyme's molecular weight in kDa"""
    min_conc: PositiveFloat | None = None
    """[Optional] The enzyme's minimal concentration in mmol⋅gDW⁻¹"""
    max_conc: PositiveFloat | None = None
    """[Optional] The enzyme's minimal concentration in mmol⋅gDW⁻¹"""
    annotation: dict[str, str | list[str]] = Field(default_factory=dict)
    """[Optional] Any annotation data for the enzyme (e.g., references). Has no effect on calculations"""
    name: str = ""
    """Colloquial name of enzyme"""


@dataclass
class ParameterReference:
    """Represents the database reference for a kinetic parameter."""

    database: str = ""
    """(If given) The database from which this parameter was read. Defaults to ''."""
    comment: str = "(no refs)"
    """Any comment given for this value (e.g. literature)? Defaults to '(no refs)'."""
    species: str = ""
    """Scientific name of the species where this value was measured. Defaults to ''."""
    substrate: str = ""
    """The metabolite (or reaction substrate) for which this value was measured. Defaults to ''."""
    pubs: list[str] = Field(default_factory=list)
    """"""
    tax_distance: int | None = None
    value: float | None = None


@dataclass
class EnzymeReactionData:
    """Represents the enzymes used by a reaction."""

    identifiers: list[str]
    """The identifiers (must be given in the associated Model enzymes instance) of the reaction's enzyme(s)"""
    k_cat: PositiveFloat
    """The reaction's k_cat (turnover numbers) in h⁻¹"""
    k_cat_references: list[ParameterReference] = Field(default_factory=list)
    """[Optional] List of references showing the source(s) of the k_cat value"""
    k_ms: dict[str, PositiveFloat] = Field(default_factory=dict)
    """[Optional] The reaction's k_ms (Michaelis-Menten constants) in M=mol⋅l⁻¹. Metabolite IDs are keys, k_ms values. Default is {}"""
    k_m_references: dict[str, list[ParameterReference]] = Field(default_factory=dict)
    """[Optional] References showing the source(s) of the k_m values. Metabolite IDs are keys, the source lists values. Default is {}"""
    k_is: dict[str, PositiveFloat] = Field(default_factory=dict)
    """[Optional] The reaction's k_is (Inhibition constants) in M=mol⋅l⁻¹. Metabolite IDs are keys, k_ms values. Default is {}"""
    k_i_references: dict[str, list[ParameterReference]] = Field(default_factory=dict)
    """[Optional] References showing the source(s) of the k_i values. Metabolite IDs are keys, the source lists values. Default is {}"""
    k_as: dict[str, PositiveFloat] = Field(default_factory=dict)
    """[Optional] The reaction's k_as (Activation constants) in M=mol⋅l⁻¹. Metabolite IDs are keys, k_ms values. Default is {}"""
    k_a_references: dict[str, list[ParameterReference]] = Field(default_factory=dict)
    """[Optional] References showing the source(s) of the k_a values. Metabolite IDs are keys, the source lists values. Default is {}"""
    special_stoichiometries: dict[str, PositiveFloat] = Field(default_factory=dict)
    """[Optional] Special (non-1) stoichiometries of polypeptides in the reaction's enzyme. Default is {}"""


@dataclass
class ExtraLinearConstraint:
    """Represents a general linear Model constraint.

    This can affect not only reactions, but also all other
    variables set in a COBRAk model. E.g., if one wants (for whatever
    reason) the following constraint:
    0.5 <= [A] - 2 * r_R1 <= 2.1
    the corresponding ExtraLinearConstraint instance would be:
    ExtraLinearConstraint(
        stoichiometries = {
            "x_A": 1.0,
            "R1": -2,
        },
        lower_value = 0.5,
        upper_value = 2.1,
    )
    lower_value or upper_value can be None if no such limit is desired.
    """

    stoichiometries: dict[str, float]
    """Keys: Model variable names; Children: Multipliers of constraint"""
    lower_value: float | None = None
    """Minimal numeric constraint value. Either this and/or upper_value must be not None. Defaults to None."""
    upper_value: float | None = None
    """Maximal numeric constraint value. Either this and/or lower_value must be not None. Defaults to None."""


@dataclass
class Metabolite:
    """Represents a Model's metabolite."""

    log_min_conc: FiniteFloat = log(1e-6)
    """Maximal logarithmic concentration (only relevant for thermodynamic constraints); Default is log(1e-6 M)"""
    log_max_conc: FiniteFloat = log(0.02)
    """Maximal logarithmic concentration (only relevant for thermodynamic constraints); Default is log(0.02 M)"""
    annotation: dict[str, str | list[str]] = Field(default_factory=dict)
    """Optional annotation (e.g., CHEBI numbers, ...); Default is {}"""
    name: str = ""
    """Colloquial name of metabolite"""
    formula: str = ""
    """Chemical formula of metabolite"""
    charge: int = 0
    """Electron charge of metabolite"""


@dataclass
class Reaction:
    """Represents a Model's reaction.

    E.g., a reaction
    A -> B [0; 1000], ΔG'°=12.1 kJ⋅mol⁻¹, catalyzed by E1 with k_cat=1000 h⁻¹
    would be
    Reaction(
        stoichiometries: {
            "A": -1,
            "B": +1,
        },
        min_flux: 0,
        max_flux: 1000,
        dG0=12.1,
        dG0_uncertainty=None,
        enzyme_reaction_data=EnzymeReactionData(
            identifiers=["E1"],
            k_cat=1000,
            k_ms=None,
            k_is=None,
            k_as=None,
        ),
        annotation={}, # Can be also ignored
    )
    """

    stoichiometries: dict[str, float]
    """Metabolite stoichiometries"""
    min_flux: float
    """Minimal flux (for COBRAk, this must be ≥ 0)"""
    max_flux: float
    """Maximal flux (must be >= min_flux)"""
    dG0: FiniteFloat | None = None
    """If given, the Gibb's free energy of the reaction (only relevant for thermodynamic constraints); Default is None"""
    dG0_uncertainty: FiniteFloat | None = None
    """If given, the Gibb's free energy's uncertainty (only relevant for thermodynamic constraints); Default is None"""
    enzyme_reaction_data: EnzymeReactionData | None = None
    """If given, enzymatic data (only relevant for enzymatic constraints); Default is None"""
    annotation: dict[str, str | list[str]] = Field(default_factory=dict)
    """Optional annotation (e.g., KEGG identifiers, ...)"""
    name: str = ""
    """Colloquial name of reaction"""


@dataclass
class Model:
    """Represents a metabolic model in COBRAk.

    This includes its Reaction instances (which define the reaction stoichiometries),
    its Metabolite instances (which are referenced in the mentioned stoichiometries),
    as well as optional enzymatic and thermodynamic data.
    """

    metabolites: dict[str, Metabolite]
    """Keys: Metabolite IDs; Children: Metabolite instances"""
    reactions: dict[str, Reaction]
    """Keys: Reaction IDs; Children: Reaction instances"""
    enzymes: dict[str, Enzyme] = Field(default_factory=dict)
    """[Only neccessary with enzymatic constraints] Keys: Enzyme IDs; Children: Enzyme instances; default is {}"""
    max_prot_pool: PositiveFloat = Field(default=1e9)
    """[Only neccessary with enzymatic constraints] Maximal usable protein pool in g/gDW; default is 1e9, i.e. basically unrestricted"""
    extra_linear_constraints: list[ExtraLinearConstraint] = Field(default_factory=list)
    """[Optional] Extra linear constraints"""
    kinetic_ignored_metabolites: list[str] = Field(default_factory=list)
    """[Optional and only works with saturation term constraints] Metabolite IDs for which no k_m is neccessary"""
    R: PositiveFloat = Field(default=STANDARD_R)
    """[Optional and only works with thermodynamic constraints] Gas constant reference for dG'° in kJ⋅K⁻¹⋅mol⁻¹; default is STANDARD_R"""
    T: PositiveFloat = Field(default=STANDARD_T)
    """[Optional and only works with thermodynamic constraints] Temperature reference for dG'° in K; default is STANDARD_T"""
    annotation: dict[str, str | list[str]] = Field(default_factory=dict)
    """[Optional] Any annotation for the model itself (e.g., its name or references). Has no effect on calculations."""
    reac_enz_separator: str = REAC_ENZ_SEPARATOR
    """[Optional] String infix that separated reaction IDs of reaction with multiple enzyme variants from their enzyme ID. Defaults to '_ENZ_'"""
    fwd_suffix: str = REAC_FWD_SUFFIX
    """[Optional] Reaction ID suffix of forward reaction variants (e.g. in a reversible reaction A→B, for the direction A→B). Default is '_FWD'"""
    rev_suffix: str = REAC_REV_SUFFIX
    """[Optional] Reaction ID suffix of reverse reaction variants (e.g. in a reversible reaction A→B, for the direction B→A). Default is '_REV'"""
    max_conc_sum: float = float("inf")
    """[Optional and only works with thermodynamic constraints] Maximal allowed sum of concentrations (for MILPs: linear approximation; for NLPs: Exact value). Inactive if set to default value of float('inf')"""
    conc_sum_ignore_prefixes: list[str] = Field(default_factory=list)
    """[Optional and only works with thermodynamic constraints] """
    conc_sum_include_suffixes: list[str] = Field(default_factory=list)
    """[Optional and only works with thermodynamic constraints] """
    conc_sum_max_rel_error: float = 0.05
    """[Optional and only works with MILPs with thermodynamic constraints] Maximal relative concentration sum approximation error"""
    conc_sum_min_abs_error: float = 1e-6
    """[Optional and only works with MILPs with thermodynamic constraints] Maximal absolute concentration sum approximation error"""

    def __enter__(self):
        """Method called when entering 'with' blocks"""
        # Return a deep copy of self
        return deepcopy(self)

    def __exit__(self, a, b, c):
        """Method called when leaving a 'with' block"""
        # Cleanup can be performed here if needed
        return  # Return None to propagate any exceptions


@dataclass
class CorrectionConfig:
    """Stores the configuration for corrections in a model (see parameter corrections chapter in documentation)."""

    error_scenario: dict[str, tuple[float, float]] = Field(default_factory=list)
    """A dictionary where keys are error scenarios and values are tuples representing the lower and upper bounds of the error. Defaults to {}."""
    add_flux_error_term: bool = False
    """Indicates whether to add flux error terms. Defaults to False."""
    add_met_logconc_error_term: bool = False
    """Indicates whether to add metabolite log concentration error terms. Defaults to False."""
    add_enzyme_conc_error_term: bool = False
    """Indicates whether to add enzyme concentration error terms. Defaults to False."""
    add_kcat_times_e_error_term: bool = False
    """Indicates whether to add k_cat ⋅ [E] error terms. Defaults to False."""
    kcat_times_e_error_cutoff: PositiveFloat = 1.0
    """The cutoff value for the k_cat ⋅ [E] error term. Defaults to 1.0."""
    max_rel_kcat_times_e_correction: PositiveFloat = QUASI_INF
    """Maximal relative correction for the k_cat ⋅ [E] error error term. Defaults to QUASI_INF."""
    add_dG0_error_term: bool = False
    """Indicates whether to add ΔG'° error terms. Defaults to False."""
    dG0_error_cutoff: PositiveFloat = 1.0
    """The cutoff value for the ΔG'° error terms. Defaults to 1.0."""
    max_abs_dG0_correction: PositiveFloat = QUASI_INF
    """Maximal absolute correction for the dG0 error term. Defaults to QUASI_INF."""
    add_km_error_term: bool = False
    """Indicates whether to add a kappa error term. Defaults to False."""
    km_error_cutoff: PositiveFloat = 1.0
    """Cutoff value for the κ error term. Defaults to 1.0."""
    max_rel_km_correction: PositiveFloat = 0.999
    """Maximal relative correction for the κ error term. Defaults to 0.999."""
    error_sum_as_qp: bool = False
    """Indicates whether to use a quadratic programming approach for the error sum. Defaults to False."""
    add_error_sum_term: bool = True
    """Whether to add an error sum term. Defaults to True."""
    use_weights: bool = False
    """Indicates whether to use weights for the corrections (otherwise, the weight is 1.0). Defaults to False."""
    weight_percentile: NonNegativeInt = 90
    """Percentile to use for weight calculation. Defaults to 90."""
    extra_weights: dict[str, float] = Field(default_factory=dict)
    """Dictionary to store extra weights for specific corrections. Defaults to {}."""
    var_lb_ub_application: Literal["", "exp", "log"] = ""
    """The application method for variable lower and upper bounds. Either '' (x=x), 'exp' or 'log'. Defaults to ''."""


@dataclass
class Solver:
    """Represents options for a pyomo-compatible solver"""

    name: str
    """The solver's name. E.g. 'scip' for SCIP and 'cplex_direct' for CPLEX."""
    solver_options: dict[str, float | int | str] = Field(default_factory=dict)
    """[Optional] Options transmitted to the solver itself."""
    solver_attrs: dict[str, float | int | str] = Field(default_factory=dict)
    """[Optional] Options set on the solver object in pyomo."""
    solve_extra_options: dict[str, Any] = Field(default_factory=dict)
    """[Optional] Options set on pyomo's solve function."""


# SHORTHAND TYPE ALIASES
ErrorScenario: TypeAlias = dict[str, tuple[float, float]]
"""A COBRAk error scenario type alias for a ConfigurationConfig; Is dict[str, tuple[float, float]]"""
OptResult: TypeAlias = dict[str, float]
"""A COBRAk variability optimization result type alias; Is dict[str, float]"""
VarResult: TypeAlias = dict[str, tuple[float | None, float | None]]
"""A COBRAk variability result type alias; Is dict[str, tuple[float | None, float | None]]"""
