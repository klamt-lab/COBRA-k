"""This module contains all COBRAk constants that are used throughout its packages.

These constants are especially used for problem constructions (to determine prefixes,
suffixes, names, ... for pyomo variables) as well as thermodynamic standard values.
"""

# CONSTANTS #
ALL_OK_KEY = "ALL_OK"
"""Shows that the result is optimal and the termination condition is ok"""

BIG_M = 10_000
"""Big M value for MILPs"""

BIGG_COMPARTMENTS = [
    "c",
    "e",
    "p",
    "m",
    "x",
    "r",
    "v",
    "n",
    "g",
    "u",
    "l",
    "h",
    "f",
    "s",
    "im",
    "cx",
    "um",
    "cm",
    "i",
    "mm",
    "w",
    "y",
]
"""List of BiGG compartment suffixes (without _) as defined in http://bigg.ucsd.edu/compartments/"""

DF_VAR_PREFIX = "f_var_"
"""Prefix for driving force problem variables"""

DG0_VAR_PREFIX = "dG0_"
"""Prefix for Gibb's free energy problem variables"""


EC_INNER_TO_OUTER_COMPARTMENTS: list[str] = [
    "c",
    "p",
    "e",
]
"""Inner to outer compartments in this order for E. coli models used here"""


EC_IONIC_STRENGTHS: dict[str, float] = {
    "c": 250,  # Source: eQuilibrator standard
    "p": 250,
    "e": 250,
}
"""Ionic strenghts (in mM) for E. coli model compartments used here"""


EC_PHS: dict[str, float] = {
    "c": 7.5,  # Source: Bionumbers ID 105980
    "p": 7.5,  # Source: Bionumbers ID 105980
    "e": 7.5,
}
"""pH values (unitless) for E. coli model compartments used here"""

EC_PMGS: dict[str, float] = {
    "c": 2.5,  # Source: eQuilibrator standard
    "p": 2.5,
    "e": 2.5,
}
"""pMg values (unitless) for E. coli model compartments used here"""

EC_POTENTIAL_DIFFERENCES: dict[tuple[str, str], float] = {
    ("c", "p"): -0.15,  # Source: eQuilibrator standard
    ("p", "e"): -0.15,
}
"""Potential differences (in V) for E. coli model compartments used here"""

ENZYME_VAR_PREFIX = "enzyme_"
"""Prefix of problem variables which stand for enzyme concentrations"""

ENZYME_VAR_INFIX = "_of_"
"""Infix for separation of enzyme name and reaction name"""

ERROR_CONSTRAINT_PREFIX = "flux_error_"
"""Prefix for the constraint that defines a scenario flux constraint"""

ERROR_BOUND_LOWER_CHANGE_PREFIX = "bound_error_change_lower_"
"""Prefix for fixed variables that show how much affected lower variable bounds have to be changed"""

ERROR_BOUND_UPPER_CHANGE_PREFIX = "bound_error_change_upper_"
"""Prefix for fixed variables that show how much affected lower variable bounds have to be changed"""

ERROR_VAR_PREFIX = "error_"
"""Prefix for error term variables for feasibility-making optimizations"""

ERROR_SUM_VAR_ID = "error_sum"
"""Name for the variable that holds the sum of all error term variables"""

FLUX_SUM_VAR_ID = "FLUX_SUM_VAR"
"""Name of optional variable that holds the sum of all reaction fluxes"""

KAPPA_PRODUCTS_VAR_PREFIX = "kappa_products_"
"""Prefix for variables representing the sum of logairthmized product concentration minus the logarithmized sum of km values"""

KAPPA_SUBSTRATES_VAR_PREFIX = "kappa_substrates_"
"""Prefix for variables representing the sum of logairthmized substrate concentration minus the logarithmized sum of km values"""

KAPPA_VAR_PREFIX = "kappa_var_"
"""Prefix for variables representing the thermodynamic restriction of a reaction (used in non-linear programs)"""

GAMMA_VAR_PREFIX = "gamma_var_"
"""Prefix for variables representing the thermodynamic restriction of a reaction (used in non-linear programs)"""

IOTA_VAR_PREFIX = "iota_var_"
"""Prefix for variables representing the inhibition of a reaction (used in non-linear programs)"""

ALPHA_VAR_PREFIX = "alpha_var_"
"""Prefix for variables representing the activation of a reaction (used in non-linear programs)"""

LNCONC_VAR_PREFIX = "x_"
"""Prefix for logarithmized concentration problem variables"""

MDF_VAR_ID = "var_B"
"""Name for minimally occuring driving force variable"""

OBJECTIVE_CONSTRAINT_NAME = "objective_constraint"
"""Name for constraint that defines the objective function's term"""

OBJECTIVE_VAR_NAME = "OBJECTIVE_VAR"
"""Name for variable that holds the objective value"""

PROTONS_PER_ELEMENT = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "Na": 11,
    "Mg": 12,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
    "Ca": 20,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Se": 34,
    "Mo": 42,
    "I": 53,
    "Bi": 83,
}
"""Protons per element (useful for eQulibrator compound searches)"""

PROT_POOL_MET_NAME = "prot_pool"
"""Identifier of the protein pool representing pseudo-metabolite"""

PROT_POOL_REAC_NAME = PROT_POOL_MET_NAME + "_delivery"
"""Identifier of the pseudo-reaction which created the protein pool pseudo-metabolite"""

QUASI_INF = 100_000
"""Big number (larger than big M) for values that would reach inf (thereby potentially causing solver problems)"""

REAC_ENZ_SEPARATOR = "_ENZ_"
"""Separator between enzyme-constrained reaction ID and attached enzyme name"""

REAC_FWD_SUFFIX = "_FWD"
"""Standard suffix for reaction IDs that represent forward directions of originally irreversible reactions"""

REAC_REV_SUFFIX = "_REV"
"""Standard suffix for reaction IDs that represent reverse directions of originally irreversible reactions"""

STANDARD_MAX_PROT_POOL = 0.25
"""Just a (for E. coli) quite high pool of metabolic enzymes on the total dry weight mass (in g⋅gDW⁻¹)."""

STANDARD_CONC_RANGES = {
    "DEFAULT": (1e-6, 0.2),
    "h_c": (1.0, 1.0),
    "h_p": (1.0, 1.0),
    "h_e": (1.0, 1.0),
    "h20_c": (1.0, 1.0),
    "h20_p": (1.0, 1.0),
    "h20_e": (1.0, 1.0),
}
"""Standard concentration ranges applicable to models with BiGG IDs; water and protons are set to one
as their effect is directly included in the ΔG'° calculation (see the eQuilibrator FAQ), while
the rest is set to wide ranges."""

SOLVER_STATUS_KEY = "SOLVER_STATUS"
"""Solver status optimization dict key"""

STANDARD_MIN_MDF = 1e-3
"""Standard minimally ocurring driving force for active reactions in kJ⋅mol⁻¹"""

STANDARD_R = 8.314e-3
"""Standard gas constant in kJ⋅K⁻1⋅mol⁻1 (Attention: Standard value is often given in J⋅K⁻1⋅mol⁻1, but we need in kJ⋅K⁻1⋅mol⁻1)"""

STANDARD_T = 298.15
"""Standard temperature in Kelvin"""

TERMINATION_CONDITION_KEY = "TERMINATION_CONDITION"
"""Solver termination condition key in optimization dict"""

USED_IDENTIFIERS_FOR_EQUILIBRATOR = [
    "inchi",
    "inchi_key",
    "metanetx.chemical",
    "bigg.metabolite",
    "kegg.compound",
    "chebi",
    "sabiork.compound",
    "metacyc.compound",
    "hmdb",
    "swisslipid",
    "reactome",
    "lipidmaps",
    "seed.compound",
]
"""Standard bunch of reaction identifier annotation names for E. coli models used here"""

Z_VAR_PREFIX = "z_var_"
"""Prefix of z variables (used with thermodynamic constraints in MI(N)LPs)"""

ZB_VAR_PREFIX = "zb_var_"
"""Extra zb variable prefix for thermodynamic bottleneck analyses"""

GENERALIZED_SUM_CONSTRAINT_NAME: str = "generalized_sum_constraint"
"""Name of constraint for generalized sum (protein pool + metabolite concentrations)"""
