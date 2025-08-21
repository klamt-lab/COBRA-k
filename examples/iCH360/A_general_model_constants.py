kinetic_ignored_metabolites: list[str] = [  # noqa: D100
    "h2o_c",
    "h2o_p",
    "h2o_e",
    "h_c",
    "h_p",
    "h_e",
    "co2_c",
    "co2_p",
    "co2_e",
    "pi_c",
    "pi_p",
    "pi_e",
]
kinetic_ignored_enzyme_ids: list[str] = ["s0001"]

conc_ranges = {
    "DEFAULT": (1e-6, 0.1),
    "h_c": (1.0, 1.0),
    "h_p": (1.0, 1.0),
    "h_e": (1.0, 1.0),
    "h2o_c": (1.0, 1.0),
    "h2o_p": (1.0, 1.0),
    "h2o_e": (1.0, 1.0),
}
