"""Runs all major analyses for the toymodel as shown in COBRA-k's initial publication"""

try:  # noqa: SIM105
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

import tempfile
from math import log

from cobrak.dataclasses import (
    ExtraLinearConstraint,
)
from cobrak.evolution import (
    perform_nlp_evolutionary_optimization,
)
from cobrak.example_models import toy_model
from cobrak.io import (
    load_annotated_sbml_model_as_cobrak_model,
    save_cobrak_model_as_annotated_sbml_model,
)
from cobrak.lps import perform_lp_variability_analysis
from cobrak.nlps import perform_nlp_reversible_optimization  # noqa: F401
from cobrak.printing import (
    print_dict,
    print_model,
    print_optimization_result,
    print_variability_result,
)
from cobrak.standard_solvers import BARON, IPOPT, SCIP  # noqa: F401


def test_toymodel_calculations() -> None:  # noqa: D103
    global toy_model  # noqa: PLW0603
    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as temp_sbml_file:
        save_cobrak_model_as_annotated_sbml_model(
            toy_model,
            filepath=temp_sbml_file.name,
        )
        toy_model = load_annotated_sbml_model_as_cobrak_model(
            filepath=temp_sbml_file.name
        )

    toy_model.extra_linear_constraints = [
        ExtraLinearConstraint(
            stoichiometries={
                "x_ATP": 1.0,
                "x_ADP": -1.0,
            },
            lower_value=log(3.0),
        )
    ]
    print_model(toy_model)

    # ecTFVA #
    variability_dict = perform_lp_variability_analysis(
        toy_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        min_flux_cutoff=1e-7,
    )
    print_variability_result(toy_model, variability_dict)

    # Evolutionary algorithm applications
    result = perform_nlp_evolutionary_optimization(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        variability_dict=variability_dict,
        with_kappa=True,
        with_gamma=True,
        with_alpha=False,
        with_iota=False,
        sampling_wished_num_feasible_starts=2,
        objvalue_json_path="",
        evolution_num_gens=10,
    )
    max_result = list(result.keys())[0]
    max_result_dict = list(result.values())[0][0]
    print_dict(max_result_dict)
    print_optimization_result(toy_model, max_result_dict)
    assert max_result > 45.430
    assert max_result < 45.432

    variability_dict["EX_S"] = (0.0, 14.0)
    result = perform_nlp_evolutionary_optimization(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        variability_dict=variability_dict,
        with_kappa=True,
        with_gamma=True,
        with_alpha=False,
        with_iota=False,
        sampling_wished_num_feasible_starts=2,
        objvalue_json_path="",
        evolution_num_gens=10,
    )
    max_result = list(result.keys())[0]
    assert max_result > 32.718
    assert max_result < 32.719
