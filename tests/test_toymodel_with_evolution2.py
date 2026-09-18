"""Runs all major analyses for the toymodel as shown in COBRA-k's initial publication"""

from cobrak.constants import OBJECTIVE_VAR_NAME

try:  # noqa: SIM105
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

import tempfile
from math import log

from cobrak._evolution2 import (
    EvolutionSettings,
    perform_nlp_evolutionary_optimization,
)
from cobrak.dataclasses import (
    ExtraLinearConstraint,
)
from cobrak.example_models import toy_model
from cobrak.io import (
    load_annotated_sbml_model_as_cobrak_model,
    save_cobrak_model_as_annotated_sbml_model,
)
from cobrak.lps import perform_lp_variability_analysis
from cobrak.nlps import (  # noqa: F401
    perform_nlp_irreversible_optimization_with_active_reacs_only,
    perform_nlp_reversible_optimization,
)
from cobrak.printing import (
    print_dict,
    print_model,
    print_optimization_result,
    print_variability_result,
)
from cobrak.standard_solvers import IPOPT, IPOPT_MA57


def test_toymodel_calculations() -> None:  # noqa: D103
    global toy_model  # noqa: PLW0603
    IPOPT.solver_options["acceptable_tol"] = 1e-10
    IPOPT.solver_options["max_iter"] = 100_000
    IPOPT.solver_options["mu_strategy"] = "adaptive"
    IPOPT.solver_options["corrector_type"] = "primal-dual"

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
    settings = EvolutionSettings(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        variability_dict=variability_dict,
        num_gens=40,
        population_size=10,
        with_kappa=True,
        with_gamma=True,
        with_alpha=False,
        with_iota=False,
        nlp_solver=IPOPT_MA57,
        use_original_ectfba_nlp_binaries_as_template=(True, True, True),
        inner_lp_objectives=("MAXZ",),
    )

    result = perform_nlp_evolutionary_optimization(
        settings=settings,
    )
    max_result = result.best_nlps[0]
    print_dict(max_result)
    print_optimization_result(toy_model, max_result)
    print(result.best_nlps[0]["ATP_Consumption"])
    assert max_result[OBJECTIVE_VAR_NAME] > 45.430
    assert max_result[OBJECTIVE_VAR_NAME] < 45.432

    settings.variability_dict["EX_S"] = (0.0, 14.0)
    result2 = perform_nlp_evolutionary_optimization(
        settings=settings,
    )
    max_result2_value = result2.best_nlps[0][OBJECTIVE_VAR_NAME]
    print(result2.best_nlps[0]["ATP_Consumption"])
    assert max_result2_value > 32.718
    assert max_result2_value < 32.719


if __name__ == "__main__":
    test_toymodel_calculations()
