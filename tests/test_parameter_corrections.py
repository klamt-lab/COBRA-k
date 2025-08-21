# Import relevant classes and functions  # noqa: D100
from math import log

from cobrak.constants import ERROR_SUM_VAR_ID, LNCONC_VAR_PREFIX
from cobrak.dataclasses import CorrectionConfig
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_optimization
from cobrak.utilities import apply_error_correction_on_model


def test_parameter_corrections() -> None:  # noqa: D103
    flux_and_concentration_error_scenario = {
        "Overflow": (1.0, 1.4),  # Overflow reaction flux between 1 and 1.4
        f"{LNCONC_VAR_PREFIX}M": (
            log(0.2),
            log(0.2),
        ),  # M concentration fixed at .2 molar
        f"{LNCONC_VAR_PREFIX}D": (
            log(0.23),
            log(0.25),
        ),  # D concentration betwewen .23 and .25 molar
    }

    # With a CorrectionConfig as optional further argument,
    # all explained extra variables for parameter corrections
    # are added to our model automatically :D
    # Now, we can just minimize the sum of correction errors.
    perform_lp_optimization(
        cobrak_model=toy_model,
        objective_target=ERROR_SUM_VAR_ID,
        objective_sense=-1,
        with_thermodynamic_constraints=True,
        correction_config=CorrectionConfig(
            error_scenario=flux_and_concentration_error_scenario,
            add_flux_error_term=True,
            add_met_logconc_error_term=True,
        ),
    )

    # k_cat*[E] correction
    # Import relevant classes and functions
    flux_and_concentration_error_scenario = {
        "Glycolysis": (40.0, 45.0),
    }

    # Again, minimize the correction error variable sum
    correction_result_2 = perform_lp_optimization(
        cobrak_model=toy_model,
        objective_target=ERROR_SUM_VAR_ID,
        objective_sense=-1,
        with_thermodynamic_constraints=True,
        with_enzyme_constraints=True,
        correction_config=CorrectionConfig(
            error_scenario=flux_and_concentration_error_scenario,
            add_kcat_times_e_error_term=True,
            add_dG0_error_term=True,
            add_km_error_term=True,
        ),
    )

    # Now, we apply the correction (i.e. set the corrected
    # parameter values to our model, overwriting the old parameter values)
    corrected_cobrak_model = apply_error_correction_on_model(
        cobrak_model=toy_model,
        correction_result=correction_result_2,
        min_abs_error_value=0.01,
        min_rel_error_value=0.01,
        verbose=True,
    )

    # Check that Glycolysis flux of at least 40 is reached
    result = perform_lp_optimization(
        cobrak_model=corrected_cobrak_model,
        objective_target="Glycolysis",
        objective_sense=+1,
        with_thermodynamic_constraints=True,
        with_enzyme_constraints=True,
    )
    assert result["Glycolysis"] >= 39.99
