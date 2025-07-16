# Parameter corrections

??? abstract "Quickstart code"
    ```py
    # Import relevant classes and functions
    from cobrak.lps import perform_lp_optimization
    from cobrak.example_models import toy_model
    from cobrak.dataclasses import CorrectionConfig
    from cobrak.constants import LNCONC_VAR_PREFIX, ERROR_SUM_VAR_ID
    from cobrak.printing import print_dict
    from math import log

    #
    flux_and_concentration_error_scenario = {
        "Overflow": (1.0, 1.4),  # Overflow reaction flux between 1 and 1.4
        f"{LNCONC_VAR_PREFIX}M": (log(0.2), log(0.2)),  # M concentration fixed at .2 molar
        f"{LNCONC_VAR_PREFIX}D": (log(0.23), log(0.25)),  # D concentration betwewen .23 and .25 molar
    }

    # With a CorrectionConfig as optional further argument,
    # all explained extra variables for parameter corrections
    # are added to our model automatically :D
    # Now, we can just minimize the sum of correction errors.
    correction_result_1 = perform_lp_optimization(
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

    print_dict(correction_result_1)

    #% k_cat*[E] correction
    # Import relevant classes and functions
    from cobrak.lps import perform_lp_optimization
    from cobrak.example_models import toy_model
    from cobrak.dataclasses import CorrectionConfig
    from cobrak.constants import LNCONC_VAR_PREFIX, ERROR_SUM_VAR_ID
    from cobrak.printing import print_dict
    from cobrak.utilities import apply_error_correction_on_model

    #
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

    print_dict(correction_result_2)

    # Now, we apply the correction (i.e. set the corrected
    # parameter values to our model, overwriting the old parameter values)
    corrected_cobrak_model = apply_error_correction_on_model(
        cobrak_model=toy_model,
        correction_result=correction_result_2,
        min_abs_error_value=0.01,
        min_rel_error_value=0.01,
        verbose=True,
    )
    ```

## Introduction

Often, when you just created a COBRA-k model, you'll find out that known *in vivo* flux/concentration/etc. measurements do not work with your model.
In other words, with your model, the given *in vivo* flux/concentration/etc. scenario is infeasible :-(

This infeasibility can be caused by too restrictive parameters in your model. It is also possible that the measurements were perfomed with an error so that their resulting values need to be corrected. In other words, the following possibilities may have caused a scenario infeasibility:

1. [If reaction fluxes caused the infeasibility] The measurement fluxes contain an error and need to be corrected
2. [If metabolite concentrations caused it] The measurement concentrations contain an error and need to be corrected
3. [If a measured flux cannot be reached] One or more $k_{cat}$ values in your model is too low.
4. [If a measured flux appears to be thermodynamically impossible] One or more $ΔG'°$ values in your model is too high, making thermodynamic term(s) $γ$ too low.
5. [If a measured flux appears to be impossible with regards to its associated saturation terms] One or more $k_M$ values in your model is either too low (for product-associated $k_M$) or too high (for substrate-associated $k_M$), making saturation term(s) $κ$ too low.

COBRA-k contains a rational way to tackle these infeasibilities and provides functions to find out the *lowest needed* measurement and/or parameter corrections to make a scenario feasible :-)

## Formulation of minimal correction calculation

### General sum of all correction types

!!! info
    The mentioned minimal correction calculation works with LPs (see LP chapter), MILPs (see MILP chapter), NLPs (see NLP chapter) and the evolutionary optimization (see evolution chapter). Therefore, we concentrate on the addition of variables to all these problems and the addition of the minimal correction objective.

As mentioned, we want to *minimize* the needed correction(s) to make a scenario feasible with our model. So in general, we have the following objective:

$$
\operatorname*{\mathbf{min}} \ \sum{corr^{fluxes} + corr^{concentrations} + corr^{k_{cat}⋅[E]}} + corr^{Δ_r G^{´°}} + corr^{k_M}  \\
s.t. \space LP, \space MILP, \space NLP \space or \space evolution \space constraints
$$

where $corr^X$ stands for a given correction error sum type. You can freely select whether you want every or just some of the correction types. If you ignore a correction type, its correction sum will be simply set to 0. Each of the correction types is explained in the following paragraphs:

### Flux scenario correction

If we have given reaction flux measurements $\mathbf{v^{measured,min}}$ and $\mathbf{v^{measured,max}}$ (i.e. lower and higher measured values, as e.g. resulting from a standard deviation), we can formulate our correction $corr^{fluxes}$ as follows for any reaction $i$ that was measured and is included in $\mathbf{v^{measured}}$:

$corr^{fluxes} = \sum_i {corr^{fluxes,+}_i + corr^{fluxes,-}_i} $

So, for each measured reaction $i$, we have an adding correction variable corr^{fluxes,+}_i and a subtracting variable corr^{fluxes,-}_i. In our optimization problem, they are introduced with the following constraints:

$v_i ≥ v^{measured,min}_i - corr^{fluxes,-}_i$
$v_i ≤ v^{measured,max}_i + corr^{fluxes,+}_i$

This way, $v_i$ is now able to be corrected such that the measured flux range can be reached :-)

### Logarithmic concentration scenario correction

!!! note
    This correction is very similar to the one with fluxes (see above)

If we have given logarithmic concentration measurements $\mathbf{\tilde{x}^{measured}}$, we can formulate our correction $corr^{concentrations}$ as follows for any metabolite $j$ that was measured and is included in $\mathbf{\tilde{x}^{measured}}$:

$corr^{concentrations} = \sum_i {corr^{concentrations,+}_i + corr^{concentrations,-}_i} $

So, for each measured metabolite $j$, we have an adding correction variable corr^{concentrations,+}_i and a subtracting variable corr^{concentrations,-}_i. In our optimization problem, they are introduced with the following constraints:

$\tilde{x}_j ≤ \tilde{x}^{measured,min}_i - corr^{concentrations,-}_i$
$\tilde{x}_j ≥ \tilde{x}^{measured,max}_i + corr^{concentrations,+}_i - corr^{concentrations,-}_i$

This way, $\tilde{x}_j$ is now able to be corrected such that the measured concentrations can be reached :-)

### $k_{cat}⋅[E]$ correction

E.g. if we have a given scenario with reaction fluxes and we see that our model cannot reach these reaction fluxes (i.e. our model's reachable fluxes are too low), it
might be a good idea to run a $V^+$ (i.e. maximal possible flux) correction. In other words, we try to make $k_{cat}⋅[E]$ so high that the fluxes can be reached.
To do so, we modify our reaction flux terms as follows:

* In LPs with enzyme constraints (see LP chapter):

$ v_i ≤ E_i ⋅ k_{cat,i} $

becomes

$ v_i ≤ E_i ⋅ k_{cat,i}^+ + corr^{k_{cat}⋅[E]}_i $

* In NLPs with enzyme constraints (see NLP chapter):

$ v_i ≤ V^{+}_i ⋅ κ_i ⋅ γ_i $

becomes

$ v_i ≤ (V^{+}_i + corr^{k_{cat}⋅[E]}) ⋅ κ_i ⋅ γ_i $

Hence, the full sum of $k_{cat}⋅[E]$ corrections is:

$corr^{k_{cat}⋅[E]} = \sum_i {corr^{k_{cat}⋅[E]}_i} $

!!! note "Why do we correct $k_{cat}⋅[E]$ and not just the $k_{cat}$?"

    If we would correct the $k_{cat}$ only, any $k_{cat,i}$ would essentially become a variable. Thus, e.g., the LP term $v_i ≤ E_i ⋅ k_{cat}^+ $ would become $v_i ≤ E_i ⋅ k_{cat} ⋅ corr^{maximal_flux}$, which would be non-linear and computationally typically too expensive. Using $k_{cat}⋅[E]$ omits this problem in LPs (see above).

### $ΔG'°$ correction

Another possibility why a scenario isn't reached by a model are too strict thermodynamic constraints. I.e. in MILPs (see MILP chapter) and NLPs, the $ΔG'°$ of a reaction
might be *too high* making a scenario either thermodynamically infeasible (in MILPs and NLPs) or lowering the thermodynamic term $γ_i$ (in NLPs) too much. To be able to correct this, we modify our thermodynamic terms as follows:

* In MILPs and NLPs with thermodynamic constraints (see MILP and NLP chapter):

$ f_i = -(Δ_r G^{´°}_i + R ⋅ T ⋅ \mathbf{N_{⋅,i}} ⋅ \mathbf{x̃}) $ becomes $ f_i = -(Δ_r G^{´°}_i + R ⋅ T ⋅ \mathbf{N_{⋅,i}} ⋅ \mathbf{x̃} - corr^{Δ_r G^{´°}}_i) $

Hence, the full sum of $Δ_r G^{´°}$ corrections is:

$corr^{Δ_r G^{´°}} = \sum_i {corr^{Δ_r G^{´°}}_i} $


### Logarithmic $k_M$ correction

Furthermore, it is possible to correct the $k_M$ (Michaelis-Menten constant) value of reactions to reach a scenario. This is only relevant if you have an NLP
and saturation terms ($κ_i$) are restricting your reactions. In these NLPs (see NLP chapter), the corrections are included as follows:

* If a metabolite is a substrate (i.e. if a *lower* $k_M$ could raise $κ_i$):

$ \tilde{s}_i = \ln ( \bar{s} ) = \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅ x̃_j ) - \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅  \ln K_{M,j,i} ) $

becomes

$ \tilde{s}_i = \ln ( \bar{s} ) = \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅ x̃_j ) - \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅  \ln K_{M,j,i} - corr^{k_{M,j,i}} ) $


* If a metabolite is a product (i.e. where a *higher* $k_M$ could raise $κ_i$):

$ \tilde{p}_i ≥ \ln ( \bar{p} ) = \sum_{k ∈ 𝖯_i} ( N_{k,i} ⋅ x̃ ) - \sum_{k ∈ 𝖯_i} ( N_{k,i} ⋅  \ln K_{M,k,i} ) $

becomes

$ \tilde{s}_i = \ln ( \bar{s} ) = \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅ x̃_j ) - \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅  \ln K_{M,j,i} + corr^{k_{M,j,i}} ) $

!!! warning
    All these $k_M$ corrections work in the logarithmic space. To get the "real" correction, you have to apply the exponential function on the respective correction value.

Hence, the full sum of $k_M$ corrections is:

$corr^{k_M} = \sum_{j,i} {k_{M,j,i}} $

#### Optional Weights

Up to now in our formulations, all $corr$ variables got the same weight, i.e. they are all multiplied with 1.
But often, we want to use other weights, especially when mixing different types of corrections, with corrections that modify parameters in different orders of magnitude.
Therefore, COBRA-k optionally allows one to set weights according to the following procedure, whereby $n$ is a user-defined percentile:

* For reaction fluxes and logarithmic concentrations: The absolute value of the given lower scenario bound
* For $k_{cat}⋅[E]$ corrections: The $n$th percentile of all possible maximal $k_{cat}⋅[E]$ values in the model, i.e. of all $Ω/W_i$ values.
* For $Δ_r G^{´°}_i$ corrections: The $n$th percentile of all absolute Δ_r G^{´°}_i in the model
* For $k_M$ corrections: The $n$th percentile of all $K_{m,i}$ in the model

#### Optional Quadratic Sum

As an alternative to the linear correction value sum (see above), one can also use a quadratic sum. Such a sum automatically penalizes larger absolute changes and prefers (potentially more) smaller changes. This changes our linear objective function into a quadratic one. E.g., the sum of flux corrections becomes:

$corr^{fluxes} = \sum_i {(corr^{fluxes,+}_i)^2 + (corr^{fluxes,-}_i)^2} $

!!! warning
    Using a quadratic instead of a linear objective function makes the correction optmization much more computationally complex.

## The CorrectionConfig dataclass

Now, in COBRA-k, we can define a corrections scenario and the corrections options using the ```CorrectionConfig``` dataclass. It is defined as follows in COBRA-k's ```dataclass``` module:

```py
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
    """Indicates whether to add a km error term. Defaults to False."""
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
```

While many of the member variables are self-explanatory in the context of the previous sub-chapters, some member variables still need to looked at in more detail:

* ```error_scenario: dict[str, tuple[float, float]]```: This member variable describes the scenario for which you run the correction. E.g., if you have a scenario where the flux of a reaction A is measured to be between 1 and 2, ```error_scenario```would be set to ```{"A": (1, 2)}```.
* ```max_rel_(...): float``` and ```max_abs_(...): float```variables: With these member variables, you can restrict the maximally possible relative (for $k_{cat}⋅[E]$ and $k_M$) or  absolute (for $Δ_r G^{´°}$) maximal correction for a parameter.
* ```var_lb_ub_application: Literal["", "exp", "log"]```: Sometimes, you may have e.g. concentration data in non-logarithmic form. Then, if you set this member variable to ```"log"```, we can directly use this data for a concentration ```error_scenario```. Conversely, with```"exp"```, you could apply the exponential function on a measurement with logarithmic values. Keep in mind that this application is applied on *all* measurements, regardless of whether they represent concentrations or something else.

## Usage examples in code

Now, after all this theory and the dataclass explanation, let's see some toy model parameter corrections in COBRA-k :-)

!!! info
    The given examples are just two of many possible combinations (or non-combinations, if you just want a single type of corrections) of parameter corrections. Feel free to use the combination or single correction type that suits your problem best.

### A flux and concentration scenario correction in a MILP

Here, we try to find the minimal changes needed to the flux and concentration scenario
so that it becomes feasible. I.e. no model parameters are changes.

```py
# Import relevant classes and functions
from cobrak.lps import perform_lp_optimization
from cobrak.example_models import toy_model
from cobrak.dataclasses import CorrectionConfig
from cobrak.constants import LNCONC_VAR_PREFIX, ERROR_SUM_VAR_ID
from cobrak.printing import print_dict
from math import log

#
flux_and_concentration_error_scenario = {
    "Overflow": (1.0, 1.4),  # Overflow reaction flux between 1 and 1.4
    f"{LNCONC_VAR_PREFIX}M": (log(0.2), log(0.2)),  # M concentration fixed at .2 molar
    f"{LNCONC_VAR_PREFIX}D": (log(0.23), log(0.25)),  # D concentration betwewen .23 and .25 molar
}

# With a CorrectionConfig as optional further argument,
# all explained extra variables for parameter corrections
# are added to our model automatically :D
# Now, we can just minimize the sum of correction errors.
correction_result_1 = perform_lp_optimization(
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

print_dict(correction_result_1)
```


### A $k_{cat}⋅[E]$, $ΔG'°$ and $k_M$ correction in an NLP, and application of correction on model

Here, we try to find the minimal changes needed to the model's $k_{cat}⋅[E]$, $ΔG'°$ and $k_M$
such that the given scenario (a gigh Glycolysis flux) can be reached. We then apply (i.e. set the corrected parameters)
using ```apply_error_correction_on_model```from COBRA-k's ```utilities``` submodule.


```py
# Import relevant classes and functions
from cobrak.lps import perform_lp_optimization
from cobrak.example_models import toy_model
from cobrak.dataclasses import CorrectionConfig
from cobrak.constants import LNCONC_VAR_PREFIX, ERROR_SUM_VAR_ID
from cobrak.printing import print_dict
from cobrak.utilities import apply_error_correction_on_model

#
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

print_dict(correction_result_2)

# Now, we apply the correction (i.e. set the corrected
# parameter values to our model, overwriting the old parameter values)
corrected_cobrak_model = apply_error_correction_on_model(
    cobrak_model=toy_model,
    correction_result=correction_result_2,
    min_abs_error_value=0.01,
    min_rel_error_value=0.01,
    verbose=True,
)
```
