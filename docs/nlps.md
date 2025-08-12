# Nonlinear COBRA-k programs

??? abstract "Quickstart code"
    ```py
    from cobrak.example_models import toy_model
    from cobrak.lps import perform_lp_variability_analysis
    from cobrak.printing import print_optimization_result, print_variability_result
    # Import MINLP functionality in our *NLP* package
    # The "reversible" means that driving forces can become negative, but
    # reactions still have to be split as irreversible ones (v_i>=0)
    from cobrak.nlps import perform_nlp_reversible_optimization
    # Import NLP functionality in our NLP package
    # The "irreversible" means that driving forces can*not* become negative
    from cobrak.nlps import perform_nlp_irreversible_optimization_with_active_reacs_only

    #%
    # Run preparatory variability analysis
    variability_dict = perform_lp_variability_analysis(
        toy_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
    )
    # Pretty-print variability result
    print_variability_result(toy_model, variability_dict)

    #%
    # Run MINLP (by default, with the SCIP solver)
    minlp_result = perform_nlp_reversible_optimization(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption", # Let's maximize ATP production
        objective_sense=+1,
        # Set the variable bounds from our preparatory variability analysis
        variability_dict=variability_dict,
        # We use the saturation term constraint (otherwise, κ is set to 1 for all reactions);
        # default is True anyway
        with_kappa=True,
        # We use the thermodynamic term constraint (otherwise, γ is set to 1 for all reactions);
        # default is True anyway
        with_gamma=True,
    )
    # Pretty-print MINLP result
    print_optimization_result(toy_model, minlp_result)

    #%
    # Run (local and fast) NLP (by default, with the IPOPT solver)
    nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
        toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        # Set the suitable set of thermodynamically active reactions
        optimization_dict=minlp_result,
        # Set the variable bounds from our preparatory variability analysis
        variability_dict=variability_dict,
        # We use the saturation term constraint (otherwise, κ is set to 1 for all reactions);
        # default is True anyway
        with_kappa=True,
        # We use the thermodynamic term constraint (otherwise, γ is set to 1 for all reactions);
        # default is True anyway
        with_gamma=True,
    )

    print_optimization_result(toy_model, nlp_result)
    ```

## Introduction

### Reversible Michaelis-Menten kinetics

In the last chapters, we looked at Linear and Mixed-Integer Linear Optimization problems. As explained, these optimization problems are able to include stoichiometric, enzyme-pool and thermodynamic constraints. But they are not able to capture full reaction kinetics. One form of reaction kinetics is the reversible Michaelis-Menten kinetics, which is [[Paper]](https://doi.org/10.1016/j.febslet.2013.07.028):

$$ v_i = V^{+}_i ⋅ κ_i ⋅ γ_i $$

$v_i$ is, again, the reaction $i$'s flux. $V^{+}_i stands for the enzyme-dependent maximal flux that we already know from the basic enzyme constraints:

$$ V^{+}_i = E_i ⋅ k_{cat} $$

$κ_i$ is the saturation term - a unitless value that lies in $[0,1]$ - which is dependent on the enzyme's Michaelis-Menten constants ($K_M$) and is:

$$ κ_i = {{\bar{s}_i} \over {1 + \bar{s}_i + \bar{p}_i}} $$

${\bar{s}_i}$ a $K_M$-dependent product of the concentrations of the reaction's *substrates*:

$$ {\bar{s}_i} = \prod_{j \in 𝖲_i} (c_j ⋅ K_{M,j,i}) ^ {|N_{i,j}|}$$

$K_{M,j,i}$ is the Michaelis-Menten constant of this substrate for this reaction, and $𝖲_i$ is the set of all indices of reaction $i$'s substrates. $c_s$ is the concentration of substrate $j$ and $|N_{i,j}|$ the absolute stoichiometry of the metabolite.

Analogously, ${\bar{p}_i}$ a value affected by the reaction's *products*:

$$ {\bar{p}_j} = \prod_{k \in 𝖯_i} (c_k ⋅ K_{M,k,i}) ^ {N_{i,k}}$$

where $k$ is the index of the product, taken from the set of reaction product indices $𝖯_i$.

The last part of the reversible Michaelis-Menten kinetics, $γ_i$, is the thermodynamic term. It is also a unitless value in $(-∞,1]$:

$$ γ_i = ( 1 - e^{-f_i / (R ⋅ T)} ) $$

As explained in the previous chapter, R is the Gas constant, T the temperature and $f_i$ the reaction's driving force.

In conclusion, $V^{+}_i$ determines the maximally possible reaction flux. Both $κ_i$ and $γ_i$ can only restrict this maximally possible flux. This is because $κ_i$ and $γ_i$ can only lie in $[0,1]$ and $(-∞,1]$, respectively. Thereby, if we only allow positive driving forces ($f_i>0$), $γ_i$ is even restricted to $(0,1]$. Also, both $κ_i$ and $γ_i$ follow the general idea of reaction kinetics: The higher the product concentrations, the lower the reaction's flux. Conversely, the higher the substrate concentrations, the higher the flux.

### (Mixed-integer) Nonlinear Programming

!!! MINLP and NLP

## Nonlinear kinetic constraints

Now that we know the formulas of the Michaelis-Menten kinetics, we want to efficiently integrate them in our framework of constriant-based modeling. From CBM, we are at least still using:

* the steady-state constraints (see LP chapter)
* extra linear constraints (see LP chapter)
* the logarithmic concentrations and their bounds (see MILP chapter)
* the reaction driving forces $f_i$ (see MILP chapter)

For our further CBM integration, "efficently integrating the kinetic formulas" means that we relax the kinetic formulas as much as possible. This means that we do not treat them as equality but as *in*equality:

$$ v_i ≤ V^{+}_i ⋅ κ_i ⋅ γ_i $$

Now, $v_i$ is allowed to become lower than what would be expected from the kinetics. However, we are also still using the central protein pool constrained introduced in the LP chapter:

$$ ∑_i W_i ⋅ E_i ≤ Ω $$

This means that, in a typical optimziation, the enzymes are still going to be used as efficiently as possible. I.e., typically, the lowest needed amount of enzymes is found, such that our inequality often becomes an equality :-) Conversely, in cases where this does not hold, the enzyme usage is not the major constraint for our optimization anyway.

Following our relaxation scheme, we also relax the maximally possible flux $V^{+}_i$ as inequality:

$$ V^{+}_i ≤ E_i ⋅ k_{cat} $$

$ κ_i $ is also now an inequality and made dependent on our *logarithmic* concentrations as follows:

$$ κ_i <= {e^{\tilde{s}_i} \over {( 1+e^{\tilde{s}_i}+e^{\tilde{p}_i} )}} $$

where $\tilde{s}_i$ and $\tilde{p}_i$ are the logarithmic variants of $\bar{s}_i$ and $\bar{p}_i$ (see above) and are using our logarithmic concentration vector $\mathbf{\tilde{x}}$ (see previous chapter):

$$ \tilde{s}_i = \ln ( \bar{s} ) = \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅ x̃_j ) - \sum_{j ∈ 𝖲_i} ( |N_{j,i}| ⋅  \ln K_{M,j,i} ) $$

$$ \tilde{p}_i ≥ \ln ( \bar{p} ) = \sum_{k ∈ 𝖯_i} ( N_{k,i} ⋅ x̃ ) - \sum_{k ∈ 𝖯_i} ( N_{k,i} ⋅  \ln K_{M,k,i} ) $$

Note that $\tilde{p}$ is even further relaxed, as lower $\tilde{p}$ could only restrict a reaction flux even further (see formula for $κ_i$).

$$ γ_i ≤ (1-e^{f_i / (R ⋅ T)})$$

whereby we also relax $f_i$, again in a direction which could only lower the flux, as (using the definitions from the previous chapter)

$$ f_i ≤ -(Δ_r G^{´°}_i + R ⋅ T ⋅ \mathbf{N_{⋅,i}} ⋅ \mathbf{x̃}) $$

And that's it :D With these additional constraints, on top of our mentioned constrained-based constrained that we introduced in earlier chapters, we could now run constraint-based analyses with full reaction kinetics :-)

### Optional concentration sum constraints

COBRA-k also provides the possibility to introduce *concentration sum constraints*. They are only activated if a Model's ```max_conc_sum``` member variable is smaller than the default value ```float("inf")```. The concentration sum constraint is simply

$$ Φ ≤ \sum{e^(x̃_j)} $$

!!! note
    We don't need any of the linear approximation tricks used for MILPs (see last chapter) here :-)

!!! warning
    Adding concentration sum constraints can cause a heavy load on a non-linear solver. Hence, if you are expereiencing very slow solving times, they might be caused by this constraint if you have set ```max_conc_sum``` to a value lower than its default ```float("inf")```.

## Preparatory Variability Analysis

One major problem with our non-linear constraints is that certain values may become very large or very small. In particular, $e^{\tilde{s}_i}$, $e^{\tilde{p}_i}$ may become way too large and $γ_i$ way too small if $f_i$ is too negative. This would lead to solver errors, making optimization impossible.

Luckily, we still use our pre-set logarithmic concentration bounds (see previous chapter) which typically restrict all our values to solver-friendly sizes. Hence, running a MILP-based ecTFVA (see previous chapter) for all linear values in our problem, namely the...

* reaction fluxes $\mathbf{v}$
* enzyme concentrations $\mathbf{E}$
* logarithmic concentrations $\mathbf{\tilde{x}}$
* driving forces $\mathbf{f}$
* logarithmic saturation values $\mathbf{\tilde{s}}$ and $\mathbf{\tilde{p}}$

...solves our problem with too high and too low values. Mathematically, the ecTFVA-based preparatory variability analysis can be written as

$$ \operatorname*{\mathbf{min}}_{\mathbf{v, E, x̃, f, z, \tilde{p}, \tilde{s}, κ, γ}, B}  \mathbf{β_i, β ∈ \{ \mathbf{v}, \mathbf{E}, \mathbf{\tilde{x}}, \mathbf{f}, \mathbf{\tilde{s}}, \mathbf{\tilde{p}} \}} \\ s.t. \space CBM \space \& \space thermodynamic \space \& \space saturation \space term \space  constraints $$

and

$$ \operatorname*{\mathbf{max}}_{\mathbf{v, E, x̃, f, z, \tilde{p}, \tilde{s}, κ, γ}, B}  \mathbf{β_i, β ∈ \{ \mathbf{v}, \mathbf{E}, \mathbf{\tilde{x}}, \mathbf{f}, \mathbf{\tilde{s}}, \mathbf{\tilde{p}} \}} \\ s.t. \space CBM \space \& \space thermodynamic \space \& \space saturation \space term \space  constraints $$


In COBRA-k, we can run (and let us pretty-print) such a preparatory variability analysis analogously as how we did it in the previous chapter:

```py
# Import our known toy model (see LP chapter) and methods
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_variability_analysis
from cobrak.printing import print_variability_result

# Run preparatory variability analysis
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
)
# Pretty-print variability result
print_variability_result(toy_model, variability_dict)
```

## Global MINLP (generally slow)

!!! warning
    While this MINLP can lead to great global result, it is typically *way* too slow and not suitable even for mid-scale metabolic models with more than a dozen reactions, even with the relaxations introduced above. For alternatives, see the next subchapter as well as the evolutionary optimization in the next chapter.

    Also, due to the non-linearity of our new constraints, only solvers which are capable of handling such constraints in pyomo can be used. E.g., CPLEX and Gurobi (at least older Gurobi versions) cannot be used with non-linear optimizations.

To *globally* optimize on a full metabolic model using our non-linear constraint, we can use COBRA-k's MINLP functionality. It is defined as

$$ \operatorname*{\mathbf{min}}_{\mathbf{v, E, x̃, f, z, \tilde{p}, \tilde{s}, κ, γ}, B}  \mathbf{g^\top} \\ s.t. \space CBM \space \& \space thermodynamic \space \& \space nonlinear \space  constraints $$

and lets us globally optimize any value introduced in our linear, mixed-integer and non-linear constraints.

In COBRA-k, we can run a MINLP on toy model (which is small enough for it to run) as follows, whereby COBRA-k's ```nlps``` subpackage is used:

```py
# Import toy model and pretty-print function and (ecT)FVA function
from cobrak.example_models import toy_model
from cobrak.printing import print_optimization_result
from cobrak.lps import perform_lp_variability_analysis
# Import MINLP functionality in our *NLP* package
# The "reversible" means that driving forces can become negative, but
# reactions still have to be split as irreversible ones (v_i>=0)
from cobrak.nlps import perform_nlp_reversible_optimization

# Run preparatory variability analysis
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
)

# Run MINLP (by default, with the SCIP solver)
minlp_result = perform_nlp_reversible_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption", # Let's maximize ATP production
    objective_sense=+1,
    # Set the variable bounds from our preparatory variability analysis
    variability_dict=variability_dict,
    # We use the saturation term constraint (otherwise, κ is set to 1 for all reactions);
    # default is True anyway
    with_kappa=True,
    # We use the thermodynamic term constraint (otherwise, γ is set to 1 for all reactions);
    # default is True anyway
    with_gamma=True,
)
# Pretty-print MINLP result
print_optimization_result(toy_model, minlp_result)
```

!!! info
    Just like all other optimization functions, ```perform_nlp_reversible_optimization``` has many other optional arguments, including the possibility to use other solvers, solver and pyomo solve funtion options. For more information, see this documentation's "API reference".


## Local NLP (fast, but restricted)

!!! warning
    The NLP shown here is fast, but *only* suitable if you are working with a selected metabolic set of reactions that can be thermodynamically feasible at the same time (not a full model). But this NLP is the basis of the evolutionary algorithm, which is shown in the next chapter and lets you quasi-globally optimize a full metabolic model with non-linear constraints.

As the MINLP can be very slow, COBRA-k also provides a fast NLP that works on a single set of reactions that can be thermodynamically feasible at the same time. I.e., the following constraint holds:

$$ f_i ≥ f^{min} $$

whereby $f^{min}$ has to be positive. This means that no full model can be used, as it typically includes (e.g. reverse direction) reactions with a negative driving force. But once you have a suitable set of reactions, you can quickly solve the following NLP (now without any binary variables, in contrast to the MINLP):

$$ \operatorname*{\mathbf{min}}_{\mathbf{v, E, x̃, f, \tilde{p}, \tilde{s}, κ, γ}, B}  \mathbf{g^\top} \\ s.t. \space CBM \space \& \space nonlinear \space constraints $$

In COBRA-k, you can find a suitable set of thermodynamically feasible reactions with an ecTFBA first, and *then* calculate the NLP:

```py
# Import toy model and pretty-print function and pathway-using ecTFBA
from cobrak.example_models import toy_model
from cobrak.printing import print_optimization_result
from cobrak.lps import perform_lp_optimization, perform_lp_variability_analysis
# Import NLP functionality in our NLP package
# The "irreversible" means that driving forces can*not* become negative
from cobrak.nlps import perform_nlp_irreversible_optimization_with_active_reacs_only

# Run preparatory variability analysis for our NLP
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
)

# Find suitable set of thermodynamically active reactions
ectfba_result = perform_lp_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    # Set enzyme constraints as they are also used in the NLP
    with_enzyme_constraints=True,
    # This following setting is important to find thermodynamically active reactions
    with_thermodynamic_constraints=True,
)

# Run (local and fast) NLP (by default, with the IPOPT solver)
nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
    toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    # Set the suitable set of thermodynamically active reactions
    optimization_dict=ectfba_result,
    # Set the variable bounds from our preparatory variability analysis
    variability_dict=variability_dict,
    # We use the saturation term constraint (otherwise, κ is set to 1 for all reactions);
    # default is True anyway
    with_kappa=True,
    # We use the thermodynamic term constraint (otherwise, γ is set to 1 for all reactions);
    # default is True anyway
    with_gamma=True,
)

print_optimization_result(toy_model, nlp_result)
```

Again, ```perform_nlp_reversible_optimization``` has many further optional arguments that you can find in this documentation's "API reference" chapter. By default, COBRA-k uses the non-linear solver IPOPT for NLPs.

### Extra non-linear flux constraints

Optionally, you can also introduce extra *non-*linear constraints (corresponding to the ```ExtraNonlinearConstraint``` dataclass, used in ```Model```) that set constrained relationships between variables. Currently, the non-linear functions "powerX" (i.e. take the X-th power of a value), "log" and "exp" are usable. "same" is the option if you do not want a function application (i.e. it's a multiplication with 1). See COBRA-k's API documentation for more. As an example, let's set a (nonsense) non-linear constraint that restricts the flux of reaction EX_P to the doubled exponential value of the concentration of C:

```py
# ...using the code imports from above
from cobrak.constants import LNCONC_VAR_PREFIX
from cobrak.dataclasses import ExtraNonlinearConstraint

# Let's define v_EX_P <= 2 * exp(x_C)
toy_model.extra_nonlinear_constraints = [
    ExtraNonlinearConstraint(
        stoichiometries={
            "EX_P": (1.0, "same"),  # first the stoichiometry, second the function application
            f"{LNCONC_VAR_PREFIX}C": (-2.0, "exp"),
        },
        upper_value=0.0,
    )
]
```

### Extra linear watches

Optionally, you can also introduce extra linear *watch variables* (corresponding to the ```ExtraNonlinearWatch``` dataclass, used in ```Model```) that add a variable with a fixed non-linear relationships to single fluxes. See COBRA-k's API documentation for more. Currently, the non-linear functions "powerX" (i.e. take the X-th power of a value), "log" and "exp" are usable. See COBRA-k's API documentation for more. Here's an example where we set a watch to the logarithm of the flux of EX_S:

```py
# ...using the code imports from above
from cobrak.dataclasses import ExtraNonlinearWatch

# Let's define v_EX_P <= 2 * exp(x_C)
toy_model.extra_nonlinear_watches = {
    "log_EX_S": ExtraNonlinearWatch(
        stoichiometries={
            "EX_S": (1.0, "log"),
        },
    )
}
```

...now, we have a variable ```log_EX_S``` that is also added to result dictionaries after LP or NLP optimizations.
