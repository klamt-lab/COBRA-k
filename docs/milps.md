# Mixed-integer linear programs

??? abstract "Quickstart code"
    ```py
    #% Thermodynamic Flux Balance Analysis
    from cobrak.example_models import toy_model
    from cobrak.lps import perform_lp_optimization
    from cobrak.printing import print_dict, print_optimization_result

    # Perform TFBA
    tfba_result = perform_lp_optimization(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        with_thermodynamic_constraints=True,
    )

    # Pretty print result as dictionary
    print_dict(tfba_result)

    # Pretty print result as tables
    print_optimization_result(toy_model, tfba_result)


    #% OptMDFpathway
    from cobrak.constants import MDF_VAR_ID

    # Perform OptMDFpathway
    optmdfpathway_result = perform_lp_optimization(
        cobrak_model=toy_model,
        objective_target=MDF_VAR_ID,
        objective_sense=+1,
        with_thermodynamic_constraints=True,
    )

    # Print MDF only
    print(optmdfpathway_result[MDF_VAR_ID])


    #% Thermodynamic Bottleneck Analysis
    from cobrak.lps import perform_lp_thermodynamic_bottleneck_analysis

    # Use model with extreme standard Gibbs free energy and enforced ATP production
    with toy_model as tba_model:
        tba_model.reactions["Glycolysis"].dG0 = 100 # A bottleneck :O

        # Perform TBA
        list_of_bottleneck_reactions, _ = perform_lp_thermodynamic_bottleneck_analysis(
            tba_model,
        ) # The second returned value is the full solution (with fluxes, concentrations, ...) which we don't need here

    # Print list of thermodynamic bottlenecks
    print(list_of_bottleneck_reactions)


    #% Thermodynamic Variability Analysis
    from cobrak.example_models import toy_model
    from cobrak.lps import perform_lp_variability_analysis
    from cobrak.printing import print_variability_result

    # Perform TFVA
    variability_dict = perform_lp_variability_analysis(
        toy_model,
        with_thermodynamic_constraints=True,
    )

    # Pretty print result as tables
    print_variability_result(toy_model, variability_dict)


    #% Perform enzyme-constrained TFBA
    # Run ecTFBA
    ectfba_result = perform_lp_optimization(
        cobrak_model=toy_model,
        objective_target="ATP_Consumption",
        objective_sense=+1,
        with_thermodynamic_constraints=True,
        with_enzyme_constraints=True,
    )

    # Pretty print result as dictionary
    print_dict(ectfba_result)

    # Pretty print result as tables,
    # now with enzyme concentrations *and* metabolite concentrations
    print_optimization_result(toy_model, ectfba_result)
    ```

## Introduction

### MILPs

We used plain Linear Programming in the last chapter, where every variable is *rational* (which includes not only whole numbers (such as $2$) but also fractions (such as $7.31$)). Now, we use concepts based on Mixed-Integer Linear Programming (MILP), where a user-selected set of binary variables which can be *either* 0 or 1, nothing in-between.

MILPs are much more difficult to solve than LPs. While the latter may contain hundreds of thousands of parameters and still be quickly solved, MILPs are restricted, at very best, to a few thousand parameters. Fortunately, that's still good enough for large genome-scale metabolic models.

??? excursion "Excursion: MILPs"
    Based on our definitions from the previous chapter, a general form of a maximizing MILP is:

    $$
    \operatorname*{\mathbf{maximize}}_{\mathbf{x}} \ \mathbf{g}ᵀ  \\
    subject \ to \ the \ constraints \\
    \mathbf{A} ⋅ \mathbf{x} ≤ \mathbf{b} \\
    x^{min}_i < x_i < x^{max}_i \\
    x_j ∈ \{0,1 \} ∀ j∈1...i^{max}
    $$

    where j stands for all indices (out of all variables indices indicated by i) whose variable is restricted to the value 0 and 1.

    Major algorithms for the solution of the selection of binary variables in MILPs are branch-and-bound [[Paper]](https://doi.org/10.2307/1910129)[[Wikipedia]](https://en.wikipedia.org/wiki/Branch_and_bound) and the cutting-plane method [[Paper]](https://doi.org/10.1287/opre.9.6.849)[[Wikipadia]](https://en.wikipedia.org/wiki/Cutting-plane_method) as well as their combination branch-and-cut [[Paper]](https://doi.org/10.1137%2F1033004)[[Wikipedia]](https://en.wikipedia.org/wiki/Branch_and_cut).

### Thermodynamic measures

Binary variables allow us to introduce *thermodynamic* constraints in our constraint-based model. Thermodynamic constraints make sure that our solution is thermodynamically feasible. This means that the concentration(s) of any active reaction's substrate are somehow high enough in comparison to the reaction's products. In exact form, whether this neccessary substrate(s)-to-product(s) concentration ratio is reached can be deduced with the driving force $f_i$ (side note: this is the negative Gibbs energy $Δ_r G^{´}$ [[Wikipedia]](https://en.wikipedia.org/wiki/Gibbs_free_energy)) with the unit kJ⋅mol⁻¹. $f_i$ is for a reaction $i$:

$$ f_i = -Δ_r G^{´°}_i + R ⋅ T ⋅ Q_i $$

R stands for the Gas constant in kJ⋅K⁻¹⋅mol⁻¹ [[Wikipedia]](https://en.wikipedia.org/wiki/Gas_constant) and T for the temperature in K; $Q_i$ is explained further below.

$Δ_r G^{´°}_i$ (with the ° at the end) is the physiologic *standard* Gibbs energy. Just like $f_i$, it has the unit kJ⋅mol⁻¹ and stands for the amount of energy (in kJ) that is used when the reaction converts 1 mole of substrates into 1 mole of products, whereby *all* substrates and products have a concentration of 1 M=1 mol⋅l⁻¹. Furthermore, to account for physiological effects, the reaction compartment's pH, pMg and ionic strength are accounted for in the complicated calculations for the determination of a $Δ_r G^{´°}_i$. This becomes even more complicated for multi-compartmental reactions. For more details, look at this fantastic paper from the authors of the fantastic eQuilibrator-API [[here]](https://doi.org/10.1093/nar/gkab1106) and the excellent eQuilibrator FAQ [[here]](https://equilibrator.weizmann.ac.il/static/classic_rxns/faq.html). This (COBRA-k-unaffiliated) API is also used by COBRA-k's functionality to retrieve $Δ_r G^{´°}_i$ of reactions, as detailed in chapter 10.

The main value of the standard Gibbs energy is its following meaning:

* If $Δ_r G^{´°}_i < 0$, it means that energy is *released* when all substrates and products have the standard concentration of 1 M. As every working reaction needs energy to be released, this means that under standard concentrations, the reaction would be thermodynamically *feasible*. In other works, it could run.

* If $Δ_r G^{´°}_i ≥ 0$, *no* energy is released at standard concentrations. This would mean that with the given standard concentrations, the reaction would be thermodynamically *infeasible*, it could not run.

!!! warning
    The driving force $f_i$ uses the *negated* sign logic as it negates the Gibbs energy. Also, as explained below with the factor $Q_i$, $f_i$ can be used for any metabolite concentrations, not only 1 M.

$f_i$ means the following:

* The *higher* the $f_i$, the more energy is released by a reaction.

* The *lower* the $f_i$, the less energy is released by a reaction.

⇒ if $f_i>0$, a reaction is thermodynamically *feasible*; if $f_i≤0$, a reaction is thermodynamically *in*feasible

$Q_i$ is the appropriately formulated thermodynamic logarithmic concentration ratio between substrates and products:

$$ Q_i = ∑_{j}( -N_{i,j} ⋅ ln(c_j) )$$

$ln$ is the natural logarithm, $c_j$ the concentration of metabolite $j$ and $N_{i,j}$ is its stoichiometry (negative for substrates, positive for products) in the stoichiometric matrix $\mathbf{N}$. The influence of substrate and product concentrations on $Q_i$, and therefore $f_i$, is as follows:

* Through the stoichiometry-negating factor $-N_{i,j}$, $Q_i$ is defined such that is becomes lower the higher the substrate concentration(s) are. In other words, the driving force $f_i$ of a reaction becomes "better" (it goes in the direction of thermodynamic feasibility, i.e. higher) the higher the the substrate concentration(s) are.

* In contrast, $Q_i$ becomes higher the higher the product concentration(s) are, bringing the reactions more towards an infeasible state as the $f_i$ also becomes lower.


??? "A toy reaction example"
    To illustrate the (often complicated looking) measures introduced in the previous paragraphs, let's look at a simple reaction called X:

    $ A → 2 B $

    Lets say that this reaction's standard Gibbs energy $Δ_r G^{´°}_{X}$ is -1 kJ⋅mol⁻¹ A's concentration ($c_A$) is 0.2 M, and $c_B§ is 0.1 M.

    Then, the driving force $f$ of reaction X is:

    $ f_{X} = Δ_r G^{´°} + -N_{X,A} ⋅ ln(c_A) + -N_{X,B} ⋅ ln(c_B) $

    $ = -1 + -1 ⋅ ln(0.2) + -(-2) ⋅ ln(0.1) ≈ -3.996 kJ⋅mol⁻¹ $

    $f_X<0$, this, the reaction is thermodynamically feasible, it can run!

## Thermodynamic constraints

*Aim:* Using the thermodynamic measured introduced above, we can now integrate thermodynamic constraints into constraint-based modeling (CBM). More precisely, we want that the driving force $f_i$ of any active reaction is positive. In other words, whenever the flux of a reaction is greater than zero ($v_i>0$), $f_i>0$ holds. As a first step to do so, we add metabolite concentrations for any metabolite $j$:

$$ \ln (c_j^{min} ) ≤ x̃_j ≤ \ln ( c_j^{max} ) $$

$c_j^{min}$ is the minimal metabolite concentration in M, $c_j^{max}$ the respective maximal concentration in M. $x̃_j$, an element of the vector $\mathbf{x̃}$, stands for a variable that holds the *logaritmic* concentration of metabolite $j$.

!!! info "How to set concentration ranges"
    Typical standard minimal and maximal concentrations are 10⁻⁶ M up to 0.02 M for intracellular metabolites, with higher maximal concentrations for extracellular metabolites. The concentrations for water ($H_2 O$) and protons ($H^+$) are often set to 0 M (i.e. their logarithm is 1). This is because both water (with a fixed "active" concentration) and proton (as pH) concentrations are integrated in the calculation of the physiologic standard Gibbs energy $Δ_r G^{´°}$. For more details of how e.g. the fantastic eQuilibrator does this, read their FAQ [here](https://equilibrator.weizmann.ac.il/static/classic_rxns/faq.html#why-can-t-i-change-the-concentration-of-water).

Now, we can formulate the driving forces $f_i$ - for any reaction $i$ for which we want thermodynamic constraints - as follows:

$$ f_i = -(Δ_r G^{´°}_i + R ⋅ T ⋅ \mathbf{N_{⋅,i}} ⋅ \mathbf{x̃}) $$

Again, $Δ_r G^{´°}_i$ is the reaction's standard Gibbs energy, R the gas constant, T the temperature. The term $\mathbf{N_{⋅,i}} ⋅ \mathbf{x̃}$ is equivalent to the formulation of $Q_i$ above and means the following: We take the stoichiometries of all (⋅) metabolites in reaction i through $\mathbf{N_{⋅,i}}$, i.e. we take the $i$-th row of the stoichiometric matrix $\mathbf{N}$. Then, we multiply these stoichiometries with the logarithmic concentration vector $\mathbf{x̃}$. This effectively gives us - as for $Q_i$ - the sum of stoichiometries multiplied with the logarithmic concentrations.

While we now have the driving force, we did not enforce it yet to be positive (i.e. to indicate feasiblity). For this, we also introduce a controlling binary variable vector $\mathbf{z}$ which holds a binary value for any reaction $i$ with thermodynamic constraints:

$$ z_i ∈ \{0,1 \} $$

Any $z_i$ must be 1 if a reaction wants to run through the constraint:

$$ v_i ≤ ub_i ⋅ z_i $$

Also, any $z_i$ must be 0 (thus making a reaction not run) if $f_i ≤ 0$:

$$ B ≤ f_i + M ⋅ ( 1-z_i ) $$

Thereby, $B$ is a variable that stands for a set *minimal* driving force that has to be reached by any active reaction:

$$ B ≥ f^{min} $$

The constant $f^{min}>0$ is the lower bound for $B$. $f^{min}>0$ must hold if we want to ensure thermodynamic feasibility.

And that's it! Through the two constraints utilizing $z_i$, we ensure that any thermodynamically infeasible reaction is inactive (its flux is 0), while any thermodynamically feasible reaction may be active.

!!! note
    Through our formulation over all reactions, we ensure the *network-wide* thermodynamic feasiblity of any solution with thermodynamic constraints as long as we set $f^{min}>0$.

### Optional concentration sum constraints

COBRA-k also provides the possibility to introduce *concentration sum constraints*. They are only activated if a Model's ```max_conc_sum``` member variable is smaller than the default value ```float("inf")```. In exact form, the concentration constraint would be

$$ Φ ≤ \sum{e^(x̃_j)} $$

whereby $Φ$ stands for the maximal concentration sum we set, and $e^(x̃_j)$ for a exponentiated logarithmic concentration. As $e^(x̃_j)$ is *non*-linear, we cannot use them directly in our MILP. Hence, we need a linearized approximation (whereby we use most of the formulation from [[this preprint](https://doi.org/10.1101/2024.03.19.585265)]).

This works as, luckily, the exponential function is monotonically rising :D This means that we can always draw a "minimum" linear constraint underneath the exponential function's curve without cutting this curve. Even more lucky, we have to set $x̃_j$ concentration bounds anyway for thermodynamic constraints (see above), so that we know for which range of logarithmic concentrations we apply the exponential function. I.e. we know the possible minimal and maximal logarithmic concentration and only have to approximate the exponential function for these values.

Now, mathematically, the exponential function's linear approximation is built as follows:

* For each metabolite $j$, we first look at its possible minimal and maximal logarithmic concentration. Then, we segment this interval into $ξ$ many values creating the value vector $\mathbf{S}$ (we will clarify how large $X$ should be, i.e. how many segments are needed). I.e., we get X many equally distant values $\mathbf{S}$ from (and including) the minimal and maximal logarithmic concentration of $j$.

* For each segment value, we approximate its delogarithmic concentration $C_j$, i.e. $e^(x̃_j)$ a linear constraint of the following form (for a visualization see [[Fig. 5 in this preprint]](https://doi.org/10.1101/2024.03.19.585265)):

$$C_j ≥ a_{s,j} ⋅ x̃_j + b_{s,j}$$

$s$ is the linear approximation segment index in $[1,S]$. $a_{s,j}$ is the slope of the exponential function at the segment's value. $b_{s,j}$ is the necessary correction such that $a_{s,j} ⋅ x̃_j + b_{s,j} = e^(x̃_j)$. It is defined as:

$e^(x̃_j) - e^(x̃_j) * S_s$

where $S_s$ is the segment $s$'s value.

* We still didn't clarify how many segments are needed, i.e. how large $ξ$ has to be. To determine $ξ$, we use the following algorithm:

1. We start with just $ξ=2$, i.e. two segments
2. For the current number of $ξ$, and define its two linear approximations as described above.
3. Then, for each segment value $S_s \space ∀ \space s<ξ$, we calculate the *Y intersect* of the linear approximation function at $S_s$ and $S_{s+1}$. I.e., this is the X value where this linear approximation and the linear approximation of the *following* segment have the same Y value. We do this because, due to the monotonically rising nature of the exponential function, this intersect is the point with the *largest possible* approximation error (again, for a visualization see [[Fig. 5 in this preprint]](https://doi.org/10.1101/2024.03.19.585265)).
4. We now calculate the relative difference between the exponential function and the linear approximation at this Y intersect. If this relative difference is larger than the set ```conc_sum_max_rel_error``` (default: ```0.05```), we increase $ξ$ by 1 and start again with step 2 of this algorithm. Otherwise, we have a linear approximation within our set maximal relative error and can stop the algorithm :D

* Finally, we can constrain our sum of the approximation of delogarithmic concentrations as follows:

$$ Φ ≤ \sum{C_j} $$

!!! warning
    For numeric reasons, linear approximations at very low concentrations may cause numeric problems :-( To mitigate this, we can set ```conc_sum_min_abs_error``` (default: ```1e-6```) in our Model instance, which sets a minimal value for $Φ$ and removing all linear constraints for any logarithmic concentration where $e^(x̃_j)$ is smaller than ```conc_sum_min_abs_error```.

## Thermodynamic Flux Balance Analysis (TFBA)

TFBA combines Flux Balance Analysis (FBA) with thermodynamic constraints to ensure that the predicted flux distribution is *both* stoichiometrically and thermodynamically feasible [[Paper]](https://doi.org/10.1529/biophysj.106.093138). The resulting optimization problem is (also using the FBA definitions from the previous chapter):

$$ \operatorname*{\mathbf{min}}_{\mathbf{v}, \mathbf{E}, \mathbf{x̃}, \mathbf{f}, B}  \mathbf{g^\top} \\ s.t. \space CBM \space \& \space thermodynamic \space constraints $$

In COBRA-k, we can run a TFBA as follows utilizing the toggle ```with_thermodynamic_constraints```:

```py
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_optimization
from cobrak.printing import print_dict, print_optimization_result

# Perform TFBA
tfba_result = perform_lp_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    with_thermodynamic_constraints=True,
)

# Pretty print result as dictionary
print_dict(tfba_result)

# Pretty print result as tables
print_optimization_result(toy_model, tfba_result)
```

Note that the printed results now also show reaction driving forces and metabolite concentrations in a solution.

!!! info Metabolite legarithmic concentration variables
    In COBRA-k result dictionaries, logarithmic metabolite concentration variables start with ```cobrak.constant.LNCONC_PREFIX```, by default ```x_```. E.g., the logarithmic metabolite concentration of metabolite ```atp_c``` (cytosolic ATP) would be called ```x_atp_c```.

!!! warning
    Just like in a non-thermodynamic FBA, there may be up to infinite many alternative solutions for the same optimal value. For an elucidation of the full solution space, you may use a TFVA (see below).

    Also, to make real use of thermodynamic constraints, your reactions need associated ```dG0``` values. For a way to automatically calculate them, see chapter 10.

## OptMDFpathway

Sometimes, it is interesting to find out the max-min driving force of a network, also called MDF. This MDF stands for the maximized minimal driving force of all active reactions. In other words, it stands for the highest possible minimal $f_i$ in a network. We can easily formulate this optimization problem, known as OptMDFpathway [[Paper]](https://doi.org/10.1371/journal.pcbi.1006492), as:

$$ \operatorname*{\mathbf{max}}_{\mathbf{v}, \mathbf{x̃}, \mathbf{f}, B}  B \\ s.t. \space CBM \space  \& \space thermodynamic \space constraints $$

In COBRA-k, the variable $B$ has the constant name ```cobrak.constants.MDF_VAR_ID```, by default ```"var_B"```. With our toy model, we can perform OptMDFpathway as follows:

```py
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_optimization
from cobrak.constants import MDF_VAR_ID

# Perform OptMDFpathway
optmdfpathway_result = perform_lp_optimization(
    cobrak_model=toy_model,
    objective_target=MDF_VAR_ID,
    objective_sense=+1,
    with_thermodynamic_constraints=True,
)

# Print MDF only
print(optmdfpathway_result[MDF_VAR_ID])
```

## Thermodynamic Bottleneck Analysis (TBA)

Often, reactions such as growth are thermodynamically infeasible in our model with a given set of $Δ_r G^{´°}_i$ and concentration ranges. This is typically caused by thermodynamic *bottleneck* reactions, i.e. reactions whose $Δ_r G^{´°}_i$ is so high that the equation

$$ f^{min} > 0 $$

cannot hold :-( To solve this problem, COBRA-k includes the thermodynamic bottleneck analysis (TBA) introduced [[here]](https://doi.org/10.1038/s41467-023-40297-8). It searches for one (not necessarily unique) *minimal* number of reactions whose $Δ_r G^{´°}_i$ has to be reduced to make a model state thermodynamically feasible. As a MILP, it introduces a new vector of binary variables $\mathbf{z^b}$ for each reaction $i$ with thermodynamic constraints whereby

$$ z^{B}_i ∈ \{0,1 \} $$

The $z^{B}_i$ are now integrated in the reformulated constraint that controls whether a reaction is thermodynamically feasible or not (see above for the formulation without $z^{B}_i$ ):

$$ f_i + M ⋅ \left( 1-z_i \right) + M ⋅ z_i^B ≥ B $$

Essentially, if $z_i^B=1$, then $z_i$ can be zero even though the reaction is thermodynamically infeasible. We can now identify a minimal set of reactions whose $Δ_r G^{´°}_i$ has to be reduced through

$$ \operatorname*{\mathbf{min}}_{\mathbf{v}, \mathbf{x̃}, \mathbf{f}, B}  ∑(z^{B}_i) \\ s.t. \space CBM \space  \& \space thermodynamic \space \& \space TBA \space constraints $$

In COBRA-k, we can perform this analysis as follows:

```py
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_thermodynamic_bottleneck_analysis

# Use model with extreme standard Gibbs free energy and enforced ATP production
with toy_model as tba_model:
    tba_model.reactions["Glycolysis"].dG0 = 100 # A bottleneck :O

    # Perform TBA
    list_of_bottleneck_reactions, _ = perform_lp_thermodynamic_bottleneck_analysis(
        tba_model,
    ) # The second returned value is the full solution (with fluxes, concentrations, ...) which we don't need here

# Print list of thermodynamic bottlenecks
print(list_of_bottleneck_reactions)
```

!!! info "Using with statement"
    In our code snippet, we used Python's ```with``` with a ```Model``` instance. This is a handy way to automatically create a deep copy of our model, which does not affect the original model.

## Thermodynamic Variability Analysis (TVA), including Concentration Variability Analysis (CVA)

Analogously to the Flux Variability Analysis shown in the previous chapter, we can now perform a general Thermodynamic Variability Analysis (TVA). Now, we may not only minimize and maximize fluxes, but also logarithmic concentrations:

$$ \operatorname*{\mathbf{min}}_{\mathbf{v},\mathbf{x̃}, \mathbf{f}, B}  β_i, β_i∈\mathbf{v}∪\mathbf{x} \\ s.t. \space CBM \space \& \space thermodynamic \space constraints $$

$$ \operatorname*{\mathbf{max}}_{\mathbf{v}, \mathbf{x̃}, \mathbf{f}, B}  β_i, β_i∈\mathbf{v}∪\mathbf{x} \\ s.t. \space CBM \space \& \space thermodynamic \space constraints $$

We can recover the actual minimal and maximal concentrations through exponentiation:

$$ c_i^{min}=e^{x̃_i^{min}} $$

$$ c_i^{max}=e^{x̃_i^{max}} $$

In COBRA-k, a TVA looks as follows:

```py
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_variability_analysis
from cobrak.printing import print_variability_result

# Perform TFVA
variability_dict = perform_lp_variability_analysis(
    toy_model,
    with_thermodynamic_constraints=True,
)

# Pretty print result as tables
print_variability_result(toy_model, variability_dict)
```

## Combining enzyme and thermodynamic constraints

Of course, you can combine thermodynamic constraints with enzyme constraints (shown in the last chapter), which e.g. makes a TFBA an *ec*TFBA. Just set the ```with_thermodynamic_constraints``` and ```with_enzyme_constraints``` to ```True```:

```py
from cobrak.example_models import toy_model
from cobrak.lps import perform_lp_optimization
from cobrak.printing import print_dict, print_optimization_result

# Perform ecTFBA
ectfba_result = perform_lp_optimization(
    cobrak_model=toy_model,
    objective_target="ATP_Consumption",
    objective_sense=+1,
    with_thermodynamic_constraints=True,
    with_enzyme_constraints=True,
)

# Pretty print result as dictionary
print_dict(ectfba_result)

# Pretty print result as tables,
# now with enzyme concentrations *and* metabolite concentrations
print_optimization_result(toy_model, ectfba_result)
```
