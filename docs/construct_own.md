# Construct Own (MI)LPs or NLPs

??? abstract "Quickstart code"
    ```py
    #%% Construct own MILP
    from cobrak.example_models import toy_model
    # Get the functionality that returns a basic (Mixed-Integer) Linear Program,
    # i.e. with steady-state constraint and variables for all reactions, and
    # if such constraints are chosen, also enzymes and metabolites.
    from cobrak.lps import get_lp_from_cobrak_model
    # Using COBRA-k's pyomo_functionality submodule, get the functions for setting
    # a custom objective and getting a pyomo solver
    from cobrak.pyomo_functionality import get_solver
    # Get a needed constants from COBRA-k
    from cobrak.constants import BIG_M
    # Get the function that converts pyomo solution states into a COBRA-k dictionary
    from cobrak.utilities import get_pyomo_solution_as_dict
    # Also get a pretty-printing function
    from cobrak.printing import print_dict
    # Get the needed classes from pyomo
    from pyomo.environ import Constraint, Binary, Var

    # Get the basic steady-state Linear Program, here with enzyme and thermodynamic
    # constraints.
    # This function has many more options that you can read in the
    # API reference.
    # The returned lp is a *pyomo* ConcreteModel object on which all operations
    # for pyomo models can be performed.
    lp = get_lp_from_cobrak_model(
        cobrak_model=toy_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        with_loop_constraints=False,
    )

    # Get all reaction IDs that start with "EX_"
    ex_reac_ids = [reac_id for reac_id in toy_model.reactions if reac_id.startswith("EX_")]

    # Now, add the EX_ reaction constraints in pyomo-style, i.e. using
    # Python's general getattr and setattr functions.
    for ex_reac_id in ex_reac_ids:
        # Get the *pyomo* LP flux variable
        flux_var = getattr(lp, ex_reac_id)

        # Now let's add a new binary variable for this EX_ reaction
        # just as we would do it in pyomo
        new_binary_var_id = f"new_binary_{ex_reac_id}"
        setattr(lp, new_binary_var_id, Var(within=Binary))

        # Finally, let's add an associated constraints through
        # which the reaction can only run if the new binary variable > 0
        # We formulate it in Big-M style (i.e. if the binary variable
        # = 1, the reaction flux can be <= Big-M, whereby Big-M is just
        # a large value)
        setattr(
            lp,
            f"new_constraint_{ex_reac_id}",
            Constraint(expr=flux_var <= getattr(lp, new_binary_var_id) * BIG_M)
        )

    # Now that we've added the constraints and variables, we can run the minimization :D
    # First, we manually choose a pyomo solver for this (here, the pre-packaged free solver SCIP)
    # (alternatively, we could also construct a Solver dataclass instance first and
    # send its data to get_solver).
    pyomo_lp_solver = get_solver(
        solver_name="scip",
        solver_options={},
        solver_attrs={},
    )

    # With this solver, we can now solve the optimization :D
    # (here, in verbose style)
    pyomo_lp_solver.solve(lp, tee=True)

    # We then retrieve the solution as a dictionary in the form of dict[str, float]
    lp_result_dict = get_pyomo_solution_as_dict(lp)

    # Finally, print the (unspecatular) result
    print_dict(lp_result_dict)


    #%% Construct own NLP
    # (See previous example for more comments)
    from cobrak.example_models import toy_model
    from cobrak.nlps import get_nlp_from_cobrak_model
    from cobrak.pyomo_functionality import get_solver
    from cobrak.utilities import get_pyomo_solution_as_dict
    from cobrak.printing import print_dict
    from cobrak.constants import LNCONC_VAR_PREFIX
    from pyomo.environ import Constraint, Reals, Var

    # Get the basic steady-state Non-Linear Program, here with enzyme and thermodynamic
    # constraints.
    # This function has many more options that you can read in the
    # API reference.
    # Again, the returned nlp is a *pyomo* ConcreteModel object on which all operations
    # for pyomo models can be performed.
    nlp = get_nlp_from_cobrak_model(
        cobrak_model=toy_model,
        ignored_reacs=[],
        with_kappa=True,
        with_gamma=True,
    )

    # Get all logarithmic metabolite IDs
    log_met_conc_ids = [f"{LNCONC_VAR_PREFIX}{met_id}" for met_id in toy_model.metabolites]

    # Now, create a sum that represents the squared sum
    # of logarithmic concentrations
    squared_log_metconc_sum = 0.0
    for log_met_conc_id in log_met_conc_ids:
        # Get the *pyomo* NLP metabolite concentration variable
        # and add its squared sum
        squared_log_metconc_sum += getattr(nlp, log_met_conc_id) ** 2

    # Now let's add a new real variable for the squared sum that
    # is always equal to it (i.e. it represents this sum, making
    # it later possible to minimize it :-)
    # Also, set it to values >=0
    squared_met_logconc_sum_var_id = "squared_met_logconc_sum"
    setattr(nlp, squared_met_logconc_sum_var_id, Var(within=Reals, bounds=(0.0, 1e6)))
    # Now fix this variable to the pyomo sum expression we created
    setattr(
        nlp,
        "squared_met_logconc_sum_constraint",
        Constraint(expr=getattr(nlp, squared_met_logconc_sum_var_id) == squared_log_metconc_sum)
    )

    # Now we get the NLP solver
    pyomo_nlp_solver = get_solver(
        solver_name="ipopt",
        solver_options={},
        solver_attrs={},
    )

    # With this solver, we can now solve the optimization :D
    # (here, in verbose style)
    pyomo_nlp_solver.solve(nlp, tee=True)

    # We then retrieve the solution as a dictionary in the form of dict[str, float]
    nlp_result_dict = get_pyomo_solution_as_dict(nlp)

    # Finally, print the (unspecatular) result
    print_dict(nlp_result_dict)
    ```

## Introduction

Up to now, we looked at the range of predefined (mixed-integer) linear programs (e.g. ecTFVA, bottleneck analyses, ...; see LP and MILP chapters) and non-linear programs (see NLP chapter) provided by COBRA-k, whereby...

- ...the general optimization functions (such as ```perform_lp_optimization```and ```perform_nlp_optimization```) allow one to optimize any objective function in the model
- ...and the special optimziation functions (such as ```perform_lp_thermodynamic_bottleneck_analysis```) provide expanded programs with additional constraints and variables.

But sometimes, for advanced optimizations, you need to add your own extra constraints and/or variables. Luckily, this is possible in COBRA-k thanks to its internal usage of pyomo [Website]()[GitHub](), as explained in the following subchapters :-)

!!! note "Alternative for simple cases: Extra (non-)linear watches and constraints"
    If you just want to restrict a (non-)linear weighted sum of any kind of model variables, you can always use the ```extra_linear_watches```, ```extra_nonlinear_watches```, ```extra_linear_constraints``` and ```extra_nonlinear_constraints``` member variable as explained in the LP and NLP chapters and the API documentation.

## Example 1: Construct own (MI)LPs

Let's say that we want to minimize the number of used exchange reactions (i.e. reactions where a metabolite is consumed or produced "out of nothing"; commonly, their ID starts with "EX_"). As they do not have an associated ΔᵣG'° (see MILP chapter), we do not have binary variables for them which control whether they are active or not. So let's introduce these binary variables and performing the subsequent minimization as follows for our toy model:

```py
# Get our toy model
from cobrak.example_models import toy_model
# Get the functionality that returns a basic (Mixed-Integer) Linear Program,
# i.e. with steady-state constraint and variables for all reactions, and
# if such constraints are chosen, also enzymes and metabolites.
from cobrak.lps import get_lp_from_cobrak_model
# Using COBRA-k's pyomo_functionality submodule, get the functions for setting
# a custom objective and getting a pyomo solver
from cobrak.pyomo_functionality import get_solver
# Get a needed constants from COBRA-k
from cobrak.constants import BIG_M
# Get the function that converts pyomo solution states into a COBRA-k dictionary
from cobrak.utilities import get_pyomo_solution_as_dict
# Also get a pretty-printing function
from cobrak.printing import print_dict
# Get the needed classes from pyomo
from pyomo.environ import Constraint, Binary, Var

# Get the basic steady-state Linear Program, here with enzyme and thermodynamic
# constraints.
# This function has many more options that you can read in the
# API reference.
# The returned lp is a *pyomo* ConcreteModel object on which all operations
# for pyomo models can be performed.
lp = get_lp_from_cobrak_model(
    cobrak_model=toy_model,
    with_enzyme_constraints=True,
    with_thermodynamic_constraints=True,
    with_loop_constraints=False,
)

# Get all reaction IDs that start with "EX_"
ex_reac_ids = [reac_id for reac_id in toy_model.reactions if reac_id.startswith("EX_")]

# Now, add the EX_ reaction constraints in pyomo-style, i.e. using
# Python's general getattr and setattr functions.
for ex_reac_id in ex_reac_ids:
    # Get the *pyomo* LP flux variable
    flux_var = getattr(lp, ex_reac_id)

    # Now let's add a new binary variable for this EX_ reaction
    # just as we would do it in pyomo
    new_binary_var_id = f"new_binary_{ex_reac_id}"
    setattr(lp, new_binary_var_id, Var(within=Binary))

    # Finally, let's add an associated constraints through
    # which the reaction can only run if the new binary variable > 0
    # We formulate it in Big-M style (i.e. if the binary variable
    # = 1, the reaction flux can be <= Big-M, whereby Big-M is just
    # a large value)
    setattr(
        lp,
        f"new_constraint_{ex_reac_id}",
        Constraint(expr=flux_var <= getattr(lp, new_binary_var_id) * BIG_M)
    )

# Now that we've added the constraints and variables, we can run the minimization :D
# First, we manually choose a pyomo solver for this (here, the pre-packaged free solver SCIP)
# (alternatively, we could also construct a Solver dataclass instance first and
# send its data to get_solver).
pyomo_lp_solver = get_solver(
    solver_name="scip",
    solver_options={},
    solver_attrs={},
)

# With this solver, we can now solve the optimization :D
# (here, in verbose style)
pyomo_lp_solver.solve(lp, tee=True)

# We then retrieve the solution as a dictionary in the form of dict[str, float]
lp_result_dict = get_pyomo_solution_as_dict(lp)

# Finally, print the (unspecatular) result
print_dict(lp_result_dict)
```

## Example 2: Construct own NLPs

Of course, we can also construct our own NLPs. Here's an example where we (for whatever reason xD) we try to minimize the
*squared* logarithmic metabolite concentrations in our model:

```py
# (See previous example for more comments)
from cobrak.example_models import toy_model
from cobrak.nlps import get_nlp_from_cobrak_model
from cobrak.pyomo_functionality import get_solver
from cobrak.utilities import get_pyomo_solution_as_dict
from cobrak.printing import print_dict
from cobrak.constants import LNCONC_VAR_PREFIX
from pyomo.environ import Constraint, Reals, Var

# Get the basic steady-state Non-Linear Program, here with enzyme and thermodynamic
# constraints.
# This function has many more options that you can read in the
# API reference.
# Again, the returned nlp is a *pyomo* ConcreteModel object on which all operations
# for pyomo models can be performed.
nlp = get_nlp_from_cobrak_model(
    cobrak_model=toy_model,
    ignored_reacs=[],
    with_kappa=True,
    with_gamma=True,
)

# Get all logarithmic metabolite IDs
log_met_conc_ids = [f"{LNCONC_VAR_PREFIX}{met_id}" for met_id in toy_model.metabolites]

# Now, create a sum that represents the squared sum
# of logarithmic concentrations
squared_log_metconc_sum = 0.0
for log_met_conc_id in log_met_conc_ids:
    # Get the *pyomo* NLP metabolite concentration variable
    # and add its squared sum
    squared_log_metconc_sum += getattr(nlp, log_met_conc_id) ** 2

# Now let's add a new real variable for the squared sum that
# is always equal to it (i.e. it represents this sum, making
# it later possible to minimize it :-)
# Also, set it to values >=0
squared_met_logconc_sum_var_id = "squared_met_logconc_sum"
setattr(nlp, squared_met_logconc_sum_var_id, Var(within=Reals, bounds=(0.0, 1e6)))
# Now fix this variable to the pyomo sum expression we created
setattr(
    nlp,
    "squared_met_logconc_sum_constraint",
    Constraint(expr=getattr(nlp, squared_met_logconc_sum_var_id) == squared_log_metconc_sum)
)

# Now we get the NLP solver
pyomo_nlp_solver = get_solver(
    solver_name="ipopt",
    solver_options={},
    solver_attrs={},
)

# With this solver, we can now solve the optimization :D
# (here, in verbose style)
pyomo_nlp_solver.solve(nlp, tee=True)

# We then retrieve the solution as a dictionary in the form of dict[str, float]
nlp_result_dict = get_pyomo_solution_as_dict(nlp)

# Finally, print the (unspecatular) result
print_dict(nlp_result_dict)
```
