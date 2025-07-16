# Evolutionary optimization

??? abstract "Quickstart code"
    ```py
    from cobrak.example_models import toy_model
    from cobrak.evolution import perform_nlp_evolutionary_optimization
    from cobrak.printing import print_dict, print_optimization_result

    # Run the evolutionary optimization, get a sorted dictionary of the form
    #  {
    #   $FOUND_OBJECTIVE_VALUE: [ALL_NLP_RESULTS_WITH_THIS_OBJECTIVE_VALUE],
    #   (...)
    #  }
    # i.e., it contains all found objective values in descending order as keys
    # and all NLP solutions (with all variable values) with this objective
    # value as values.
    complete_result = perform_nlp_evolutionary_optimization(
        cobrak_model=toy_model,
        objective_target="EX_ATP",
        objective_sense=+1,
        variability_dict={}, # No variability dict given -> An ecTFVA is automatically run for us
    )
    print_dict(complete_result)
    best_result = complete_result[list(complete_result.keys())[0]][0] # 0->The first element is the best
    print_optimization_result(toy_model, best_result)

    # %% Postprocesing
    from cobrak.evolution import postprocess

    postprocess_results, best_postprocess_result = postprocess(
        cobrak_model=toy_model,
        opt_dict=best_result,
        objective_target="EX_ATP",
        objective_sense=+1,
        variability_data={}, # No variability dict given -> An ecTFVA is automatically run for us
    )

    print("All postprocessing tries")
    print_dict(postprocess_results)
    print("Best postprocess result (i.e., with best objective value :-)")
    print_optimization_result(toy_model, best_postprocess_result)
    ```

## Introduction

!!! note
    Also check out another toy model evolution example in COBRA-k's repository which you may find
    under ```examples/toymodel/run_toymodel_calculations.py```.

As mentioned in the previous chapter, both the global MINLP and the local NLP have major disadvantages:

* The MINLP quickly becomes way too slow with increasing model size
* The NLP is quick, but only works for one selected set of thermodynamically feasible reactions (i.e. reactions with $f_i>0$)

To overcome these disadvantages, COBRA-k provides an efficient and often quite quick evolutionary algorithm [[Wikipedia]](https://en.wikipedia.org/wiki/Evolutionary_algorithm). This evolutionary algorithm allows one to quite quickly optmize a COBRA-k Model with a full integration of non-linear thermodynamic and kinetic constraints.

!!! warning "Global solutions"
    The evolutionary algorithm cannot *guarantee* that the found solution is global. Hence, it is advisable to let it run until no improvement can be found over a set amount of rounds (see below). Nevertheless, any solution found by the evolutionary algorithm is *valid* and a solution to the given constraints, just maybe not one that is globally optimal.

    For a higher confidence in a found result, it is advisable to run multiple evolutionary optimizations with the same objective function. If the returned results are the same or at least nearly identical, the confidence in the result may rise.

# The evolutionary algorithm

### Input

As *obligatory* input, COBRA-k's evolutionary algorithm takes:

- A *COBRA-k Model*, including all its kinetic and thermodynamic parameters
- An *objective target*, which can be - as for linear optmizations (see LP chapter) - either a single model variable ID (```str```) or a linear term of multiple model variables (```dict[str, float]```). As with all other optmizations, the objective target can be any kind of model variable, including fluxes, concentrations, κ and γ values etc.
- An *objective sense*, i.e. either a minimization (-1) or maximization (+1) of a value.
- Furthermore, COBRA-k's evolutionary algorithm provides many further *optional* inputs. These include e.g. whether or not κ and/or γ constraints shall be regarded and solver settings. For more, see this documentation's API reference for the submodule ```evolution```.

COBRA-k's evolutionary algorithm itself is composed of the following main steps:

1. Preparatory Variability Analysis - to find variable bounds
2. Initialization - to find some feasbile solutions as starting points
3. Evolutionary Run - to find a (hopefully) optimal solution

The next sub-paragraphs explain these steps in more detail.

### Preparatory Variability Analysis

Again, as we use the NLP, we have to run the preparatory variability analysis as explained in the previous chapter. In COBRA-k, you do not have to call this variability analysis separately in your code as such a variablity analysis is automatically calculated if you do not provide already existing variability data.

!!! note
    As for the NLP (see NLP chapter), this preparatory variability analysis is essential to avoid non-linear solver errors.


### Initialization

Using the mentioned input values, the following routines are run for our initializing sampling:

1. We have a given (e.g. randomly selected) set of enforced deactivated reactions (i.e. $v_i=0$).
2. A first ecTFBA (see MILP chapter) is run with the COBRA-k Model and the given objective function. From this first ecTFBA, we store the objective value.
3. Using the objective value from the previous step, we now maximize the number of thermodynamically active reactions, i.e. $max \sum{z_i}$ with an ecTFBA. We store the resulting flux pattern (i.e. all resulting active reactions).
4. Using the flux pattern from the previous step, we create a reduced COBRA-k model where no inactive reactions (i.e. $v_i=0$ in the flux pattern) occur.
5. Using this reduced COBRA-k model, we now run our fast NLP with the general input objective function. We return the objective value.
6. Exception step: If any of the previous optmizations results in infeasibility (i.e. no solution can be found with the given set of deavtivated reactions), we return an appropriate signal to let us know that.

If enough feasible solutions are found in a given maximal number of initializing sampling rounds (all which can be set in COBRA-k's evolution method), we now have an initializing set of feasible flux distributions and associated objective values.


### Evolutionary run

Now that we have variable bounds (from the preparatory variability analysis) and some feasible starting points (from the initialization), we can call the actual evolutionary run. The evolutionary run consists of an outer and an inner optimization, as well as a small preparatory calculations:

#### Identification of stoichiometric couples

First, we identify all stoichiometric couples in our model's stoichiometric matrix. Such couples are reactions which *always* have to be active or inactive together. E.g. if we have a linear pathway such as ->A->B->C->, then all reactions have to be active together in order to achieve a steady-state.

#### Outer optimization

The *outer* optimization is affecting a binary (or $\[0,1\]$) vector which we'll call $\mathbf{β}$. This vector is as long as the number of the previously identified stoichiometric couples. If an element $i$ of $\mathbf{β}$ is 0 (or close to 0), then the flux of all reactions in the respective $i$th stoichiometric couple is set to 0. Otherwise, the reactions are allowed to run.

$\mathbf{β}$ is optmized through a genetic algorithm [[Wikipedia]](https://en.wikipedia.org/wiki/Genetic_algorithm). This algorithm works in *rounds* and has a given *population (or candidate) size*:

1. In each round, as many inner optimizations (see next subparagraph) as there are population members are run. Thus, we actually have as many $\mathbf{β}$ as there are population members.
2. Then, the fitnesses (here, the final NLP optimization value from the inner optimizations) are collected for each population member.
3. Using these fitness values, the $\mathbf{β}$ of the population members are mutated according to the genetic algorithm, with the intention and hope that better $\mathbf{β}$ are found. Then, we start again with step 1 unless the maximal number of rounds or the maximal number of rounds without an increase in the best optimal value is reached.

#### Inner optimization

The *inner* optimization consists of two ecTFBAs (see MILP chapter) and a subsequent NLP (see NLP chapter):

1. The first ecTFBA maximizies our objective.
2. The second ecTFBA maximizied the number of thermodynamically feasible active reactions under our previously calculated objective value. This allows to find solutions which need more than a minimal or near-minimal set of active reactions.
3. The NLP is then run on all active reactions of the last ecTFBA and returns a fitness value :-)

If any of these three steps shows that there is an infeasiblity, an infeasibility signal is sent to the outer optimization.

## Running an evolutionary optimization in COBRA-k

Now let's run an evolutionary optimization with our toy model, which is much simpler to do than to describe the many steps of the algorithm :-) COBRA-k has all its evolutionary algorithm functionality in the ``evolution```subpackage.

```py
from cobrak.example_models import toy_model
from cobrak.evolution import perform_nlp_evolutionary_optimization
from cobrak.printing import print_dict, print_optimization_result

# Run the evolutionary optimization, get a sorted dictionary of the form
#  {
#   $FOUND_OBJECTIVE_VALUE: [ALL_NLP_RESULTS_WITH_THIS_OBJECTIVE_VALUE],
#   (...)
#  }
# i.e., it contains all found objective values in descending order as keys
# and all NLP solutions (with all variable values) with this objective
# value as values.
complete_result = perform_nlp_evolutionary_optimization(
    cobrak_model=toy_model,
    objective_target="EX_ATP",
    objective_sense=+1,
    variability_dict={}, # No variability dict given -> An ecTFVA is automatically run for us
)
print_dict(complete_result)
best_result = complete_result[list(complete_result.keys())[0]][0] # 0->The first element is the best
print_optimization_result(toy_model, best_result)
```

COBRA-k's ```perform_nlp_evolutionary_optimization```function has many more possible parameters (see also the following info box), which you can look up in this documentation's API overview for the ```evolution```module.

!!! info "Working with larger models"
    The evolutionary algorithm can also work efficiently with much larger models. E.g. in COBRA-k's initial publication [[Paper]](), the mid-scale model iCH360 [[Paper]]() was succesfully used. iCH360 consists of hundreds of reactions and metabolites.

    If you have trouble finding a good objective value from the evolutionary algorithm for your model, consider trying out different population sizes (through the ```perform_nlp_evolutionary_optimization``` parameter ```pop_size```, which defaults to your computer's number of CPU cores) and the ```evolution_num_gens```total evolutionary rounds parameter. Usually, the more population members and the more rounds, the better the optimization should become, at the expense of more needed computational time.

## Postprocessing

Oftentimes, the genetic algorithm may fail to identify better sets of active reactions that only *slightly* differ from the best found solution. To mitigate this problem, COBRA-k also provides a postprocessing routine, which simply looks for single (and a low number of) reaction inactivations and activations, trying to identify better solutions. Here's how to use its functionality in COBRA-k:

```py
# ...continuing from the previous code block
from cobrak.evolution import postprocess

postprocess_results, best_postprocess_result = postprocess(
    cobrak_model=toy_model,
    opt_dict=best_result,
    objective_target="EX_ATP",
    objective_sense=+1,
    variability_data={}, # No variability dict given -> An ecTFVA is automatically run for us
)

print("All postprocessing tries")
print_dict(postprocess_results)
print("Best postprocess result (i.e., with best objective value :-)")
print_optimization_result(toy_model, best_postprocess_result)
```

!!! note
    In our toy model, the postprocessing does not lead to a better solution as the genetic algorithm
    already found the optimal solution. However, in larger models, it is strongly advised to double-check
    your genetic algorithm results through running the postprocessing.
