# Installation

COBRA-k is a Python package [hosted on PyPI](https://pypi.org/project/cobrak/), compatible with Python version 3.10 or later. If you have already installed an appropriate Python version on your system, you can easily install COBRA-k with pip:

```sh
pip install cobrak
```

??? info "Optional alternative installation for conda or mamba users"
    If you're using conda or mamba and want to have a clean COBRA-k environment, create an empty environment with Python and pip, activate it and install COBRA-k afterwards. E.g., on a plain bash console, the steps are as follows (again, any Python≥3.10 should work):

    ```sh
    (base) conda create --name cobrak python=3.10 pip -c conda-forge
    (base) conda activate cobrak
    (cobrak) pip install cobrak
    ```
    If you use mamba, just switch "conda" to "mamba".


??? info "Alternative installation for COBRA-k developers"
    If you want to directly play with COBRA-k's code base, follow the [developer instructions in COBRA-k's README]()

## Installation of third-party solvers

COBRA-k always comes pre-packaged with the quite capable open source mixed-integer (non-)linear program solver [SCIP](https://scipopt.org/) and the non-linear program solver [IPOPT](https://github.com/coin-or/Ipopt). With these solvers, all optimizations provided by COBRA-k can be run.

!!! note
    IPOPT comes pre-packages with the free linear subsolver MA27. However, for larger models, the much faster linear subsolver MA57 or one of its alternatives are highly recommended. See here on how to obtain them (e.g. a free academic license is provided):

    [https://licences.stfc.ac.uk/product/coin-hsl](https://licences.stfc.ac.uk/product/coin-hsl)

    In the ```cobrak.standard_solvers```subpackage (which is further explained in later chapters), you can find ```IPOPT_MA57```, an example on how to configure a ```Solver```instance for MA57. The standard solver ```IPOPT```uses the slow MA27.

Additionally, all solvers supported by [pyomo](https://github.com/Pyomo/pyomo) can be used. These include the very fast commercial mixed-integer linear solvers CPLEX and Gurobi. However, as both of these solvers are not automatically bundled with COBRA-k, these solvers *and* their Python bindings have to be installed seperately. Thereby, the Python bindings have to be installed in your COBRA-k Python environment.

For CPLEX, you can find installation instructions here:

* [CPLEX itself](https://www.ibm.com/de-de/products/ilog-cplex-optimization-studio)

* [CPLEX Python bindings](https://www.ibm.com/docs/en/icos/22.1.1?topic=cplex-setting-up-python-api)

For Gurobi, you can find installation instructions here:

* [Gurobi itself](https://www.gurobi.com/)

* [Gurobi Python bindings](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python)
