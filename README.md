# COBRA-k Package

<table style="border: none; border-collapse: collapse;">
    <tr>
        <td>
            <div style="background-color: #2E4053; padding: 4px 8px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);">
            <a style="color: #FFFFFF; text-decoration: none; font-weight: bold;" href="https://pypi.org/project/cobrak/">
            pip &#124; version 0.0.1
            </a>
            </div>
        </td>
        <td>
            <div style="background-color: #964B00; padding: 4px 8px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);">
            <a style="color: #FFFFFF; text-decoration: none; font-weight: bold;" href="https://klamt-lab.github.io/COBRA-k/">
            Documentation
            </a>
            </div>
        </td>
        <td>
            <div style="background-color: #3E8E41; padding: 4px 8px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);">
            <a style="color: #FFFFFF; text-decoration: none; font-weight: bold;" href="https://github.com/klamt-lab/COBRA-k/issues">
            GitHub issues
            </a>
            </div>
        </td>
        <td>
            <div style="background-color: #666666; padding: 4px 8px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);">
            <a style="color: #FFFFFF; text-decoration: none; font-weight: bold;" href="#troubleshooting-and-contact">
            Contact
            </a>
            </div>
        </td>
        <td>
            <div style="background-color: #444444; padding: 4px 8px; border-radius: 4px; cursor: pointer; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);">
            <a style="color: #FFFFFF; text-decoration: none; font-weight: bold;" href="#publication">
            Publication
            </a>
            </div>
        </td>
    </tr>
</table>

<details>
 <summary>Index</summary>

 * [Introduction](#introduction)
 * [Installation](#installation)
 * [Documentation](#documentation)
 * [Troubleshooting and contact](#troubleshooting-and-contact)
 * [Contributing as developer](#contributing-as-developer)
 * [License](#license)
</details>

## Graphical Abstract

<img src="docs/img/COBRAk_figure.png" alt="Graphical COBRA-k abstract" style="max-height:450px;">

## Introduction

COBRA-k stands for "Constraint-Based Reconstruction and Analysis (COBRA) **with kinetics**" [[Paper]](#publication), a generalized variant of the metabolic modeling framework known as COBRA. Both frameworks allow one to analyze constraint-based metabolic models with the help of optimization approaches. While COBRA uses exclusively linear constraints to describe steady state flux distributions in a given metabolic network, COBRA-k expands COBRA by allowing also the integration of non-linear kinetic rate laws (herein we use reversible Michaelis-Menten kinetics. In particular, this enables also a consistent integration of (steady state) fluxes, metabolite concentration and enzyme abundances in COBRA models. While classical COBRA model are usually formulated as linear or mixed-integer linear program (LP / MILP), COBRA-k requires solving mixed-integer **non-liner** programs (MINLP).

The COBRA-k package is a general COBRA/COBRA-k suite written as a Python module, while also being the reference implementation of COBRA-k. Some of COBRA-k's major features are:

1. "Classical" COBRA methods such as:
    * Flux Balance Analysis (FBA) [[Review]](https://doi.org/10.1038/nbt.1614)
    * Parsimonious Flux Balance Analysis (pFBA) [[Paper]](https://doi.org/10.1038/msb.2010.47)
    * Flux Variability Analysis (FVA) [[Paper]](https://doi.org/10.1016/j.ymben.2003.09.002) (written, as all other provided variability analyses, in a multi-core-parallelized way)

    These standard methods can also optionally be treated in combination with general enzyme constraints *and/or* thermodynamic constraints, enabling methods such as:

    * Enzyme-constrained FBA based on MOMENT [[Paper]](https://doi.org/10.1371/journal.pcbi.1002575), GECKO [[Paper]](https://doi.org/10.15252/msb.20167411) or sMOMENT [[Paper]](https://doi.org/10.1186/s12859-019-3329-9)
    * Thermodynamic Flux Balance Analysis (TFBA) [[Paper]](https://doi.org/10.1529/biophysj.106.093138)
    * OptMDFpathway [[Paper]](https://doi.org/10.1371/journal.pcbi.1006492)
    * Concentration Variability Analysis (CVA) [[Paper]](https://doi.org/10.1371/journal.pcbi.1009093)
    * ...and more!

You may also freely choose your constraints and objectives, allowing you to construct new COBRA methods.


2. COBRA-k methods operating on extended COBRA models with *reversible Michaelis-Menten kinetics* [[Paper]](https://www.sciencedirect.com/science/article/pii/S0014579313005577):

    * Our developed COBRA-k solver (combining a genetic algorithm with enzyme-constrained FBA + non-linear optimization) to solve also larger COBRA-k problems.
    * A MINLP formulation that can be used in combination with global MINLP solvers (works only with small models).

    Again, you can freely choose whether you want to enable enzymatic *and/or* thermodynamic constraints, even on the reaction level.

3. Automatic thermodynamic and/or enzyme-kinetic data retrieval for your model:

    * $k_{cat}$ (turnover numbers) and $K_M$ (Michaelis constants) using the databases SABIO-RK [[Website]](http://sabio.h-its.org/) and BRENDA [[Website]](https://www.brenda-enzymes.org/) as well as taxonomic information from NCBI TAXONOMY [[Website]](https://www.ncbi.nlm.nih.gov/taxonomy); this retrieval algorithm is similar to the one used by AutoPACMEN [[GitHub]](https://github.com/klamt-lab/autopacmen)
    * $Δ_rG'°$ (standard reaction-based Gibb's free energies) using the fantastic eQuilibrator-API [[GitLab]](https://gitlab.com/equilibrator/equilibrator-api)
    * $W_i$ (molecular weights of enzymes) using Uniprot [[Website]](https://www.uniprot.org/)
    * Of course, you can also manually change or set any parameters at your will

4. Pretty-printing and exporting results:

    * With the help of rich [[GitHub]](https://github.com/Textualize/rich), COBRA-k allows one to display colored tables of any COBRA(k) analysis results, as well as models themselves.
    * As storage options, results can be saved and loaded as JSON text files, and also exported as XLSX spreadsheets.
    * Results can also be directly exported as CNApy [[GitHub]](https://github.com/cnapy-org/CNApy) scenarios, making it possible to visualize them on *interactive* metabolic maps.

5. Interoperable loading, editing and saving of metabolic models:

    * Metabolic models and their components (reactions, metabolites, enzymes and related parameters) can be freely constructed, changed and/or augmented using COBRA-k.
    * Models can be imported and exported in the widely used SBML format [[Site]](https://sbml.org/) with the help of COBRApy [[GitHub]](https://github.com/opencobra/cobrapy). Thereby, COBRA-k can store and reload any additional kinetic and/or thermodynamic information in SBML annotations. The SBML format is also used by other popular constraint-based packages for various programming languages.
    * Additionally, models, model components and results can be quickly saved and loaded in COBRA-k's custom JSON format.
    * Also, models with an associated COBRA-k solution can be directly exported as kinetic models in the Antimony and/or SBML format, making it possible to use them with popular kinetic systems biology packages such as Tellurium [[GitHub]](https://github.com/sys-bio/tellurium) and COPASI [[GitHub]](https://github.com/copasi/COPASI).

Programmatically, COBRA-k relies on the optimization framework pyomo [[GitHub]](https://github.com/Pyomo/pyomo), making it compatible with any solver supported by it (see [here in pyomo's documentation](https://pyomo.readthedocs.io/en/6.8.0/solving_pyomo_models.html#supported-solvers) for more). For the conversion of COBRA-k models into kinetic models, Tellurium [[GitHub]](https://github.com/sys-bio/tellurium) is used.
<br>

## Installation

COBRA-k is provided as a PyPI package (follows) and can be installed through the following command:

```sh
pip install cobrak
```

*Note:* If you encounter any trouble due to a missing SCIP installation (which may happen on some systems), you have to install the free and open-source mixed-integer linear solver [SCIP](https://scipopt.org/) on your system. To do so, follow the SCIP download instructions [here](https://scipopt.org/index.php#download) or or, if no matching download is provided for your system, the compilation instructions [here](https://scipopt.org/doc/html/INSTALL.php).

After the installation, which should just take a few minutes, the package can be imported as ```cobrak```.

Also note that only the free and open-source solver [IPOPT](https://github.com/coin-or/Ipopt) (and on some systems, also automatically the solver SCIP) come pre-installed. For more on how to install other solvers, [visit the respective documentation page here](https://klamt-lab.github.io/COBRA-k/installation.html#installation-of-third-party-solvers).

COBRA-k requires Python≥3.10, making it work on Python-3.10-compatible OS such as Windows, MacOS and Linux. Thereby, COBRA-k is mainly tested on Ubuntu on an x86-64 computer and on a MacBook with an ARM processor.

<br>

## Documentation

An introduction to using COBRA-k, including a tutorial with integrated quickstart sections, is provided through its documentation:

https://www.klamt-lab.github.io/COBRA-k/

Also, you may check out COBRA-k's standard examples, provided in the "examples" subfolder of this repository. They also contain all the files to reproduce the results of COBRA-k's publication, as also further explained in the respective chapter in COBRA-k's documentation.

<br>

## Troubleshooting and Contact

If you encounter any issues with COBRA-k, don't hesitate to open a GitHub issue here:

[https://github.com/klamt-lab/COBRA-k/issues](https://github.com/klamt-lab/COBRA-k/issues)

Or write an e-mail, which is provided on the COBRA-k package's main maintainer [employer page profile](https://www.mpi-magdeburg.mpg.de/person/98416).

<br>

## Contributing as developer

You're welcome to contribute to COBRA-k :-), feel free to open pull requests. Keep in mind that you have to be ok with COBRA-k's open-source license (see next paragraph).

### Code structure

* Main folder: The PyPI descriptions and environment settings can be found in the main folder's "pyproject.toml" and "setup.py". The "mkdocs.yml" contains COBRA-k's documentation index and module settings.
* Subfolder "cobrak": Here, you can find all COBRA-k package Python functionalities
* Subfolder "docs": Contains the mkdocs Markdown text for COBRA-k's documentation
* Subfolder "examples": Contains all main COBRA-k examples, as mentioned in COBRA-k's publication; "examples/common_needed_external_resources" contains files that are used by multiple examples, and are created by the .py files in "examples" itself. For a short demo that usually takes only up to a minute, run and look at "examples/toymodel/run_toymodel_calculations.py".
* Subfolder "tests": Contains pytest tests for many COBRA-k functions; also imports the toymodel generation and analyses as test in "tests/test_toymodel.py"


### Python code environment

In any code contribution, Python's mypy typing capabilities should be used [(see this mypy cheat sheet)](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html). Also, ruff [[GitHub]](https://github.com/astral-sh/ruff) format is used for code formatting, together with the ruff linter import sorter. As testing framework for testable functions, pytest is used.

To set up a development environment, you may e.g. use uv [[GitHub]](https://github.com/astral-sh/uv) as follows:

1. Make sure you've installed uv and ruff:

```
pip install uv
uv pip install ruff
```

2. Checkout the latest COBRA-k development version using git

```
git clone https://github.com/klamt-lan/COBRA-k.git
```

3. (optional) Run the toy model example

```
cd COBRA-k
uv run examples/toymodel/run_toymodel_calculations.py
```

uv will automatically install a correct Python version ≥ 3.10 and COBRA-k dependencies.

4. After code changes: Apply ruff

To maintain a consistent style, please apply ruff on you code before you create a pull request as follows (you may run these commands in COBRA-k's main folder):

```
ruff check --select I --fix
ruff format
ruff check
# All three commands can also be found in this reposiroty's format_code.sh
```

The first command sorts imports, the second one fixes code formatting. The third command tells you remaining issues that you have to change manually. If you want to deviate from ruff's rules, you can e.g. use a ```# noqa: $RULE_NAME``` command at the end of each affected line. This will make ruff ignore the found issue.

5. After code changes: Test with pytest

If applicable, you are invited to write new or changed tests in COBRA-k's code. Tests are located in the "tests" subfolder and use pytest [[GitHub]](https://github.com/pytest-dev/pytest) which you may install as follows:

```
uv pip install pytest
```

You can then run the tests in the "test" subfolder as follows:

```
pytest tests
```

Please make sure that your code changes do not unintendently break any tests.

### Documentation development

To work on COBRA-k's documentation, you need mkdocs [[GitHub]](https://www.mkdocs.org/) and some dependencies. To install them, you may run the following commands using uv:

```
uv pip install mkdocs
uv pip install mkdocs-material
uv pip install pymdown-extensions
uv pip install mkdocstrings[python]
```

To now self-run the documentation on a local server, you can use mkdocs itself with the following command in COBRA-k's main folder:

```
mkdocs serve
```

<br>

## License

COBRA-k is free and open-source and licensed under the Apache License, Version 2.0

Note that COBRA-k's documentation integrated the open-source JavaScript packages KaTeX [[Link to its open-source license]](docs/javascript/katex/LICENSE.txt), iframe-worker [[Link to its open-source license]](docs/javascript/.cache/assets/external/unpkg.com/iframe-worker/LICENSE.txt) and mermaid [[Link to its open-source license]](docs/javascript/.cache/assets/external/unpkg.com/mermaid@11/LICENSE.txt).

Some of COBRA-k's code is also adapted and modified from AutoPACMEN [[GitHub]](https://github.com/klamt-lab/autopacmen) and TCOSA [[GitHub]](https://github.com/klamt-lab/TCOSA). For both repositories, their main author is the same as the one of COBRA-k and both packages use the same license type as COBRA-k.

<br>

## Publication

If you use COBRA-k in your academic work, please cite its publication:

* Bekiaris PS & Klamt S (2025). COBRA-k: a powerful framework bridging constraint-based and kinetic metabolic modeling. *Submitted*.
