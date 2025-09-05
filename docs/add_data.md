# Automatically add kinetic and thermodynamic data

If you want to perform COBRA-k's thermodynamic and/or enzyme kinetic ("thermokinetic") analyses, you need appropriate data. When loading an existing model (see previous chapter) or create a new Model instance from scratch (see second to last chapter), such data is often missing or at least not directly included in the model instance. To collect such data, COBRA-k provides two major ways, a fully automated one and one where you can add data manually thanks to its ```model_instantiation``` submodule:

## Automatic way

The automatic way of adding thermokinetic data to a newly loaded SBML file uses the function ```get_cobrak_model_with_kinetic_data_from_sbml_model_alone```. Regarding the SBML model, it needs the following identifier usage for reections and metabolites:

* for the reactions: An EC number annotation using the annotation key "ec-code"
* for the metabolites: BiGG IDs as identifiers

Furthermore, some extra non-optional settings have to be provided, whereby the latino-greek organism name is very important.

In addition, the following databases have to be downloaded manually beforehand:

* The BRENDA .json.tar.gz from <https://www.brenda-enzymes.org/download.php>
* The BiGG metabolites txt from <http://bigg.ucsd.edu/data_access>
* taxdmp.zip from <https://ftp.ncbi.nih.gov/pub/taxonomy/>

SABIO-RK and UniProt data are downloaded automatically into a file in the given folder of the function (see below). Keep in mind that this download may take several dozens of minutes! Once the database is downloaded, it is cached, and no new download is triggered.

Using all this information the automatic procedure collects the following information and adds it to the Model:

* $k_{cat}$, $K_M$ and $K_I$ data: From SABIO-RK and BRENDA
* Hill coefficients: From SABIO-RK
* Molecular enzyme weights: From UniProt
* Taxonomic distances (used to collect taxonomically nearer enzyme kinetic data): From NCBI TAXONOMY
* $Δ_r G^{'°}$: Using the eQuilibrator API

Here's a usage example:

```py
from cobrak.model_instantiation import get_cobrak_model_with_kinetic_data_from_sbml_model_alone

cobrak_model = get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
    sbml_model_path="/path/to/sbml.xml",
    path_to_external_resources="/path/where/the/manually/downloaded/datafiles/are",
    folder_of_sabio_database: "/path/where/the/sabiork/database/shall/be/downloaded",
    brenda_version="$CURRENT_BRENDA_VERSION", # E.g. 2023_1
    prefer_brenda=True, # Whether or not BRENDA k_cat or k_M, ... values shall be used if SABIO-RK data is available
    base_species="$MODEL_SPECIES", # E.g. Escherichia coli
    max_prot_pool=0.5,
    conc_ranges={
        # E.g., for all metabolites without a given identifier,
        # we can use the key "DEFAULT":
        "DEFAULT": (1e-6, 0.02),
        # (...)
    },
    inner_to_outer_compartments=["INNERMOST_COMPARTMENT", "NEXT_TO_INNERMOST_", ], # E.g., ["c", "p", "e"], used for dG0 calculation
    phs={"c": 7.0, } # dict[str, float], shows ph of each compartment in the model, used for dG0 calculation
    pmgs={"c": 2.5, } # dict[str, float], shows pMg of each compartment in the model, used for dG0 calculation
    ionic_strenghts={"c": 250, }, # dict[str, float], shows ionic strength in mM of each compartment in the model, used for dG0 calculation
    potential_differences={("c", "p"): 0.15, }, # dict[tuple[str, str], float], shows potential difference from first to second given compartment in mV, used for dG0 calculation
    kinetic_ignored_enzymes=["IDS", "OF", "IGNORED", "ENZYMES", ], # Enzymes for which no kinetic shall be found
    custom_kms_and_kcats={}, # Can be dict[str, EnzymeReactionData | None] if you want to overwrite some kms or kcats
    kinetic_ignored_metabolites=["IDS", "OF", "IGNORED", "METABOLITES",], # IDs of metabolites for which no enzyme kinetic value (e.g., K_M) shall be found
    do_model_fullsplit = True, # Explained below
    do_delete_enzymatically_suboptimal_reactions = True, # Explained below
    ignore_dG0_uncertainty=True, # Whether or not eQuilibrator-calculated dG0 uncertainties shall be simply set to 0
    enzyme_conc_ranges={}, # Is dict[str, tuple[float, float] | None]
    dG0_exclusion_prefixes=[], # Prefixes (first parts of IDs) for which no dG0 shall be set, a common one would be "EX_"; is list[str]
    dG0_exclusion_inner_parts=[], # Infixes (inner parts of IDs) for which no dG0 shall be set, is list[str]
    extra_flux_constraints=[], # list[ExtraFluxConstraint]
    extra_conc_ratios=[], # Is list[ExtraConcRatios]
    data_cache_folder="/path/to/folder/for/uniprot/cache",
    R=$GAS_CONSTANT, # Default is STANDARD_R
    T=$TEMPERATURE, # Default is STANDARD_T
    add_hill_coefficients=True,  # Default is True, if False, no Hill coefficeints are loaded
)
```

Two of the arguments have the following non-obvious meanings:

* ```do_model_fullsplit```: This means that each reaction is going to be split i) for forward & reverse directions and ii) for each enzyme (complex) catalyzing it. E.g., a reversible reaction
```R1: A → B``` catalyzed by the enzyme $E_1$ and the enzyme complex $E_{2,sub1} \space and \space E_{2,sub2}$ is going to be split into the four reactions ```R1_ENZ_E1_FWD: A → B$,  $R1_ENZ_E2SUB1_AND_E2SUB2_FWD: A → B``` and ```R1_ENZ_E1_REV: B → A```,  ```R1_ENZ_E2SUB1_AND_E2SUB2_REV: B → A```. This fullsplit is neccessary in order to perform thermodynamic and enzymatic calculations later on.
* ```do_delete_enzymatically_suboptimal_reactions```: Akin to the enzyme constraint method sMOMENT [!], all (fullsplit) variants of a reaction which do not have the lowest $k_{cat}/MW$ ratio (i.e., which have higher enzyme costs à flux) are *deleted*. Keep in mind that, while this can drastically reduce a model's size, this also means that any $K_M$, $K_I$ etc. variants of reactions are not considered.

## Manual and semi-automatic way

If you want to automatically create only a select amount of data, look up COBRA-k's submodules
```equilibrator_functionality``` (for $Δ_r G^{'°}$), ```uniprot_functionality``` (for molecular enzyme weights),
```sabio_rk_functionality``` (for enzyme kinetic data from SABIO-RK), ```brenda_functionality``` (for enzyme kinetic
data from BRENDA) and ```ncbi_taxonomy_functionality``` (for taxonomy distance data). Their functions are also described in this documentation's API reference.

Finally, if you already have some data and want to add it to an SBML file Model generation, you can use ```get_cobrak_model_from_sbml_and_thermokinetic_data``` as follows:

```py
from cobrak.model_instantiation import get_cobrak_model_from_sbml_and_thermokinetic_data
from cobrak.dataclasses import EnzymeReactionData

cobrak_model = get_cobrak_model_from_sbml_and_thermokinetic_data(
    sbml_path="path/to/sbml.xml",
    extra_flux_constraints: list[ExtraFluxConstraint],
    dG0s={ # Is dict[str, float], unit is kJ/mol
        "$REAC_ID": dG0_of_reaction,
        # (...)
    },
    dG0_uncertainties={ # Is dict[str, float], unit is kJ/mol
        "$REAC_ID": dG0_of_reaction,
        # (...)
    },
    conc_ranges={ # Is dict[str, tuple[float, flaot]], these are the concentrations in M
        "$MET_ID": (min_conc_of_met, max_conc_of_met),
        # (...)
    },
    extra_conc_ratios=[], # Is a list[ExtraConcRatios]
    enzyme_molecular_weights={ # Is dict[str, float], MW in g/mmol
        "$ENZYME_ID": molecular_weight,
        # (...)
    },
    enzyme_reaction_data: dict[str, EnzymeReactionData | None],  # Contains k_cats, k_ms, k_is, k_as and Hill coefficients
    max_prot_pool=0.5, # In g/gDW
    kinetic_ignored_metabolites=["h2_c", "h2_p",], # Is list[str]
    enzyme_conc_ranges = { # Is dict[str, tuple[float, float] | None]
        "$ENZYME_ID": (min_enzyme_conc, max_enzyme_conc),
        # (...)
    },
    do_model_fullsplit: bool = False, # Explained below
    do_delete_enzymatically_suboptimal_reactions: bool = True, # Explained below
    R: float = STANDARD_R, # Standard gas constant
    T: float = STANDARD_T, # Standard temperature
)
```

!!! note
    For a real-life example of a semi-automated way of collecting enzymatic data, look into the file ["examples/iCH360/B_prepare_external_data_for_cobrak.py" in COBRA-k's "examples/iCH360" repository subfolder](https://github.com/klamt-lab/COBRA-k/blob/main/examples/iCH360/B_prepare_external_data_for_cobrak.py). It also contains separate calls to the BRENDA and SABIO-RK functionalities.
