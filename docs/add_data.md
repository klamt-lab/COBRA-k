# Automatically add kinetic and thermodynamic data

If you want to perform COBRA-k's thermodynamic and/or enzyme kinetic (in short, "thermokinetic") analyses, you need appropriate data, namely:

* If you want to use enzyme constraints (see LP chapter): $k_{cat}$ values for reactions in h⁻¹ and molecular weights ($MW$) for enzymes in kDa
* If you want to use thermodynamic constraints (see MILP & NLP chapters): ΔG'° values in kJ⋅mol⁻¹
* If you want to use kinetic saturation term κ constraints (see NLP chapter): $k_M$ values in M for metabolites of reactions
* If you want to use advanced kinetics: $k_I$ for inhibition, $k_A$ for activation and Hill coefficients for binding effects

!!! note
    It's not necessary to collect all of this data for COBRA-k. E.g., if you want to calculate with thermodynamic constraints only, you only need ΔG'° values. On the other hand, it also doesn't hurt to collect more data than necessary for your model as long as you don't use it in your calculations. As explained in the previous chapters, the COBRA-k package allows you to control which constraints are active in a flexible way.

If you already have such data, follow the subchapter ["manually adding data"](#manually-adding-data) at the bottom of this page. If you only have some of this data and want to use COBRA-k's automatic collection functions for the rest, read the [second-to-last chapter](#automatically-collecting-some-of-the-data). But if you don't have any such data, COBRA-k provides automatic data retrieval functions that are explained in the [second-to-next subchapter](#full-automatic-data-collection).

Before you can automatically retrieve some or all of the thermokinetic data for your model, make sure that your model must be correctly split and have the right usage of identifiers and annotations, as explained in the following subchapters :-)

## Model requirements for any COBRA-k model

As always when using the COBRA-k package (see e.g. the chapter about Linear Programming), models must be "fullsplit" which you can automatically do for SBML models using `load_annotated_sbml_model_as_cobrak_model` in `cobrak.io`.

"Fullsplit" means that each original reaction is split i) for forward & reverse directions and ii) for each enzyme (complex) catalyzing it. E.g., a reversible reaction `R1: A → B` catalyzed by the enzyme $E_1$ and the enzyme complex $E_{2,sub1} \space and \space E_{2,sub2}$ is going to be split into the four reactions `R1_ENZ_E1_FWD: A → B,  R1_ENZ_E2SUB1_AND_E2SUB2_FWD: A → B` and `R1_ENZ_E1_REV: B → A`,  `R1_ENZ_E2SUB1_AND_E2SUB2_REV: B → A`. This fullsplit is necessary in order to perform thermodynamic and enzymatic calculations later on.

## Additional model requirements for automatic data collection

!!! note
    This model preparation is not necessary if you just want to add *already existing* thermokinetic data manually, as explained in the ["manually adding data"](#manually-adding-data) subchapter.

Before we explain the naming requirements, we'll illustrate these requirements with a very small toy model that we'll call ```data_toy_model```. This small model represents the forward glucose-6-phosphate isomerase (BiGG ID PGI_fw) with just one substrate (cytosolic glucose-6-phosphate; BiGG ID g6p_c) and one product (cytosolic fructose-6-phosphate; BiGG ID f6p_c):

```py
from cobrak.dataclasses import Model, Reaction, EnzymeReactionData, Metabolite, Enzyme

# Let's define the COBRA-k Model instance (for more about the dataclass Model
# and the other dataclasses, read chapter "Create model from scratch")
data_toy_model = Model(
    metabolites={
        # Most parameter searches need the metabolites to use BiGG
        # database IDs, which are already used in many existing
        # metabolite networks. If not, you can look them up (and
        # download the database itself) from https://bigg.ucsd.edu/
        "g6p_c": Metabolite(
            # ...additionally, for the search of ΔG'° values, you
            # have to use eQuilibrator-API-compatible
            # annotations, i.e. metabolite IDs from a multitude of databases
            # Let's define three examples (for a full list of supported
            # identifiers, check out the USED_IDENTIFIERS_FOR_EQUILIBRATOR
            # constant in cobrak.constants; note: Sometimes, INCHI strings
            # and keys cannot be read out correctly)
            annotation={
                "bigg.metabolite": "g6p",  # From the BiGG database (https://bigg.ucsd.edu/)
                "kegg.compound": "C00092",  # From the KEGG database (https://www.genome.jp/kegg/)
                "metanetx.chemical": "MNXM160",  # From MetaNetX (https://www.metanetx.org/)
                # ...again, this extra annotation is fully optional as long
                # as the eQuilibrator-API can read out your metabolite ID
            },
        ),
        # Also for f6p, we have to define an eQuilibrator-API-compatible annotation explicitly
        "f6p_c": Metabolite(
            annotation={
                "bigg.metabolite": "f6p",  # From the BiGG database (https://bigg.ucsd.edu/)
            }
        ),
    },
    reactions={
        # While this reaction uses a BiGG ID, the reaction ID can actually
        # be of any format (unlike many other IDs, as explained here)
        "PGI_fw": Reaction(
            stoichiometries={
                # Note again that the metabolites use BiGG IDs
                "g6p_c": -1.0,
                "f6p_c": 1.0,
            },
            annotation={
                # For k_cat, k_M and other kinetic reaction parameters,
                # it is *necessary* to give the reaction a valid Enzyme Commission
                # (EC) number through such an 'ec-code' annotation, which is already
                # included in many published metabolic models.
                # If you do not know the EC code of your reaction, you can try to
                # look it up through databases such as, amongst others,
                # EXPASY ENZYME (https://enzyme.expasy.org/) or also BiGG (https://bigg.ucsd.edu/)
                "ec-code": "5.3.1.9",
            },
            # As always with COBRA-k, reactions have to be *irreversible*,
            # reversible reactions have to be split up beforehand,
            # which you can e.g. automatically do for SBML models with the COBRA-k function
            # `load_annotated_sbml_model_as_cobrak_model` in `cobrak.io`
            # while keeping the `do_model_fullsplit` argument at `True`.
            min_flux=0.0,
            max_flux=1000.0,
            # And, again, as always with COBRA-k, not only reversible reactions have to be split
            # into separate ones, but also reactions that are catalyzed by multiple enzymes (isozymes).
            # In this case, the reaction is split into as many variants as there are enzymes that catalyze it.
            # Again, you can do this automatically for SBML models with the COBRA-k function
            # `load_annotated_sbml_model_as_cobrak_model` in `cobrak.io`
            # while keeping the `do_model_fullsplit` argument at `True`.
            enzyme_reaction_data=EnzymeReactionData(
                identifiers=["b4025"],
            ),
        ),
        # To make this model work, we'll also add pseudo-reactions that deliver
        # g6p_c and take up f6p_c into the environment. As these reactions are not
        # mass-balanced (mass of substrates ≠ mass of products), thermokinetic parameters
        # and constraints do not make any sense here, so that we can omit any annotations.
        "EX_g6p_c": Reaction(
            stoichiometries={"g6p_c": +1.0}
        ),  # Produce substrate glucose-6-phosphate
        "EX_f6p_c": Reaction(
            stoichiometries={"f6p_c": -1.0}
        ),  # Take up product fructose-6-phosphate
    },
    enzymes={
        # For the automated collection of molecular weights, enzymes need an ID and/or
        # name that can be found in Uniprot (https://www.uniprot.org/). Make sure
        # that the ID and/or name is not ambiguous in your modeled organism (the
        # automatic routine only searches for enzymes of the modeled organism), so that
        # the right enzyme can be chosen.
        "b4025": Enzyme(
            name="pgi",
        ),
    },
)

# ...we can also find this model as cobrak.example_models.data_toy_model
```

!!! note "Requirements for SBML models"
    While the annotation and ID requirements were explained for a COBRA-k Model instance, the same requirements exist for SBML models. I.e., metabolites should use BiGG IDs, reactions should have an `ec-code` annotation and so on...

To summarize the key points for automatic data collection (in addition to the "fullsplit" requirement explained above):

* **Additionally for $k_{cat}$ and $k_M$ values**: Metabolite IDs must be BiGG IDs, enzyme IDs and/or names must be readable in Uniprot,
and reactions need an ```ec-code``` EC number annotation
* **Additionally for ΔG'° values**: Metabolite annotations must contain any eQuilibrator-API compatible values and annotation keys


## Full automatic data collection

!!! warning
    Before you can effectively use these functions, make sure that you performed the necessary manual downloads (if you want to collect $k_{cat}$ or $k_M$ values) *and* that your model uses identifiers/annotations as explained in the [previous subchapter](#additional-model-requirements-for-automatic-data-collection).

If you do not have any thermokinetic data for your model, COBRA-k provides functions for adding it automatically. Thereby, COBRA-k uses the following databases:

* $k_{cat}$ and $k_M$: From [SABIO-RK](https://www.sabio.h-its.org/) (using its web API) and [BRENDA](https://www.brenda-enzymes.org/) (using its downloadable JSON file)
* Molecular enzyme weights: From [UniProt](https://www.uniprot.org/) (using its web API)
* Taxonomic distances (used to collect the taxonomically nearest enzyme kinetic data): From [NCBI TAXONOMY](https://www.ncbi.nlm.nih.gov/taxonomy)
* $Δ_r G^{'°}$: Using the [eQuilibrator-API](https://equilibrator.readthedocs.io/en/latest/) (which downloads a huge database on first usage; this only happens once!)
* In addition, for translating BiGG IDs: The [BiGG database](https://bigg.ucsd.edu/)
* Optionally, for finding EC number transfers (e.g. if an EC number became obsolete and was superseded by a new one): [EXPASY ENZYME](https://enzyme.expasy.org/)

Now, to make the automatic data collection work, you have to manually download some of the databases:
In addition, also for legal reasons, the following databases have to be downloaded manually beforehand *and put into a single folder*:

* The BRENDA .json.tar.gz from <https://www.brenda-enzymes.org/download.php>
* The BiGG metabolites txt from <http://bigg.ucsd.edu/data_access>
* The NCBI TAXONOMY data from taxdmp.zip from <https://ftp.ncbi.nih.gov/pub/taxonomy/>
* Optionally (but recommended), the EC number transfer enzyme.rdf from <https://ftp.expasy.org/databases/enzyme/>

SABIO-RK and UniProt data are downloaded automatically into a file in the given folder of the function (see below). Keep in mind that the SABIO-RK download may take several dozens of minutes! Once the database is downloaded, it is cached, and no new download is triggered.

!!! info
    If you encounter any server problems with SABIO-RK, you can also download an existing (but older) COBRA-k-compatible SABIO-RK cache from here: https://github.com/klamt-lab/COBRA-k/blob/main/examples/common_needed_external_resources/sabio_single_tsvs.zip. Again, put this file into the same folder as the other mentioned database files.
    
    Using this cache, you omit the long SABIO-RK download at the price of missing newer data that was added after late 2024.

Using all this information the automatic procedure collects the following information and adds it to the Model:

Phew 😸! Now that you've come so far, adding data automatically is easy :D Depending on whether you want to i) create a COBRA-k Model directly from an existing SBML, or ii) add the data to an existing COBRA-k Model, use one of the following methods:

### i) Create model with full data directly from an existing SBML

For SBML-based full data usage, we can use COBRA-k's `model_instantiation` (i.e. functions that create Model instances from other types) submodule as follows:

!!! note
    You can safely ignore any `'' is not a valid SBML 'SId'. No objective coefficients in model. Unclear what should be optimized`
    messages :-)

```py
# First, we'll simulate creating an SBML model out of our data toy model
# in a temporary directory (just ignore this part if you already have an SBML,
# then just set sbml_path to where your SBML is located)
from cobrak.example_models import data_toy_model
from cobrak.io import save_cobrak_model_as_annotated_sbml_model
from tempfile import TemporaryDirectory

temp_directory = TemporaryDirectory()
sbml_path = temp_directory.name + "temp.xml"
save_cobrak_model_as_annotated_sbml_model(
    data_toy_model,
    sbml_path,
)

# Now, we can create a COBRA-k Model instance with full data out of this SBML :-)
# *Note*: We assume that you put the downloaded needed files (see above) into a subfolder
# called 'database_data'
from cobrak.model_instantiation import (
    get_cobrak_model_with_kinetic_data_from_sbml_model_alone,
)

cobrak_model_with_full_data = get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
    sbml_path=sbml_path,
    database_data_folder="./database_data",
    brenda_version="2025_1",  # or 2024_1 if you use the older BRENDA database version
    base_species="Escherichia coli",
)
```

!!! warning "Cache files"
    You may notice that in the `database_data_folder`, several new files starting with "_cache" were generated. These contain the found database data. But if you change your model or run this routine with a different model, you have to delete these cache files first! They exist so that, if you do not change your model, the expensive database searches do not have to be performed again.

!!! note "Important setting"
    `get_cobrak_model_with_kinetic_data_from_sbml_model_alone` has (amongst other) the important argument: `do_delete_enzymatically_suboptimal_reactions`: Akin to the enzyme constraint method sMOMENT [!], all (fullsplit) variants of a reaction which do not have the lowest $k_{cat}/MW$ ratio (i.e., which have higher enzyme costs à flux) are *deleted*. Keep in mind that, while this can drastically reduce a model's size, this also means that any $K_M$, $K_I$ etc. variants of reactions are not considered.

Now, you have a COBRA-k model with k_cats, k_Ms and ΔG'°. Let's have a look at its data, which reveal its $MW$, ΔG'°, $k_{cat}$ and $k_{M}$ values:

```py
# ...following the last code block
from cobrak.printing import print_model

print_model(cobrak_model=cobrak_model_with_full_data)
```

!!! info "cobra_global_settings" reactions
    You may notice an empty cobra_global_settings reaction. This is an artifact from the SBML ex- and import, but you can safely ignore it.

You can also e.g. save the newly generated COBRA-k model  using

```py
# ...following the second-to-last code block
from cobrak.io import json_write

json_write("cobrak_model.json", cobrak_model_with_full_data)
```
and/or perform any other calculations with it :-) You can also save it again as an (this time with COBRA-k-compatible
annotations for the newly generated thermokinetic data) SBML model (for more, see chapter "Load/Save Models").

### ii) Add full data to an existing COBRA-k Model

If you already have a COBRA-k Model object, you don't need the SBML-based functionality. Instead, you can use COBRA-k's ```thermokinetic_data_retrieval``` submodule. For our toy model, it looks like this:

```py
from cobrak.example_models import data_toy_model
from cobrak.thermokinetic_data_retrieval import (
    automatically_add_database_thermokinetic_data_to_cobrak_model,
)

cobrak_model = automatically_add_database_thermokinetic_data_to_cobrak_model(
    data_toy_model,
    database_data_path="./database_data",
    brenda_version="2025_1",  # or 2024_1 if you use the older BRENDA database version
    base_species="Escherichia coli",
)
```

Again, we can look at the added data (ΔG'°, $k_{cat}$, $k_M$ and $MW$) in our model:

```py
# ...following the last code block
from cobrak.printing import print_model

print_model(cobrak_model=cobrak_model)
```

## Automatically collecting some of the data

!!! warning
    Before you can effectively use these functions, make sure that your model uses identifiers/annotations as explained in the [second-to-last subchapter](#additional-model-requirements-for-automatic-data-collection).

### ΔG'° values

To add ΔG'° to a COBRA-k model, use the function ```cobrak.thermokinetic_data_retrieval.get_database_dG0s_for_cobrak_model```. It returns a dictionary with the ΔG'° values that you can then add as explained in subchapter ["Manually adding data"](#manually-adding-data)

### $k_{cat}$ and $k_{M}$ values

To add $k_{cat}$ and $k_{M}$ values to a COBRA-k model, use the function ```cobrak.thermokineticdata_retrieval.get_database_kcats_kms_kis_and_kas_for_cobrak_model```. It returns a dictionary with the kinetic values that you can then add as explained in subchapter ["Manually adding data"](#manually-adding-data)

### $MW$ values

To add $MW$ to a COBRA-k model, use the function `cobrak.thermokinetic_data_retrieval.get_database_mws_for_cobrak_model`. It returns a dictionary with the $MW$ values that you can then add as explained in subchapter ["Manually adding data"](#manually-adding-data)

### Alternative functions for an SBML

If you want to automatically create only a select amount of data directly for an SBML, look up COBRA-k's submodules
`equilibrator_functionality` (for $Δ_r G^{'°}$), `uniprot_functionality` (for molecular enzyme weights),
`sabio_rk_functionality` (for enzyme kinetic data from SABIO-RK), `brenda_functionality` (for enzyme kinetic
data from BRENDA) and `ncbi_taxonomy_functionality` (for taxonomy distance data). Their functions are also described in this documentation's API reference.

## Manually adding data

!!! warning
    Before using thermokinetic data, make sure that your model is "fullsplit" as explained in subsection [Model requirements for any COBRA-k model](#model-requirements-for-any-cobra-k-model).

If you have some or all thermokinetic data already ready, make sure that you have it in the following type form:

* ΔG'° as ```dict[str, float]``` with reaction IDs as keys and ΔG'° as values in kJ⋅mol⁻¹
* ΔG'° uncertainties as ```dict[str, float]```  with reaction IDs as keys and ΔG'° uncertainties as values in kJ⋅mol⁻¹
* $MW$ as ```dict[str, float]```  with enzyme IDs as keys and the weights in kDa
* $k_{cat}$, $k_{M}$ and so on as ```dict[str, EnzymeReactionData]```  with reaction IDs as keys and the kinetic data as `EnzymeReactionData` instances (for more about this dataclass, see "Building metabolic models from scratch" chapter)

!!! info
    The methods presented in [Automatically collecting some of the data](#automatically-collecting-some-of-the-data) already
    adhere to these types :D

Now that you have such data, we can use the major `add_thermokinetic_data_to_cobrak_model` function in `cobrak.thermokinetic_data_retrieval`, which is defined as follows:

```py
def add_thermokinetic_data_to_cobrak_model(
    cobrak_model: Model,
    mws: dict[str, float] = {},
    kcats: dict[str, float] = {},
    kms: dict[str, dict[str, float]] = {},
    kis: dict[str, dict[str, float]] = {},
    kas: dict[str, dict[str, float]] = {},
    dG0s: dict[str, float] = {},
    dG0_uncertainties: dict[str, float] = {},
    conc_ranges: dict[str, tuple[float, float]] = {},
    delete_old_dG0s: bool = False,
    overwrite_existing_dG0s: bool = True,
    overwrite_existing_enzyme_reaction_data: bool = True,
```

Just add the data that you want to add to your COBRA-k model as corresponding argument :-)
