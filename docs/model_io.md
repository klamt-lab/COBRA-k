# Load/Save Models

All functionality regarding loading and saving a Model instance or parts of it can be found in the ```cobrak.io``` submodule. For specific information on its functions, you can also see the ```io``` section in this documentation's "API reference" chapter.

!!! note "Pydantic validation"
    As mentioned in the previous chapter, COBRA-k uses pydantic's [[GitHub]](https://github.com/pydantic/pydantic) validation features for its dataclasses, including ```Model```. If your model file (in whatever format) contains an invalid value (e.g., a negative turnover number $k_cat$ for a Reaction's enzyme reaction data), you'll get a ```ValidationError```from pydantic. The error message shows more details about why the validation failed.

    If this happens, going through pydantic's validation error messages helps you to find these strange values in your model. Once fixed, you can then load the model into COBRA-k with more confidence in the integritiy of its data :-)

There are two formats, JSON and SBML, in which you can interchangeably load and save a COBRA-k model, both of which store it in an human readable way:

## 1. JSON

Any of COBRA-k's dataclasses (see previous chapter), including Model instances, can be stored and loaded as JSON [[Wikipedia]](https://en.wikipedia.org/wiki/JSON) files. JSONs can be quickly saved and loaded tend to be easier to read and have a smaller file size than SBMLs. However, they are not interoperable with other program packages.

!!! note
    While some other constraint-based modeling packages also provide JSON as input and output format, COBRA-k's JSON definition is incompatible with their definitions. To use an interoperable format that works with virtually all packages, use the SBML format (see below).

### Save Model as JSON

Use ```json_write```:

```py
from cobrak.io import json_write

# Assuming that we already have a cobrak_model variable :-)
json_write(
    path="/wished/save/path/filename.json",
    cobrak_model=cobrak_model,
)
```

### Load Model as JSON

Use ```json_load``` using ```Model```as type argument:

```py
from cobrak.io import json_load
from cobrak.dataclasses import Model

# The ": Model" addition and import of the Model dataclass
# are not neccessary, however, they help if you utilize
# Python's typing functionality, which provides e.g.
# automatic completion for model instances.
cobrak_model: Model = json_load(
    path="/path/to/json/model/filename.json",
    Model,
)
```

The type argument of ```json_load```is optional to load a JSON. Using this type, COBRA-k utilizes pydantic to automatically verify the correctness of the given JSON (i.e. that it matches the given type or dataclass). If you're unsure about your JSON's type (e.g. if it is not a COBRA-k dataclass), you can also just use the ```Any```type from Python's ```typing``` package.

## 2. **(Annotated) SBML**

The Systems Biology Markup Language [[Paper]](https://doi.org/10.1093/bioinformatics/btg015) is a widely used format for storing metabolic models. COBRA-k can directly load such models and convert them to Model instances. However, as typical SBML files lack thermodynamic and kinetic extra information, this information has to be added later, as detailed in the next chapter about thermokinetic data retrieval.

To overcome the limitation of lacking thermodynamic and kinetic data in SBML, COBRA-k stores this relevant information in extra annotation fields of the reactions, metabolites and genes (enzymes are interpreted as genes), making it an *annotated SBML*. COBRA-k can then load such an annotated SBML to create a Model with full thermodynamic and kinetic data again. Other program packages such as COBRApy [[Paper]](https://doi.org/10.1186/1752-0509-7-74) can load the annotated SBML. However, they typically cannot directly use most of the extra thermodynamic and kinetic information. This means that only basic analyses (as the ones described in the "Linear Programs" chapter) can be conducted with the other program packages.

### Save as annotated SBML

To save a Model instance as annotated SBML, use ```save_cobrak_model_as_annotated_sbml_model```:

```py
from cobrak.io import save_cobrak_model_as_annotated_sbml_model

# Assuming that we already have a cobrak_model variable :-)
save_cobrak_model_as_annotated_sbml_model(
    cobrak_model,
    filepath="/wished/path/filename.xml",
    combine_base_reactions=False, # Default is False
    add_enzyme_constraints=False, # Default is False
)
```

The ```combine_base_reactions``` argument controls whether split forward and reverse as well as split enzyme (for reactions which are catalyzed by multiple enzymes) reactions are to be merged in the SBML as a single reaction or not. This merging is then automatically reversed when loading an annotated SBML again (see following subchapter).

The ```add_enzyme_constraints``` parameter controls whether or not the exported metabolic model shall include an expansion of the stoichiometric matrix (see chapter about Linear Programs) with enzyme constraints akin to the method GECKO ([Paper](https://doi.org/10.15252/msb.20167411)). If True, pseudo-metabolites representing enzymes (if given, with their respective enzyme concentration bounds) and pseudo-reactions representing the protein pool and the delivery of these enzymes are added to the model, and enzymatically catalyzed reactions must consume the new pseudo-metabolites. For more information, read GECKO's paper. This expansion allows other program packages, including COBRApy, to directly use the enzyme constraints defined by the Model instance.

When loading such an annotated XML again with COBRA-k, the enzyme constraints are automatically detected and the "normal" COBRA-k form of enzyme constraints is used again (as detailed in the chapter about Linear Programs).

### Load (annotated) SBML

Loading an (annotated) SBML as a Model is simple:

```py
from cobrak.io import load_annotated_sbml_model_as_cobrak_model

cobrak_model = load_annotated_sbml_model_as_cobrak_model(
    "/path/to/sbml/filename.xml"
)
```

Keep in mind that only with COBRA-k annotations, COBRA-k can successfully directly read the thermodynamic and kinetic data of the model. Otherwise, t has to be included later to the Model instance. Options to do so are detailed in the next chapter. Without COBRA-k annotations in the SBML, a standard metabolite concentration range of $10^{-6}$ up to 0.02 M is assumed. Also, the protein pool is set to 1_000 g⋅gDW⁻¹, essentially deactivating any possible enzyme constraints.
