# IMPORTS SECTION #
import z_add_path  # noqa: F401

from cobrak.dataclasses import Model
from cobrak.io import json_load, json_write

cobrak_model: Model = json_load(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak_prestepC_uncalibrated.json",
    Model,
)
cobrak_model.max_prot_pool = 0.224  # new protpool

json_write(
    "examples/iCH360/prepared_external_resources/iCH360_cobrak.json",
    cobrak_model,
)
