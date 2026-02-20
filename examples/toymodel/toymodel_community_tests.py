"""Runs all analyses for the toymodel as shown in COBRA-k's initial publication"""
from cobrak.lps import perform_lp_optimization

try:  # noqa: SIM105
    import z_add_path  # noqa: F401
except ModuleNotFoundError:
    pass

from cobrak.example_models import toy_model
from cobrak.io import json_write
from cobrak._community import create_multiplied_model, create_community_model_with_fixed_growth

multiplied_toy_model = create_multiplied_model(
    {"ONE": toy_model, "TWO": toy_model},
)
fixed_growth_model = create_community_model_with_fixed_growth(
    {
        "ONE": (toy_model, "ATP_Consumption"),
        "TWO": (toy_model, "ATP_Consumption"),
    },
    growth_rate=0.5,
)
json_write("examples/toymodel/community_multiplied_toy_model.json", multiplied_toy_model)
result = perform_lp_optimization(multiplied_toy_model, {"ATP_Consumption_ONE" : 1.0, "ATP_Consumption_TWO": -0.1}, +1)
json_write("examples/toymodel/xxx.json", result)
