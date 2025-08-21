from cobrak.expasy_functionality import get_ec_number_transfers  # noqa: D100
from cobrak.io import json_write

ec_number_transfers = get_ec_number_transfers(
    "examples/common_needed_external_resources/enzyme.rdf"
)
json_write(
    "examples/common_needed_external_resources/ec_number_transfers.json",
    ec_number_transfers,
)
