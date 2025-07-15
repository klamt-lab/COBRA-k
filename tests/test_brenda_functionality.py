from cobrak.brenda_functionality import _brenda_parse_full_json, _is_fitting_ec_numbers


def test_brenda_parse_full_json():
    _brenda_parse_full_json(
        brenda_json_targz_file_path="examples/common_needed_external_resources/brenda_2024_1.json.tar.gz",
        bigg_metabolites_json_path="examples/common_needed_external_resources/bigg_models_metabolites.json",
        brenda_version="2024_1",
        min_ph=0.0,
        max_ph=150.0,
        accept_nan_ph=True,
        min_temperature=0.0,
        max_temperature=100.0,
        accept_nan_temperature=True,
    )


def test_is_fitting_ec_numbers():
    assert not _is_fitting_ec_numbers("1.1.1.1", "1.2.2.2", 2)
    assert _is_fitting_ec_numbers("1.1.1.1", "1.2.2.2", 3)
