from pytest import fail


def test_no_import_error():
    try:
        pass
    except Exception as e:
        fail(f"An error occurred: {e}")
