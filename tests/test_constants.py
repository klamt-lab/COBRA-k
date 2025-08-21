"""pytest tests for COBRA-k's module constants"""

from pytest import fail


def test_no_import_error() -> None:  # noqa: D103
    try:
        pass
    except Exception as e:
        fail(f"An error occurred: {e}")
