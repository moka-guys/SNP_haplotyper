import pytest
from helper_functions import custom_order_generator


@pytest.mark.parametrize(
    "input_value, expected",
    [
        (1, ["0-1MB_from_start", "within_gene", "0-1MB_from_end"]),
        (
            2,
            [
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
            ],
        ),
        (
            3,
            [
                "2-3MB_from_start",
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
                "2-3MB_from_end",
            ],
        ),
        (
            4,
            [
                "3-4MB_from_start",
                "2-3MB_from_start",
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
                "2-3MB_from_end",
                "3-4MB_from_end",
            ],
        ),
        (
            5,
            [
                "4-5MB_from_start",
                "3-4MB_from_start",
                "2-3MB_from_start",
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
                "2-3MB_from_end",
                "3-4MB_from_end",
                "4-5MB_from_end",
            ],
        ),
    ],
)
def test_custom_order_generator(input_value, expected):
    result = custom_order_generator(input_value)
    assert (
        result == expected
    ), f"Expected {expected} for input {input_value}, but got {result}"


def test_custom_order_generator_with_zero():
    expected = ["within_gene"]
    result = custom_order_generator(0)
    assert result == expected, f"Expected {expected} for input 0, but got {result}"
