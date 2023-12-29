from helper_functions import custom_order_generator
from SNPDataClass import SNPData
import pandas as pd


def test_annotate_distances():
    # Define a sample series
    distances = pd.Series([-4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

    # Run the function on the sample series
    result = SNPData.annotate_distances(distances)

    # Define expected result
    expected = pd.Categorical(
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
            "4-5MB_from_end",
        ],
        categories=custom_order_generator(5),
        ordered=True,
    )

    # Assert that the result matches the expected categories
    assert all(result == expected), f"Expected {expected}, but got {result}"
