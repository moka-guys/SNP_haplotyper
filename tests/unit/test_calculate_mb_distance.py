import pandas as pd
import pytest

from SNPDataClass import SNPData


def test_calculate_mb_distance():
    # Sample gene positions for the test
    gene_positions = pd.Series(
        [
            999_999,
            1_000_000,
            1_000_001,
            1_999_999,
            2_000_000,
            2_000_001,
            2_999_999,
            3_000_000,
            3_000_001,
            3_999_999,
            4_000_000,
            4_000_001,
            4_999_999,
            5_000_000,
            5_000_001,
            5_999_999,
            6_000_000,
            6_000_001,
            6_999_999,
            7_000_000,
            7_000_001,
            7_999_999,
            8_000_000,
            8_000_001,
            8_999_999,
            9_000_000,
            9_000_001,
        ]
    )
    gene_start = 4_000_000
    gene_end = 5_000_000

    expected = pd.Series(
        [
            -4,
            -3,
            -3,
            -3,
            -2,
            -2,
            -2,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            2,
            2,
            2,
            3,
            3,
            3,
            4,
            4,
            4,
            5,
        ]
    )

    # Calling the function
    result = SNPData.calculate_mb_distance(gene_positions, gene_start, gene_end)

    # Asserting that the result matches the expected values
    assert all(
        result == expected
    ), f"Expected {expected.tolist()}, but got {result.tolist()}"
