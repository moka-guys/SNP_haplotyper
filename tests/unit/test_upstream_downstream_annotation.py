import pandas as pd
import pytest
from SNPDataClass import SNPData


def test_categorize_gene_regions():
    # Define sample data for the test
    gene_positions = pd.Series(
        [
            0,
            999_999,
            1_000_000,
            1_000_001,
            2_999_999,
            3_000_000,
            3_000_001,
            3_999_999,
            4_000_000,
            4_000_001,
            5_999_999,
            6_000_000,
            6_000_001,
            7_000_000,
        ]
    )
    gene_start = 3_000_000
    gene_end = 4_000_000
    max_region_mb = 2  # This is equivalent to 1 million base pairs

    # Call the function
    result = SNPData.categorize_SNP_by_position(
        gene_positions, gene_start, gene_end, max_region_mb
    )
    # Expected categories based on the sample data
    expected = pd.Categorical(
        [
            "outside_ROI",
            "outside_ROI",
            "upstream",
            "upstream",
            "upstream",
            "within_gene",
            "within_gene",
            "within_gene",
            "within_gene",
            "downstream",
            "downstream",
            "downstream",
            "outside_ROI",
            "outside_ROI",
        ],
        categories=["upstream", "within_gene", "downstream", "outside_ROI"],
        ordered=True,
    )

    # Assert that the result matches the expected values
    assert all(result == expected), f"Expected {expected}, but got {result}"
