import pandas as pd
import pytest
from exceptions import ArgumentInputError


def test_ArgumentInputError():
    d = {
        "male_partner": [
            "AB",
        ],
        "female_partner": ["AA"],
        "reference": [
            "AA",
        ],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "unaffected"
    consanguineous = True

    # Should throw exception if reference_status = "unaffected" and consanguineous = True
    with pytest.raises(ArgumentInputError):
        results_df = autosomal_recessive_analysis(
            test_df,
            "male_partner",
            "female_partner",
            "reference",
            reference_status,
            consanguineous,
        )
