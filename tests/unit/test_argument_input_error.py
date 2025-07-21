import pandas as pd
import pytest
from EnumDataClasses import Relationship, Sex, Status
from exceptions import ArgumentInputError
from inheritance_logic import AutosomalRecessiveLogic


def test_ArgumentInputError(capsys):
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
    reference_status = Status.UNAFFECTED
    reference_relationship = Relationship.GRANDPARENT
    reference_sex = Sex.FEMALE
    consanguineous = True

    # Should throw exception if reference_status = "unaffected" and consanguineous = True
    captured = AutosomalRecessiveLogic(
            test_df,
            "male_partner",
            "female_partner",
            "reference",
            reference_status,
            reference_relationship,
            reference_sex,
            consanguineous,
        )
    captured = capsys.readouterr()
    assert "Cannot have unaffected reference in consanguineous analysis" in captured.out
