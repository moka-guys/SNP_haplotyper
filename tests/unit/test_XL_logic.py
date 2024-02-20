from io import StringIO

import pandas as pd
import pytest
from EnumDataClasses import InheritanceMode, Relationship, Sex, Status
from helper_functions import set_inherited_from_category_dtype, set_risk_category_dtype
from inheritance_logic import XLinkedLogic
from pandas import testing as tm


@pytest.fixture
def setup_all_combination_of_inputs_XL():
    # Dictionary of all possible combinations of haplotypes for the trio
    d = {
        "carrier_female_partner": [
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
        ],
        "unaffected_male_partner": [
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
        ],
        "reference": [
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
        ],
    }
    return d


@pytest.mark.x_linked_logic
@pytest.mark.ref_affected
@pytest.mark.ref_child
def test_ref_affected_male_ref_XL(setup_all_combination_of_inputs_XL):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_XL)

    x_link_logic = XLinkedLogic(
        df=test_df,
        unaffected_male_partner="unaffected_male_partner",
        carrier_female_partner="carrier_female_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.CHILD,
        reference_sex=Sex.MALE,
        consanguineous=False,
    )
    results_df = x_link_logic.get_results()

    data = """
    reference,affected_partner,unaffected_partner,snp_risk_category_AB,snp_risk_category_AA,snp_risk_category_BB
    AA,AA,AA,uninformative,uninformative,uninformative
    AA,AA,BB,uninformative,uninformative,uninformative
    AA,AA,AB,uninformative,uninformative,uninformative
    AA,BB,AA,uninformative,uninformative,uninformative
    AA,BB,BB,uninformative,uninformative,uninformative
    AA,BB,AB,uninformative,uninformative,uninformative
    AA,AB,AA,low_risk,high_risk,low_risk
    AA,AB,BB,high_risk,high_risk,low_risk
    AA,AB,AB,uninformative,uninformative,uninformative
    BB,AA,AA,uninformative,uninformative,uninformative
    BB,AA,BB,uninformative,uninformative,uninformative
    BB,AA,AB,uninformative,uninformative,uninformative
    BB,BB,AA,uninformative,uninformative,uninformative
    BB,BB,BB,uninformative,uninformative,uninformative
    BB,BB,AB,uninformative,uninformative,uninformative
    BB,AB,AA,high_risk,low_risk,high_risk
    BB,AB,BB,low_risk,low_risk,high_risk
    BB,AB,AB,uninformative,uninformative,uninformative
    AB,AA,AA,uninformative,uninformative,uninformative
    AB,AA,BB,uninformative,uninformative,uninformative
    AB,AA,AB,uninformative,uninformative,uninformative
    AB,BB,AA,uninformative,uninformative,uninformative
    AB,BB,BB,uninformative,uninformative,uninformative
    AB,BB,AB,uninformative,uninformative,uninformative
    AB,AB,AA,uninformative,uninformative,uninformative
    AB,AB,BB,uninformative,uninformative,uninformative
    AB,AB,AB,uninformative,uninformative,uninformative
    """

    # Using StringIO to simulate a file object
    data_io = StringIO(data)

    # Reading the data into a pandas DataFrame
    expected_results_df = pd.read_csv(data_io)
    expected_results_df = set_risk_category_dtype(expected_results_df)

    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_AA"], expected_results_df["snp_risk_category_AA"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_BB"], expected_results_df["snp_risk_category_BB"]
    )


@pytest.mark.x_linked_logic
@pytest.mark.ref_affected
@pytest.mark.ref_child
def test_ref_affected_female_ref_XL(setup_all_combination_of_inputs_XL):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_XL)

    x_link_logic = XLinkedLogic(
        df=test_df,
        unaffected_male_partner="unaffected_male_partner",
        carrier_female_partner="carrier_female_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.CHILD,
        reference_sex=Sex.FEMALE,
        consanguineous=False,
    )
    results_df = x_link_logic.get_results()

    data = """
    reference,affected_partner,unaffected_partner,snp_risk_category_AB,snp_risk_category_AA,snp_risk_category_BB
    AA,AA,AA,uninformative,uninformative,uninformative
    AA,AA,BB,uninformative,uninformative,uninformative
    AA,AA,AB,uninformative,uninformative,uninformative
    AA,BB,AA,uninformative,uninformative,uninformative
    AA,BB,BB,uninformative,uninformative,uninformative
    AA,BB,AB,uninformative,uninformative,uninformative
    AA,AB,AA,low_risk,high_risk,low_risk
    AA,AB,BB,high_risk,high_risk,low_risk
    AA,AB,AB,uninformative,uninformative,uninformative
    BB,AA,AA,uninformative,uninformative,uninformative
    BB,AA,BB,uninformative,uninformative,uninformative
    BB,AA,AB,uninformative,uninformative,uninformative
    BB,BB,AA,uninformative,uninformative,uninformative
    BB,BB,BB,uninformative,uninformative,uninformative
    BB,BB,AB,uninformative,uninformative,uninformative
    BB,AB,AA,high_risk,low_risk,high_risk
    BB,AB,BB,low_risk,low_risk,high_risk
    BB,AB,AB,uninformative,uninformative,uninformative
    AB,AA,AA,uninformative,uninformative,uninformative
    AB,AA,BB,uninformative,uninformative,uninformative
    AB,AA,AB,uninformative,uninformative,uninformative
    AB,BB,AA,uninformative,uninformative,uninformative
    AB,BB,BB,uninformative,uninformative,uninformative
    AB,BB,AB,uninformative,uninformative,uninformative
    AB,AB,AA,high_risk,low_risk,high_risk
    AB,AB,BB,high_risk,high_risk,low_risk
    AB,AB,AB,uninformative,uninformative,uninformative
    """

    # Using StringIO to simulate a file object
    data_io = StringIO(data)
    # Reading the data into a pandas DataFrame
    expected_results_df = pd.read_csv(data_io)
    expected_results_df = set_risk_category_dtype(expected_results_df)

    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_AA"], expected_results_df["snp_risk_category_AA"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_BB"], expected_results_df["snp_risk_category_BB"]
    )
