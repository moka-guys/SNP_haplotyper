from io import StringIO

import pandas as pd
import pytest
from EnumDataClasses import InheritanceMode, Relationship, Sex, Status
from helper_functions import set_risk_category_dtype
from inheritance_logic import AutosomalDominantLogic
from pandas import testing as tm


@pytest.fixture
def setup_all_combination_of_inputs_AD():
    # Dictionary of all possible combinations of haplotypes for the trio
    d = {
        "affected_partner": [
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
        "unaffected_partner": [
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


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_affected
@pytest.mark.ref_grandparent
def test_ref_affected_grandparent_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    autosomal_dominant_logic = AutosomalDominantLogic(
        df=test_df,
        affected_partner="affected_partner",
        unaffected_partner="unaffected_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.GRANDPARENT,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_dominant_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "high_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
            ]
        }
    )

    expected_results_df = set_risk_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_unaffected
@pytest.mark.ref_grandparent
def test_ref_unaffected_grandparent_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)

    autosomal_dominant_logic = AutosomalDominantLogic(
        df=test_df,
        affected_partner="affected_partner",
        unaffected_partner="unaffected_partner",
        reference="reference",
        reference_status=Status.UNAFFECTED,
        reference_relationship=Relationship.GRANDPARENT,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_dominant_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "high_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
            ],
        }
    )
    expected_results_df = set_risk_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_affected
@pytest.mark.ref_child
def test_ref_affected_child_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    autosomal_dominant_logic = AutosomalDominantLogic(
        df=test_df,
        affected_partner="affected_partner",
        unaffected_partner="unaffected_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.CHILD,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_dominant_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "high_risk",
                "uninformative",
            ],
        }
    )
    expected_results_df = set_risk_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_unaffected
@pytest.mark.ref_child
def test_ref_unaffected_child_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    autosomal_dominant_logic = AutosomalDominantLogic(
        df=test_df,
        affected_partner="affected_partner",
        unaffected_partner="unaffected_partner",
        reference="reference",
        reference_status=Status.UNAFFECTED,
        reference_relationship=Relationship.CHILD,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_dominant_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "low_risk",
                "uninformative",
            ],
        }
    )
    expected_results_df = set_risk_category_dtype(expected_results_df)

    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
