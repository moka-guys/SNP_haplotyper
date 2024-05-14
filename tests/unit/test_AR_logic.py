from io import StringIO

import pandas as pd
import pytest
from EnumDataClasses import InheritanceMode, Relationship, Sex, Status
from helper_functions import set_inherited_from_category_dtype, set_risk_category_dtype
from inheritance_logic import AutosomalRecessiveLogic
from pandas import testing as tm


@pytest.fixture
def setup_all_combination_of_inputs_AR():
    # Dictionary of all possible combinations of haplotypes for the trio
    d = {
        "female_partner": [
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
        "male_partner": [
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


@pytest.mark.autosomal_recessive_logic
@pytest.mark.ref_affected
def test_ref_affected_AR(setup_all_combination_of_inputs_AR):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AR)
    autosomal_recessive_logic = AutosomalRecessiveLogic(
        df=test_df,
        male_partner="male_partner",
        female_partner="female_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.GRANDPARENT,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_recessive_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "low_risk",
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
                "low_risk",
                "uninformative",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "uninformative",
                "uninformative",
                "high_risk",
                "high_risk",
                "high_risk",
                "uninformative",
            ],
            "snp_inherited_from": [
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "male_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "male_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "female_partner",
                "male_partner",
                "male_partner",
                "unassigned",
            ],
        }
    )

    expected_results_df = set_risk_category_dtype(expected_results_df)
    expected_results_df = set_inherited_from_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
    tm.assert_series_equal(
        results_df["snp_inherited_from"], expected_results_df["snp_inherited_from"]
    )


@pytest.mark.autosomal_recessive_logic
@pytest.mark.ref_unaffected
def test_ref_unaffected_AR(setup_all_combination_of_inputs_AR):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AR)

    autosomal_recessive_logic = AutosomalRecessiveLogic(
        df=test_df,
        male_partner="male_partner",
        female_partner="female_partner",
        reference="reference",
        reference_status=Status.UNAFFECTED,
        reference_relationship=Relationship.GRANDPARENT,
        reference_sex=Sex.UNKNOWN,
        consanguineous=False,
    )
    results_df = autosomal_recessive_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "high_risk",
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
                "high_risk",
                "uninformative",
                "high_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "low_risk",
                "uninformative",
                "uninformative",
                "low_risk",
                "low_risk",
                "low_risk",
                "uninformative",
            ],
            "snp_inherited_from": [
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "male_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "male_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "female_partner",
                "male_partner",
                "male_partner",
                "unassigned",
            ],
        }
    )

    expected_results_df = set_risk_category_dtype(expected_results_df)
    expected_results_df = set_inherited_from_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
    tm.assert_series_equal(
        results_df["snp_inherited_from"], expected_results_df["snp_inherited_from"]
    )


@pytest.mark.autosomal_recessive_logic
@pytest.mark.ref_affected
@pytest.mark.consanguineous
def test_ref_affected_consang_AR(setup_all_combination_of_inputs_AR):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AR)
    autosomal_recessive_logic = AutosomalRecessiveLogic(
        df=test_df,
        male_partner="male_partner",
        female_partner="female_partner",
        reference="reference",
        reference_status=Status.AFFECTED,
        reference_relationship=Relationship.GRANDPARENT,
        reference_sex=Sex.UNKNOWN,
        consanguineous=True,
    )
    results_df = autosomal_recessive_logic.get_results()

    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category_AB": [
                "uninformative",
                "uninformative",
                "low_risk",
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
                "low_risk",
                "uninformative",
                "low_risk",
                "uninformative",
                "uninformative",
                "uninformative",
                "high_risk",
                "uninformative",
                "uninformative",
                "high_risk",
                "high_risk",
                "high_risk",
                "uninformative",
            ],
            "snp_risk_category_AA": [
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
            ],
            "snp_risk_category_BB": [
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
            ],
            "snp_inherited_from": [
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "unassigned",
                "male_partner",
                "unassigned",
                "both_partners",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "male_partner",
                "both_partners",
                "unassigned",
                "unassigned",
                "female_partner",
                "unassigned",
                "unassigned",
                "female_partner",
                "male_partner",
                "male_partner",
                "unassigned",
            ],
        }
    )

    expected_results_df = set_risk_category_dtype(expected_results_df)
    expected_results_df = set_inherited_from_category_dtype(expected_results_df)
    tm.assert_series_equal(
        results_df["snp_risk_category_AB"], expected_results_df["snp_risk_category_AB"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_AA"], expected_results_df["snp_risk_category_AA"]
    )
    tm.assert_series_equal(
        results_df["snp_risk_category_BB"], expected_results_df["snp_risk_category_BB"]
    )

    tm.assert_series_equal(
        results_df["snp_inherited_from"], expected_results_df["snp_inherited_from"]
    )
