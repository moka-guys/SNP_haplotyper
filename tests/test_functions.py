import pytest
import pandas as pd
from autosomal_recessive_logic import autosomal_recessive_analysis


def test_criteria_AR1():
    d = {
        "male_partner": ["AB", "AB", "AB", "AB"],
        "female_partner": ["AA", "AA", "BB", "BB"],
        "reference": ["AA", "BB", "AA", "BB"],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "affected"
    consanguineous = True
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        "ref_status",
        consanguineous,
    )
    assert results_df["snp_risk_category"].unique().tolist() == [
        "male_partner_low_risk"
    ]


def test_criteria_AR2():
    d = {
        "male_partner": ["AA", "AA", "BB", "BB"],
        "female_partner": ["AB", "AB", "AB", "AB"],
        "reference": ["AA", "BB", "AA", "BB"],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "affected"
    consanguineous = True
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        "ref_status",
        consanguineous,
    )
    assert results_df["snp_risk_category"].unique().tolist() == [
        "female_partner_low_risk"
    ]


def test_criteria_AR3():
    d = {
        "male_partner": [
            "AB",
            "AB",
        ],
        "female_partner": [
            "AB",
            "AB",
        ],
        "reference": [
            "AA",
            "BB",
        ],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "affected"
    consanguineous = True
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        "ref_status",
        consanguineous,
    )
    assert results_df["snp_risk_category"].unique().tolist() == ["shared_high_risk"]


def test_criteria_AR4():
    d = {
        "male_partner": [
            "AA",
            "BB",
        ],
        "female_partner": [
            "AB",
            "AB",
        ],
        "reference": [
            "AB",
            "AB",
        ],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "affected"
    consanguineous = True
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        "ref_status",
        consanguineous,
    )
    assert results_df["snp_risk_category"].unique().tolist() == [
        "female_partner_unique_high_risk"
    ]


def test_criteria_AR5():
    d = {
        "male_partner": [
            "AB",
            "AB",
        ],
        "female_partner": [
            "AA",
            "BB",
        ],
        "reference": [
            "AB",
            "AB",
        ],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "affected"
    consanguineous = True
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        "ref_status",
        consanguineous,
    )
    assert results_df["snp_risk_category"].unique().tolist() == [
        "male_partner_unique_high_risk"
    ]
