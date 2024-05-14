import pandas as pd
from EmbryoDataClass import update_embryo_risk_column
from EnumDataClasses import InheritanceMode
import pytest


def load_test_data():
    return pd.read_csv("test_data/ado_or_miscall_test_data.csv")


def test_ado_or_miscall_AD_logic():
    test_data_df = load_test_data()
    test_data_df["reference"] = test_data_df[
        "female_partner"
    ]  # Duplicate female partner to pass as reference (Reference column is only checked for 'NoCall' status)
    AD_basher_results = update_embryo_risk_column(
        df=test_data_df,
        risk_col_name="calculated_result",
        male_partner_col="male_partner",
        female_partner_col="female_partner",
        reference_col="reference",
        embryo_col="embryo",
        mode_of_inheritance=InheritanceMode.AUTOSOMAL_DOMINANT,
    )
    assert AD_basher_results["calculated_result"].equals(
        AD_basher_results["AD_expected_result"]
    ), f"Mismatch between calculated and expected results in AD logic."


def test_ado_or_miscall_AR_logic():
    test_data_df = load_test_data()
    test_data_df["reference"] = test_data_df[
        "female_partner"
    ]  # Duplicate female partner to pass as reference (Reference column is only checked for 'NoCall' status)
    AR_basher_results = update_embryo_risk_column(
        df=test_data_df,
        risk_col_name="calculated_result",
        male_partner_col="male_partner",
        female_partner_col="female_partner",
        reference_col="reference",
        embryo_col="embryo",
        mode_of_inheritance=InheritanceMode.AUTOSOMAL_RECESSIVE,
    )
    assert AR_basher_results["calculated_result"].equals(
        AR_basher_results["AR_expected_result"]
    ), f"Mismatch between calculated and expected results in AR logic."


def test_nocall_reference_AD_logic():
    test_data_df = load_test_data()
    test_data_df["reference"] = "NoCall"
    test_data_df["AD_expected_result"] = "NoCall"
    AD_basher_results = update_embryo_risk_column(
        df=test_data_df,
        risk_col_name="calculated_result",
        male_partner_col="male_partner",
        female_partner_col="female_partner",
        reference_col="reference",
        embryo_col="embryo",
        mode_of_inheritance=InheritanceMode.AUTOSOMAL_DOMINANT,
    )
    assert AD_basher_results["calculated_result"].equals(
        AD_basher_results["AD_expected_result"]
    ), f"If reference is a NoCall then result should be NoCall."


def test_nocall_reference_AR_logic():
    test_data_df = load_test_data()
    test_data_df["reference"] = "NoCall"
    test_data_df["AR_expected_result"] = "NoCall"
    AR_basher_results = update_embryo_risk_column(
        df=test_data_df,
        risk_col_name="calculated_result",
        male_partner_col="male_partner",
        female_partner_col="female_partner",
        reference_col="reference",
        embryo_col="embryo",
        mode_of_inheritance=InheritanceMode.AUTOSOMAL_RECESSIVE,
    )
    assert AR_basher_results["calculated_result"].equals(
        AR_basher_results["AR_expected_result"]
    ), f"If reference is a NoCall then result should be NoCall."
