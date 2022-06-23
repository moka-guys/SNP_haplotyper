import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'


def autosomal_dominant_analysis(
    df,
    affected_partner,
    affected_partner_sex,
    unaffected_partner,
    unaffected_partner_sex,
    reference,
    reference_status,
    reference_relationship,
):
    """
    TODO: Add docstring
    """
    if reference_relationship in [
        "grandparent",  # TODO populate with all appropriate relationships
    ]:
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label high Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AD1
            (df[reference] == "AA")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria AD2
            (df[reference] == "BB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria AD3
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria AD4
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria AD1
                "high_risk",  # Criteria AD2
                "low_risk",  # Criteria AD3
                "low_risk",  # Criteria AD4
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria AD1
                "low_risk",  # Criteria AD2
                "high_risk",  # Criteria AD3
                "high_risk",  # Criteria AD4
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
    elif reference_relationship in [
        "child",  # TODO populate with all appropriate relationships
    ]:
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label high Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AD1
            (df[reference] == "AB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria AD2
            (df[reference] == "AB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria AD3
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria AD4
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria AD1
                "high_risk",  # Criteria AD2
                "low_risk",  # Criteria AD3
                "low_risk",  # Criteria AD4
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria AD1
                "low_risk",  # Criteria AD2
                "high_risk",  # Criteria AD3
                "high_risk",  # Criteria AD4
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")

    return df
