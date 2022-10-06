import numpy as np
import pandas as pd
from exceptions import ArgumentInputError


def autosomal_recessive_analysis(
    df,
    male_partner,
    female_partner,
    reference,
    reference_status,
    consanguineous,
):
    """
    TODO: Add docstring
    """
    # Cosanguineous samples should always have an affected reference
    if consanguineous and reference_status == "unaffected":
        raise ArgumentInputError(
            "Unexpected Input: unaffected reference status should not be used if Cosanguineous = true, check input parameters"
        )

    if reference_status == "affected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria ARto label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AR1
            (df[reference] == "AA")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria AR2
            (df[reference] == "AA")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria AR3 - consanguineous SNPs only
            (df[reference] == "AA")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB")
            & consanguineous,
            # Criteria AR4
            (df[reference] == "BB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria AR5
            (df[reference] == "BB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
            # Criteria AR6 - consanguineous SNPs only
            (df[reference] == "BB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB")
            & consanguineous,
            # Criteria AR7
            (df[reference] == "AB")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria AR8
            (df[reference] == "AB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria AR9
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria AR10
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "low_risk",  # Criteria AR1
                "low_risk",  # Criteria AR2
                "high_risk",  # Criteria AR3
                "low_risk",  # Criteria AR4
                "low_risk",  # Criteria AR5
                "high_risk",  # Criteria AR6
                "high_risk",  # Criteria AR7
                "high_risk",  # Criteria AR8
                "high_risk",  # Criteria AR9
                "high_risk",  # Criteria AR10
            ]

        elif reference_status == "unaffected":
            values = [
                "high_risk",  # Criteria AR1
                "high_risk",  # Criteria AR2
                "high_risk",  # Criteria AR3
                "high_risk",  # Criteria AR4
                "high_risk",  # Criteria AR5
                "low_risk",  # Criteria AR6
                "low_risk",  # Criteria AR7
                "low_risk",  # Criteria AR8
                "low_risk",  # Criteria AR9
                "low_risk",  # Criteria AR10
            ]

        snp_inherited_from = [
            "female_partner",  # Criteria AR1
            "male_partner",  # Criteria AR2
            "both_partners",  # Criteria AR3
            "female_partner",  # Criteria AR4
            "male_partner",  # Criteria AR5
            "both_partners",  # Criteria AR6
            "female_partner",  # Criteria AR7
            "female_partner",  # Criteria AR8
            "male_partner",  # Criteria AR9
            "male_partner",  # Criteria AR10
        ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
        df["snp_inherited_from"] = np.select(
            conditions, snp_inherited_from, default="uninformative"
        )

    return df
