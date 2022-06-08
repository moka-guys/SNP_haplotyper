import numpy as np
import pandas as pd
from exceptions import ArgumentInputError


def autosomal_recessive_analysis(
    df,
    male_partner,
    female_partner,
    reference,
    reference_status,
    reference_relationship,
    consanguineous,
):
    """
    TODO: Add docstring
    """
    # Process consanguineus with an affected reference (would not be unaffected) TODO add condition to catch that and issue warning
    if consanguineous and reference_status == "affected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria ARto label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AR1
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria AR2
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & ((df[male_partner] == "AA") | (df[male_partner] == "BB"))
            & (df[female_partner] == "AB"),
            # Criteria AR3
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB"),
            # Criteria AR4
            (df[reference] == "AB")
            & ((df[male_partner] == "AA") | (df[male_partner] == "BB"))
            & (df[female_partner] == "AB"),
            # Criteria AR5
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
        ]
        # Assign the correct labels depending upon reference status
        if consanguineous and reference_status == "affected":
            values = [
                "low_risk",  # Criteria AR1
                "low_risk",  # Criteria AR2
                "high_risk",  # Criteria AR3
                "unique_high_risk",  # Criteria AR4
                "unique_high_risk",  # Criteria AR5
            ]

            snp_inherited_from = [
                "male_partner",  # Criteria AR1
                "female_partner",  # Criteria AR2
                "shared",  # Criteria AR3
                "female_partner",  # Criteria AR4
                "male_partner",  # Criteria AR5
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
        df["snp_inherited_from"] = np.select(
            conditions, snp_inherited_from, default="uninformative"
        )

    elif consanguineous and reference_status == "unaffected":
        raise ArgumentInputError(
            "Unexpected Input: unaffected reference status should not be used if Cosanguineous = true, check input parameters"
        )

    # process non-consanguineous with affected reference
    elif reference_status == "affected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria ARto label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AR7
            (df[reference] == "AA")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria AR8
            (df[reference] == "BB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
            # Criteria AR9
            (df[reference] == "AA")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria AR10
            (df[reference] == "BB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria AR11
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria AR12
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
            # Criteria AR13
            (df[reference] == "AB")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria AR14
            (df[reference] == "AB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria AR15
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "low_risk",  # Criteria AR7
                "low_risk",  # Criteria AR8
                "low_risk",  # Criteria AR9
                "low_risk",  # Criteria AR10
                "high_risk",  # Criteria AR11
                "high_risk",  # Criteria AR12
                "high_risk",  # Criteria AR13
                "high_risk",  # Criteria AR14
                "low_&_high_risk",  # Criteria AR15
            ]

            snp_inherited_from = [
                "male_partner",  # Criteria AR7
                "male_partner",  # Criteria AR8
                "female_partner",  # Criteria AR9
                "female_partner",  # Criteria AR10
                "male_partner",  # Criteria AR11
                "male_partner",  # Criteria AR12
                "female_partner",  # Criteria AR13
                "female_partner",  # Criteria AR14
                "shared",  # Criteria AR15 TODO check that shared is correct
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
        df["snp_inherited_from"] = np.select(
            conditions, snp_inherited_from, default="uninformative"
        )
    # process non-consanguineous with unaffected reference
    elif reference_status == "unaffected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria ARto label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria AR16
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria AR17
            ((df[reference] == "AB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria AR18
            ((df[reference] == "AA"))
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria AR19
            ((df[reference] == "BB"))
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria AR20
            (
                (df[reference] == "AB")
                | (df[reference] == "AA")
                | (df[reference] == "BB")
            )
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "AA")),
            # Criteria AR21
            (
                (df[reference] == "AB")
                | (df[reference] == "AA")
                | (df[reference] == "BB")
            )
            & ((df[male_partner] == "AA") | (df[reference] == "BB"))
            & (df[female_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "unaffected":
            values = [
                "high_risk",  # Criteria AR16
                "low_risk",  # Criteria AR17
                "high_risk",  # Criteria AR18
                "low_risk",  # Criteria AR19
                "low_and_high_risk",  # Criteria AR20 TODO check if this should be split into low/high high/low?
                "low_and_high_risk",  # Criteria AR21 TODO check if this should be split into low/high high/low?
            ]

            snp_inherited_from = [
                "male_parther",  # Criteria AR16
                "male_parther",  # Criteria AR17
                "female_parther",  # Criteria AR18
                "female_parther",  # Criteria AR19
                "shared",  # Criteria AR20 TODO check if this is correct
                "shared",  # Criteria AR21 TODO check if this is correct
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
        df["snp_inherited_from"] = np.select(
            conditions, snp_inherited_from, default="uninformative"
        )
    return df
