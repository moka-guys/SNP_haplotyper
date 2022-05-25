import numpy as np
import pandas as pd


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
    if consanguineous and reference_status == "affected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria 1
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria 2
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & ((df[male_partner] == "AA") | (df[male_partner] == "BB"))
            & (df[female_partner] == "AB"),
            # Criteria 3
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if consanguineous and reference_status == "affected":
            values = [
                "male_partner_low_risk",  # Criteria 1
                "female_partner_low_risk",  # Criteria 2
                "shared_high_risk",  # Criteria 3
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")

    elif reference_status == "affected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria 1
            (df[reference] == "AA")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria 2
            (df[reference] == "BB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
            # Criteria 3
            (df[reference] == "AA")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria 4
            (df[reference] == "BB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria 5
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AA"),
            # Criteria 6
            (df[reference] == "AB")
            & (df[male_partner] == "AB")
            & (df[female_partner] == "BB"),
            # Criteria 7
            (df[reference] == "AB")
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria 8
            (df[reference] == "AB")
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria 9
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & (df[female_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "male_partner_low_risk",  # Criteria 1
                "male_partner_low_risk",  # Criteria 2
                "female_partner_low_risk",  # Criteria 3
                "female_partner_low_risk",  # Criteria 4
                "male_partner_high_risk",  # Criteria 5
                "male_partner_high_risk",  # Criteria 6
                "female_partner_high_risk",  # Criteria 7
                "female_partner_high_risk",  # Criteria 8
                "low_&_high_risk",  # Criteria 9 TODO check if this should be split into low/high high/low?
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")

    elif reference_status == "unaffected":
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label low Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria 1
            ((df[reference] == "AA") | (df[reference] == "BB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria 2
            ((df[reference] == "AB"))
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "BB")),
            # Criteria 3
            ((df[reference] == "AA"))
            & (df[male_partner] == "AA")
            & (df[female_partner] == "AB"),
            # Criteria 4
            ((df[reference] == "BB"))
            & (df[male_partner] == "BB")
            & (df[female_partner] == "AB"),
            # Criteria 5
            (
                (df[reference] == "AB")
                | (df[reference] == "AA")
                | (df[reference] == "BB")
            )
            & (df[male_partner] == "AB")
            & ((df[female_partner] == "AA") | (df[female_partner] == "AA")),
            # Criteria 6
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
                "male_parther_high_risk",  # Criteria 1
                "male_parther_low_risk",  # Criteria 2
                "female_parther_high_risk",  # Criteria 3
                "female_parther_low_risk",  # Criteria 4
                "low_and_high_risk",  # Criteria 5 TODO check if this should be split into low/high high/low?
                "low_and_high_risk",  # Criteria 6 TODO check if this should be split into low/high high/low?
            ]

        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")

    return df
