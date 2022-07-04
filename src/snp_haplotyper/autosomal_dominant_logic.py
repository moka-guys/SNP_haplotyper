import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'


def autosomal_dominant_analysis(
    df,
    affected_partner,
    unaffected_partner,
    reference,
    reference_status,
    reference_relationship,
):
    """Identifies any SNPs which are "informative" and then categorises them as "high_risk" or "low_risk"
    # TODO change google doc link to link to read the docs.
    The full logic behind the function is described here https://docs.google.com/document/d/1ZsdSTQ_oliSDM-1EO65XbKfifaPK2-lIPwK1Ksl8L6s/edit?usp=sharing
    This function asks the question - For the combination of haplotypes present in the reference trio (reference,
    unaffected_partner and affected_partner) what information would this SNP provide us with if the SNP was
    AB in the embryo - does the SNP indicate "high_risk", "low_risk", or is it "uninformative"?


    Args:
        df (dataframe): A dataframe which requires rendering as HTML for inclusion in the HTML report
            affected_partner (string): Column name in dataframe refering to affected_partner's data
            unaffected_partner (string): Column name in dataframe refering to unaffected_partner's data
            reference (string) : Column name in dataframe refering to reference's data
            reference_status (string) : "affected" or "unaffected"
            reference_relationship (string) : "grandparent" or "child"

    Returns:
        dataframe: Dataframe with a "snp_risk_category" column added, used to categorise the SNPs as
        "high_risk", "low_risk", and "uninformative"

    """
    # NOTE: Lines marked with a Critera_AD# and Option_AD# ID  refer to passages in the specification for
    # this project. During code review this allows the programs logic to be easily compared to the spec.

    if reference_relationship in [
        "grandparent",  # TODO populate with all appropriate relationships
    ]:
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label high Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria_AD1
            (df[reference] == "AA")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria_AD2
            (df[reference] == "BB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria_AD3
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria_AD4
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria_AD1, Option_AD1
                "high_risk",  # Criteria_AD2, Option_AD2
                "low_risk",  # Criteria_AD3, Option_AD3
                "low_risk",  # Criteria_AD4, Option_AD4
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria_AD1, Option_AD5
                "low_risk",  # Criteria_AD2, Option_AD6
                "high_risk",  # Criteria_AD3, Option_AD7
                "high_risk",  # Criteria_AD4, Option_AD8
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
    elif reference_relationship in [
        "child",
    ]:
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label high Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria_AD5
            (df[reference] == "AB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria_AD6
            (df[reference] == "AB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria_AD7
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria_AD8
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria_AD5, Option_AD9
                "high_risk",  # Criteria_AD6, Option_AD10
                "low_risk",  # Criteria_AD7, Option_AD11
                "low_risk",  # Criteria_AD8, Option_AD12
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria_AD5, Option_AD13
                "low_risk",  # Criteria_AD6, Option_AD14
                "high_risk",  # Criteria_AD7, Option_AD15
                "high_risk",  # Criteria_AD8, Option_AD16
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")

    return df
