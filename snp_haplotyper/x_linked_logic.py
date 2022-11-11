import numpy as np
import pandas as pd

# Logic is the same whether using affected son of carrier or grandmother as reference
def x_linked_analysis(
    df,
    carrier_female_partner,  # Always female partner for X-linked
    unaffected_male_partner,
    reference,  # Either carrier or affected
):
    """Identifies any SNP site which could be used to inform a decision regarding inheriting an X-linked condition ("informative" SNPs) and categorizes
    the site as indicating "high_risk" or "low_risk" of inheriting an X-linked condition for a known haplotype in the embryo (See below for details).
    Identifies any SNPs which are "informative" and then categorises them as "high_risk" or "low_risk" for any cases which are x-linked inheritance.
    # TODO change google doc link to link to read the docs.
    The full logic behind the function is described here https://docs.google.com/document/d/1QKdQ-XpD8TFaxhdP41X8-BrUc0j9gENGrIBnYhPtJiE/edit?usp=sharing
    This function asks the question - For the combination of haplotypes present in the reference trio (reference,
    unaffected_partner and affected_partner) what information would this SNP provide us with if the SNP was:

    AB in a female embryo
    AA in a male embryo
    BB in a male embryo

    i.e. does this SNP indicate "high_risk", "low_risk", or is it "uninformative" in regard to inheriting the X-linked condition?


    Args:
        df (dataframe): A dataframe containing SNP input data for both partners, reference, & embryos
        carrier_female_partner (string): Column name in dataframe refering to carrier_female_partner's data
        unaffected_male_partner (string): Column name in dataframe refering to unaffected_male_partner's data
        reference (string) : Column name in dataframe refering to reference's data (Always child)"

    Returns:
        Dataframe: Dataframe containing 3 new "snp_risk_category" columns, used to categorise the SNPs as
        "high_risk", "low_risk", and "uninformative" for the three different embryo catergories - female_AB_snp_risk_category,
        male_AA_snp_risk_category, male_BB_snp_risk_category
    """

    # Calculate & classify the informative for female AB embryos
    conditions = [
        # Criteria_XL1
        (df[reference] == "AA")
        & (df[unaffected_male_partner] == "AA")
        & (df[carrier_female_partner] == "AB"),
        # Criteria_XL2
        (df[reference] == "AA")
        & (df[unaffected_male_partner] == "BB")
        & (df[carrier_female_partner] == "AB"),
        # Criteria_XL3
        (df[reference] == "BB")
        & (df[unaffected_male_partner] == "AA")
        & (df[carrier_female_partner] == "AB"),
        # Criteria_XL4
        (df[reference] == "BB")
        & (df[unaffected_male_partner] == "BB")
        & (df[carrier_female_partner] == "AB"),
    ]
    # Assign the correct labels

    values = [
        "low_risk",  # Criteria_XL1
        "high_risk",  # Criteria_XL2
        "high_risk",  # Criteria_XL3
        "low_risk",  # Criteria_XL4
    ]

    df["female_AB_snp_risk_category"] = np.select(
        conditions, values, default="uninformative"
    )

    # Calculate & classify the informative for male embryos
    conditions = [
        # Criteria_XL5, Criteria_XL7
        (df[reference] == "AA")
        & (df[unaffected_male_partner] == "BB")
        & (df[carrier_female_partner] == "AB"),
        # Criteria_XL6, Criteria_XL8
        (df[reference] == "BB")
        & (df[unaffected_male_partner] == "AA")
        & (df[carrier_female_partner] == "AB"),
    ]
    # Assign the correct labels for male AA

    values = [
        "high_risk",  # Criteria_XL5
        "low_risk",  # Criteria_XL6
    ]

    df["male_AA_snp_risk_category"] = np.select(
        conditions, values, default="uninformative"
    )

    values = [
        "low_risk",  # Criteria_XL7
        "high_risk",  # Criteria_XL8
    ]

    # Assign the correct labels for male BB

    df["male_BB_snp_risk_category"] = np.select(
        conditions, values, default="uninformative"
    )

    return df
