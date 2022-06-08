import numpy as np
import pandas as pd

# Logic is the same whether using affected son of carrier grandmother as reference
def x_linked_analysis(
    df,
    affected_partner,
    unaffected_partner,
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
            # Criteria 1
            (df[reference] == "AA")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria 2
            (df[reference] == "BB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria 3
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria 4
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria 1
                "high_risk",  # Criteria 2
                "low_risk",  # Criteria 3
                "low_risk",  # Criteria 4
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria 1
                "low_risk",  # Criteria 2
                "high_risk",  # Criteria 3
                "high_risk",  # Criteria 4
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
    elif reference_relationship in [
        "child",  # TODO populate with all appropriate relationships
    ]:
        # Label alleles as high or low risk
        conditions = [
            # Criteria to label high Risk SNPs if reference affected, or low risk SNPs if reference unaffected
            # Criteria 1
            (df[reference] == "AB")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria 2
            (df[reference] == "AB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
            # Criteria to label low Risk SNPs if reference affected, or high risk SNPs if reference unaffected
            # Criteria 3
            (df[reference] == "AA")
            & (df[unaffected_partner] == "AA")
            & (df[affected_partner] == "AB"),
            # Criteria 4
            (df[reference] == "BB")
            & (df[unaffected_partner] == "BB")
            & (df[affected_partner] == "AB"),
        ]
        # Assign the correct labels depending upon reference status
        if reference_status == "affected":
            values = [
                "high_risk",  # Criteria 1
                "high_risk",  # Criteria 2
                "low_risk",  # Criteria 3
                "low_risk",  # Criteria 4
            ]
        elif reference_status == "unaffected":
            values = [
                "low_risk",  # Criteria 1
                "low_risk",  # Criteria 2
                "high_risk",  # Criteria 3
                "high_risk",  # Criteria 4
            ]
        df["snp_risk_category"] = np.select(conditions, values, default="uninformative")
    return df
