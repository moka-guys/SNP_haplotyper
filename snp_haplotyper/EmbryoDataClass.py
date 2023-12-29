import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from FamilyDataClass import FamilyData
from pydantic import BaseModel, root_validator
from snp_plot import plot_results
from SNPDataClass import SNPData


class EmbryoAlleleCategorizer(ABC):
    def __init__(self, embryo_id, embryo_sex, family_genetic_data, snp_data_df):
        self.embryo_id = embryo_id
        self.embryo_sex = embryo_sex
        self.family_genetic_data = family_genetic_data
        self.embryo_category_df = snp_data_df[
            [
                "probeset_id",
                "rsID",
                "Position",
                "gene_distance",
                "snp_position",
                "snp_risk_category_AB",
                "snp_inherited_from",
                "snp_risk_category_AA",
                "snp_risk_category_BB",
                self.family_genetic_data.male_partner,
                self.family_genetic_data.female_partner,
                self.family_genetic_data.reference,
            ]
            + [self.embryo_id]
        ].copy()

    @abstractmethod
    def categorize(self, *args, **kwargs):
        pass


class AutosomalDominantCategorizer(EmbryoAlleleCategorizer):
    def categorize(self) -> pd.DataFrame:
        """
        For Autosomal Dominant mode of inheritance, calculate the embryo risk category for each SNP.
        """
        embryo_risk_col = f"embryo_risk_category"

        conditions = [
            (self.embryo_category_df["snp_risk_category_AB"] == "high_risk")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            (self.embryo_category_df["snp_risk_category_AB"] == "low_risk")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            (self.embryo_category_df["snp_risk_category_AB"] != "uninformative")
            & (self.embryo_category_df[self.embryo_id] == "NoCall"),
        ]
        values = ["high_risk", "low_risk", "NoCall"]

        self.embryo_category_df[embryo_risk_col] = np.select(
            conditions, values, default="uninformative"
        )
        self.embryo_category_df[embryo_risk_col] = pd.Categorical(
            self.embryo_category_df[embryo_risk_col],
            categories=[
                "high_risk",
                "low_risk",
                "uninformative",
                "ADO",
                "miscall",
                "NoCall",
                "NoCall_in_trio",
            ],
            ordered=True,
        )

        return self.embryo_category_df


class AutosomalRecessiveCategorizer(EmbryoAlleleCategorizer):
    def categorize(self):
        """
        For Autosomal Recessive mode of inheritance, calculate the embryo risk category for each SNP.
        """
        embryo_risk_col = f"embryo_risk_category"
        consanguineous = self.family_genetic_data.consanguineous

        conditions = [
            # male_partner
            (self.embryo_category_df["snp_risk_category_AB"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "male_partner")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            (self.embryo_category_df["snp_risk_category_AB"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "male_partner")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            # female_partner
            (self.embryo_category_df["snp_risk_category_AB"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "female_partner")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            (self.embryo_category_df["snp_risk_category_AB"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "female_partner")
            & (self.embryo_category_df[self.embryo_id] == "AB"),
            # both_partners
            (self.embryo_category_df["snp_risk_category_AA"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "AA"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_BB"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "BB"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_AA"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "AA"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_BB"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "BB"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_AA"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "AA"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_BB"] == "high_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "BB"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_AA"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "AA"))
            & consanguineous,
            (self.embryo_category_df["snp_risk_category_BB"] == "low_risk")
            & (self.embryo_category_df["snp_inherited_from"] == "both_partners")
            & ((self.embryo_category_df[self.embryo_id] == "BB"))
            & consanguineous,
            # NoCall
            (self.embryo_category_df["snp_risk_category_AB"] != "uninformative")
            & (self.embryo_category_df["snp_risk_category_AB"] == "NoCall"),
        ]
        values = [
            "high_risk",
            "low_risk",
            "high_risk",
            "low_risk",
            "high_risk",
            "high_risk",
            "low_risk",
            "low_risk",
            "high_risk",
            "high_risk",
            "low_risk",
            "low_risk",
            "NoCall",
        ]
        self.embryo_category_df[embryo_risk_col] = np.select(
            conditions, values, default="uninformative"
        )
        self.embryo_category_df[embryo_risk_col] = pd.Categorical(
            self.embryo_category_df[embryo_risk_col],
            categories=[
                "high_risk",
                "low_risk",
                "uninformative",
                "ADO",
                "miscall",
                "NoCall",
                "NoCall_in_trio",
            ],
            ordered=True,
        )
        return self.embryo_category_df


class XLinkedCategorizer(EmbryoAlleleCategorizer):
    def categorize(self) -> pd.DataFrame:
        """
        For X-linked mode of inheritance, calculate the embryo risk category for each SNP.
        """
        embryo_risk_col = f"embryo_risk_category"

        if self.embryo_sex == Sex.UNKNOWN:
            raise ArgumentInputError(
                f"'unknown' embryo sex not allowed for x-linked mode of inheritance. Check that correct mode of inheritance has been entered for {self.embryo_id}, or enter correct sex for {self.embryo_id}"
            )

        if self.embryo_sex == Sex.FEMALE:
            conditions = [
                (self.embryo_category_df["snp_risk_category_AB"] == "high_risk")
                & (self.embryo_category_df[self.embryo_id] == "AB"),
                (self.embryo_category_df["snp_risk_category_AB"] == "low_risk")
                & (self.embryo_category_df[self.embryo_id] == "AB"),
                (self.embryo_category_df["snp_risk_category_AB"] != "uninformative")
                & (self.embryo_category_df[self.embryo_id] == "NoCall"),
                (self.embryo_category_df["snp_risk_category_AB"] == "NoCall_in_trio"),
            ]
            values = [
                "high_risk",
                "low_risk",
                "NoCall",
                "NoCall_in_trio",
            ]

            self.embryo_category_df[embryo_risk_col] = np.select(
                conditions, values, default="uninformative"
            )

        elif self.embryo_sex == Sex.MALE:
            conditions = [
                (self.embryo_category_df["snp_risk_category_AA"] == "high_risk")
                & (self.embryo_category_df[self.embryo_id] == "AA"),
                (self.embryo_category_df["snp_risk_category_AA"] == "low_risk")
                & (self.embryo_category_df[self.embryo_id] == "AA"),
                (self.embryo_category_df["snp_risk_category_BB"] == "high_risk")
                & (self.embryo_category_df[self.embryo_id] == "BB"),
                (self.embryo_category_df["snp_risk_category_BB"] == "low_risk")
                & (self.embryo_category_df[self.embryo_id] == "BB"),
                (self.embryo_category_df["snp_risk_category_AA"] != "uninformative")
                & (self.embryo_category_df["snp_risk_category_BB"] != "uninformative")
                & (self.embryo_category_df[self.embryo_id] == "NoCall"),
                (self.embryo_category_df["snp_risk_category_AA"] == "NoCall_in_trio"),
                (self.embryo_category_df["snp_risk_category_BB"] == "NoCall_in_trio"),
            ]
            values = [
                "high_risk",
                "low_risk",
                "high_risk",
                "low_risk",
                "NoCall",
                "NoCall_in_trio",
                "NoCall_in_trio",
            ]

            self.embryo_category_df[embryo_risk_col] = np.select(
                conditions, values, default="uninformative"
            )

        else:
            raise ArgumentInputError(
                f"Invalid embryo sex for {self.embryo_id} in x-linked categorization. Check your input data."
            )
        self.embryo_category_df[embryo_risk_col] = pd.Categorical(
            self.embryo_category_df[embryo_risk_col],
            categories=[
                "high_risk",
                "low_risk",
                "uninformative",
                "ADO",
                "miscall",
                "NoCall",
                "NoCall_in_trio",
            ],
            ordered=True,
        )
        return self.embryo_category_df


def get_allele_categorizer(mode_of_inheritance: InheritanceMode):
    """
    Select the correct categorizer based on the mode of inheritance.
    """
    categorizers = {
        InheritanceMode.AUTOSOMAL_DOMINANT: AutosomalDominantCategorizer,
        InheritanceMode.AUTOSOMAL_RECESSIVE: AutosomalRecessiveCategorizer,
        InheritanceMode.X_LINKED: XLinkedCategorizer,
    }
    return categorizers.get(mode_of_inheritance, None)


class EmbryoRiskSummariser(ABC):
    """
    Abstract class for summarising embryo risk data into summary tables grouping the number of SNPs in each category.
    """

    def __init__(self, embryo_id, embryo_category_df):
        self.embryo_id = embryo_id
        self.embryo_category_df = embryo_category_df.copy()


class AutosomalDominantRiskSummariser(EmbryoRiskSummariser):
    """
    Class for Autosomal Dominant cases, summarising embryo risk data into summary tables grouping the number of SNPs in each category.
    """

    def risk_summary(self, embryo_id, embryo_category_df: pd.DataFrame) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category.
        """
        risk_summary_df = (
            embryo_category_df.groupby("embryo_risk_category")
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_df

    def risk_summary_per_region(
        self, embryo_id, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category and region.
        """
        risk_summary_per_region = (
            embryo_category_df.groupby(["embryo_risk_category", "gene_distance"])
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_per_region

    def risk_summary_for_testing(
        self, embryo_id: str, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data in format expected by pytests.

        Parameters:
        - embryo_id: ID of the embryo
        - embryo_category_df: DataFrame containing embryo risk categories and SNP positions

        Returns:
        - pd.DataFrame: Filtered DataFrame containing only "high_risk" and "low_risk" categories
        """

        # Group by embryo_risk_category and snp_position and count occurrences
        grouped_df = (
            embryo_category_df.groupby(["embryo_risk_category", "snp_position"])
            .size()
            .reset_index(name=embryo_id)
        )

        # Filter rows where embryo_risk_category is either 'high_risk' or 'low_risk'
        filtered_df = grouped_df[
            grouped_df["embryo_risk_category"].isin(["high_risk", "low_risk"])
        ]

        return filtered_df


class AutosomalRecessiveRiskSummariser(EmbryoRiskSummariser):
    """
    Class for Autosomal Recessive cases, summarising embryo risk data into summary tables grouping the number of SNPs in each category.
    """

    def risk_summary(self, embryo_id, embryo_category_df: pd.DataFrame) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category.
        """
        risk_summary_df = (
            embryo_category_df.groupby(
                [
                    "snp_inherited_from",
                    "embryo_risk_category",
                ]
            )
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_df

    def risk_summary_per_region(
        self, embryo_id, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category and region.
        """
        risk_summary_per_region = (
            embryo_category_df.groupby(
                [
                    "snp_inherited_from",
                    "gene_distance",
                    "embryo_risk_category",
                ]
            )
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_per_region

    def risk_summary_for_testing(
        self, embryo_id: str, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data in format expected by pytests.

        Parameters:
        - embryo_id: ID of the embryo
        - embryo_category_df: DataFrame containing embryo risk categories and SNP positions

        Returns:
        - pd.DataFrame: Filtered DataFrame containing only "high_risk" and "low_risk" categories
        """

        # Group by embryo_risk_category and snp_position and count occurrences
        grouped_df = (
            embryo_category_df.groupby(
                [
                    "snp_inherited_from",
                    "snp_position",
                    "embryo_risk_category",
                ]
            )
            .size()
            .reset_index(name=embryo_id)
        )

        # Filter rows where embryo_risk_category is either 'high_risk' or 'low_risk'
        filtered_df = grouped_df[
            grouped_df["embryo_risk_category"].isin(["high_risk", "low_risk"])
        ]

        filtered_df = filtered_df[
            grouped_df["snp_inherited_from"].isin(
                ["male_partner", "female_partner", "both_partners"]
            )
        ]

        return filtered_df


class XLinkedRiskSummariser(EmbryoRiskSummariser):
    """
    Class for X-linked cases, summarising embryo risk data into summary tables grouping the number of SNPs in each category.
    """

    def risk_summary(self, embryo_id, embryo_category_df: pd.DataFrame) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category.
        """
        risk_summary_df = (
            embryo_category_df.groupby("embryo_risk_category")
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_df

    def risk_summary_per_region(
        self, embryo_id, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data by SNP risk category and region.
        """
        risk_summary_per_region = (
            embryo_category_df.groupby(["embryo_risk_category", "gene_distance"])
            .size()
            .reset_index(name=embryo_id)
        )
        return risk_summary_per_region

    def risk_summary_for_testing(
        self, embryo_id: str, embryo_category_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarises embryo risk data in format expected by pytests.

        Parameters:
        - embryo_id: ID of the embryo
        - embryo_category_df: DataFrame containing embryo risk categories and SNP positions

        Returns:
        - pd.DataFrame: Filtered DataFrame containing only "high_risk" and "low_risk" categories
        """

        # Group by embryo_risk_category and snp_position and count occurrences
        grouped_df = (
            embryo_category_df.groupby(["embryo_risk_category", "snp_position"])
            .size()
            .reset_index(name=embryo_id)
        )

        # Filter rows where embryo_risk_category is either 'high_risk' or 'low_risk'
        filtered_df = grouped_df[
            grouped_df["embryo_risk_category"].isin(["high_risk", "low_risk"])
        ]

        return filtered_df


def get_embryo_risk_summariser(mode_of_inheritance):
    """
    Select the correct summariser based on the mode of inheritance.
    """
    summarisers = {
        InheritanceMode.AUTOSOMAL_DOMINANT: AutosomalDominantRiskSummariser,
        InheritanceMode.AUTOSOMAL_RECESSIVE: AutosomalRecessiveRiskSummariser,
        InheritanceMode.X_LINKED: XLinkedRiskSummariser,
    }
    return summarisers.get(mode_of_inheritance, None)


def detect_miscall_or_ado(
    male_partner_haplotype: str,
    female_partner_haplotype: str,
    reference_haplotype: str,
    embryo_haplotype: str,
    mode_of_inheritance: InheritanceMode,
):
    """QC identify miscalls or ADOs (Allele Drop Outs)
    Takes the haplotypes for the male partner, female partner and embryo and calculates whether
    it indicates a miscall or ADO (Allele dropout) in the embryo for that SNP.
    The definition of a miscall is any haplotype in the embryo which is inconsistent
    with the haplotype of the parents i.e. Parents "AA", "BB" and an embryo "AA". This is due
    to technical error in the measurement.  NOTE: that the miscall could be in any one of the
    trio even though it is recorded under the embryo.
    The definition of ADO (Allele dropout) is used when there is a suspected biological origin for the
    mismatch in haplotypes, due to uniparental inheritance of the allele i.e Parents AA, BB and an embryo AA,
    the B allele has dropped out.  NOTE: that the ADO could have occured in any of the trio even though it is
    recorded under the embryo.  It is expected that the Genomic Scientist will look at the SNP plots and
    use their judgement as to whether allele dropout is observed.
    Args:
        male_partner_haplotype (string): Either "AA", "BB", "AB", or "NoCall"
        female_partner_haplotype (string): Either "AA", "BB", "AB", or "NoCall"
        embryo_haplotype (string): Either "AA", "BB", "AB", or "NoCall"
    Returns:
        string: 'call', 'miscall', 'ADO', or an Error message
    """

    # Validate input data or throw exception
    illegal_args = set(
        [
            male_partner_haplotype,
            female_partner_haplotype,
            reference_haplotype,
            embryo_haplotype,
        ]
    ) - set(["AA", "BB", "AB", "NoCall"])
    if len(illegal_args) != 0:
        raise ArgumentInputError(
            f"Function detect_miscall_or_ado() only excepts 'AA','BB', 'AB', NoCall' as arguments, recieved {str(illegal_args)}"
        )

    alleles = [
        male_partner_haplotype,
        female_partner_haplotype,
        reference_haplotype,
        embryo_haplotype,
        mode_of_inheritance,
    ]

    # Check if mode of inheritance is null
    if mode_of_inheritance == None:
        raise InvalidParameterSelectedError(
            f"Invalid mode of inheritance: {mode_of_inheritance}"
        )

    match alleles:
        # For autosomal dominant and autosomal recessive mark as NoCall if any column NoCall
        case [
            "NoCall",
            _,
            _,
            _,
            InheritanceMode.AUTOSOMAL_DOMINANT | InheritanceMode.AUTOSOMAL_RECESSIVE,
        ] | [
            _,
            "NoCall",
            _,
            _,
            InheritanceMode.AUTOSOMAL_DOMINANT | InheritanceMode.AUTOSOMAL_RECESSIVE,
        ] | [
            _,
            _,
            "NoCall",
            _,
            InheritanceMode.AUTOSOMAL_DOMINANT | InheritanceMode.AUTOSOMAL_RECESSIVE,
        ] | [
            _,
            _,
            _,
            "NoCall",
            InheritanceMode.AUTOSOMAL_DOMINANT | InheritanceMode.AUTOSOMAL_RECESSIVE,
        ]:
            result = "NoCall"
        # For x-linked mark as NoCall if any column other than male_partner is NoCall
        case [_, "NoCall", _, _, InheritanceMode.X_LINKED] | [
            _,
            _,
            "NoCall",
            _,
            InheritanceMode.X_LINKED,
        ] | [_, _, _, "NoCall", InheritanceMode.X_LINKED]:
            result = "NoCall"
        # Compare embryo haplotype to parent haplotypes
        case ["AA", "AA", _, _, _]:
            result = "miscall" if embryo_haplotype != "AA" else "uninformative"
        case ["BB", "BB", _, _, _]:
            result = "miscall" if embryo_haplotype != "BB" else "uninformative"
        case ["AA", "BB", _, _, _]:
            result = "ADO" if embryo_haplotype != "AB" else "uninformative"
        case ["BB", "AA", _, _, _]:
            result = "ADO" if embryo_haplotype != "AB" else "uninformative"
        case ["AA", "AB", _, _, _]:
            result = (
                "ADO"
                if embryo_haplotype
                not in [
                    "AA",
                    "AB",
                ]
                else "uninformative"
            )
        case ["AB", "AA", _, _, _]:
            result = (
                "ADO"
                if embryo_haplotype
                not in [
                    "AA",
                    "AB",
                ]
                else "uninformative"
            )
        case ["BB", "AB", _, _, _]:
            result = (
                "ADO"
                if embryo_haplotype
                not in [
                    "BB",
                    "AB",
                ]
                else "uninformative"
            )
        case ["AB", "BB", _, _, _]:
            result = (
                "ADO"
                if embryo_haplotype
                not in [
                    "BB",
                    "AB",
                ]
                else "uninformative"
            )
        case ["AB", "AB", _, _, _]:
            result = (
                "uninformative"
                if embryo_haplotype
                not in [
                    "AA",
                    "BB",
                    "AB",
                ]
                else "uninformative"
            )
        # Match anything not caught by the above cases
        case _:
            result = "Error: Unexpected combination of parent alleles."

    return result


def update_embryo_risk_column(
    df: pd.DataFrame,
    risk_col_name: str,
    male_partner_col: str,
    female_partner_col: str,
    reference_col: str,
    embryo_col: str,
    mode_of_inheritance: InheritanceMode,
):
    """
    Update the risk column in the DataFrame with miscalls and ADOs information.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing embryo data
    - risk_col_name (str): The column name to update
    - male_partner_col (str): The column name for male partner data
    - female_partner_col (str): The column name for female partner data
    - embryo_col (str): The column name for embryo data

    Returns:
    - pd.DataFrame: Updated DataFrame with modified risk column
    """
    df[risk_col_name] = df.apply(
        lambda row: row[risk_col_name]
        if row[risk_col_name] != "uninformative"
        else detect_miscall_or_ado(
            row[male_partner_col],
            row[female_partner_col],
            row[reference_col],
            row[embryo_col],
            mode_of_inheritance,
        ),
        axis=1,
    )
    return df


@dataclass
class EmbryoData:
    embryo_id: str
    embryo_sex: Sex
    family_genetic_data: FamilyData
    embryo_category_df: pd.DataFrame

    def __post_init__(self):
        self.mode_of_inheritance: InheritanceMode = (
            self.family_genetic_data.mode_of_inheritance
        )
        self.male_partner: str = self.family_genetic_data.male_partner
        self.female_partner: str = self.family_genetic_data.female_partner
        self.consanguineous: bool = self.family_genetic_data.consanguineous
        self.reference: str = self.family_genetic_data.reference

        self.rsid: "pd.Series" = self.embryo_category_df["rsID"]

        CategorizerClass = get_allele_categorizer(self.mode_of_inheritance)

        if CategorizerClass:
            categorizer_instance = CategorizerClass(
                self.embryo_id,
                self.embryo_sex,
                self.family_genetic_data,
                self.embryo_category_df,
            )
            # Store the result from the categorize method
            self.embryo_category_df = categorizer_instance.categorize()
        else:
            raise ValueError(f"Invalid mode of inheritance: {self.mode_of_inheritance}")
        self.embryo_category_df = update_embryo_risk_column(
            self.embryo_category_df,
            "embryo_risk_category",
            self.male_partner,
            self.female_partner,
            self.reference,
            self.embryo_id,
            self.mode_of_inheritance,
        )

        # Define the categorical type with specific categories
        cat_type = pd.CategoricalDtype(
            categories=[
                "high_risk",
                "low_risk",
                "uninformative",
                "ADO",
                "miscall",
                "NoCall",
                "NoCall_in_trio",
            ],
            ordered=True,
        )

        self.embryo_category_df["embryo_risk_category"] = self.embryo_category_df[
            "embryo_risk_category"
        ].astype(cat_type)

        SummarizerClass = get_embryo_risk_summariser(self.mode_of_inheritance)

        summariser_instance = SummarizerClass(
            self.embryo_id,
            self.embryo_category_df,
        )

        self.risk_summary_df = summariser_instance.risk_summary(
            self.embryo_id, self.embryo_category_df
        )
        self.risk_summary_per_region_df = summariser_instance.risk_summary_per_region(
            self.embryo_id, self.embryo_category_df
        )
        self.risk_summary_for_testing_df = summariser_instance.risk_summary_for_testing(
            self.embryo_id, self.embryo_category_df
        )
        self.results_plot = plot_results(
            self.embryo_category_df,
            self.embryo_id,
            self.embryo_sex,
            self.family_genetic_data.gene_start,
            self.family_genetic_data.gene_end,
            self.mode_of_inheritance,
            self.risk_summary_for_testing_df,
            self.family_genetic_data.flanking_region_size,
            self.consanguineous,
        )

    def get_embryo_results(self):
        return self.embryo_category_df

    def get_pytest_benchmark(self):
        return self.risk_summary_for_testing_df
