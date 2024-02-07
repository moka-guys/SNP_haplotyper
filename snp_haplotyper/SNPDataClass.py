import logging
import os
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum, auto
from io import IOBase
from pathlib import Path
from typing import Any, Dict, List, Optional

import config as config  # TODO add code to use this dependency
import numpy as np
import pandas as pd
import pandera as pa
import pdfkit
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from helper_functions import custom_order_generator
from inheritance_logic import (
    AutosomalDominantLogic,
    AutosomalRecessiveLogic,
    XLinkedLogic,
)
from jinja2 import Environment, PackageLoader
from pandas import DataFrame
from pandas.api.types import CategoricalDtype
from pydantic import (
    BaseModel,
    ValidationError,
    conlist,
    constr,
    fields,
    root_validator,
    validator,
)
from pydantic.dataclasses import dataclass

logger = logging.getLogger("BASHer_logger")

# Add the directory containing this script to the PYTHONPATH
sys.path.append(os.path.dirname(__file__))
mod_path = Path(__file__).parent

from dataclasses import dataclass

import config as config
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


def single_value(s: pd.Series) -> bool:
    # checks if a series has a unique value
    return s.nunique() == 1


def no_duplicates(s: pd.Series) -> bool:
    # checks if a series has duplicate values
    return not s.duplicated().any()


snp_df_schema = pa.DataFrameSchema(
    {
        "probeset_id": pa.Column(
            str,
            checks=[
                pa.Check.str_startswith("AX-"),
                pa.Check(
                    no_duplicates,
                    error="Duplicate values found in 'probeset_id' column",
                ),
            ],
            required=True,
        ),
        "Chr": pa.Column(
            pa.Category,
            checks=[
                pa.Check.isin(
                    [
                        "X",
                        "Y",
                        "1",
                        "2",
                        "3",
                        "4",
                        "5",
                        "6",
                        "7",
                        "8",
                        "9",
                        "10",
                        "11",
                        "12",
                        "13",
                        "14",
                        "15",
                        "16",
                        "17",
                        "18",
                        "19",
                        "20",
                        "21",
                        "22",
                    ]
                ),
                pa.Check(
                    single_value,
                    error="More than one chromosome specified in 'Chr' column",
                ),
            ],
            required=True,
        ),
        "Position": pa.Column(
            int,
            checks=[
                pa.Check(
                    lambda x: x >= 0, error="Negative values in 'Position' column"
                ),
            ],
            required=True,
        ),
        # Check all columns with SNP calls matching the regex
        ".*rhchp$": pa.Column(
            pa.Category,
            checks=pa.Check.isin(["AA", "BB", "AB", "NoCall"]),
            required=True,
            regex=True,
        ),
    },
)


@dataclass
class SNPData:
    family_data: FamilyData
    snp_df: pd.DataFrame
    affy_2_rs_ids_df: Optional[pd.DataFrame] = None
    """
    SNPData dataclass
    =================

    This class represents and manipulates Single Nucleotide Polymorphism (SNP) data imported from an Affymetrix array. It further annotates this P data
    based on genomic coordinates and other associated data.

    """

    def __post_init__(self):
        self.number_snps_imported = self.snp_df.shape[0]
        self.validate_snp_df()
        self.bp_flanking_region_size: int = (
            self.family_data.flanking_region_size * 10**6
        )
        if self.affy_2_rs_ids_df is None:
            self.affy_2_rs_ids_df = self._load_affy_2_rs_ids()
        self._add_rsid_column()
        self._add_duplicated_probeset_id_column()
        self._annotate_snp_data(self.family_data.flanking_region_size)
        self.inheritance_logic = InheritanceLogicFactory.create_logic(
            self.family_data.mode_of_inheritance
        )
        self.remove_snps_outside_roi()
        self.categorise_snp_risk_category()
        self.add_risk_summary_column()
        self.summary_logic = SummaryFactory.create_summary(
            self.family_data.mode_of_inheritance, self.family_data.flanking_region_size
        )
        self.informative_snps_summary = self.group_informative_snps_by_region()
        self.summarise_test_data()

    @staticmethod
    def _load_affy_2_rs_ids():
        """
        Loads the AffyID2rsid data containing the mapping for the probes_set IDs to dbSNP rsIDs.

        Returns:
            pd.DataFrame: A DataFrame containing the mapping for Affy probes_set IDs to dbSNP rsIDs.
        """
        # TODO add this to config
        mod_path = Path(__file__).parent
        rsid_data_path = (
            mod_path / "../test_data/AffyID2rsid.txt"
        ).resolve()  #  TODO add this to config
        if not rsid_data_path.exists():
            raise FileNotFoundError(f"File not found at {rsid_data_path}")
        df = pd.read_csv(rsid_data_path, delimiter=",", low_memory=False)
        if df is None or df.empty:
            raise ValueError(
                f"Failed to load the required DataFrame from the provided path."
            )
        return df

    def _add_rsid_column(self) -> None:
        """Provides dbsnp rsIDs
        New column created in the dataframe, self.snp_df, matching the probes_set IDs to dbSNP rsIDs.
        """
        if self.affy_2_rs_ids_df is None:
            raise ValueError("affy_2_rs_ids_df is None. Expected a DataFrame.")
        self.snp_df = pd.merge(
            self.snp_df,
            self.affy_2_rs_ids_df[["probeset_id", "affy_snp_id", "rsID"]],
            on="probeset_id",
            how="left",
        )

        # Rearrange columns so that rsID is next to Affy Id
        self.snp_df.insert(0, "probeset_id", self.snp_df.pop("probeset_id"))
        self.snp_df.insert(1, "rsID", self.snp_df.pop("rsID"))

    def _add_duplicated_probeset_id_column(self):
        """
        Adds a boolean column to the dataframe indicating if 'probeset_id' is duplicated.

        Parameters:
        snp_df (DataFrame): The dataframe to which the new column will be added.

        Returns:
        DataFrame: The dataframe with the new boolean column.
        """
        if not isinstance(self.snp_df, pd.DataFrame):
            raise ValueError("Expected a DataFrame as input.")

        # Checking for duplicate 'probeset_id' and adding a new column
        self.snp_df["probe_is_duplicated"] = self.snp_df.duplicated(
            "probeset_id", keep=False
        )

    def validate_snp_df(self) -> None:
        """
        Validates the SNP dataframe against a predefined schema.
        """
        snp_df_schema.validate(self.snp_df, lazy=True)

    @staticmethod
    def categorize_SNP_by_position(
        position_series: pd.Series, gene_start: int, gene_end: int, max_region_mb: int
    ) -> pd.Series:
        """
        Categorizes SNPs as 'upstream', 'within_gene', 'downstream', or 'outside_ROI' based on positions.

        Args:
            position_series (pd.Series): Series containing SNP positions.
            gene_start (int): Starting position of the SNP in base pairs.
            gene_end (int): Ending position of the SNP in base pairs.
            max_region_mb (int): Maximum region distance in megabases from the gene.

        Returns:
            pd.Series: Categorized series.
        """

        # Convert max_region_mb to base pairs and compute region_start and region_end
        max_region_bp = max_region_mb * 10**6
        region_start = gene_start - max_region_bp
        region_end = gene_end + max_region_bp

        # Define conditions based on the regions
        conditions = [
            (position_series >= region_start) & (position_series < gene_start),
            (position_series >= gene_start) & (position_series <= gene_end),
            (position_series > gene_end) & (position_series <= region_end),
        ]

        # Corresponding labels for the regions
        labels = ["upstream", "within_gene", "downstream"]

        # Use numpy select to categorize based on conditions and labels
        categorized_series = pd.Series(
            pd.Categorical(
                np.select(conditions, labels, default="outside_ROI"),
                categories=["upstream", "within_gene", "downstream", "outside_ROI"],
                ordered=True,
            )
        )

        return categorized_series

    @staticmethod
    def calculate_mb_distance(
        gene_positions: pd.Series, gene_start: int, gene_end: int
    ) -> pd.Series:
        # Convert gene_positions to numpy array for vectorized operations
        positions = gene_positions.to_numpy()

        # Calculate distances
        distances_to_start = (positions - gene_start) / 1_000_000
        distances_to_end = (positions - gene_end) / 1_000_000

        # Use np.where to decide which distance to use
        # If within gene range, set to 0
        distances = np.where(
            (positions >= gene_start) & (positions <= gene_end),
            0,
            np.where(positions < gene_start, distances_to_start, distances_to_end),
        )

        # Round away from zero using np.copysign and np.ceil
        rounded_distances = np.copysign(np.ceil(np.abs(distances)), distances).astype(
            int
        )

        return pd.Series(rounded_distances, index=gene_positions.index)

    @staticmethod
    def annotate_distances(distances: pd.Series) -> pd.Series:
        max_range_mb = max(abs(distances.max()), abs(distances.min()))

        def annotate(value):
            if value < 0:
                return f"{abs(value) - 1}-{abs(value)}MB_from_start"
            elif value == 0:
                return "within_gene"
            else:
                return f"{value-1}-{value}MB_from_end"

        category_order = custom_order_generator(max_range_mb)
        return pd.Series(
            pd.Categorical(
                distances.apply(annotate), categories=category_order, ordered=True
            )
        )

    def _annotate_snp_data(self, distance_mb: int):
        """
        Annotates the SNP based on the provided genomic coordinates by adding two new columns.
        "gene_distance" column is used to group the SNPs based on their distance from the gene i.e.
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
        "snp_position" column to group the gene distance as "upstream", "within_gene", or "downstream".

        Args:
            distance_mb (int): The distance in megabases to categorize SNPs relative to genes.
        """

        # Ensure gene_start and gene_end are integers
        start = int(self.family_data.gene_start)
        end = int(self.family_data.gene_end)

        # Use static method to categorize SNP by position
        self.snp_df["snp_position"] = self.categorize_SNP_by_position(
            self.snp_df["Position"], start, end, distance_mb
        )

        # Use static method to calculate distances in MB
        distances_in_mb = self.calculate_mb_distance(
            self.snp_df["Position"], start, end
        )
        self.snp_df["gene_distance"] = self.annotate_distances(distances_in_mb)

        # Convert gene_distance to Categorical type with custom order
        self.snp_df["gene_distance"] = pd.Categorical(
            self.snp_df["gene_distance"],
            categories=custom_order_generator(distance_mb),
            ordered=True,
        )

        # Convert snp_position to Categorical type with predefined order
        categories_list = [
            "upstream",
            "within_gene",
            "downstream",
        ]
        self.snp_df["snp_position"] = pd.Categorical(
            self.snp_df["snp_position"], categories=categories_list, ordered=True
        )

    def remove_snps_outside_roi(self):
        """
        Filters out SNPs outside the region of interest.
        """
        # Check if snp_df attribute exists and is a DataFrame
        if hasattr(self, "snp_df") and isinstance(self.snp_df, pd.DataFrame):
            # Check if 'snp_position' column exists in snp_df
            if "snp_position" in self.snp_df.columns:
                categories_list = [
                    "upstream",
                    "within_gene",
                    "downstream",
                    "outside_ROI",
                ]
                self.snp_df["snp_position"] = pd.Categorical(
                    self.snp_df["snp_position"],
                    categories=categories_list,
                    ordered=True,
                )
                # Filter out SNPs outside ROI and reset index
                self.outside_roi = self.snp_df[
                    self.snp_df["snp_position"] == "outside_ROI"
                ].reset_index(drop=True)
                # Filter in SNPs inside ROI and reset index
                self.snp_df = self.snp_df[
                    self.snp_df["snp_position"] != "outside_ROI"
                ].reset_index(drop=True)
            else:
                print("Error: 'snp_position' column not found in snp_df.")
        else:
            print("Error: snp_df attribute is missing or not a DataFrame.")

    def categorise_snp_risk_category(self):
        if self.inheritance_logic:
            InheritanceSpecificAnalysis = self.inheritance_logic(
                self.snp_df,
                self.family_data.get_partner1(),
                self.family_data.get_partner2(),
                reference=self.family_data.reference,
                reference_status=self.family_data.reference_status,
                reference_relationship=self.family_data.reference_relationship,
                reference_sex=self.family_data.reference_sex,
                consanguineous=self.family_data.consanguineous,
            )
            self.snp_df = InheritanceSpecificAnalysis.df
        else:
            raise ValueError(
                "Invalid inheritance logic or mode of inheritance not set!"
            )
        categories_list = [
            "upstream",
            "within_gene",
            "downstream",
        ]
        self.snp_df["snp_position"] = pd.Categorical(
            self.snp_df["snp_position"], categories=categories_list, ordered=True
        )

    def add_risk_summary_column(self):
        def categorize(row):
            values = [
                row["snp_risk_category_AB"],
                row["snp_risk_category_AA"],
                row["snp_risk_category_BB"],
            ]

            if values == ["uninformative", "uninformative", "uninformative"]:
                return "uninformative"
            elif "high_risk" in values and "low_risk" in values:
                return "high_or_low_risk"
            elif "high_risk" in values:
                return "high_risk"
            elif "low_risk" in values:
                return "low_risk"
            else:
                return "uninformative"

        # Apply the function to each row and assign to new column
        self.snp_df["snp_risk_category_summary"] = self.snp_df.apply(categorize, axis=1)

        # Define the order of categories
        ordered_categories = [
            "high_risk",
            "high_or_low_risk",
            "low_risk",
            "uninformative",
            "ADO",
            "miscall",
            "NoCall",
            "NoCall_in_trio",
        ]

        # Convert the new column to an ordered category
        self.snp_df["snp_risk_category_summary"] = pd.Categorical(
            self.snp_df["snp_risk_category_summary"],
            categories=ordered_categories,
            ordered=True,
        )

    def get_snp_df(self):
        return self.snp_df

    def group_informative_snps_by_region(self) -> pd.DataFrame:
        """Calculates a summary of informative SNPs and sets the respective attribute."""
        if self.inheritance_logic:
            informative_snps_summary = self.summary_logic.summarise(
                self.snp_df, self.family_data.consanguineous
            )
        else:
            raise ValueError(
                "Invalid inheritance logic or mode of inheritance not set!"
            )
        if informative_snps_summary is None:
            raise ValueError(
                "Failed to calculate informative SNPs summary. Expected a DataFrame."
            )

        if not isinstance(informative_snps_summary, pd.DataFrame):
            raise ValueError(
                f"Failed to calculate informative SNPs summary. Expected a DataFrame, but got {type(informative_snps_summary)}."
            )
        return informative_snps_summary

    def summarise_test_data(self):
        """Calculates a summary of informative SNPs and sets the respective attribute."""
        if self.inheritance_logic:
            self.pytest_format_snp_df = self.summary_logic.get_test_data(
                self.snp_df, self.family_data.consanguineous
            )
        else:
            raise ValueError(
                "Invalid inheritance logic or mode of inheritance not set!"
            )

    def get_test_df(self):
        return self.pytest_format_snp_df


class InheritanceLogicFactory:
    @staticmethod
    def create_logic(mode_of_inheritance):
        """
        Select the appropriate logic class based on the mode of inheritance.
        """
        logic_map = {
            InheritanceMode.AUTOSOMAL_DOMINANT: AutosomalDominantLogic,
            InheritanceMode.AUTOSOMAL_RECESSIVE: AutosomalRecessiveLogic,
            InheritanceMode.X_LINKED: XLinkedLogic,
        }
        return logic_map.get(mode_of_inheritance)


class SummaryStrategy(ABC):
    def __init__(self, custom_order):
        self.custom_order = custom_order

    @abstractmethod
    def summarise(self, df, consanguineous):
        pass


class BaseSummary(SummaryStrategy):
    def filter_snps(self, df, group_by_cols, snp_risk_categories):
        snps_by_region = df.value_counts(group_by_cols).reset_index(name="snp_count")
        summarised_snps = (
            snps_by_region.groupby(by=group_by_cols).sum(numeric_only=True).fillna(0)
        )
        return summarised_snps[
            np.in1d(
                summarised_snps.index.get_level_values("snp_risk_category_summary"),
                snp_risk_categories,
            )
        ]


class AutosomalDominantSummary(BaseSummary):
    def summarise(self, df: pd.DataFrame, consanguineous: bool):
        concise_df = self.filter_snps(
            df,
            ["snp_risk_category_summary", "gene_distance"],
            ["high_risk", "low_risk"],
        )
        snp_total = concise_df["snp_count"].sum()
        totals = pd.DataFrame(
            {"snp_count": [snp_total]},
            index=pd.MultiIndex.from_tuples(
                [("Total_SNPs", "")],
                names=["gene_distance", "snp_risk_category_summary"],
            ),
        )
        return pd.concat([concise_df, totals])

    def get_test_data(self, df: pd.DataFrame, consanguineous: bool):
        return self.filter_snps(
            df,
            ["snp_position", "snp_risk_category_summary"],
            ["high_risk", "low_risk"],
        )


class AutosomalRecessiveSummary(BaseSummary):
    def summarise(self, df: pd.DataFrame, consanguineous: bool):
        if consanguineous:
            valid_combinations = [
                ("male_partner", "high_risk"),
                ("male_partner", "low_risk"),
                ("both_partners", "high_or_low_risk"),
                ("female_partner", "high_risk"),
                ("female_partner", "low_risk"),
            ]
        else:
            valid_combinations = [
                ("male_partner", "high_risk"),
                ("male_partner", "low_risk"),
                ("female_partner", "high_risk"),
                ("female_partner", "low_risk"),
            ]

        # Remove 'uninformative' and 'unassigned' categories
        filtered_snps = df[
            (df["snp_risk_category_summary"] != "uninformative")
            & (df["snp_inherited_from"] != "unassigned")
        ]

        concise_df = (
            filtered_snps.groupby(
                [
                    "snp_inherited_from",
                    "snp_risk_category_summary",
                    "gene_distance",
                ]
            )
            .size()
            .reset_index(name="snp_count")
        )

        # Filter concise_df with valid combinations
        concise_df = concise_df[
            concise_df[["snp_inherited_from", "snp_risk_category_summary"]]
            .apply(tuple, axis=1)
            .isin(valid_combinations)
        ]

        # Set index back if needed
        concise_df.set_index(
            ["snp_inherited_from", "snp_risk_category_summary", "gene_distance"],
            inplace=True,
        )

        snp_total = concise_df["snp_count"].sum()
        totals = pd.DataFrame(
            {"snp_count": [snp_total]},
            index=pd.MultiIndex.from_tuples(
                [("Total_SNPs", "", "")],
                names=[
                    "snp_inherited_from",
                    "snp_risk_category_summary",
                    "gene_distance",
                ],
            ),
        )
        return pd.concat([concise_df, totals])

    def get_test_data(self, df: pd.DataFrame, consanguineous: bool):
        if consanguineous:
            valid_combinations = [
                ("male_partner", "high_risk"),
                ("male_partner", "low_risk"),
                ("both_partners", "high_or_low_risk"),
                ("female_partner", "high_risk"),
                ("female_partner", "low_risk"),
            ]
        else:
            valid_combinations = [
                ("male_partner", "high_risk"),
                ("male_partner", "low_risk"),
                ("female_partner", "high_risk"),
                ("female_partner", "low_risk"),
            ]

        # Remove 'uninformative' and 'unassigned' categories
        filtered_snps = df[
            (df["snp_risk_category_summary"] != "uninformative")
            & (df["snp_inherited_from"] != "unassigned")
        ]

        if consanguineous:
            ordered_categories = [
                "high_risk",
                "high_or_low_risk",
                "low_risk",
            ]
            inherited_from_categories = [
                "male_partner",
                "both_partners",
                "female_partner",
            ]
        else:
            ordered_categories = [
                "high_risk",
                "low_risk",
            ]
            inherited_from_categories = [
                "male_partner",
                "female_partner",
            ]

        filtered_snps["snp_risk_category_summary"] = pd.Categorical(
            filtered_snps["snp_risk_category_summary"],
            categories=ordered_categories,
            ordered=True,
        )
        filtered_snps = filtered_snps.dropna(subset=["snp_risk_category_summary"])

        filtered_snps["snp_inherited_from"] = pd.Categorical(
            filtered_snps["snp_inherited_from"],
            categories=inherited_from_categories,
            ordered=True,
        )
        filtered_snps = filtered_snps.dropna(subset=["snp_inherited_from"])

        concise_df = (
            filtered_snps.groupby(
                [
                    "snp_inherited_from",
                    "snp_risk_category_summary",
                    "snp_position",
                ]
            )
            .size()
            .reset_index(name="snp_count")
        )

        # Filter concise_df with valid combinations
        concise_df = concise_df[
            concise_df[["snp_inherited_from", "snp_risk_category_summary"]]
            .apply(tuple, axis=1)
            .isin(valid_combinations)
        ]

        # Set index back if needed
        concise_df.set_index(
            ["snp_inherited_from", "snp_risk_category_summary", "snp_position"],
            inplace=True,
        )

        return concise_df


class XLinkedSummary(BaseSummary):
    # TODO Filter out ADO and MisCalls if needed
    def summarise(self, df: DataFrame, consanguineous: bool):
        filtered_snps = df[
            (df["snp_risk_category_AB"] != "uninformative")
            & (df["snp_risk_category_AA"] != "uninformative")
            & (df["snp_risk_category_BB"] != "uninformative")
        ]

        AB_df = filtered_snps.groupby(
            [
                "snp_risk_category_AB",
                "gene_distance",
            ]
        ).size()

        BB_df = filtered_snps.groupby(
            [
                "snp_risk_category_BB",
                "gene_distance",
            ]
        ).size()

        AA_df = filtered_snps.groupby(
            [
                "snp_risk_category_AA",
                "gene_distance",
            ]
        ).size()

        AA_and_BB_df = AA_df + BB_df

        summarised_snps_by_region = pd.concat([AB_df, AA_and_BB_df], axis=1, sort=False)
        summarised_snps_by_region.columns = [
            "snp_count_female_AB",
            "snp_count_male_AA_and_BB",
        ]
        summarised_snps_by_region.index.names = [
            "snp_risk_category_summary",
            "gene_distance",
        ]

        # Filter rows where summarised_snps_by_region is either 'high_risk' or 'low_risk'
        summarised_snps_by_region = summarised_snps_by_region[
            summarised_snps_by_region.index.isin(["high_risk", "low_risk"], level=0)
        ]

        snp_count_female_AB_total = summarised_snps_by_region[
            "snp_count_female_AB"
        ].sum()
        snp_count_male_AA_and_BB_total = (
            summarised_snps_by_region["snp_count_male_AA_and_BB"].sum() / 2
        )  # Each SNP will be counted twice as it is present in both AA and BB but can
        # can only be one or the other

        totals = pd.DataFrame(
            {
                "snp_count_female_AB": [snp_count_female_AB_total],
                "snp_count_male_AA_and_BB": [snp_count_male_AA_and_BB_total],
            },
            index=pd.MultiIndex.from_tuples(
                [("Total_SNPs", "")],
                names=[
                    "snp_risk_category_summary",
                    "gene_distance",
                ],
            ),
        )
        return pd.concat([summarised_snps_by_region, totals])

    def get_test_data(self, df: pd.DataFrame, consanguineous: bool):
        filtered_snps = df[
            (df["snp_risk_category_AB"] != "uninformative")
            & (df["snp_risk_category_AA"] != "uninformative")
            & (df["snp_risk_category_BB"] != "uninformative")
        ]

        AB_df = filtered_snps.groupby(
            [
                "snp_position",
                "snp_risk_category_AB",
            ]
        ).size()

        BB_df = filtered_snps.groupby(
            [
                "snp_position",
                "snp_risk_category_BB",
            ]
        ).size()

        AA_df = filtered_snps.groupby(
            [
                "snp_position",
                "snp_risk_category_AA",
            ]
        ).size()

        test_df = pd.concat([AB_df, AA_df, BB_df], axis=1, sort=False)
        test_df.columns = [
            "snp_count_female_AB",
            "snp_count_male_AA",
            "snp_count_male_BB",
        ]
        test_df.index.names = [
            "snp_position",
            "snp_risk_category_summary",
        ]

        test_df = test_df[test_df.index.isin(["low_risk", "high_risk"], level=1)]

        return test_df


class SummaryFactory:
    @staticmethod
    def create_summary(mode_of_inheritance: InheritanceMode, max_range_mb: int):
        custom_order = custom_order_generator(max_range_mb)
        if mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT:
            return AutosomalDominantSummary(custom_order)
        elif mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE:
            return AutosomalRecessiveSummary(custom_order)
        elif mode_of_inheritance == InheritanceMode.X_LINKED:
            return XLinkedSummary(custom_order)
        else:
            raise ValueError("Invalid mode_of_inheritance")
