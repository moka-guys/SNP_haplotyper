import argparse
import json
import logging
import os
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum, auto
from io import IOBase
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import config as config  # TODO add code to use this dependency
import numpy as np
import pandas as pd
import pdfkit
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from helper_functions import (
    create_human_readable_heading,
    custom_order_generator,
    dict2html,
    format_plot_html_str,
    generate_html_plot,
    generate_pdf_plot,
    generate_plots,
    produce_html_table,
    replace_column_names,
)
from jinja2 import Environment, PackageLoader
from pandas.io.formats.style import Styler
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

# Import mode of inheritance specific code


logger = logging.getLogger("BASHer_logger")

# Add the directory containing this script to the PYTHOPATH
sys.path.append(os.path.dirname(__file__))
mod_path = Path(__file__).parent

import config as config
from EmbryoDataClass import EmbryoData
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)

# Import mode of inheritance specific code
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from FamilyDataClass import FamilyData
from ReportDataClass import ReportData, ReportGenerator
from SNPDataClass import SNPData

# Imports variables for genome_build, allow_autosomal_dominant_cases, allow_autosomal_recessive_cases,
# allow_x_linked_cases,allow_consanguineous_cases, basher_version, released_to_production


class EmbryoTableFormatter(ABC):
    """
    Abstract class for collating summary tables from multiple EmbryoData objects and formatting for display in the HTML report.
    """

    def __init__(self):
        pass

    def set_consanguinity_flag(self, consanguinity_flag: bool) -> None:
        self.consanguinity_flag = consanguinity_flag

    def add_total_row_to_df(
        self,
        df: pd.DataFrame,
        index_columns: List[str],
        filter_values_1: Optional[List[str]] = None,
        filter_column_1: Optional[str] = None,
        filter_values_2: Optional[List[str]] = None,
        filter_column_2: Optional[str] = None,
    ) -> pd.DataFrame:
        # Group by the specified columns and sum the data
        grouped_df = df.groupby(index_columns).sum()

        # Initialize a mask that selects all rows
        mask = pd.Series([True] * len(grouped_df), index=grouped_df.index)

        # Update the mask based on filter criteria, if provided
        if filter_values_1 and filter_column_1:
            mask &= grouped_df.index.get_level_values(filter_column_1).isin(
                filter_values_1
            )
        if filter_values_2 and filter_column_2:
            mask &= grouped_df.index.get_level_values(filter_column_2).isin(
                filter_values_2
            )

        # Apply the mask to the grouped DataFrame
        df = grouped_df[mask]

        # Calculate the sum for each column
        sum_series = df.sum()

        # Convert the sum to a DataFrame and transpose it
        total_df = sum_series.to_frame().T

        # Create and set index for the total row
        if len(index_columns) > 1:
            index_values = [("Total_SNPs",) + ("",) * (len(index_columns) - 1)]
            total_df.index = pd.MultiIndex.from_tuples(
                index_values, names=index_columns
            )
        else:
            total_df.index = pd.Index(["Total_SNPs"], name=index_columns[0])

        # Concatenate the total row with the original DataFrame
        df = pd.concat([df, total_df], ignore_index=False)

        return df

    def create_table_html(
        self,
        summary_df: pd.DataFrame,
        table_identifier: str,
        include_index: bool = False,
        classes: str = "table table-striped",
    ) -> str:
        # This method currently just returns the table_name for demonstration purposes.
        # Ideally, you would have the logic to produce the HTML table here.
        summary_table_html = produce_html_table(
            summary_df,
            table_identifier,
            classes,
            include_index,
        )
        return summary_table_html  # or the actual function produce_html_table


class AutosomalDominantFormatter(EmbryoTableFormatter):
    def embryo_snp_table(
        self,
        embryo_risk_summary_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        embryo_risk_summary_df = self.add_total_row_to_df(
            embryo_risk_summary_df, ["embryo_risk_category"]
        )

        embryo_risk_summary_df = replace_column_names(
            embryo_risk_summary_df, human_readable_headings
        )

        summary_snps_table = self.create_table_html(
            summary_df=embryo_risk_summary_df,
            table_identifier="summary_embryo_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_snps_table

    def embryo_snp_table_by_region(
        self,
        summary_embryo_by_region_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        summary_embryo_by_region_df = self.add_total_row_to_df(
            summary_embryo_by_region_df,
            ["embryo_risk_category", "gene_distance"],
            filter_values_1=["high_risk", "low_risk"],
            filter_column_1="embryo_risk_category",
        )

        summary_embryo_by_region_df = replace_column_names(
            summary_embryo_by_region_df, human_readable_headings
        )

        summary_embryo_by_region_table = self.create_table_html(
            summary_df=summary_embryo_by_region_df,
            table_identifier="risk_summary_per_region_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_embryo_by_region_table


class AutosomalRecessiveFormatter(EmbryoTableFormatter):
    # get consanguinity flag from family data
    def embryo_snp_table(
        self,
        summary_snps_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        if self.consanguinity_flag:
            filter_values_2 = [
                "male_partner",
                "both_partners",
                "female_partner",
            ]
        else:
            filter_values_2 = [
                "male_partner",
                "female_partner",
            ]

        # Adding a total row without additional filtering
        summary_snps_df = self.add_total_row_to_df(
            summary_snps_df,
            index_columns=["snp_inherited_from", "embryo_risk_category"],
            filter_values_1=["low_risk", "high_risk"],
            filter_column_1="embryo_risk_category",
            filter_values_2=filter_values_2,
            filter_column_2="snp_inherited_from",
        )

        summary_snps_df = replace_column_names(summary_snps_df, human_readable_headings)

        summary_snps_table = self.create_table_html(
            summary_df=summary_snps_df,
            table_identifier="summary_embryo_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_snps_table

    def embryo_snp_table_by_region(
        self,
        summary_embryo_by_region_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        if self.consanguinity_flag:
            filter_values_1 = [
                "male_partner",
                "both_partners",
                "female_partner",
            ]
        else:
            filter_values_1 = [
                "male_partner",
                "female_partner",
            ]

        # Adding a total row with filtering based on the specified values
        summary_embryo_by_region_df = self.add_total_row_to_df(
            summary_embryo_by_region_df,
            ["snp_inherited_from", "embryo_risk_category", "gene_distance"],
            filter_values_1=filter_values_1,
            filter_column_1="snp_inherited_from",
            filter_values_2=["high_risk", "low_risk"],
            filter_column_2="embryo_risk_category",
        )

        summary_embryo_by_region_df = replace_column_names(
            summary_embryo_by_region_df, human_readable_headings
        )

        summary_embryo_by_region_table = self.create_table_html(
            summary_df=summary_embryo_by_region_df,
            table_identifier="risk_summary_per_region_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_embryo_by_region_table


class XLinkedFormatter(EmbryoTableFormatter):
    def embryo_snp_table(
        self,
        embryo_risk_summary_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        embryo_risk_summary_df = self.add_total_row_to_df(
            embryo_risk_summary_df, ["embryo_risk_category"]
        )

        embryo_risk_summary_df = replace_column_names(
            embryo_risk_summary_df, human_readable_headings
        )

        summary_snps_table = self.create_table_html(
            summary_df=embryo_risk_summary_df,
            table_identifier="summary_embryo_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_snps_table

    def embryo_snp_table_by_region(
        self,
        summary_embryo_by_region_df: pd.DataFrame,
        human_readable_headings: Dict[str, str],
    ) -> str:
        summary_embryo_by_region_df = self.add_total_row_to_df(
            summary_embryo_by_region_df,
            ["embryo_risk_category", "gene_distance"],
            filter_values_1=["high_risk", "low_risk"],
            filter_column_1="embryo_risk_category",
        )

        summary_embryo_by_region_df = replace_column_names(
            summary_embryo_by_region_df, human_readable_headings
        )

        summary_embryo_by_region_table = self.create_table_html(
            summary_df=summary_embryo_by_region_df,
            table_identifier="risk_summary_per_region_table",
            include_index=True,
            classes="table table-striped",
        )

        return summary_embryo_by_region_table


def format_tables_for_html_report(self):
    # Format summary tables for display in HTML report
    SummaryFormatterClass = self.get_embryo_summary_formatter(
        self.family_data.mode_of_inheritance
    )
    summary_formatter_instance = SummaryFormatterClass()

    summary_formatter_instance.set_consanguinity_flag(self.family_data.consanguineous)

    if self.family_data.trio_only == False and self.family_data.embryo_ids:
        summary_snps_table = summary_formatter_instance.embryo_snp_table(
            self.embryo_risk_summary_df,
            self.human_readable_headings,
        )
        summary_embryo_by_region_table = (
            summary_formatter_instance.embryo_snp_table_by_region(
                self.embryo_risk_summary_by_region_df,
                self.human_readable_headings,
            )
        )
    else:
        summary_snps_table = ""
        summary_embryo_by_region_table = ""

    return summary_snps_table, summary_embryo_by_region_table


class SNPAnalysis:
    family_data: FamilyData
    snp_data_df: pd.DataFrame
    """
    Pipeline for SNP Analysis using family data and external files.
    """

    def __init__(self, family_data: FamilyData, snp_data_df: pd.DataFrame):
        """
        Initializes the SNPAnalysisPipeline with a FamilyData object and a dataframe produced from a ChAS csv output file.

        :param family_data: A FamilyData object containing relevant genetic and family information.
        :param file_path: Path to an external file for further analysis.
        """
        self.family_data = family_data
        self.create_embryo_sex_lookup()
        self.human_readable_headings = create_human_readable_heading(
            embryo_sex_lookup=self.embryo_sex_lookup,
            male_partner=self.family_data.male_partner,
            female_partner=self.family_data.female_partner,
            reference=self.family_data.reference,
        )
        self.initialize_snp_data(snp_data_df)
        self.calculate_and_set_qc_metrics()
        if self.family_data.trio_only == False and self.family_data.embryo_ids:
            self.initialize_embryos_dict()
            self.embryo_summary_df = SNPAnalysis.collate_embryo_results(
                self.snp_data.snp_df, self.embryos
            )
            # Create summary embryo results table across all embryos
            self.embryo_risk_summary_df = self.summarise_embryo_results(
                self.embryos, "risk_summary_df"
            )
            self.embryo_risk_summary_by_region_df = self.summarise_embryo_results(
                self.embryos, "risk_summary_per_region_df"
            )

            self.pytest_format_embryo_df = self.summarise_embryo_results(
                self.embryos, "risk_summary_for_testing_df"
            )

        (
            self.summary_snps_table,
            self.summary_embryo_by_region_table,
        ) = format_tables_for_html_report(self)
        self.initialise_report_data()

    def create_embryo_sex_lookup(self) -> None:
        """Creates a lookup dictionary for embryo sex based on provided family data."""
        if self.family_data.embryo_sex and self.family_data.embryo_ids:
            self.embryo_sex_lookup = dict(
                zip(self.family_data.embryo_ids, self.family_data.embryo_sex)
            )
        else:
            self.embryo_sex_lookup = {}

    def initialize_snp_data(self, snp_df: pd.DataFrame) -> None:
        """Initializes the SNP data attribute."""
        self.snp_data = SNPData(family_data=self.family_data, snp_df=snp_df)

    def calculate_and_set_qc_metrics(self) -> None:
        """Calculates and sets the QC metrics attribute."""
        self.qc_metrics = self.calculate_qc_metrics()

    def initialize_embryos_dict(self) -> None:
        """Initializes the embryos dictionary and creates an EmbryoData object for each embryo ID."""
        self.embryos: Dict[str, "EmbryoData"] = {}
        for embryo_id in self.family_data.embryo_ids:
            self.add_embryo(embryo_id)

    def filter_out_nocalls(self) -> None:
        """
        Adds a 'filter_out_nocalls' column in snp_data DataFrame.

        This method creates a new column 'filter_out_nocalls' in the snp_data DataFrame. If the male partner, female partner, or reference
        has "NoCall" for a probeset then this probeset is marked as False in the 'filter_out_nocalls' column. For the x-linked conditions,
        NoCalls in male partner data are not considered, while for the autosomal and recessive cases, NoCalls in any of the male partner,
        female partner or reference data are considered.

        The 'filter_out_nocalls' column can be used later for filtering out these probesets. If a probeset has 'False' in this column, it
        should be filtered out, while if it has 'True', it should be retained.

        Note:
            This method does not return anything. It modifies the snp_data DataFrame in-place.
        """
        # TODO Check the logic in this function
        if self.family_data.mode_of_inheritance == InheritanceMode.X_LINKED:
            # For x-link cases we do not care if the male sample is a NoCall as only female or reference NoCalls should be filtered out
            self.snp_data_df["filter_out_nocalls"] = (
                self.snp_data_df[self.family_data.female_partner] != "NoCall"
            ) | (self.snp_data_df[self.family_data.reference] != "NoCall")
        else:
            # TODO Check the logic in this function
            # For autosomal and recessive cases all NoCalls should be filtered out
            self.snp_data_df["filter_out_nocalls"] = (
                (self.snp_data_df[self.family_data.male_partner] != "NoCall")
                & (self.snp_data_df[self.family_data.female_partner] != "NoCall")
                & (self.snp_data_df[self.family_data.reference] != "NoCall")
            )

    def calculate_qc_metrics(self) -> Styler:
        """
        Calculates QC metrics based on the number of NoCalls per sample (measure of DNA quality)
        and also the percentage of NoCalls per sample.

        This method calculates QC metrics based on the number of NoCalls per sample which can be used as a metric of DNA quality,
        and also calculates the percentage of NoCalls per sample.

        Returns:
            A DataFrame containing the calculated QC metrics.
        """
        samples = [
            self.family_data.male_partner,
            self.family_data.female_partner,
            self.family_data.reference,
            *self.family_data.embryo_ids,
        ]

        qc_metrics = (
            pd.DataFrame(
                {
                    sample: self.snp_data.snp_df[sample].value_counts()
                    for sample in samples
                }
            )
            .reindex(["AA", "BB", "AB", "NoCall"])
            .fillna(0)
            .astype(int)
        )

        qc_metrics = replace_column_names(qc_metrics, self.human_readable_headings)

        # Calculate the percentage of NoCalls
        summary_styler = (
            qc_metrics.agg(
                [
                    "sum",
                ]
            )
            .style.format(precision=2)
            .relabel_index(
                [
                    "Total",
                ]
            )
        )

        percentage_styler = (
            (qc_metrics.loc["NoCall"] / qc_metrics.sum())
            .to_frame()
            .transpose()
            .style.format("{:.2%}")
            .relabel_index(
                [
                    "NoCall %",
                ]
            )
        )

        qc_metrics = (
            qc_metrics.style.format(precision=1)
            .concat(summary_styler)
            .concat(percentage_styler)
        )

        return qc_metrics

    def add_embryo(self, embryo_id):
        embryo_category_df = self.snp_data.snp_df.copy()

        # Remove current embryo_id from list of embryo IDs to create a list of irrelevant columns for this embryo
        irrelevant_columns = list(
            filter(lambda x: x != embryo_id, self.family_data.embryo_ids)
        )

        # Remove irrelevant columns from embryo_category_df
        embryo_category_df.drop(columns=irrelevant_columns, inplace=True)

        self.embryos[embryo_id] = EmbryoData(
            embryo_id=embryo_id,
            embryo_sex=self.embryo_sex_lookup[embryo_id],
            family_genetic_data=self.family_data,
            embryo_category_df=embryo_category_df,
        )

    @staticmethod
    def collate_embryo_results(snp_data_df, embryos):
        """
        This static method takes a DataFrame snp_data_df and a dictionary embryos,
        and returns a DataFrame with added columns based on embryo data.

        :param snp_data_df: DataFrame with SNP data
        :param embryos: Dictionary with embryo data
        :return: Updated DataFrame
        """
        embryo_results_df = snp_data_df[
            ["probeset_id", "rsID", "Position", "gene_distance", "snp_position"]
        ].copy()

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

        for embryo_id, embryo_data in embryos.items():
            # Assign the new column and convert it to the categorical type
            embryo_results_df[embryo_id] = embryo_data.embryo_category_df[
                "embryo_risk_category"
            ].astype(cat_type)

        return embryo_results_df

    @staticmethod
    def summarise_embryo_results(
        embryos: Dict[str, EmbryoData], attribute_name: str
    ) -> pd.DataFrame:
        """
        Summarise embryo results for each embryo in embryo_ids
        """
        df_list = []

        for embryo_id in embryos:
            embryo_data = embryos[embryo_id]

            # Check if the specified attribute exists in the EmbryoData object
            if (
                hasattr(embryo_data, attribute_name)
                and getattr(embryo_data, attribute_name) is not None
            ):
                df_list.append(getattr(embryo_data, attribute_name))

        # Initialize merged_df as None
        merged_df = None

        # Iterate through the list of DataFrames and merge them one by one
        for df in df_list:
            if merged_df is None:
                # If merged_df is None, set it to the first DataFrame
                merged_df = df
            else:
                # Merge the current DataFrame with the existing merged_df
                merged_df = merged_df.merge(df)

        if merged_df is None:
            # If merged_df is still None, return an empty DataFrame
            return pd.DataFrame()
        return merged_df

    def collate_figures(self, embryos: Dict[str, EmbryoData]) -> Dict[str, str]:
        # Reiterates over all the EmbryoData objects and collates the figures
        figures_dict = {}
        for embryo_id, embryo_data in embryos.items():
            embryo_data.results_plot
            figures_dict[embryo_id] = embryo_data.results_plot
        return figures_dict

    def get_qc_metrics(self, format="dataframe"):
        # If user wants the data in HTML format
        # TODO Maybe delete this method, standardise return type
        if format == "html":
            return self.qc_metrics.to_html(index=False)
        return self.qc_metrics

    def get_embryo(self, embryo_id=None):
        """
        Retrieve embryos from the dictionary.

        Parameters:
        - embryo_id: If provided, returns the embryo with the matching ID.

        Returns:
        - A single EmbryoData object if embryo_id is provided.
        - The whole embryos dictionary otherwise.
        """
        if embryo_id:
            return self.embryos.get(embryo_id, None)

        return self.embryos

    def get_table(
        self,
        summary_df,
        table_name,
        index=False,
        classes="table table-striped",
    ):
        table = produce_html_table(summary_df, table_name, classes, index)
        return table

    @staticmethod
    def get_embryo_summary_formatter(mode_of_inheritance: InheritanceMode):
        """
        Formats summary reports for display in the HTML report.  Different modes of inheritance require different formatting.
        """
        formatter = {
            InheritanceMode.AUTOSOMAL_DOMINANT: AutosomalDominantFormatter,
            InheritanceMode.AUTOSOMAL_RECESSIVE: AutosomalRecessiveFormatter,
            InheritanceMode.X_LINKED: XLinkedFormatter,
        }
        return formatter.get(mode_of_inheritance, None)

    def initialise_report_data(self):
        self.report_data = ReportData(
            header_html=dict2html(self.family_data.report_header_info),
            mode_of_inheritance=str(self.family_data.mode_of_inheritance.value),
            gene_symbol=self.family_data.gene_symbol,
            chromosome=self.family_data.chr.value.upper(),  # corrected typo here
            gene_start=f"{int(self.family_data.gene_start):,}",
            gene_end=f"{int(self.family_data.gene_end):,}",
            genome_build=config.genome_build,
            basher_version=config.basher_version,
            input_file=self.family_data.input_ChAS_filepath.name
            if isinstance(self.family_data.input_ChAS_filepath, Path)
            else self.family_data.input_ChAS_filepath,
            male_partner=self.family_data.male_partner,
            male_partner_status=self.family_data.male_partner_status.value,
            female_partner=self.family_data.female_partner,
            female_partner_status=self.family_data.female_partner_status.value,
            reference=self.family_data.reference,
            reference_status=self.family_data.reference_status.value,
            reference_relationship=self.family_data.reference_relationship.value,
            results_table_1=produce_html_table(self.snp_data.snp_df, "results_table_1"),
            qc_table=self.qc_metrics.to_html(
                table_uuid="qc_table",
                index=True,
            ),
            report_date=datetime.today().strftime("%Y-%m-%d %H:%M:%S"),
            summary_snps_table=produce_html_table(
                summary_df=self.snp_data.informative_snps_summary,
                table_identifier="summary_snps_table",
                include_index=True,
                classes="table table-striped",
            ),
            summary_embryo_table=self.summary_snps_table,
            summary_embryo_by_region_table=self.summary_embryo_by_region_table,
            # If self.family_data.trio_only is True, assign an empty string, else assign the variable
            html_text_for_plots=""
            if self.family_data.trio_only
            else format_plot_html_str(
                generate_plots(self.collate_figures(self.embryos), static_plots=False),
                add_dropdown_selection=False,
            ),
            # If self.family_data.trio_only is True, assign an empty string, else assign the variable
            pdf_text_for_plots=""
            if self.family_data.trio_only
            else format_plot_html_str(
                generate_plots(self.collate_figures(self.embryos), static_plots=True),
                add_dropdown_selection=True,
            ),
            text_for_plots="",  # Dynamically updated when rendered selecting either html or pdf
            warning=""
            if config.released_to_production
            else f"<h3 style='color:red'> This is a pre-release version of the BASHer tool. Please contact the BASHer team if you have any questions.</h3>",
            consanguinity_flag="non-consanguineous"
            if self.family_data.consanguineous == False
            else "consanguineous",
        )
