import argparse
from io import IOBase
from jinja2 import Environment, PackageLoader
import json
import os
import pandas as pd
from pathlib import Path
import pdfkit
import numpy as np
from datetime import datetime

import sys

import logging

logger = logging.getLogger("BASHer_logger")

# Add the directory containing this script to the PYTHOPATH
sys.path.append(os.path.dirname(__file__))
mod_path = Path(__file__).parent

# Import mode of inheritance specific code
from autosomal_dominant_logic import autosomal_dominant_analysis
from autosomal_recessive_logic import autosomal_recessive_analysis
import config as config

# Imports variables for genome_build, allow_autosomal_dominant_cases, allow_autosomal_recessive_cases,
# allow_x_linked_cases,allow_consanguineous_cases, basher_version, released_to_production

from x_linked_logic import x_linked_analysis
from snp_plot import plot_results

from exceptions import ArgumentInputError, InvalidParameterSelectedError

# TODO Copy rsID from hover tap
# TODO Check telomeric/centromeric genes work with 2mb window (FHSD1 - D4Z4 repeat, PKD1)
# TODO Add support for no embryos (just TRIOs being run to check if enough informative SNPs)
# TODO Add ADO % to table

# Import environment variables set by docker-compose
UPLOAD_FOLDER = os.getenv("UPLOAD_FOLDER")

# Import command line arguments (these can be automatically generated from the sample sheet using sample_sheet_reader.py)
parser = argparse.ArgumentParser(description="SNP Haplotying from SNP Array data")

# File input/output data
parser.add_argument(
    "-i",
    "--input_file",
    type=str,
    help="Input txt file containing SNP Array output",
)

parser.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    help="Output filename prefix",
)

parser.add_argument(
    "-f",
    "--output_folder",
    type=str,
    help="Output folder path",
)

# Patient data
parser.add_argument(
    "-m",
    "--mode_of_inheritance",
    type=str,
    choices=["autosomal_dominant", "autosomal_recessive", "x_linked"],
    help="The mode of inheritance",
)

parser.add_argument(
    "-mp",
    "--male_partner",
    type=str,
    help="ID in input table for male_partner",
)

parser.add_argument(
    "-mps",
    "--male_partner_status",
    type=str,
    choices=["affected", "unaffected", "carrier"],
    help="Status of male_partner",
)

parser.add_argument(
    "-fp",
    "--female_partner",
    type=str,
    help="ID in input table for male_partner",
)

parser.add_argument(
    "-fps",
    "--female_partner_status",
    choices=["affected", "unaffected", "carrier"],
    type=str,
    help="ID in input table for female_partner",
)

parser.add_argument(
    "-consang",
    "--consanguineous",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to indicate that partners are consanguineous",
)

parser.add_argument(
    "-r",
    "--reference",
    type=str,
    help="ID in input table for reference sample",
)

parser.add_argument(
    "-rs",
    "--reference_status",
    type=str,
    choices=["affected", "unaffected", "carrier"],
    help="Status of Reference",
)

parser.add_argument(
    "-rr",
    "--reference_relationship",
    type=str,
    choices=[
        "grandparent",
        "child",
    ],
    help="Reference relationship to pro-band",
)

parser.add_argument(
    "-e",
    "--embryo_ids",
    nargs="+",
    type=str,
    help="IDs of embryos in the input table",
)

parser.add_argument(
    "-es",
    "--embryo_sex",
    nargs="+",
    type=str,
    choices=[
        "male",
        "female",
        "unknown",
    ],
    help="Embryo sex - must be in same order as embryo_ids (cannot be unknown for X-linked diseases)",
)

# Gene/ROI data
parser.add_argument(
    "-g",
    "--gene_symbol",
    type=str,
    help="Gene Symbol",
)

parser.add_argument(
    "-gs",
    "--gene_start",
    type=int,
    help="Gene Start genomic co-ordinate (1-based referencing)",
)

parser.add_argument(
    "-ge",
    "--gene_end",
    type=int,
    help="Gene End genomic co-ordinate (1-based referencing)",
)

parser.add_argument(
    "-c",
    "--chr",
    type=str,
    choices=[
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
        "x",
        "y",
    ],
    help="Chromosome of ROI/gene",
)

parser.add_argument(
    "--flanking_region_size",
    type=str,
    nargs="?",
    choices=["2mb", "3mb", "4mb", "5mb", "6mb", "7mb", "8mb", "9mb", "10mb"],
    const="2mb",
    help="Size of the flanking region either side of the gene",
)

parser.add_argument(
    "--trio_only",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to produce a preliminary report without looking at embryos, must be used if not embryo data is provided.",
)

parser.add_argument(
    "--header_info",
    type=str,
    required=True,
    help="Pass a string to populate the report header. A field will be created for each entry field_title=field_value separated by ';', for example 'PRU=1234;Hospital No=1234;Biopsy No=111' will produce 3 fields in the header with the titles PRU, Hospital No, and Biopsy No.",
)

# If no arguments are provided, print the help message
# sys.argv includes a list of elements starting with the program
# if len(sys.argv) < 2:
#     parser.print_help()
#     parser.exit()


def header_to_dict(header_str):  # TODO keep this function
    """
    Converts a string of header_info into a dictionary
    Args:
        header_info (str): A string in the key=value pairs like "PRU=1234;Hospital No=1234;Biopsy No:111", where the keys will be the titles of the fields in the header
    Returns:
        dict: A dictionary of the header info with field titles as keys and values as values
    """
    if header_str is None:
        return None
    else:
        d = dict(x.split("=") for x in header_str.split(";"))
        return d


def add_rsid_column(df, affy_2_rs_ids_df):  # TODO remove this function
    """Provides dbsnp rsIDs
    New column created in the dataframe, df, matching the probes_set IDs to dbSNP rsIDs.
    Args:
        df (dataframe): A dataframe with a "probeset_id" column
        affy_2_rs_ids_df (dataframe): A dataframe with columns "probeset_id" & "rsID" used to map between the 2 identifiers
    Returns:
        dataframe: Original dataframe, df, with columns for "rsID" added next to the "probeset_id" column (these columns are now the 1st columns of the dataframe)
    """
    df = pd.merge(
        df, affy_2_rs_ids_df[["probeset_id", "rsID"]], on="probeset_id", how="left"
    )

    # Rearrange columns so that rsID is next to Affy Id
    df.insert(0, "probeset_id", df.pop("probeset_id"))
    df.insert(1, "rsID", df.pop("rsID"))
    return df


def export_json_data_as_csv(input_json, output_csv):
    """Import a JSON file and save the data as a CSV
    Imports a simple JSON file and exports it as a CSV file.
    Used to import test data from JSON files (informative_snp_validation.json, embryo_validation_data.json,
    launch.json) and export it as human readable CSV.  These CSV can be shared with Genomic Scientists during
    the validation process.
    Args:
        input_json (string): The path to a JSON file
        output_csv (string): The path and filename for the output csv file
    """
    p = Path(input_json)
    with p.open("r", encoding="utf-8") as f:
        data = json.loads(f.read())
    df = pd.json_normalize(data)
    df.to_csv(output_csv, index=False, encoding="utf-8")


# filter dataframe on region of interest
def filter_dataframe(
    df, gene_start, gene_end, flanking_region_size
):  # TODO remove this function
    """
    Filters a dataframe to only include rows where the SNP is within the region of interest.
    Args:
        df (pandas dataframe): Dataframe containing SNP data
        int(args.gene_start) (int): Start position of gene of interest
        args.gene_end (int): End position of gene of interest
        args.flanking_region_size (str): Size of flanking region either side of gene of interest "2mb" or "3mb"
    Returns:
        df (pandas dataframe): Dataframe containing only SNPs within the region of interest
    """
    if flanking_region_size == "2mb":
        region_start = int(gene_start) - 2000000
        region_end = int(gene_end) + 2000000
    elif flanking_region_size == "3mb":
        region_start = int(gene_start) - 3000000
        region_end = int(gene_end) + 3000000
    df = df[df["Position"] >= region_start]
    df = df[df["Position"] <= region_end]
    return df


def annotate_distance_from_gene(df, chr, start, end):  # TODO remove this function
    """Annotates the probeset based on the provided genomic co-ordinates
    New column created, "gene_distance", in the dataframe, df, annotating the region the SNP is in. SNPs allocated to "within_gene", "0-1MB_from_start", "1-2MB_from_start", "0-1MB_from_end", and "1-2MB_from_end",
    Args:
        df (dataframe): A dataframe with a "probeset_id" column and the feature's genomic co-ordinates,  "Position"
        chr (string):  The chromsome the gene of interest is on
        start (int): The start coordinate of the gene (1-based)
        end (int): The end coordinate of the gene (1-based)
    Returns:
        dataframe: Original dataframe, df, with "gene_distance" column added characterising the probeset in relation to the gene of interest
    """
    # insure inputs are integers
    start = int(start)
    end = int(end)

    conditions = [
        (df["Position"] > start) & (df["Position"] <= end),
        (df["Position"] <= start) & (df["Position"] > start - 1000000),
        (df["Position"] <= start - 1000000) & (df["Position"] > start - 2000000),
        (df["Position"] > end) & (df["Position"] <= end + 1000000),
        (df["Position"] > end + 1000000) & (df["Position"] <= end + 2000000),
    ]

    values = [
        "within_gene",
        "0-1MB_from_start",
        "1-2MB_from_start",
        "0-1MB_from_end",
        "1-2MB_from_end",
    ]

    # # TODO check correct chr has been given and boundaries are correct
    df["gene_distance"] = np.select(conditions, values, default="outside_range")

    df["gene_distance"] = df["gene_distance"].astype("category")

    df["gene_distance"] = df["gene_distance"].cat.set_categories(
        [
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
        ],
    )
    return df


def filter_out_nocalls(
    df, male_partner, female_partner, reference, filter_male_nocalls=True
):
    """Filters out no calls
    If the male partner, female partner, or reference has "NoCall" for a probeset then this probeset should be filtered out.
    Args:
        df (dataframe): A dataframe with the SNP array data
        male_partner (string):  Column name representing the data for the male partner
        female_partner (string):  Column name representing the data for the female partner
        reference (string):  Column name representing the data for the reference
        filter_male_nocalls (boolean): Should male partner NoCalls be filtered for the analysis (for male embryos in x-linked conditions they should be retained)
    Returns:
        dataframe: Original dataframe, df, with any rows where the male partner, female partner or reference has a "NoCall" filtered out
    """
    if filter_male_nocalls == True:
        filtered_df = df[
            (
                (df[male_partner] != "NoCall")
                & (df[female_partner] != "NoCall")
                & (df[reference] != "NoCall")
            )
        ]
    elif filter_male_nocalls == False:
        filtered_df = df[
            ((df[female_partner] != "NoCall") | (df[reference] != "NoCall"))
        ]

    # TODO add logger - how many NoCalls filtered
    return filtered_df


# TODO standardise the order of fp/mp args across functions


def calculate_qc_metrics(df, male_partner, female_partner, reference, embryo_ids=None):
    """Calculate QC metrics based on the number of NoCalls per sample (measure of DNA quality)
    Calculate QC metrics based on the number of NoCalls per sample which can be used as a metric of DNA quality.
    Args:
        df (dataframe): A dataframe with the SNP array data
        male_partner (string):  Column name representing the data for the male partner
        female_partner (string):  Column name representing the data for the female partner
        reference (string):  Column name representing the data for the reference
        embryo_ids (list): List of column names representing the data for 1>n embryo samples
    Returns:
        dataframe: Dataframe summarising the number of NoCalls per sample
    """
    # Initiate dataframe
    qc_df = pd.DataFrame(index=["AA", "BB", "AB", "NoCall"])
    # Populate dataframe
    qc_df[female_partner] = [
        df[df[female_partner] == "AA"].shape[0],
        df[df[female_partner] == "BB"].shape[0],
        df[df[female_partner] == "AB"].shape[0],
        df[df[female_partner] == "NoCall"].shape[0],
    ]
    qc_df[male_partner] = [
        df[df[male_partner] == "AA"].shape[0],
        df[df[male_partner] == "BB"].shape[0],
        df[df[male_partner] == "AB"].shape[0],
        df[df[male_partner] == "NoCall"].shape[0],
    ]
    qc_df[reference] = [
        df[df[reference] == "AA"].shape[0],
        df[df[reference] == "BB"].shape[0],
        df[df[reference] == "AB"].shape[0],
        df[df[reference] == "NoCall"].shape[0],
    ]
    if embryo_ids is not None:
        for embryo in embryo_ids:
            qc_df[embryo] = [
                df[df[embryo] == "AA"].shape[0],
                df[df[embryo] == "BB"].shape[0],
                df[df[embryo] == "AB"].shape[0],
                df[df[embryo] == "NoCall"].shape[0],
            ]
    # Clean up dataframe
    qc_df = qc_df.reset_index()
    qc_df = qc_df.rename(
        columns={"index": "call_type"},
    )
    return qc_df


def calculate_nocall_percentages(df):
    """Calculate the percentage of nocalls
    Takes a dataframe produced from calculate_qc_metrics() and calculates
    the % of nocalls per sample.
    Args:
        df (dataframe): A dataframe produced by calculate_qc_metrics()
    Returns:
        dataframe: Dataframe summarising the % of NoCalls per sample,
        with the same format as the dataframe produced by calculate_qc_metrics()
        Essentially provides a row for the % of NoCalls per sample
    """
    nocall = df[df["call_type"] == "NoCall"]
    if nocall.shape[0] > 0:
        # Remove non-numeric field 'call_type' from the dataframe
        numeric_nocall_df = nocall.drop("call_type", axis=1)
        numeric_df = df.drop("call_type", axis=1)
        # Calculate the percentage of NoCalls
        nocall_percentage = (
            numeric_nocall_df / numeric_df.sum(axis=0, numeric_only=True) * 100
        )
        # Add descriptive column
        nocall_percentage.insert(0, "call_type", "NoCall")
        return nocall_percentage
    else:
        # TODO add logger message
        pass


def detect_miscall_or_ado(
    male_partner_haplotype, female_partner_haplotype, embryo_haplotype
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
        [male_partner_haplotype, female_partner_haplotype, embryo_haplotype]
    ) - set(["AA", "BB", "AB", "NoCall"])
    if len(illegal_args) != 0:
        raise ArgumentInputError(
            f"Function detect_miscall_or_ado() only excepts 'AA','BB', 'AB', NoCall' as arguments, recieved {str(illegal_args)}"
        )

    parent_alleles = [male_partner_haplotype, female_partner_haplotype]
    if embryo_haplotype == "NoCall":
        result = "NoCall"
    else:
        match parent_alleles:
            case ["AA", "AA"]:
                result = "miscall" if embryo_haplotype != "AA" else "call"
            case ["BB", "BB"]:
                result = "miscall" if embryo_haplotype != "BB" else "call"
            case ["AA", "BB"]:
                result = "ADO" if embryo_haplotype != "AB" else "call"
            case ["BB", "AA"]:
                result = "ADO" if embryo_haplotype != "AB" else "call"
            case ["AA", "AB"]:
                result = (
                    "ADO"
                    if embryo_haplotype
                    not in [
                        "AA",
                        "AB",
                    ]
                    else "call"
                )
            case ["AB", "AA"]:
                result = (
                    "ADO"
                    if embryo_haplotype
                    not in [
                        "AA",
                        "AB",
                    ]
                    else "call"
                )
            case ["BB", "AB"]:
                result = (
                    "ADO"
                    if embryo_haplotype
                    not in [
                        "BB",
                        "AB",
                    ]
                    else "call"
                )
            case ["AB", "BB"]:
                result = (
                    "ADO"
                    if embryo_haplotype
                    not in [
                        "BB",
                        "AB",
                    ]
                    else "call"
                )
            case ["AB", "AB"]:
                result = (
                    "Error! Parent haplotypes cannot be AB & AB these shoulds have been filtered out"
                    if embryo_haplotype
                    not in [
                        "AA",
                        "BB",
                        "AB",
                    ]
                    else "call"
                )
    return result


def snps_by_region(df, mode_of_inheritance):
    """Summarise the number of SNPs by regions around the gene of interest
    Takes a results_df dataframe produced from either autosomal_dominant_analysis(),
    autosomal_recessive_analysis(), or x_linked_analysis() and counts
    the SNPs per "gene_distance":
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
     and "snp_risk_category":
            "low_risk",
            "high_risk",
    For autosomal dominant one "snp_risk_category" column is produced for where the embryo SNP is AB,
    for x-linked three columns are produced for where the embryo SNP is female_AB, male_AA, or male_BB,
    for autosomal recessive cases an "snp_inherited_from" is also added to show which partner the SNP
    is inherited from .
    Args:
        df (dataframe): A dataframe produced by either autosomal_dominant_analysis(),
    autosomal_dominant_analysis(), or x_linked_analysis()
    Returns:
        dataframe: Dataframe summarising the SNPs per genome region with additional columns for each relevant haplotype in the embryo.
    """

    if mode_of_inheritance == "autosomal_dominant":
        snps_by_region = df.value_counts(
            ["gene_distance", "snp_risk_category"]
        ).to_frame()
        # Extract snps_by_region data from index into columns
        snps_by_region = snps_by_region.reset_index()
        # Rename columns
        snps_by_region.columns = [
            "gene_distance",
            "snp_risk_category",
            "snp_count",
        ]
    elif mode_of_inheritance == "autosomal_recessive":
        snps_by_region = df.value_counts(
            ["gene_distance", "snp_risk_category", "snp_inherited_from"]
        ).to_frame()
        # Extract snps_by_region data from index into columns
        snps_by_region = snps_by_region.reset_index()
        # Rename columns
        snps_by_region.columns = [
            "gene_distance",
            "snp_risk_category",
            "snp_inherited_from",
            "snp_count",
        ]
    elif mode_of_inheritance == "x_linked":
        snps_by_region_female_AB = df.value_counts(
            ["gene_distance", "female_AB_snp_risk_category"]
        ).reset_index()
        snps_by_region_male_AA = df.value_counts(
            ["gene_distance", "male_AA_snp_risk_category"]
        ).reset_index()
        snps_by_region_male_BB = df.value_counts(
            ["gene_distance", "male_BB_snp_risk_category"]
        ).reset_index()

        # Merge dataframes into single dataframe
        snps_by_region = pd.merge(
            snps_by_region_female_AB,
            snps_by_region_male_AA,
            how="outer",
            left_on=["gene_distance", "female_AB_snp_risk_category"],
            right_on=["gene_distance", "male_AA_snp_risk_category"],
        )

        snps_by_region = pd.merge(
            snps_by_region,
            snps_by_region_male_BB,
            how="outer",
            left_on=["gene_distance", "female_AB_snp_risk_category"],
            right_on=["gene_distance", "male_BB_snp_risk_category"],
        )
        # Set column names
        snps_by_region.columns = [
            "gene_distance",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "male_AA_snp_risk_category",
            "male_AA_snp_count",
            "male_BB_snp_risk_category",
            "male_BB_snp_count",
        ]
        # Replace any NaNs in numerical columns with 0s
        snps_by_region["female_AB_snp_count"] = (
            snps_by_region["female_AB_snp_count"].fillna(0).astype("Int64")
        )
        snps_by_region["male_AA_snp_count"] = (
            snps_by_region["male_AA_snp_count"].fillna(0).astype("Int64")
        )
        snps_by_region["male_BB_snp_count"] = (
            snps_by_region["male_BB_snp_count"].fillna(0).astype("Int64")
        )

        # Fill any gaps in the risk categories if no SNPs are for a particular category
        snps_by_region.male_AA_snp_risk_category.fillna(
            snps_by_region.female_AB_snp_risk_category, inplace=True
        )
        snps_by_region.male_AA_snp_risk_category.fillna(
            snps_by_region.male_BB_snp_risk_category, inplace=True
        )
        snps_by_region.female_AB_snp_risk_category.fillna(
            snps_by_region.male_AA_snp_risk_category, inplace=True
        )
        snps_by_region.male_BB_snp_risk_category.fillna(
            snps_by_region.male_AA_snp_risk_category, inplace=True
        )

    else:
        pass  # TODO exception
    return snps_by_region


def summarised_snps_by_region(df, mode_of_inheritance):
    """
    Summarises the SNPs per genome region for a given mode of inheritance.
    """
    if mode_of_inheritance == "autosomal_dominant":
        # Filter out "uninformative" from summary
        categorised_snps_by_region = df[df["snp_risk_category"] != "uninformative"]
        # Group informative 'low_risk' and 'high_risk' SNPs together per region
        summary_categorised_snps_by_region = categorised_snps_by_region.groupby(
            by=["gene_distance"]
        ).sum(numeric_only=True)
        # Ensure rows are in logical order
        summary_categorised_snps_by_region = summary_categorised_snps_by_region.reindex(
            [
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
            ]
        )
        # Replace any NaN with 0 and convert any floats to ints (caused when NaN included in column)
        summary_categorised_snps_by_region = summary_categorised_snps_by_region.fillna(
            0
        )
        summary_categorised_snps_by_region[
            "snp_count"
        ] = summary_categorised_snps_by_region["snp_count"].astype(int)
        # Calculate totals
        summary_categorised_snps_by_region.loc[
            "total_snps"
        ] = summary_categorised_snps_by_region.sum(numeric_only=True, axis=0)
        summary_categorised_snps_by_region = (
            summary_categorised_snps_by_region.reset_index()
        )
    elif mode_of_inheritance == "autosomal_recessive":
        # Filter out "uninformative" from summary
        categorised_snps_by_region = df[df["snp_risk_category"] != "uninformative"]
        # Group informative 'low_risk' and 'high_risk' SNPs together per region
        summary_categorised_snps_by_region = (
            categorised_snps_by_region.groupby(
                by=["gene_distance", "snp_risk_category", "snp_inherited_from"]
            )
            .sum(numeric_only=True)
            .reset_index()
        )

        # Ensure rows are in logical order
        summary_categorised_snps_by_region["snp_inherited_from"] = pd.Categorical(
            summary_categorised_snps_by_region["snp_inherited_from"],
            ["male_partner", "female_partner"],
        )

        summary_categorised_snps_by_region["gene_distance"] = pd.Categorical(
            summary_categorised_snps_by_region["gene_distance"],
            [
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
            ],
        )

        summary_categorised_snps_by_region = (
            summary_categorised_snps_by_region.sort_values(
                by=["snp_inherited_from", "snp_risk_category", "gene_distance"],
                ascending=[True, True, True],
            )
        )

        # Replace any NaN with 0 and convert any floats to ints (caused when NaN included in column)

        # TODO Add total column in HTML

    elif mode_of_inheritance == "x_linked":
        # Filter out "uninformative" from summary
        categorised_snps_by_region = df[
            (df["female_AB_snp_risk_category"] != "uninformative")
            & (df["male_AA_snp_risk_category"] != "uninformative")
            & (df["male_BB_snp_risk_category"] != "uninformative")
        ]
        # Group informative 'low_risk' and 'high_risk' SNPs together per region
        summary_categorised_snps_by_region = categorised_snps_by_region.groupby(
            by=["gene_distance"]
        ).sum(numeric_only=True)
        # Ensure rows are in logical order
        summary_categorised_snps_by_region = summary_categorised_snps_by_region.reindex(
            [
                "1-2MB_from_start",
                "0-1MB_from_start",
                "within_gene",
                "0-1MB_from_end",
                "1-2MB_from_end",
            ]
        )
        # Replace any NaN with 0 and convert any floats to ints (caused when NaN included in column)
        summary_categorised_snps_by_region = summary_categorised_snps_by_region.fillna(
            0
        )
        summary_categorised_snps_by_region[
            "female_AB_snp_count"
        ] = summary_categorised_snps_by_region["female_AB_snp_count"].astype(int)
        summary_categorised_snps_by_region[
            "male_AA_snp_count"
        ] = summary_categorised_snps_by_region["male_AA_snp_count"].astype(int)
        summary_categorised_snps_by_region[
            "male_BB_snp_count"
        ] = summary_categorised_snps_by_region["male_BB_snp_count"].astype(int)
        # Calculate totals
        summary_categorised_snps_by_region.loc[
            "total_snps"
        ] = summary_categorised_snps_by_region.sum(numeric_only=True, axis=0)
        summary_categorised_snps_by_region = (
            summary_categorised_snps_by_region.reset_index()
        )
    return summary_categorised_snps_by_region


def categorise_embryo_alleles(
    df,
    male_partner,
    female_partner,
    embryo_ids,
    embryo_sex,
    mode_of_inheritance,
    consanguineous,
):
    """For each embryo this fuction categorises their SNPs
    Note the usable/informative genotypes for each mode of inheritance are hardcoded into this function.
    These have been defined with the PGD team. Note genotypes are defined as usable based on whether we
    can trace the inheritance from a parent AND if it's shared by an affected/unaffected reference and NOT
    if it matches with the required genotype to perform a SNV analysis of that site, eg In AR condition,
    if performing SNV analysis would be interested in homozygous sites but these sites are not informative
    in this application - only heterozygous sites allow us to determine who an allele was inherited from and
    if it's shared by the reference.
    Args:
        df (dataframe): A results_df dataframe produce
        male_partner (string):
        female_partner (string):
        embryo_ids (list): List of embryo_ids matching the columns names in df
        embryo_sex (list) : List of sexes in the same order as embryo_ids (required for x-linked cases)
        mode_of_inheritance (string): "autosomal_dominant", "autosomal_recessive", "x_linked"
        consanguineous (boolean): Boolean value indicating if the parents are consanguineous
    Returns:
        dataframe: Dataframe with new column for each embryo annotate with a risk_category.
    """

    embryo_sex_lookup = dict(
        zip(embryo_ids, embryo_sex)
    )  # TODO This has been moved to the class initialization

    # Initiate dataframe for results
    embryo_category_df = df[
        [
            "probeset_id",
            "rsID",
            "Position",
            "gene_distance",
            male_partner,
            female_partner,
        ]
        + embryo_ids
    ].copy()
    # For autosomal_resessive also include the "snp_inherited_from" column
    if mode_of_inheritance == "autosomal_recessive":
        embryo_category_df = embryo_category_df.join(df["snp_inherited_from"].copy())

    for embryo in embryo_ids:
        embryo_risk_col = f"{embryo}_risk_category"
        # categorise risk category for each SNP in the embryo
        if mode_of_inheritance == "autosomal_dominant":
            conditions = [
                (df["snp_risk_category"] == "high_risk") & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "low_risk") & (df[embryo] == "AB"),
                (df["snp_risk_category"] != "uninformative") & (df[embryo] == "NoCall"),
            ]
            values = [
                "high_risk",
                "low_risk",
                "NoCall",
            ]
            embryo_category_df[embryo_risk_col] = np.select(
                conditions, values, default="uninformative"
            )
        elif mode_of_inheritance == "autosomal_recessive":
            conditions = [
                # male_partner
                (df["snp_risk_category"] == "high_risk")
                & (df["snp_inherited_from"] == "male_partner")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "high_risk")
                & (df["snp_inherited_from"] == "male_partner")
                & ((df[embryo] == "AA") | (df[embryo] == "BB"))
                & consanguineous,
                (df["snp_risk_category"] == "low_risk")
                & (df["snp_inherited_from"] == "male_partner")
                & (df[embryo] == "AB"),
                # female_partner
                (df["snp_risk_category"] == "high_risk")
                & (df["snp_inherited_from"] == "female_partner")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "high_risk")
                & (df["snp_inherited_from"] == "female_partner")
                & ((df[embryo] == "AA") | (df[embryo] == "BB"))
                & consanguineous,
                (df["snp_risk_category"] == "low_risk")
                & (df["snp_inherited_from"] == "female_partner")
                & (df[embryo] == "AB"),
                # NoCall
                (df["snp_risk_category"] != "uninformative") & (df[embryo] == "NoCall"),
            ]
            values = [
                "high_risk",
                "high_risk",
                "low_risk",
                "high_risk",
                "high_risk",
                "low_risk",
                "NoCall",
            ]

            embryo_category_df[embryo_risk_col] = np.select(
                conditions, values, default="uninformative"
            )

        elif mode_of_inheritance == "x_linked":
            if embryo_sex_lookup[embryo] == "female":
                conditions = [
                    (df["female_AB_snp_risk_category"] == "high_risk")
                    & (df[embryo] == "AB"),
                    (df["female_AB_snp_risk_category"] == "low_risk")
                    & (df[embryo] == "AB"),
                    (df["female_AB_snp_risk_category"] != "uninformative")
                    & (df[embryo] == "NoCall"),
                ]
                values = [
                    "high_risk",
                    "low_risk",
                    "NoCall",
                ]
                embryo_category_df[embryo_risk_col] = np.select(
                    conditions, values, default="uninformative"
                )
            elif embryo_sex_lookup[embryo] == "male":
                conditions = [
                    (df["male_AA_snp_risk_category"] == "high_risk")
                    & (df[embryo] == "AA"),
                    (df["male_AA_snp_risk_category"] == "low_risk")
                    & (df[embryo] == "AA"),
                    (df["male_BB_snp_risk_category"] == "high_risk")
                    & (df[embryo] == "BB"),
                    (df["male_BB_snp_risk_category"] == "low_risk")
                    & (df[embryo] == "BB"),
                    (df["male_AA_snp_risk_category"] != "uninformative")
                    & (df["male_BB_snp_risk_category"] != "uninformative")
                    & (df[embryo] == "NoCall"),
                ]
                values = [
                    "high_risk",
                    "low_risk",
                    "high_risk",
                    "low_risk",
                    "NoCall",
                ]
                embryo_category_df[embryo_risk_col] = np.select(
                    conditions, values, default="uninformative"
                )
            elif embryo_sex_lookup[embryo] == "unknown":
                raise ArgumentInputError(
                    f"'unknown' embryo sex not allowed for x-linked mode of inheritance.  Check that correct mode of inheritance has been entered for {embryo}, or enter correct sex for {embryo}"
                )
            else:
                raise ArgumentInputError(
                    f"Check x_linked code in categorise_embryo_alleles() function and input data for {embryo} - SNP not assigned risk category"
                )

        # Populate embryo_category_df with data regarding miscalls and ADOs
        embryo_category_df[embryo_risk_col] = embryo_category_df.apply(
            lambda row: row[embryo_risk_col]
            if row[embryo_risk_col] != "uninformative"
            else detect_miscall_or_ado(
                row[male_partner], row[female_partner], row[embryo]
            ),
            axis=1,
        )
        # Rename any uninformative calls from "call" to "uninformative" so they are consistent
        # with plotting functions
        embryo_category_df[embryo_risk_col] = embryo_category_df[
            embryo_risk_col
        ].replace("call", "uninformative")

        embryo_category_df[embryo_risk_col] = embryo_category_df[
            embryo_risk_col
        ].astype("category")
        embryo_category_df[embryo_risk_col] = embryo_category_df[
            embryo_risk_col
        ].cat.set_categories(
            ["high_risk", "low_risk", "uninformative", "NoCall", "miscall", "ADO"]
        )
    return embryo_category_df


def annotate_snp_position(df):
    """For a dataframe with a "gene_distance" column this adds a "snp_position" column.  This is useful for summarising
    data in the column,
    Args:
        df (dataframe): A dataframe with "gene_distance" column with category values in the range:
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
    Returns:
        dataframe: Dataframe with new column "snp_position", with the category values "upstream", "within_gene", and "downstream".
    """

    df["snp_position"] = df["gene_distance"].str.replace(
        r"\d-\dMB_from_",
        "",
        regex=True,
    )
    df["snp_position"] = (
        df["snp_position"].replace("start", "upstream").astype("category")
    )
    df["snp_position"] = (
        df["snp_position"].replace("end", "downstream").astype("category")
    )
    df["snp_position"] = df["snp_position"].cat.set_categories(
        [
            "upstream",
            "within_gene",
            "downstream",
        ],
    )
    return df


def summarise_snps_per_embryo_pretty(
    df,
    embryo_ids,
):
    """
    This function groups a results data frame by gene_distance and risk category, and then sums the number of SNPs in each category.
    It then adds a new column to the dataframe, "snp_position", which is either "upstream", "downstream", or "within_gene".
    Where upstream is 0-2MB from the start of the gene (5' direction) and downstream is 0-2MB from the end of the gene in the 3' direction.
    Args:
        df (dataframe): A dataframe with "gene_distance" column with category values in the range:
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
            and risk category columns for each embryo in embryo_ids, with category values in the range:
            "high_risk",
            "low_risk",
            "uninformative",
            "NoCall",
            "miscall",
            "ADO"
        embryo_ids (list): A list of embryo columns in the dataframe to be summarised.
    """
    counter = 0
    for embryo in embryo_ids:
        # Code below acts like a pivot table producing a breakdown of the number in each category of
        # gene_distance and risk category. size() tells the function to count the number of occurrences
        # rather than say sum(). reset_index() converts it from a groupby object with two indexes,
        # gene_distance and risk_category into one index of "gene_distance, risk_category".
        if "snp_inherited_from" in df:
            if counter == 0:
                output_df = (
                    df.groupby(
                        [
                            "gene_distance",
                            "snp_inherited_from",
                            f"{embryo}_risk_category",
                        ],
                    )
                    .size()
                    .reset_index()
                )
                output_df.columns = [
                    "gene_distance",
                    "snp_inherited_from",
                    "risk_category",
                    embryo,
                ]
                counter = 1
            elif counter > 0:
                temp_df = (
                    df.groupby(
                        [
                            "gene_distance",
                            "snp_inherited_from",
                            f"{embryo}_risk_category",
                        ],
                    )
                    .size()
                    .reset_index()
                )
                temp_df.columns = [
                    "gene_distance",
                    "risk_category",
                    "snp_inherited_from",
                    embryo,
                ]
                output_df[embryo] = temp_df[embryo].values
        else:
            if counter == 0:
                output_df = (
                    df.groupby(["gene_distance", f"{embryo}_risk_category"])
                    .size()
                    .reset_index()
                )
                output_df.columns = ["gene_distance", "risk_category", embryo]
                counter = 1
            elif counter > 0:
                temp_df = (
                    df.groupby(["gene_distance", f"{embryo}_risk_category"])
                    .size()
                    .reset_index()
                )
                temp_df.columns = ["gene_distance", "risk_category", embryo]
                output_df[embryo] = temp_df[embryo].values

    # Add new column- 'upstream', 'downstream', or 'within_gene'
    output_df = annotate_snp_position(output_df)

    return output_df


def summarise_embryo_results(df, embryo_ids):
    """
    Summarise embryo results for each embryo in embryo_ids
    """
    summary_embryo_results = pd.DataFrame()
    for embryo in embryo_ids:
        new_column = df[f"{embryo}_risk_category"].value_counts()
        summary_embryo_results = pd.concat([summary_embryo_results, new_column], axis=1)
    summary_embryo_results = (
        summary_embryo_results.fillna("0").astype(int).reset_index()
    )
    # Ensure same ordering of table accross samples
    summary_embryo_results = summary_embryo_results.sort_values(
        [
            "index",
        ],
        ascending=False,
    )
    return summary_embryo_results


def produce_html_table(
    df,
    table_identifier,
    include_index=False,
    include_total=False,
):
    """HTML table for pandas dataframe
    Converts a pandas dataframe into an HTML table ready for inclusion in an HTML report
    Args:
        df (dataframe): A dataframe which requires rendering as HTML for inclusion in the HTML report
        table_identifier (string): Sets id attribute for the table in the HTML
    Returns:
        String: HTML formated table with the provide table_id used to set the HTML table id attribute.
    """

    # styled_df = df.style.highlight_null(null_color='red').hide_columns(['hap1_risk_category','hap2_risk_category'])
    # styled_df = df.style.apply(highlight_risk)
    html_table = df.to_html(
        table_id=table_identifier, index=include_index, classes="display"
    )
    return html_table


def add_embryo_sex_to_column_name(html_string, embryo_ids, embryo_sex):
    """
    Annotated any table with with embryo data with the sex of the embryos
    Args:
        html_string (string): A HTML formated table with embryo ID column names
        embryo_ids (list): A list of embryo IDs
        embryo_sex (list): A list of embryo sexes coressponding to the embryo ids
    Returns:
        String: HTML formated table with the column headings annotated with the embryo sex.
    """
    embryo_sex_lookup = dict(zip(embryo_ids, embryo_sex))
    for embryo_id in embryo_sex_lookup:
        html_string = html_string.replace(
            f"{embryo_id}",
            f"{embryo_id} (Sex:{embryo_sex_lookup[embryo_id].title()})",
        )
    return html_string


# convert header dictionary into html
def dict2html(header_dictionary):
    """
    Converts a dictionary of header into an html table for displaying in the report.
    Flexible way of allowing the user to add any information they want to the report header.
    Args:
        header_dictionary (dict): Dictionary of header information
        For example: {"Analysis Name": "test", "Analysis Date": "2020-01-01"}
    Returns:
        header_html (str): HTML table of header information
    """
    header_html = f'<h2> Analysis Details </h2> <table style="width:100%"><tr>'
    for key in header_dictionary:
        header_html = header_html + f"<td><b>{key}:</b> {header_dictionary[key]}</td>"
    header_html = header_html + f"</tr> </table>"
    return header_html


def main(args):
    # Check config.py file to see which paramters are currently supported.
    # Typically this is used when the script has been validated for some modes of inheritance
    # and we want to ensure that the script is not run for other, unvalidated, modes of inheritance.

    logger.info(f"snp_haplotyper version: called successfully.")

    if (
        config.allow_autosomal_dominant_cases == False
        and args.mode_of_inheritance == "autosomal_dominant"
    ):
        raise InvalidParameterSelectedError(
            "Please check the config.py file to see whether autosomal dominant samples are supported in the current release"
        )
    elif (
        config.allow_autosomal_recessive_cases == False
        and args.mode_of_inheritance == "autosomal_recessive"
    ):
        raise InvalidParameterSelectedError(
            "As per the config.py file autosomal recessive samples are not supported in this release"
        )
    elif (
        config.allow_x_linked_cases == False and args.mode_of_inheritance == "x_linked"
    ):
        raise InvalidParameterSelectedError(
            "As per the config.py file x-linked samples are not supported in this release"
        )

    elif config.allow_consanguineous_cases == False and args.consanguineous == True:
        raise InvalidParameterSelectedError(
            "As per the config.py file consanguineous samples are not supported in this release"
        )
    elif config.allow_trio_only_analysis == False and args.trio_only == True:
        raise InvalidParameterSelectedError(
            "As per the config.py file trio only analysis are not supported in this release"
        )

    # Check input variables are valid
    if (
        (args.mode_of_inheritance == "autosomal_dominant")
        | (args.mode_of_inheritance == "autosomal_recessive")
    ) and (args.chr == "x"):
        raise ArgumentInputError(
            "Chromosome X is not a valid chromosome for autosomal dominant or recessive samples, please check input"
        )
    elif (args.mode_of_inheritance == "x_linked") and (args.chr != "x"):
        raise ArgumentInputError(
            "Chromosome X is the only valid chromosome for x-linked samples, please check input"
        )

    # import haplotype data from text file
    df = pd.read_csv(
        args.input_file,
        delimiter="\t",
    )
    # Remove space from column titles and make lower case
    df = df.rename(
        columns={
            "Probeset ID": "probeset_id",
        }
    )

    number_snps_imported = df.shape[0]

    logger.info(
        f"Number of SNPs imported from SNP Array File = {number_snps_imported}."
    )

    # Import mapping of Affy IDs to dbSNP rs IDs
    mod_path = Path(__file__).parent
    rsid_data_path = (mod_path / "../test_data/AffyID2rsid.txt").resolve()
    affy_2_rs_ids_df = pd.read_csv(rsid_data_path, delimiter="\t", low_memory=False)

    # Assign the correct partner to 'affected' and 'unaffected'
    if args.mode_of_inheritance == "autosomal_dominant":
        if args.male_partner_status == "affected":
            affected_partner = args.male_partner
            affected_partner_sex = "male_partner"
            unaffected_partner = args.female_partner
            unaffected_partner_sex = "female_partner"
        elif args.female_partner_status == "affected":
            affected_partner = args.female_partner
            affected_partner_sex = "female_partner"
            unaffected_partner = args.male_partner
            unaffected_partner_sex = "male_partner"

    # Filter out any rows not in the region of interest #TODO Now marked in imported data
    df = filter_dataframe(
        df, int(args.gene_start), int(args.gene_end), args.flanking_region_size
    )

    # Add column describing how far the SNP is from the gene of interest #TODO Now done in object
    df = annotate_distance_from_gene(
        df, args.chr, int(args.gene_start), int(args.gene_end)
    )

    # Add column of dbSNP rsIDs  #TODO Now done in object
    df = add_rsid_column(df, affy_2_rs_ids_df)

    # Calculate qc metrics before filtering out Nocalls #TODO Now marked in imported data
    if args.trio_only == True:
        qc_df = calculate_qc_metrics(
            df, args.male_partner, args.female_partner, args.reference, None
        )
    else:
        qc_df = calculate_qc_metrics(
            df, args.male_partner, args.female_partner, args.reference, args.embryo_ids
        )

    # Calculate NoCall percentages

    nocall_percentages = calculate_nocall_percentages(
        qc_df
    )  # TODO Now done in QC calculation

    # Filter out any rows where the partners or reference have a NoCall as these cannot be used in the analysis
    filtered_df = filter_out_nocalls(
        df, args.male_partner, args.female_partner, args.reference
    )

    if args.mode_of_inheritance == "autosomal_dominant":
        results_df = autosomal_dominant_analysis(
            filtered_df,
            affected_partner,
            unaffected_partner,
            args.reference,
            args.reference_status,
            args.reference_relationship,
        )
    elif args.mode_of_inheritance == "autosomal_recessive":
        results_df = autosomal_recessive_analysis(
            filtered_df,
            args.male_partner,
            args.female_partner,
            args.reference,
            args.reference_status,
            args.consanguineous,
        )
    elif args.mode_of_inheritance == "x_linked":
        results_df = x_linked_analysis(
            filtered_df,
            args.female_partner,
            args.male_partner,
            args.reference,
        )

    # Informative SNPs
    informative_snps_by_region = snps_by_region(results_df, args.mode_of_inheritance)

    # Get total of informative SNPs
    summary_snps_by_region = summarised_snps_by_region(
        informative_snps_by_region,
        args.mode_of_inheritance,
    )

    # Do not calculate embryo results for pre-cases and trio_only analysis is required
    if args.trio_only == False:
        # Categorise embryo alleles
        embryo_category_df = categorise_embryo_alleles(
            results_df,
            args.male_partner,
            args.female_partner,
            args.embryo_ids,
            args.embryo_sex,
            args.mode_of_inheritance,
            args.consanguineous,
        )

        embryo_count_data_df = summarise_snps_per_embryo_pretty(
            embryo_category_df,
            args.embryo_ids,
        )

        summary_embryo_df = embryo_count_data_df.groupby(by=["risk_category"]).sum(
            numeric_only=True
        )

        # Group by risk category and SNP position (and snp_inherited_from for AR) and sum the counts
        if args.mode_of_inheritance == "autosomal_recessive":
            summary_embryo_by_region_df = embryo_count_data_df.groupby(
                by=["snp_inherited_from", "risk_category", "gene_distance"]
            ).sum(numeric_only=True)
        else:
            summary_embryo_by_region_df = embryo_count_data_df.groupby(
                by=["risk_category", "gene_distance"]
            ).sum(numeric_only=True)
    ##############################################################################

    # Produce report
    results_table_1 = produce_html_table(
        results_df,
        "results_table_1",
    )

    nocall_table = produce_html_table(
        qc_df,
        "nocall_table",
    )

    nocall_percentages_table = produce_html_table(
        nocall_percentages,
        "nocall_percentages_table",
    )

    if args.mode_of_inheritance == "autosomal_dominant":
        summary_snps_table = produce_html_table(
            summary_snps_by_region,
            "summary_snps_table",
        )
    elif args.mode_of_inheritance == "autosomal_recessive":
        temp_df = annotate_snp_position(summary_snps_by_region)
        temp_df = temp_df.groupby(
            by=["snp_inherited_from", "snp_risk_category", "gene_distance"]
        ).sum(numeric_only=True)
        summary_snps_table = produce_html_table(
            temp_df,
            "summary_snps_table",
            True,
        )
    elif args.mode_of_inheritance == "x_linked":
        temp_df = pd.DataFrame().assign(
            gene_distance=summary_snps_by_region["gene_distance"],
            female_embryo_snp_count=summary_snps_by_region["female_AB_snp_count"],
            male_snp_count=summary_snps_by_region["male_AA_snp_count"]
            + summary_snps_by_region["male_BB_snp_count"],
        )
        summary_snps_table = produce_html_table(
            temp_df,
            "summary_snps_table",
        )
    else:
        pass  # TODO raise exception

    # Do not produce plots for embryos if only a trio is being run
    if (
        args.trio_only == False
    ):  # If only a trio is being run do not produce tables/plots for embryos
        summary_embryo_table = produce_html_table(
            summary_embryo_df,
            "summary_embryo_table",
            True,
        )
        # Annotate column names with the sex of the embryo
        summary_embryo_table = add_embryo_sex_to_column_name(
            summary_embryo_table, args.embryo_ids, args.embryo_sex
        )

        # Filter out any rows for NoCall, MisCall, ADO, or uninformative SNPs so as not to clutter the report tables with unnecessary detail as per user feedback
        concise_embryo_df = summary_embryo_by_region_df[
            np.in1d(
                summary_embryo_by_region_df.index.get_level_values("risk_category"),
                ["high_risk", "low_risk"],
            )
        ]

        if args.mode_of_inheritance == "autosomal_recessive":
            concise_embryo_df = concise_embryo_df[
                np.in1d(
                    concise_embryo_df.index.get_level_values("snp_inherited_from"),
                    [
                        "male_partner",
                        "female_partner",
                    ],  # Filters out uninformative SNPs which may have been allocated to ADO
                )
            ]
        summary_embryo_by_region_table = produce_html_table(
            concise_embryo_df,
            "summary_embryo_by_region_table",
            True,
        )
        # Annotate column names with the sex of the embryo
        summary_embryo_by_region_table = add_embryo_sex_to_column_name(
            summary_embryo_by_region_table, args.embryo_ids, args.embryo_sex
        )

        html_list_of_dynamic_plots, html_list_of_static_plots = plot_results(
            embryo_category_df,
            args.embryo_ids,
            args.embryo_sex,
            int(args.gene_start),
            int(args.gene_end),
            args.mode_of_inheritance,
            embryo_count_data_df,
            args.flanking_region_size,
        )

        html_text_for_plots = "<br><hr><br>" + "<br><hr><br>".join(
            html_list_of_dynamic_plots
        )
        pdf_text_for_plots = "<br><hr><br>" + "<br><hr><br>".join(
            html_list_of_static_plots
        )

    elif args.trio_only == True:
        html_text_for_plots = ""
        pdf_text_for_plots = ""
        embryo_count_data_df = None

    env = Environment(loader=PackageLoader("snp_haplotype", "templates"))

    if type(args.header_info) is dict:
        header_html = args.header_info
    else:
        header_html = dict2html(header_to_dict(args.header_info))

    if config.released_to_production == True:
        warning_text = ""
    elif config.released_to_production == False:
        warning_text = (
            f"<h3 style='color:red'> This is a pre-release version of the BASHer tool. "
            f"Please contact the BASHer team if you have any questions.</h3>"
        )

    template = env.get_template("report_template.html")
    place_holder_values = {
        "header_html": header_html,
        "mode_of_inheritance": args.mode_of_inheritance,
        "gene_symbol": args.gene_symbol,
        "chromsome": args.chr.upper(),
        "gene_start": f"{int(args.gene_start):,}",  # Format with 1000s comma separator
        "gene_end": f"{int(args.gene_end):,}",  # Format with 1000s comma separator
        "genome_build": config.genome_build,  # Imported from config.py file
        "basher_version": config.basher_version,  # Imported from config.py file
        "input_file": args.input_file.name
        if isinstance(args.input_file, IOBase)
        else args.input_file,  # Check if input file is a file object or a string
        "male_partner": args.male_partner,
        "male_partner_status": args.male_partner_status,
        "female_partner": args.female_partner,
        "female_partner_status": args.female_partner_status,
        "reference": args.reference,
        "reference_status": args.reference_status,
        "reference_relationship": args.reference_relationship,
        "results_table_1": results_table_1,
        "nocall_table": nocall_table,
        "nocall_percentages_table": nocall_percentages_table,
        "report_date": datetime.today().strftime("%Y-%m-%d %H:%M:%S"),
        "summary_snps_table": summary_snps_table,
        "summary_embryo_table": summary_embryo_table if args.trio_only == False else "",
        "summary_embryo_by_region_table": summary_embryo_by_region_table
        if args.trio_only == False
        else "",
        "html_text_for_plots": html_text_for_plots,
        "warning": warning_text,  # Warning text, for example if the tool is not released to production
    }

    for file_type in ["html", "pdf"]:
        if file_type == "html":
            place_holder_values["html_text_for_plots"] = html_text_for_plots
            html_string = template.render(place_holder_values)
        elif file_type == "pdf":
            place_holder_values["html_text_for_plots"] = pdf_text_for_plots
            pdf_string = template.render(place_holder_values)

    return (
        args.mode_of_inheritance,
        args.output_prefix,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    )


# Code when running as a script
if __name__ == "__main__":
    args = parser.parse_args()
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        summary_embryo_df,
        html_string,
        pdf_string,
    ) = main(args)

    # Save HTML report to file in output folder, including timestamp in filename

    timestr = datetime.now().strftime("%Y%m%d-%H%M%S")
    with open(
        os.path.join(args.output_folder, args.output_prefix + "_" + timestr + ".html"),
        "w",
    ) as f:
        f.write(html_string)

    # Convert HTML report to PDF
    pdfkit.from_string(
        pdf_string,
        os.path.join(args.output_folder, args.output_prefix + "_" + timestr + ".pdf"),
    )
