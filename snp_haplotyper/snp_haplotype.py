import argparse
from jinja2 import Environment, PackageLoader
import json
import pandas as pd
from pathlib import Path
import numpy as np
import re
import sys

# Import mode of inheritance specific code
from autosomal_dominant_logic import autosomal_dominant_analysis
from autosomal_recessive_logic import autosomal_recessive_analysis
from config import (
    genome_build,
    allow_autosomal_dominant_cases,
    allow_autosomal_recessive_cases,
    allow_x_linked_cases,
    allow_cosanguineous_cases,
)
from x_linked_logic import x_linked_analysis
from snp_plot import plot_results, summarise_snps_per_embryo
from stream_output import (
    stream_autosomal_dominant_output,
    stream_autosomal_recessive_output,
    stream_x_linked_output,
)

from exceptions import ArgumentInputError, InvalidParameterSelectedError

# TODO Copy rsID from hover tap
# TODO Check telomeric/centromeric genes work with 2mb window (FHSD1 - D4Z4 repeat, PKD1)
# TODO Add support for no embryos (just TRIOs being run to check if enough informative SNPs)
# TODO Add ADO % to table

# Import command line arguments (these can be automatically generated from the sample sheet using sample_sheet_reader.py)
parser = argparse.ArgumentParser(description="SNP Haplotying from SNP Array data")

# File input/output data
parser.add_argument(
    "-i",
    "--input_file",
    type=argparse.FileType("r"),
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
    "--testing",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to produce JSON output easily parsed by pytest and prevent HTML reports being produced",
)

parser.add_argument(
    "--trio_only",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to produce a preliminary report without looking at embryos, must be used if not embryo data is provided.",
)

parser.add_argument(
    "--header_info",
    type=json.loads,
    help='Pass a string in the format of a dictionary to populate the report header. A field will be created for each key provided, for example \'{"PRU":"1234", "Hospital No":"1234", "Biopsy No":"111"}\' will produce 3 fields in the header.',
)


def add_rsid_column(df, affy_2_rs_ids_df):

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


def annotate_distance_from_gene(df, chr, start, end):
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


def filter_out_nocalls(df, male_partner, female_partner, reference):
    """Filters out no calls

    If the male partner, female partner, or reference has "NoCall" for a probeset then this probeset should be filtered out.

    Args:
        df (dataframe): A dataframe with the SNP array data
        male_partner (string):  Column name representing the data for the male partner
        female_partner (string):  Column name representing the data for the female partner
        reference (string):  Column name representing the data for the reference

    Returns:
        dataframe: Original dataframe, df, with any rows where the male partner, female partner or reference has a "NoCall" filtered out

    """
    filtered_df = df[
        (df[male_partner] != "NoCall")
        & (df[female_partner] != "NoCall")
        & (df[reference] != "NoCall")
    ]
    # TODO add logger - how many NoCalls filtered
    return filtered_df


# TODO standardise the order of fp/mp args across functions


def calculate_qc_metrics(df, male_partner, female_partner, reference, embryo_ids):
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
        dataframe: Dataframe summarising the % of NoCalls per sample
    """
    nocall = df[df["call_type"] == "NoCall"]
    if nocall.shape[0] > 0:
        trimmed_nocall = nocall.iloc[:, 1:]  # Trim first column
        trimmed_df = df.iloc[:, 1:]  # Trim first column
        nocall_percentage = trimmed_nocall / trimmed_df.sum(axis=0)
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
    with the haplotype of the parents i.e. Parents AA, AA and an embryo AB. This is due
    to technical error in the measurement.  NOTE: that the miscall could be in any one of the
    trio even though it is recorded under the embryo.

    The definition of ADO (Allele dropout) is used when there is a suspected biological origin for the
    mismatch in haplotypes, due to uniparental inheritance of the allele i.e Parents AA, AB and an embryo AA,
    the B allele has dropped out.  NOTE: that the ADO could have occured in any of the trio even though it is
    recorded under the embryo.  It is expected that the Genomic Scientist will look at the SNP plots and
    use their judgement as to whether allele dropout is observed.

    Args:
        male_partner_haplotype (string): Either "AA", "BB", or "AB" (NoCalls will throw an exception))
        female_partner_haplotype (string): Either "AA", "BB", or "AB" (NoCalls will throw an exception))
        embryo_haplotype (string): Either "AA", "BB", or "AB" (NoCalls will throw an exception))
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
                    "Error! Haplotypes most be AA, BB, or AB"
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
    TODO - check what this function does compared to snps_by_region()
    """
    if mode_of_inheritance == "autosomal_dominant":
        # Filter out "uninformative" from summary
        categorised_snps_by_region = df[df["snp_risk_category"] != "uninformative"]
        # Group informative 'low_risk' and 'high_risk' SNPs together per region
        summary_categorised_snps_by_region = categorised_snps_by_region.groupby(
            by=["gene_distance"]
        ).sum()
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
            .sum()
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
        ).sum()
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

    Takes a dataframe produced

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

    embryo_sex_lookup = dict(zip(embryo_ids, embryo_sex))

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
                pass  # TODO raise exception
            else:
                pass  # TODO raise exception

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
            ["uninformative", "NoCall", "miscall", "ADO", "high_risk", "low_risk"]
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
    mode_of_inheritance,
):
    counter = 0
    for embryo in embryo_ids:
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
    # TODO Docstring
    """
    summary_embryo_results = pd.DataFrame()
    for embryo in embryo_ids:
        new_column = df[f"{embryo}_risk_category"].value_counts()
        summary_embryo_results = pd.concat([summary_embryo_results, new_column], axis=1)
    summary_embryo_results = (
        summary_embryo_results.fillna("0").astype(int).reset_index()
    )
    # Ensure same ordering of table accross samples TODO check it doesn't break if index not present
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
            f"{embryo_id} (Sex:{embryo_sex_lookup[embryo_id]})",
        )
    return html_string


def main(args=None):  # default argument allows pytest to override argparse for testing
    if args is None:
        args = parser.parse_args()

    # Check config.py file to see whether the script should run for the parameters provided.
    # Typically this is used when the script has been validated for one or more modes of inheritance
    # and we want to ensure that the script is not run for other, unvalidated, modes of inheritance.
    # This allows the software to be ushed in production without the risk of running the script on
    # parameters which have not been validated.

    if (
        allow_autosomal_dominant_cases == False
        and args.mode_of_inheritance == "autosomal_dominant"
    ):
        raise InvalidParameterSelectedError(
            "Please check the config.py file to see whether autosomal dominant samples are supported in the current release"
        )
    elif (
        allow_autosomal_recessive_cases == False
        and args.mode_of_inheritance == "autosomal_recessive"
    ):
        raise InvalidParameterSelectedError(
            "Please check the config.py file to see whether autosomal recessive samples are supported in the current release"
        )
    elif allow_x_linked_cases == False and args.mode_of_inheritance == "x_linked":
        raise InvalidParameterSelectedError(
            "Please check the config.py file to see whether x-linked samples are supported in the current release"
        )

    elif allow_cosanguineous_cases == False and args.consanguieous == True:
        raise InvalidParameterSelectedError(
            "Please check the config.py file to see whether consanguineous samples are supported in the current release"
        )
    else:
        pass

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

    # Import mapping of Affy IDs to dbSNP rs IDs
    affy_2_rs_ids_df = pd.read_csv(
        "test_data/AffyID2rsid.txt", delimiter="\t", low_memory=False
    )

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

    # Add column describing how far the SNP is from the gene of interest
    df = annotate_distance_from_gene(df, args.chr, args.gene_start, args.gene_end)

    # Add column of dbSNP rsIDs
    df = add_rsid_column(df, affy_2_rs_ids_df)

    # Calculate qc metrics before filtering out Nocalls
    qc_df = calculate_qc_metrics(
        df, args.male_partner, args.female_partner, args.reference, args.embryo_ids
    )

    # Calculate NoCall percentages

    nocall_percentages = calculate_nocall_percentages(qc_df)

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

    # Summarise embryo results
    embryo_snps_summary_df = summarise_snps_per_embryo(
        embryo_category_df,
        args.embryo_ids,
        args.mode_of_inheritance,
    )

    # summary_embryo_df = summarise_embryo_results(embryo_category_df, args.embryo_ids)

    ##############################################################################
    embryo_count_data_df = summarise_snps_per_embryo_pretty(
        embryo_category_df, args.embryo_ids, args.mode_of_inheritance
    )

    summary_embryo_df = embryo_count_data_df.groupby(by=["risk_category"]).sum()
    summary_embryo_by_region_df = embryo_count_data_df.groupby(
        by=["risk_category", "gene_distance"]
    ).sum()
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
        ).sum()
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
            + summary_snps_by_region["male_AA_snp_count"],
        )
        summary_snps_table = produce_html_table(
            temp_df,
            "summary_snps_table",
        )
    else:
        pass  # TODO raise exception

    # Initiate list to hold HTML for each plot produced below
    html_list_of_plots = []

    # Produce plot for trios

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

        # summary_embryo_by_region_df = embryo_snps_summary_df.set_index("embryo_id")
        # summary_embryo_by_region_df = summary_embryo_by_region_df.transpose()
        summary_embryo_by_region_table = produce_html_table(
            summary_embryo_by_region_df,
            "summary_embryo_by_region_table",
            True,
        )
        # Annotate column names with the sex of the embryo
        summary_embryo_by_region_table = add_embryo_sex_to_column_name(
            summary_embryo_by_region_table, args.embryo_ids, args.embryo_sex
        )

        if args.mode_of_inheritance == "autosomal_dominant":
            html_list_of_plots = html_list_of_plots + plot_results(
                embryo_category_df,
                embryo_snps_summary_df,
                args.embryo_ids,
                args.embryo_sex,
                args.gene_start,
                args.gene_end,
                args.mode_of_inheritance,
            )
        elif args.mode_of_inheritance == "autosomal_recessive":
            html_list_of_plots = html_list_of_plots + plot_results(
                embryo_category_df,
                embryo_snps_summary_df,
                args.embryo_ids,
                args.embryo_sex,
                args.gene_start,
                args.gene_end,
                args.mode_of_inheritance,
            )
        elif args.mode_of_inheritance == "x_linked":
            html_list_of_plots = html_list_of_plots + plot_results(
                embryo_category_df,
                embryo_snps_summary_df,
                args.embryo_ids,
                args.embryo_sex,
                args.gene_start,
                args.gene_end,
                args.mode_of_inheritance,
            )

    html_text_for_plots = "<br>".join(html_list_of_plots)

    env = Environment(loader=PackageLoader("snp_haplotype", "templates"))

    # convert header dictionary into html
    def dict2html(header_dictionary):
        """ """
        header_html = f'<h2> Analysis Details </h2> <table style="width:100%"><tr>'
        for key in header_dictionary:
            header_html = (
                header_html + f"<td><b>{key}:</b> {header_dictionary[key]}</td>"
            )
        header_html = header_html + f"</tr> </table>"
        return header_html

    header_html = dict2html(args.header_info)

    template = env.get_template("report_template.html")
    place_holder_values = {
        "header_html": header_html,
        "mode_of_inheritance": args.mode_of_inheritance,
        "gene_symbol": args.gene_symbol,
        "chromsome": args.chr.upper(),
        "gene_start": args.gene_start,
        "gene_end": args.gene_end,
        "genome_build": genome_build,  # Imported from config.py file
        "input_file": args.input_file.name,
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
        "summary_snps_table": summary_snps_table,
        "summary_embryo_table": summary_embryo_table,
        "summary_embryo_by_region_table": summary_embryo_by_region_table,
        "html_text_for_plots": html_text_for_plots,
    }

    html_string = template.render(place_holder_values)

    # Stream machine readable JSON output to stdout for testing purposes
    # args.testing = True
    if args.testing:

        export_json_data_as_csv(
            "test_data/embryo_validation_data.json",
            "test_data/embryo_validation_data.csv",
        )
        export_json_data_as_csv(
            "test_data/informative_snp_validation.json",
            "test_data/informative_snp_validation.csv",
        )

        if args.mode_of_inheritance == "autosomal_dominant":
            informative_snp_data, embryo_cat_data = stream_autosomal_dominant_output(
                args.mode_of_inheritance,
                informative_snps_by_region,
                embryo_snps_summary_df,
                number_snps_imported,
                args.output_prefix,
            )
        elif args.mode_of_inheritance == "autosomal_recessive":
            informative_snp_data, embryo_cat_data = stream_autosomal_recessive_output(
                args.mode_of_inheritance,
                informative_snps_by_region,
                embryo_snps_summary_df,
                number_snps_imported,
                args.output_prefix,
            )
        elif args.mode_of_inheritance == "x_linked":
            informative_snp_data, embryo_cat_data = stream_x_linked_output(
                args.mode_of_inheritance,
                informative_snps_by_region,
                embryo_snps_summary_df,
                number_snps_imported,
                args.output_prefix,
            )

        json.dump(
            {
                "informative_snp_data": informative_snp_data,
                "embryo_cat_json": embryo_cat_data,
            },
            sys.stdout,
            indent=4,
        )

    else:
        # Produce human readable HTML report
        with open(f"{args.output_folder}{args.output_prefix}.html", "w") as f:
            f.write(html_string)


if __name__ == "__main__":
    main()
