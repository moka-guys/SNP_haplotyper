import argparse
import json
import pandas as pd
import numpy as np
import sys

# Import mode of inheritance specific code
from autosomal_dominant_logic import autosomal_dominant_analysis
from autosomal_recessive_logic import autosomal_recessive_analysis
from x_linked_logic import x_linked_analysis
from snp_plot import plot_results, summarise_snps_per_embryo


# TODO Copy rsID from hover tap
# TODO Add in gene count to plot
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
    nargs="*",
    type=str,
    choices=[
        "male",
        "female",
        "unknown",
    ],
    help="Embryo sex - must be in same order as embryo_ids (required for X-linked diseases)",
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
    help="Flag to produce JSON output easily parsed by pytest and prevent HTML reports being produced",
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
        (df["Position"] > start) & (df["Position"] < end),
        (df["Position"] < start) & (df["Position"] > start - 1000000),
        (df["Position"] < start - 1000000) & (df["Position"] > start - 2000000),
        (df["Position"] > end) & (df["Position"] < end + 1000000),
        (df["Position"] > end + 1000000) & (df["Position"] < end + 2000000),
    ]

    values = [
        "within_gene",
        "0-1MB_from_start",
        "1-2MB_from_start",
        "0-1MB_from_end",
        "1-2MB_from_end",
    ]

    # TODO check correct chr has been given and boundaries are correct
    df["gene_distance"] = np.select(conditions, values, default="outside_range")
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


def calculate_miscalls(df, male_partner, female_partner, embryo_ids):
    """
    QC identify miscalls
    """
    miscall_df = df[
        [
            "probeset_id",
            "rsID",
            "Position",
        ]
    ].copy()
    for embryo in embryo_ids:
        mis_list = []
        for row in df.iterrows():
            parent_alleles = [row[1][male_partner], row[1][female_partner]]
            if row[1][embryo] == "NoCall":
                mis_list.append("NoCall")
            else:
                match parent_alleles:
                    case ["AA", "AA"]:
                        mis_list.append("miscall") if row[1][
                            embryo
                        ] != "AA" else mis_list.append("call")
                    case ["BB", "BB"]:
                        mis_list.append("miscall") if row[1][
                            embryo
                        ] != "BB" else mis_list.append("call")
                    case ["AA", "BB"]:
                        mis_list.append("ADO") if row[1][
                            embryo
                        ] != "AB" else mis_list.append("call")
                    case ["BB", "AA"]:
                        mis_list.append("ADO") if row[1][
                            embryo
                        ] != "AB" else mis_list.append("call")
                    case ["AA", "AB"]:
                        mis_list.append("ADO") if row[1][embryo] not in [
                            "AA",
                            "AB",
                        ] else mis_list.append("ADO")
                    case ["AB", "AA"]:
                        mis_list.append("ADO") if row[1][embryo] not in [
                            "AA",
                            "AB",
                        ] else mis_list.append("call")
                    case ["BB", "AB"]:
                        mis_list.append("ADO") if row[1][embryo] not in [
                            "BB",
                            "AB",
                        ] else mis_list.append("call")
                    case ["AB", "BB"]:
                        mis_list.append("ADO") if row[1][embryo] not in [
                            "BB",
                            "AB",
                        ] else mis_list.append("call")
                    # TODO chcek that AB AB are all filtered out
                    case ["AB", "AB"]:
                        mis_list.append("placeholder") if row[1][embryo] not in [
                            "AA",
                            "BB",
                            "AB",
                        ] else mis_list.append("call")
        miscall_df[embryo] = mis_list
    return miscall_df


def summarise_miscalls():
    # TODO
    pass


def snps_by_region(df):
    snps_by_region = df.value_counts(["gene_distance", "snp_risk_category"]).to_frame()
    # Extract snps_by_region data from index into columns
    snps_by_region = snps_by_region.reset_index()
    # Rename columns
    snps_by_region.columns = [
        "gene_distance",
        "snp_risk_category",
        "snp_count",
    ]
    return snps_by_region


def summarised_snps_by_region(df):
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
    summary_categorised_snps_by_region = summary_categorised_snps_by_region.fillna(0)
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
    return summary_categorised_snps_by_region


def categorise_embryo_alleles(df, miscall_df, embryo_ids, mode_of_inheritance):
    embryo_category_df = df[
        ["probeset_id", "rsID", "Position", "snp_inherited_from", "gene_distance"]
    ].copy()
    for embryo in embryo_ids:
        # Initiate empty database for results
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
            embryo_category_df[f"{embryo}_risk_category"] = np.select(
                conditions, values, default="uninformative"
            )
        elif mode_of_inheritance == "autosomal_recessive":
            conditions = [
                (df["snp_risk_category"] == "male_partner_low_risk")
                & (df[embryo] == "AA"),
                (df["snp_risk_category"] == "male_partner_low_risk")
                & (df[embryo] == "BB"),
                (df["snp_risk_category"] == "male_partner_high_risk")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "male_partner_high_risk")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "male_partner_low_and_high_risk")
                & ((df[embryo] == "AA") | (df[embryo] == "BB")),
                (df["snp_risk_category"] == "female_partner_low_risk")
                & (df[embryo] == "AA"),
                (df["snp_risk_category"] == "female_partner_low_risk")
                & (df[embryo] == "BB"),
                (df["snp_risk_category"] == "female_partner_high_risk")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "female_partner_high_risk")
                & (df[embryo] == "AB"),
                (df["snp_risk_category"] == "female_partner_low_and_high_risk")
                & ((df[embryo] == "AA") | (df[embryo] == "BB")),
                (df["snp_risk_category"] != "uninformative") & (df[embryo] == "NoCall"),
            ]
            values = [
                "male_partner_low_risk",
                "male_partner_low_risk",
                "male_partner_high_risk",
                "male_partner_high_risk",
                "male_partner_low_and_high_risk",
                "female_partner_low_risk",
                "female_partner_low_risk",
                "female_partner_high_risk",
                "female_partner_high_risk",
                "female_partner_low_and_high_risk",
                "uninformative",
            ]
            embryo_category_df[f"{embryo}_risk_category"] = np.select(
                conditions, values, default="uninformative"
            )
        elif mode_of_inheritance == "x_linked":
            pass

    # Populate embryo_category_df with data from miscall_df
    for embryo in embryo_ids:
        embryo_category_df.loc[
            miscall_df[embryo] == "miscall", f"{embryo}_risk_category"
        ] = "miscall"
        embryo_category_df.loc[
            miscall_df[embryo] == "ADO", f"{embryo}_risk_category"
        ] = "ADO"
    return embryo_category_df


def summarise_embryo_results(df, embryo_ids):
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


def highlight_risk(category):
    criteria = ["low_risk_allele", "high_risk_allele"]
    return [
        "background-color: red" if i == "high_risk_allele" else "" for i in category
    ]


def produce_html_table(
    df,
    table_identifier,
):
    """HTML table for pandas dataframe

    Converts a pandas dataframe into an HTML table ready for inclusion in an HTML report

    Args:
        df (dataframe): A dataframe which requires rendering as HTML for inclusion in the HTML report
        table_identifier (): A

    Returns:
        String: HTML formated table with the provide table_id used to set the HTML table id attribute.

    """

    # styled_df = df.style.highlight_null(null_color='red').hide_columns(['hap1_risk_category','hap2_risk_category'])
    # styled_df = df.style.apply(highlight_risk)
    html_table = df.to_html(table_id=table_identifier, index=False, classes="display")
    return html_table


# TODO fill out functions below

# def autosomal_recessive_analysis(carrier_female, carrier_male, reference):
#     pass

# def x_chromosome_linked_analysis(carrier_female, unaffected_male_partner, reference):
#     pass


def main(args=None):  # default argument allows pytest to override argparse for testing
    if args is None:
        args = parser.parse_args()

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
            affected_partner_sex,
            unaffected_partner,
            unaffected_partner_sex,
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
            args.reference_relationship,
            args.consanguineous,
        )
    elif args.mode_of_inheritance == "x_linked":
        results_df = x_linked_analysis(
            filtered_df,
            affected_partner,
            unaffected_partner,
            args.reference,
            args.reference_status,
            args.reference_relationship,
        )

    # Informative SNPs
    informative_snps_by_region = snps_by_region(results_df)

    # Get total of informative SNPs
    summary_snps_by_region = summarised_snps_by_region(informative_snps_by_region)

    # Produce summary of miscalled SNPs
    miscall_df = calculate_miscalls(
        results_df, args.male_partner, args.female_partner, args.embryo_ids
    )

    # Categorise embryo alleles
    embryo_category_df = categorise_embryo_alleles(
        results_df,
        miscall_df,
        args.embryo_ids,
        args.mode_of_inheritance,
    )

    # Summarise embryo results
    embryo_snps_summary_df = summarise_snps_per_embryo(
        embryo_category_df, args.embryo_ids
    )

    summary_embryo_df = summarise_embryo_results(embryo_category_df, args.embryo_ids)

    # Produce report
    results_table_1 = produce_html_table(
        results_df,
        "results_table_1",
    )

    summary_snps_table = produce_html_table(
        summary_snps_by_region,
        "summary_snps_table",
    )

    nocall_table = produce_html_table(
        qc_df,
        "nocall_table",
    )

    nocall_percentages_table = produce_html_table(
        nocall_percentages,
        "nocall_percentages_table",
    )

    summary_embryo_table = produce_html_table(
        summary_embryo_df,
        "summary_embryo_table",
    )

    # TODO remove flow control
    if args.mode_of_inheritance == "autosomal_dominant":
        html_list_of_plots = plot_results(
            embryo_category_df,
            embryo_snps_summary_df,
            args.embryo_ids,
            args.gene_start,
            args.gene_end,
        )
    elif args.mode_of_inheritance == "autosomal_recessive":
        html_list_of_plots = plot_results(
            embryo_category_df,
            embryo_snps_summary_df,
            args.embryo_ids,
            args.gene_start,
            args.gene_end,
        )
    elif args.mode_of_inheritance == "x_linked":
        html_list_of_plots = plot_results(
            embryo_category_df,
            embryo_snps_summary_df,
            args.embryo_ids,
            args.gene_start,
            args.gene_end,
        )

    html_text_for_plots = "<br>".join(html_list_of_plots)

    # TODO Place html into separate folder
    html_string = (
        """
    <html>
    <head>
        <meta charset="utf-8" />
        <title></title>
        <script src="https://code.jquery.com/jquery-3.5.1.js"></script>
        <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.11.5/js/dataTables.bootstrap5.min.js"></script>
        <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css">
    </head>

        <body>
            <h1>Placeholder 1</h1>

            <!-- *** Section 1 *** --->
            <h2>Probe Classification Table</h2>
                <tbody><tr>
                <td>Minimum Genomic Coordinate:</td>
                <td><input type="text" id="min" name="min"></td>
            </tr>
            <tr>
                <td>Maximum Genomic Coordinate:</td>
                <td><input type="text" id="max" name="max"></td>
            </tr>
            """
        + results_table_1
        + """
            <h2>NoCalls per Sample</h2>
        """
        + nocall_table
        + """
        <h2>NoCalls Proportion per Sample</h2>
        """
        + nocall_percentages_table
        + """
        <h2>Informative SNPs by region</h2>
        """
        + summary_snps_table
        + """
            <h2>Embryo Alleles</h2>
        """
        + summary_embryo_table
        + """
            <h2>Plot Embryo Results</h2>
        """
        + html_text_for_plots
        + """

        </body>

        <script>
        $(document).ready( function () {
        $('#results_table_1 thead tr')
            .clone(true)
            .addClass('filters')
            .appendTo('#results_table_1 thead');
            var table = $('#results_table_1').DataTable({
            orderCellsTop: true,
            fixedHeader: true,
            initComplete: function () {
                var api = this.api();
    
                // For each column
                api
                    .columns()
                    .eq(0)
                    .each(function (colIdx) {
                        // Set the header cell to contain the input element
                        var cell = $('.filters th').eq(
                            $(api.column(colIdx).header()).index()
                        );
                        var title = $(cell).text();
                        $(cell).html('<input type="text" placeholder="' + title + '" />');
    
                        // On every keypress in this input
                        $(
                            'input',
                            $('.filters th').eq($(api.column(colIdx).header()).index())
                        )
                            .off('keyup change')
                            .on('keyup change', function (e) {
                                e.stopPropagation();
    
                                // Get the search value
                                $(this).attr('title', $(this).val());
                                var regexr = '({search})'; //$(this).parents('th').find('select').val();
    
                                var cursorPosition = this.selectionStart;
                                // Search the column for that value
                                api
                                    .column(colIdx)
                                    .search(
                                        this.value != ''
                                            ? regexr.replace('{search}', '(((' + this.value + ')))')
                                            : '',
                                        this.value != '',
                                        this.value == ''
                                    )
                                    .draw();
    
                                $(this)
                                    .focus()[0]
                                    .setSelectionRange(cursorPosition, cursorPosition);
                            });
                    });
            },
        });
    });
        </script>

        <script>
        $(document).ready( function () {
        $('#summary_snps_table').DataTable({
            "paging":   false,
            "ordering": false,
            "info":     false
        } );
        } );
        </script>

        <script>
        $(document).ready( function () {
        $('#nocall_table').DataTable({
            "paging":   false,
            "ordering": false,
            "info":     false
        } );
        } );
        </script>

        <script>
        $(document).ready( function () {
        $('#nocall_percentages_table').DataTable({
            "paging":   false,
            "ordering": false,
            "info":     false
        } );
        } );
        </script>

        <script>
        $(document).ready( function () {
        $('#summary_embryo_table').DataTable({
            "paging":   false,
            "ordering": false,
            "info":     false
        } );
        } );
        </script>

        <script>
            /* Custom filtering function which will search data in column four between two values */
        $.fn.dataTable.ext.search.push(
            function( settings, data, dataIndex ) {
                var min = parseInt( $('#min').val(), 10 );
                var max = parseInt( $('#max').val(), 10 );
                var Position = parseFloat( data[3] ) || 0; // use data for the Position column
        
                if ( ( isNaN( min ) && isNaN( max ) ) ||
                    ( isNaN( min ) && Position <= max ) ||
                    ( min <= Position   && isNaN( max ) ) ||
                    ( min <= Position   && Position <= max ) )
                {
                    return true;
                }
                return false;
            }
        );
        
        $(document).ready(function() {
            var table = $('#results_table_1').DataTable();
            
            // Event listener to the two range filtering inputs to redraw on input
            $('#min, #max').keyup( function() {
                table.draw();
            } );
        } );
        </script>

    </html>"""
    )
    if args.testing:
        # Stream machine readable output to stdout for testing purposes
        informative_snp_data = {
            "mode": args.mode_of_inheritance,
            "sample_id": args.output_prefix,
            "num_snps": number_snps_imported,
            "info_snps_upstream_2mb": int(
                informative_snps_by_region[
                    (informative_snps_by_region["gene_distance"].str.endswith("_start"))
                    & (
                        informative_snps_by_region["snp_risk_category"]
                        != "uninformative"
                    )
                ].snp_count.sum()
            ),
            "info_snps_in_gene": int(
                informative_snps_by_region[
                    (informative_snps_by_region["gene_distance"] == "within_gene")
                    & (
                        informative_snps_by_region["snp_risk_category"]
                        != "uninformative"
                    )
                ].snp_count.sum()
            ),
            "info_snps_downstream_2mb": int(
                informative_snps_by_region[
                    (informative_snps_by_region["gene_distance"].str.endswith("_end"))
                    & (
                        informative_snps_by_region["snp_risk_category"]
                        != "uninformative"
                    )
                ].snp_count.sum()
            ),
            "total_info_snps": int(
                informative_snps_by_region[
                    informative_snps_by_region["snp_risk_category"] != "uninformative"
                ].snp_count.sum()
            ),
            # TODO included within gene
            "high_risk_snps_upstream_2mb": int(
                informative_snps_by_region[
                    (
                        informative_snps_by_region["gene_distance"].str.endswith(
                            "_start"
                        )
                        | (informative_snps_by_region["gene_distance"] == "within_gene")
                    )
                    & (informative_snps_by_region["snp_risk_category"] == "high_risk")
                ].snp_count.sum()
            ),
            # TODO included within gene
            "high_risk_snps_downstream_2mb": int(
                informative_snps_by_region[
                    (
                        informative_snps_by_region["gene_distance"].str.endswith("_end")
                        | (informative_snps_by_region["gene_distance"] == "within_gene")
                    )
                    & (informative_snps_by_region["snp_risk_category"] == "high_risk")
                ].snp_count.sum()
            ),
            "low_risk_snps_upstream_2mb": int(
                informative_snps_by_region[
                    (informative_snps_by_region["gene_distance"].str.endswith("_start"))
                    & (informative_snps_by_region["snp_risk_category"] == "low_risk")
                ].snp_count.sum()
            ),
            "low_risk_snps_downstream_2mb": int(
                informative_snps_by_region[
                    (informative_snps_by_region["gene_distance"].str.endswith("_end"))
                    & (informative_snps_by_region["snp_risk_category"] == "low_risk")
                ].snp_count.sum()
            ),
        }
        json.dump(informative_snp_data, sys.stdout, indent=4)

    else:
        # Produce human readable HTML report
        with open(f"{args.output_prefix}.html", "w") as f:
            f.write(html_string)


if __name__ == "__main__":
    main()
