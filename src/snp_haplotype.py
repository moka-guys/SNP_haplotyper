import argparse
import pandas as pd
import numpy as np
import plotly.express as px

# TODO Add split of calls by region to plot & table
# TODO Add support for ADOs in plots
# TODO Check telomeri/centromeric genes work with 2mb window (FHSD1 - D4Z4 repeat, PKD1)
# TODO Add support for no embryos (just TRIOs being run to check if enough informative SNPs)
# TODO Autosomal Recessive

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


# Import mapping of Affy IDs to dbSNP rs IDs
affy_2_rs_ids_df = pd.read_csv(
    "test_data/AffyID2rsid.txt",
    delimiter="\t",
)


def add_rsid_column(df, affy_2_rs_ids_df):
    # Add column with dbSNP rsIDs
    df = pd.merge(
        df, affy_2_rs_ids_df[["probeset_id", "rsID"]], on="probeset_id", how="left"
    )
    # Rearrange columns so that rsID is next to Affy Id
    df.insert(1, "rsID", df.pop("rsID"))
    return df


def annotate_distance_from_gene(df, chr, start, end):
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
    filtered_df = df[
        (df[male_partner] != "NoCall")
        & (df[female_partner] != "NoCall")
        & (df[reference] != "NoCall")
    ]
    # TODO add logger - how many NoCalls filtered
    return filtered_df


number_snps_imported = df.shape[0]


def autosomal_dominant_analysis(
    df,
    affected_partner,
    unaffected_partner,
    reference,
    reference_status,
):

    if args.reference_relationship in [
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
    elif args.reference_relationship in [
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


def calculate_qc_metrics(df, female_partner, male_partner, reference, embryo_ids):
    """
    Calculate QC metrics such as number of NoCalls per sample (measure of DNA quality)
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
        nocall_percentage = trimmed_nocall / trimmed_df.sum(axis=0) * 100
        # Add descriptive column
        nocall_percentage.insert(0, "call_type", "NoCall(%)")
        # TODO add formatting for %
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
    # Extract data from index into columns
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


def categorise_embryo_alleles(df, miscall_df, embryo_ids):
    embryo_category_df = df[
        [
            "probeset_id",
            "rsID",
            "Position",
        ]
    ].copy()
    for embryo in embryo_ids:
        # Initiate empty database for results
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
    # Populate embryo_category_df with data from miscall_df
    for embryo in embryo_ids:
        embryo_category_df.loc[
            miscall_df[embryo] == "miscall", f"{embryo}_risk_category"
        ] = "miscall"
    return embryo_category_df


def summarise_embryo_results(df, embryo_ids):
    summary_embryo_results = pd.DataFrame()
    for embryo in embryo_ids:
        new_column = df[f"{embryo}_risk_category"].value_counts()
        summary_embryo_results = pd.concat([summary_embryo_results, new_column], axis=1)
    summary_embryo_results = (
        summary_embryo_results.fillna("0").astype(int).reset_index()
    )
    # Ensure same ordering of table accross samples TODO check it doesnt break if index not present
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
    # styled_df = df.style.highlight_null(null_color='red').hide_columns(['hap1_risk_category','hap2_risk_category'])
    # styled_df = df.style.apply(highlight_risk)
    html_table = df.to_html(table_id=table_identifier, index=False, classes="display")
    return html_table


def plot_results(
    df,
    embryo_ids,
    gene_start,
    gene_end,
):
    plots_as_html = []
    for embryo in embryo_ids:
        fig = px.scatter(
            df,
            x="Position",
            y=df[f"{embryo}_risk_category"]
            .str.replace("high_risk", "2")
            .str.replace("low_risk", "-2")
            .str.replace("NoCall", "-1")
            .str.replace("uninformative", "0")
            .str.replace("miscall", "1")
            .astype(int),
            color=f"{embryo}_risk_category",
            color_discrete_map={
                "high_risk": "#e60e0e",
                "low_risk": "#0ee60e",
                "NoCall": "#0818a6",
                "miscall": "#f0690a",
                "uninformative": "#52555e",
            },
            symbol=f"{embryo}_risk_category",
            symbol_map={
                "high_risk": "line-ns-open",
                "low_risk": "line-ns-open",
                "NoCall": "x",
                "miscall": "diamond-tall-open",
                "uninformative": "line-ns-open",
            },
            category_orders={
                f"{embryo}_risk_category": [
                    "high_risk",
                    "low_risk",
                    "miscall",
                    "NoCall",
                    "uniformative",
                ],
            },
            labels={
                "y": "SNP Category",
                "Position": "Genomic coordinates",
            },
            size=df[f"{embryo}_risk_category"]
            .str.replace("high_risk", "8")
            .str.replace("low_risk", "8")
            .str.replace("NoCall", "8")
            .str.replace("uninformative", "1")
            .str.replace("miscall", "8")
            .astype(int),
            hover_data=[
                "Position",
                "probeset_id",
                "rsID",
            ],
        )
        fig.add_vrect(
            x0=gene_start,
            x1=gene_end,
            annotation_text="Gene",
            annotation_position="outside top",
        )
        fig.add_vline(
            x=gene_start - 2000000,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene Start",
            annotation_position="left top",
            annotation_textangle=90,
        )
        fig.add_vline(
            x=gene_end + 2000000,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene End",
            annotation_position="right top",
            annotation_textangle=90,
        )
        fig.update_xaxes(
            range=[gene_start - 2000000, gene_end + 2000000], exponentformat="none"
        )
        fig.update_yaxes(range=[-2.2, 2.2], showticklabels=False)
        fig.update_layout(height=250, width=1800, title_text=f"Results for {embryo}")

        plots_as_html.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))
    return plots_as_html


# TODO fill out functions below

# def autosomal_recessive_analysis(carrier_female, carrier_male, reference):
#     pass

# def x_chromosome_linked_analysis(carrier_female, unaffected_male_partner, reference):
#     pass

# Assign the correct partner to 'affected' and 'unaffected'
if args.male_partner_status == "affected":
    affected_partner = args.male_partner
    unaffected_partner = args.female_partner
elif args.female_partner_status == "affected":
    affected_partner = args.female_partner
    unaffected_partner = args.male_partner

# Add column describing how far the SNP is from the gene of interest
df = annotate_distance_from_gene(df, args.chr, args.gene_start, args.gene_end)

# Add column of dbSNP rsIDs
df = add_rsid_column(df, affy_2_rs_ids_df)

# Calculate qc metrics before filtering out Nocalls
qc_df = calculate_qc_metrics(
    df, args.female_partner, args.male_partner, args.reference, args.embryo_ids
)

# Calculate NoCall percentages

nocall_percentages = calculate_nocall_percentages(qc_df)

# Filter out any rows where the partners or reference have a NoCall as these cannot be used in the analysis
filtered_df = filter_out_nocalls(
    df, args.male_partner, args.female_partner, args.reference
)

# TODO add if statement for mode of inheritance
results_df = autosomal_dominant_analysis(
    filtered_df,
    affected_partner,
    unaffected_partner,
    args.reference,
    args.reference_status,
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
)


# Summarise embryo results
summary_embryo_df = summarise_embryo_results(embryo_category_df, args.embryo_ids)


# Filter out any rows which are noninformative across all embryos (produces less cluttered plots)
"_risk_category"


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

# Filter out uninformative
# embryo_category_df[embryo_category_df[""]]

embryo_columns = list(
    embryo_category_df.columns[
        embryo_category_df.columns.str.endswith("_risk_category")
    ]
)
embryo_category_df[embryo_columns] != "uninformative"
mask = ~(embryo_category_df[embryo_columns] == "uninformative").all(axis=1)
uncluttered_df_for_plotting = embryo_category_df[mask]

mask = ~(embryo_category_df[embryo_columns] == "uninformative").all(axis=1)

uncluttered_df_for_plotting = embryo_category_df

html_list_of_plots = plot_results(
    uncluttered_df_for_plotting,
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
   <link rel="stylesheet" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css">
  </head>

    <body>
        <h1>Placeholder 1</h1>

        <!-- *** Section 1 *** --->
        <h2>Probe Classification Table</h2>
        """
    + results_table_1
    + """
        <h2>NoCalls per Sample</h2>
      """
    + nocall_table
    + """
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

</html>"""
)

with open(f"{args.output_prefix}.html", "w") as f:
    f.write(html_string)
