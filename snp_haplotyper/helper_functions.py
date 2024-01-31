import math
import re
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import plotly as plt
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from jinja2 import BaseLoader, DictLoader, Environment


def custom_order_generator(max_range_mb: int) -> List[str]:
    """
    Return a custom order for gene distance to ensure consistent ordering in the report.
    """
    start_range = [f"{i}-{i+1}MB_from_start" for i in reversed(range(max_range_mb))]
    end_range = [f"{i}-{i+1}MB_from_end" for i in range(max_range_mb)]
    return start_range + ["within_gene"] + end_range


def import_haplotype_data(input_file_path: str) -> pd.DataFrame:
    """
    Imports haplotype data from a given text file and processes the columns.
    """

    # Import haplotype data from text file
    snp_data = pd.read_csv(
        input_file_path,
        delimiter="\t",
        dtype={"Chr": str, "Pos": float, "Probeset ID": str},
    )

    # Rename column
    snp_data = snp_data.rename(columns={"Probeset ID": "probeset_id"})

    # Convert "Chr" to category
    snp_data["Chr"] = snp_data["Chr"].astype("category")

    # Convert columns matching *rhchp$ to category
    rhchp_cols = snp_data.filter(like="rhchp", axis=1).columns
    snp_data[rhchp_cols] = snp_data[rhchp_cols].astype("category")

    return snp_data


def dict2html(header_dictionary: dict) -> str:
    """
    Converts a dictionary of header info into an html table for displaying in the report.
    Flexible way of allowing the user to add any information they want to the report header.
    """

    template_string = """
    <h2> Analysis Details </h2>
    <table style="width:100%">
        <tr>
            {% for key, value in header_dictionary.items() %}
                <td><b>{{ key }}:</b> {{ value }}</td>
            {% endfor %}
        </tr>
    </table>
    """

    loader = DictLoader({"header_template": template_string})
    environment = Environment(loader=loader)
    template = environment.get_template("header_template")
    return template.render(header_dictionary=header_dictionary)


def create_human_readable_heading(
    embryo_sex_lookup: dict,
    male_partner: str,
    female_partner: str,
    reference: str,
) -> dict:
    def process_id(input_id: str, suffix: str) -> str:
        # Strip the .rhchp suffix if present and add the corresponding suffix
        return (
            input_id[:-6] + suffix if input_id.endswith(".rhchp") else input_id + suffix
        )

    # Part 1: Process the embryo_sex_lookup dictionary
    human_readable_embryo_lookup = {}
    for embryo_id, sex in embryo_sex_lookup.items():
        human_readable_id = process_id(
            embryo_id, ""
        )  # No additional suffix for embryo IDs
        human_readable_embryo_lookup[embryo_id] = f"{human_readable_id} {sex}"

    # Part 2: Process individual IDs
    human_readable_trio_lookup = {
        male_partner: process_id(male_partner, " (MP)"),
        female_partner: process_id(female_partner, " (FP)"),
        reference: process_id(reference, " (Ref)"),
    }

    # Combine the two dictionaries
    combined_result = {**human_readable_embryo_lookup, **human_readable_trio_lookup}

    return combined_result


def replace_column_names(
    df: pd.DataFrame, human_readable_headings: dict
) -> pd.DataFrame:
    # Replace the column names
    df.rename(columns=human_readable_headings, inplace=True)

    return df


def produce_html_table(
    summary_df: pd.DataFrame,
    table_identifier: str,
    classes: str = "display",
    include_index: bool = False,
) -> str:
    """HTML table for pandas dataframe
    Converts a pandas dataframe into an HTML table ready for inclusion in an HTML report
    Args:
        df (dataframe): A dataframe which requires rendering as HTML for inclusion in the HTML report
        table_identifier (string): Sets id attribute for the table in the HTML
    Returns:
        String: HTML formated table with the provide table_id used to set the HTML table id attribute.
    """
    html_table = summary_df.to_html(
        table_id=table_identifier, index=include_index, classes=classes
    )
    return html_table


def check_affy_duplicates(rsid_data_path):
    """
    Checks for duplicate Affymetrix probesets in the input data.  Some SNPs in the genome
    may have multiple probesets for the same SNP.  This function will produce of summary of
    SNPs with multiple probesets and the number of probesets per SNP.  This is useful as
    we want to count the number of informative SNPs per region, but we don't want to count
    the same SNP multiple times if it has multiple probesets.

    Args:
        str: The path to the input data.
    Returns:
        pd.DataFrame: The dataframe of duplicate probesets.
    """

    df = pd.read_csv(rsid_data_path, delimiter="\t", low_memory=False)

    # Create a single column with the position data
    df["position"] = (
        df["Chr"].astype(str).str.strip() + "_" + df["Min"].astype(str).str.strip()
    )

    # Using keep=False to mark all duplicates
    duplicates_mask = df["position"].duplicated(keep=False)
    duplicate_rows = df[duplicates_mask]

    # Group by 'position' and count the size of each group
    multiple_probeset_summary = duplicate_rows.groupby("position").size()
    multiple_probeset_summary_df = multiple_probeset_summary.reset_index(name="count")

    # Optionally, sort the result by 'count' in descending order
    sorted_result_df = multiple_probeset_summary_df.sort_values(
        by="count", ascending=False
    )

    lookup_df = df[df["position"].isin(sorted_result_df["position"])]
    lookup_df = lookup_df[["position", "rsID"]]
    lookup_df = lookup_df.drop_duplicates(subset=["position"])
    lookup_df = lookup_df.dropna()
    # Merge sorted_result_df with lookup_df

    sorted_result_df = pd.merge(sorted_result_df, lookup_df, on="position", how="left")

    # Function to reformat the 'position' values
    def reformat_position(pos):
        chromosome, coordinate = pos.split("_")
        return f"chr{chromosome}:{coordinate}-{coordinate}"

    sorted_result_df["position"] = sorted_result_df["position"].apply(reformat_position)

    return sorted_result_df


def generate_pdf_plot(fig) -> str:
    """
    Generate static plot for use in PDF report
    """
    return plt.io.to_image(fig, format="svg").decode("utf-8")


def generate_html_plot(fig, embryo_id) -> str:
    """
    Generate dynamic plot for use in HTML report
    """
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id=embryo_id)


def generate_plots(fig_dict: Dict[str, str], static_plots=False) -> Dict[str, str]:
    """
    Generate multiple dynamic plots for use in HTML report
    """
    plots_as_html = {}

    # Iterate over the keys and values of figures_dict
    for key, value in fig_dict.items():
        # Generate HTML for the figure and store it in the dictionary
        if static_plots:
            plots_as_html[key] = generate_pdf_plot(value)
        else:
            plots_as_html[key] = generate_html_plot(value, key)

    return plots_as_html


def format_plot_html_str(
    plots_as_html: Dict[str, str], add_dropdown_selection=False
) -> str:
    """
    Concatenate multiple dynamic plots into a single HTML string
    """
    # TODO Add dropdown functionality and update this function and docstring
    # Initialize an empty list to hold divs
    div_list = []

    if add_dropdown_selection:
        # Initialize an empty list to hold dropdown options
        dropdown_options = []

        # Iterate through each HTML plot and wrap it in a div
        for plot_id, plot_html in plots_as_html.items():
            # Create a div with a unique id using the plot_id
            div = f'<div id="plot-{plot_id}" class="plot-container" style="display:none;">\n{plot_html}\n</div>'
            div_list.append(div)

            # Create a dropdown option for the plot
            dropdown_options.append(
                f'<option value="plot-{plot_id}">Plot {plot_id}</option>'
            )

        # Combine all the divs into a single HTML string
        all_divs = "\n".join(div_list)

        # Create a dropdown HTML element
        dropdown = (
            """<select id="plot-dropdown">
            <option value="" disabled selected>Select a Plot</option>
            """
            + "\n".join(dropdown_options)
            + """
        </select>"""
        )

        # JavaScript code to handle dropdown selection and div visibility
        script = """
<script>
    document.addEventListener('DOMContentLoaded', function () {
        var dropdown = document.getElementById("plot-dropdown");

        if (dropdown && dropdown.options.length > 1) {
            // Set the first real option as the default
            dropdown.selectedIndex = 1;
            var defaultPlotId = dropdown.options[1].value;
            var defaultPlot = document.getElementById(defaultPlotId);
            if (defaultPlot) {
                defaultPlot.style.display = "block"; // Show the default plot
            }

            dropdown.addEventListener("change", function() {
                var plots = document.getElementsByClassName("plot-container");
                for (var i = 0; i < plots.length; i++) {
                    plots[i].style.display = 'none'; // Hide all plots
                }

                var selectedValue = dropdown.options[dropdown.selectedIndex].value;
                var selectedPlot = document.getElementById(selectedValue);
                if (selectedPlot) {
                    selectedPlot.style.display = "block"; // Show the selected plot
                }
            });
        }
    });
</script>
        """

        # Combine dropdown, divs, and script
        final_html = f"{dropdown}\n{all_divs}\n{script}"
    else:
        # Iterate through each HTML plot and wrap it in a div
        for plot_id, plot_html in plots_as_html.items():
            # Create a div with a unique id using the plot_id
            # div = f'<div id="plot-{plot_id}" class="plot-container" style="display:none;">\n{plot_html}\n</div>'
            div = plot_html
            div_list.append(div)

        # Combine all the divs into a single HTML string
        all_divs = "\n".join(div_list)
        final_html = all_divs

    return final_html


def set_risk_category_dtype(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sets columns starting with "snp_risk_category" in a DataFrame to ordered Categorical dtype.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.

    Returns:
    - pd.DataFrame: The modified DataFrame.
    """

    # Define the categories for the categorical dtype
    categories = [
        "high_risk",
        "high_or_low_risk",
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    # Iterate through each column and check if it starts with "snp_risk_category"
    for col in df.columns:
        if col.startswith("snp_risk_category"):
            df[col] = pd.Categorical(df[col], categories=categories, ordered=True)

    return df


def set_inherited_from_category_dtype(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sets column "snp_inherited_from" in a DataFrame to ordered Categorical dtype.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.

    Returns:
    - pd.DataFrame: The modified DataFrame.
    """

    # Define the categories for the categorical dtype
    categories = [
        "male_partner",
        "both_partners",
        "female_partner",
        "unassigned",
    ]

    df["snp_inherited_from"] = pd.Categorical(
        df["snp_inherited_from"], categories=categories, ordered=True
    )

    return df


def get_clean_filename(path):
    # Split the path on either \ or /
    parts = re.split(r"[\\/]", path)
    # Return the last element
    return parts[-1]
