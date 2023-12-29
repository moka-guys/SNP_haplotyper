import logging

import pandas as pd
import plotly as plt
import plotly.express as px
import plotly.graph_objects as go
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)

logger = logging.getLogger("BASHer_logger")


def plot_results(
    df,
    embryo_id,
    embryo_sex,
    gene_start,
    gene_end,
    mode_of_inheritance,
    embryo_count_data_df,
    flanking_region_size,
    consanguineous,
):
    """Plots SNP data

    For AD and XL produces a single plotly plot of the gene + mb flanking region with SNP information and summaries.
    For AR a faceted plot is produced further splitting the info by partner on the SNP inherited from.

    Args:
        df:
        embryo_ids:
        mode_of_inheritance:
    Returns:

    """

    bp_flanking_region_size = int(flanking_region_size * 10**6)

    # Ensure that the gene start and end are integers
    gene_start = int(gene_start)
    gene_end = int(gene_end)

    # Replace "unassigned" in "snp_inherited_from" with "uninformative"
    df["snp_inherited_from"] = df["snp_inherited_from"].replace(
        {"unassigned": "uninformative"}
    )
    if consanguineous:
        df["snp_inherited_from"] = df["snp_inherited_from"].replace(
            {"both_partners": "uninformative"}
        )
    else:
        # Ignore both_partners category for non-consanguineous families
        df["snp_inherited_from"] = pd.Categorical(
            df["snp_inherited_from"],
            categories=["male_partner", "female_partner", "uninformative"],
            ordered=True,
        )

    # Mapping for y-axis values
    y_mapped = df["embryo_risk_category"].map(
        {
            "high_risk": 2,
            "low_risk": -2,
            "NoCall": -1,
            "uninformative": 0,
            "miscall": 1,
            "ADO": 1,
        }
    )

    # Determine if faceting is needed based on the mode of inheritance
    facet_col = (
        "snp_inherited_from"
        if mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE
        else None
    )
    facet_col_wrap = (
        1 if mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE else 0
    )

    # Color mapping
    color_discrete_map = {
        "high_risk": "#e60e0e",
        "low_risk": "#0ee60e",
        "NoCall": "#0818a6",
        "miscall": "#f0690a",
        "uninformative": "#52555e",
        "ADO": "#00ccff",
    }

    # Symbol mapping
    symbol_map = {
        "high_risk": "line-ns-open",
        "low_risk": "line-ns-open",
        "NoCall": "x",
        "miscall": "line-ns-open",
        "ADO": "line-ns-open",
        "uninformative": "line-ns-open",
    }

    # Size mapping
    size_mapped = df["embryo_risk_category"].map(
        {
            "high_risk": 8,
            "low_risk": 8,
            "NoCall": 4,
            "NoCall_in_trio": 4,
            "miscall": 1,
            "ADO": 1,
            "uninformative": 1,
        }
    )

    # Plotting

    fig = px.scatter(
        df,
        x="Position",
        y=y_mapped,
        facet_col=facet_col,
        facet_col_wrap=facet_col_wrap,
        color="embryo_risk_category",
        color_discrete_map=color_discrete_map,
        symbol="embryo_risk_category",
        symbol_map=symbol_map,
        category_orders={
            "embryo_risk_category": [
                "high_risk",
                "low_risk",
                "miscall",
                "ADO",
                "NoCall",
                "NoCall_in_trio",
                "uniformative",
            ],
            "snp_inherited_from": ["male_partner", "uninformative", "female_partner"],
        },
        labels={"y": "SNP Category", "Position": "Genomic coordinates"},
        size=size_mapped,
        hover_data={"Position": ":.0f", "probeset_id": True, "rsID": True},
    )

    # Format facet plot labels
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.for_each_annotation(lambda a: a.update(xshift=-680))

    # Highlight gene region
    fig.add_vrect(
        x0=gene_start,
        x1=gene_end,
        annotation_text="Gene",
        annotation_position="outside top",
        fillcolor="blue",
        opacity=0.25,
        line_width=0,
    )
    # add downstream line for flanking region lines
    fig.add_vline(
        x=gene_start - bp_flanking_region_size,
        line_width=3,
        line_dash="dash",
        line_color="green",
        annotation_text=f"{flanking_region_size}Mb from Gene Start",
        annotation_position="left top",
        annotation_textangle=90,
    )
    # add upstream line for flanking region lines
    fig.add_vline(
        x=gene_end + bp_flanking_region_size,
        line_width=3,
        line_dash="dash",
        line_color="green",
        annotation_text=f"{flanking_region_size}Mb from Gene End",
        annotation_position="right top",
        annotation_textangle=90,
    )
    # Set reasonable axis size
    fig.update_xaxes(
        range=[
            gene_start
            - (
                bp_flanking_region_size + 100000
            ),  # nicely place the annotation text within the plot
            gene_end
            + (
                bp_flanking_region_size + 100000
            ),  # nicely place the annotation text within the plot
        ],
        exponentformat="none",
    )

    # Functions to add SNP count annotations to plot
    # add annotations to three different types of plots (and AR has three faceted plots)
    def add_snp_count_annotation(
        facet_row,
        annotation_name_high_risk,
        annotation_name_low_risk,
        upstream_high_sum,
        within_gene_high_sum,
        downstream_high_sum,
        upstream_low_sum,
        within_gene_low_sum,
        downstream_low_sum,
        bp_flanking_region_size,
    ):
        fig.add_trace(
            go.Scatter(
                name=annotation_name_high_risk,
                x=[
                    gene_start
                    - (
                        bp_flanking_region_size - 100000
                    ),  # nicely place the annotation text within the plot
                    (gene_start + gene_end) / 2,
                    gene_end
                    + (
                        bp_flanking_region_size - 100000
                    ),  # nicely place the annotation text within the plot
                ],
                y=[
                    3,
                    3,
                    3,
                ],
                mode="text",
                textfont_color="#e60e0e",
                text=[
                    upstream_high_sum,
                    within_gene_high_sum,
                    downstream_high_sum,
                ],
                textposition="top center",
            ),
            row=facet_row,  # The facet plot to anotate (1=bottom, 2=middle, 3=top)
            col=1,
        )

        fig.add_trace(
            go.Scatter(
                name=annotation_name_low_risk,
                x=[
                    gene_start
                    - (
                        bp_flanking_region_size - 100000
                    ),  # nicely place the annotation text within the plot
                    (gene_start + gene_end) / 2,
                    gene_end
                    + (
                        bp_flanking_region_size - 100000
                    ),  # nicely place the annotation text within the plot
                ],
                y=[-3, -3, -3],
                mode="text",
                textfont_color="#0ee60e",
                text=[
                    upstream_low_sum,
                    within_gene_low_sum,
                    downstream_low_sum,
                ],
                textposition="bottom center",
            ),
            row=facet_row,  # The facet plot to anotate (1=bottom, 2=middle, 3=top)
            col=1,
        )

    if (
        mode_of_inheritance == InheritanceMode.X_LINKED
        or mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT
    ):
        summary_df = embryo_count_data_df.groupby(
            ["embryo_risk_category", "snp_position"]
        ).sum()
        add_snp_count_annotation(
            1,
            "High_risk count",
            "Low_risk count",
            summary_df.at[("high_risk", "upstream"), embryo_id],
            summary_df.at[("high_risk", "within_gene"), embryo_id],
            summary_df.at[("high_risk", "downstream"), embryo_id],
            summary_df.at[("low_risk", "upstream"), embryo_id],
            summary_df.at[("low_risk", "within_gene"), embryo_id],
            summary_df.at[("low_risk", "downstream"), embryo_id],
            bp_flanking_region_size,
        )
    # A faceted 3 plot figure is created for AR so that male and female SNPs can be separated out.
    elif mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE:
        summary_df = embryo_count_data_df.groupby(
            ["embryo_risk_category", "snp_position", "snp_inherited_from"]
        ).sum()
        add_snp_count_annotation(  # Male plot
            3,
            "High_risk count male",
            "Low_risk count male",
            summary_df.at[("high_risk", "upstream", "male_partner"), embryo_id],
            summary_df.at[("high_risk", "within_gene", "male_partner"), embryo_id],
            summary_df.at[("high_risk", "downstream", "male_partner"), embryo_id],
            summary_df.at[("low_risk", "upstream", "male_partner"), embryo_id],
            summary_df.at[("low_risk", "within_gene", "male_partner"), embryo_id],
            summary_df.at[("low_risk", "downstream", "male_partner"), embryo_id],
            bp_flanking_region_size,
        )

        add_snp_count_annotation(  # female plot
            1,
            "High_risk count female",
            "Low_risk count female",
            summary_df.at[("high_risk", "upstream", "female_partner"), embryo_id],
            summary_df.at[("high_risk", "within_gene", "female_partner"), embryo_id],
            summary_df.at[("high_risk", "downstream", "female_partner"), embryo_id],
            summary_df.at[("low_risk", "upstream", "female_partner"), embryo_id],
            summary_df.at[("low_risk", "within_gene", "female_partner"), embryo_id],
            summary_df.at[("low_risk", "downstream", "female_partner"), embryo_id],
            bp_flanking_region_size,
        )

    fig.update_yaxes(range=[-4, 4], showticklabels=False, automargin=True)
    fig.update_layout(
        height=540,
        width=1700,
        title_text=f"Results for {embryo_id} (Embryo: {embryo_sex})",
    )

    return fig
