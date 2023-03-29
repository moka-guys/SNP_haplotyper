import plotly.express as px
import plotly.graph_objects as go
import plotly as plt
import pandas as pd

import logging

logger = logging.getLogger("BASHer_logger")


def plot_results(
    df,
    embryo_ids,
    embryo_sex,
    gene_start,
    gene_end,
    mode_of_inheritance,
    embryo_count_data_df,
    flanking_region_size,
):
    """Plots SNP data

    For AD and XL produces a single plotly plot of the gene + 2mb flanking region with SNP information and summaries.
    For AR a faceted plot is produced further splitting the info by partner theSNP inherited from.

    Args:
        df:
        embryo_ids:
        mode_of_inheritance:
    Returns:

    """
    if flanking_region_size == "2mb":
        flanking_region_size = 2000000
    elif flanking_region_size == "3mb":
        flanking_region_size = 3000000

    # Create lists to store plots as html and pdf
    plots_as_html = []
    plots_as_pdf = []

    # Create lookup dictionary for embryo sex
    embryo_dict = dict(zip(embryo_ids, embryo_sex))

    # preprocess dataframe for input to plotting function
    for embryo in embryo_ids:

        fig = px.scatter(
            df,
            x="Position",
            # Replace risk category with numerical values for plotting
            y=df[f"{embryo}_risk_category"].map(
                {
                    "high_risk": 2,
                    "low_risk": -2,
                    "NoCall": -1,
                    "uninformative": 0,
                    "miscall": 1,
                    "ADO": 1,
                }
            ),
            # If AR mode of inheritance, facet by partner the SNP was inherited from
            # (i.e. produce an additional plot for male_partner & female_partner in addition
            # to the main plot - 3 plots in total)
            facet_col="snp_inherited_from"
            if mode_of_inheritance == "autosomal_recessive"
            else None,
            facet_col_wrap=1 if mode_of_inheritance == "autosomal_recessive" else 0,
            color=f"{embryo}_risk_category",
            #  Set color scheme for risk categories
            color_discrete_map={
                "high_risk": "#e60e0e",
                "low_risk": "#0ee60e",
                "NoCall": "#0818a6",
                "miscall": "#f0690a",
                "uninformative": "#52555e",
                "ADO": "#00ccff",
            },
            symbol=f"{embryo}_risk_category",
            # Set symbol for risk categories
            symbol_map={
                "high_risk": "line-ns-open",
                "low_risk": "line-ns-open",
                "NoCall": "x",
                "miscall": "line-ns-open",
                "ADO": "line-ns-open",
                "uninformative": "line-ns-open",
            },
            # Set order for risk categories for consistently ordered plotting
            category_orders={
                f"{embryo}_risk_category": [
                    "high_risk",
                    "low_risk",
                    "miscall",
                    "ADO",
                    "NoCall",
                    "uniformative",
                ],
                "snp_inherited_from": [
                    "male_partner",
                    "uninformative",
                    "female_partner",
                ],
            },
            labels={
                "y": "SNP Category",
                "Position": "Genomic coordinates",
            },
            # Set size for risk categories
            size=df[f"{embryo}_risk_category"].map(
                {
                    "high_risk": 8,
                    "low_risk": 8,
                    "NoCall": 4,
                    "miscall": 1,
                    "ADO": 1,
                    "uninformative": 1,
                }
            ),
            # Turn on hover data
            hover_data={
                "Position": ":.0f",
                "probeset_id": True,
                "rsID": True,
            },
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
        # add downstream line for 2mb flanking region lines
        fig.add_vline(
            x=gene_start - flanking_region_size,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene Start",
            annotation_position="left top",
            annotation_textangle=90,
        )
        # add upstream line for 2mb flanking region lines
        fig.add_vline(
            x=gene_end + flanking_region_size,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene End",
            annotation_position="right top",
            annotation_textangle=90,
        )
        # Set reasonable axis size
        fig.update_xaxes(
            range=[
                gene_start - (flanking_region_size + 100000),
                gene_end + (flanking_region_size + 100000),
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
            flanking_region_size,
        ):
            fig.add_trace(
                go.Scatter(
                    name=annotation_name_high_risk,
                    x=[
                        gene_start - (flanking_region_size - 100000),
                        (gene_start + gene_end) / 2,
                        gene_end + (flanking_region_size - 100000),
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
                        gene_start - (flanking_region_size - 100000),
                        (gene_start + gene_end) / 2,
                        gene_end + (flanking_region_size - 100000),
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
            mode_of_inheritance == "x_linked"
            or mode_of_inheritance == "autosomal_dominant"
        ):
            summary_df = embryo_count_data_df.groupby(
                ["risk_category", "snp_position"]
            ).sum()
            add_snp_count_annotation(
                1,
                "High_risk count",
                "Low_risk count",
                summary_df.at[("high_risk", "upstream"), embryo],
                summary_df.at[("high_risk", "within_gene"), embryo],
                summary_df.at[("high_risk", "downstream"), embryo],
                summary_df.at[("low_risk", "upstream"), embryo],
                summary_df.at[("low_risk", "within_gene"), embryo],
                summary_df.at[("low_risk", "downstream"), embryo],
                flanking_region_size,
            )
        # A faceted 3 plot figure is created for AR so that male and female SNPs can be separated out.
        elif mode_of_inheritance == "autosomal_recessive":
            summary_df = embryo_count_data_df.groupby(
                ["risk_category", "snp_position", "snp_inherited_from"]
            ).sum()
            add_snp_count_annotation(  # Male plot
                3,
                "High_risk count male",
                "Low_risk count male",
                summary_df.at[("high_risk", "upstream", "male_partner"), embryo],
                summary_df.at[("high_risk", "within_gene", "male_partner"), embryo],
                summary_df.at[("high_risk", "downstream", "male_partner"), embryo],
                summary_df.at[("low_risk", "upstream", "male_partner"), embryo],
                summary_df.at[("low_risk", "within_gene", "male_partner"), embryo],
                summary_df.at[("low_risk", "downstream", "male_partner"), embryo],
                flanking_region_size,
            )

            add_snp_count_annotation(  # female plot
                1,
                "High_risk count female",
                "Low_risk count female",
                summary_df.at[("high_risk", "upstream", "female_partner"), embryo],
                summary_df.at[("high_risk", "within_gene", "female_partner"), embryo],
                summary_df.at[("high_risk", "downstream", "female_partner"), embryo],
                summary_df.at[("low_risk", "upstream", "female_partner"), embryo],
                summary_df.at[("low_risk", "within_gene", "female_partner"), embryo],
                summary_df.at[("low_risk", "downstream", "female_partner"), embryo],
                flanking_region_size,
            )

        fig.update_yaxes(range=[-4, 4], showticklabels=False)
        fig.update_layout(
            height=540,
            width=1700,
            title_text=f"Results for {embryo} (Embryo Sex: {embryo_dict[embryo]})",
        )

        plots_as_pdf.append(plt.io.to_image(fig, format="svg").decode("utf-8"))
        plots_as_html.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))

    return plots_as_html, plots_as_pdf
