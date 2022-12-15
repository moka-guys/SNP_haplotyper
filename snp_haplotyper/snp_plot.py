import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


def filter_snps_by_partner_sex(
    df,
    partner_sex="all",
):
    # Check that df is not empty
    if partner_sex == "all":
        filtered_df = df
    else:
        filtered_df = df[df["snp_inherited_from"] == partner_sex]

    return filtered_df


def filter_snps_by_catergory(
    df,
    lookup_category,
    embryo,
):
    # Check that df is not empty
    filtered_df = df[df[f"{embryo}_risk_category"] == lookup_category]
    return filtered_df


def filter_snps_by_region(
    df,
    required_region,
):
    # Check that df is not empty
    filtered_df = df[df["gene_distance"].isin(required_region)]
    return filtered_df


def summarise_snps_per_embryo(
    df,
    embryo_ids,
    mode_of_inheritance,
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
    if mode_of_inheritance == "autosomal_dominant" or mode_of_inheritance == "x_linked":
        summary_snps_for_embryos_df = pd.DataFrame(
            columns=[
                "embryo_id",
                "total_high_risk_snps",
                "downstream_2mb_high_risk_snps",
                "within_gene_high_risk_snps",
                "upstream_2mb_high_risk_snps",
                "total_low_risk_snps",
                "downstream_2mb_low_risk_snps",
                "within_gene_low_risk_snps",
                "upstream_2mb_low_risk_snps",
            ]
        )
    elif mode_of_inheritance == "autosomal_recessive":
        # For AR cases create dataframe with "snp_inherited_from" info
        summary_snps_for_embryos_df = pd.DataFrame(
            columns=[
                "embryo_id",
                "total_high_risk_snps",
                "downstream_2mb_high_risk_snps",
                "within_gene_high_risk_snps",
                "upstream_2mb_high_risk_snps",
                "total_low_risk_snps",
                "downstream_2mb_low_risk_snps",
                "within_gene_low_risk_snps",
                "upstream_2mb_low_risk_snps",
                "total_male_high_risk_snps",
                "downstream_2mb_male_high_risk_snps",
                "within_gene_male_high_risk_snps",
                "upstream_2mb_male_high_risk_snps",
                "total_female_high_risk_snps",
                "downstream_2mb_female_high_risk_snps",
                "within_gene_female_high_risk_snps",
                "upstream_2mb_female_high_risk_snps",
                "downstream_2mb_male_low_risk_snps",
                "within_gene_male_low_risk_snps",
                "upstream_2mb_male_low_risk_snps",
                "total_female_low_risk_snps",
                "downstream_2mb_female_low_risk_snps",
                "within_gene_female_low_risk_snps",
                "upstream_2mb_female_low_risk_snps",
            ]
        )

    for embryo in embryo_ids:
        # Filter dataframes required to calculate summary
        total_high_risk_snps_df = filter_snps_by_catergory(
            df,
            "high_risk",
            embryo,
        )
        upstream_2mb_high_risk_snps_df = filter_snps_by_region(
            total_high_risk_snps_df,
            [
                "0-1MB_from_start",
                "1-2MB_from_start",
            ],
        )
        within_gene_high_risk_snps_df = filter_snps_by_region(
            total_high_risk_snps_df,
            [
                "within_gene",
            ],
        )
        downstream_2mb_high_risk_snps_df = filter_snps_by_region(
            total_high_risk_snps_df,
            [
                "0-1MB_from_end",
                "1-2MB_from_end",
            ],
        )
        total_low_risk_snps_df = filter_snps_by_catergory(
            df,
            "low_risk",
            embryo,
        )
        upstream_2mb_low_risk_snps_df = filter_snps_by_region(
            total_low_risk_snps_df,
            [
                "0-1MB_from_start",
                "1-2MB_from_start",
            ],
        )
        within_gene_low_risk_snps_df = filter_snps_by_region(
            total_low_risk_snps_df,
            [
                "within_gene",
            ],
        )
        downstream_2mb_low_risk_snps_df = filter_snps_by_region(
            total_low_risk_snps_df,
            [
                "0-1MB_from_end",
                "1-2MB_from_end",
            ],
        )

        # Count SNPs for each condition
        total_high_risk_snps = len(total_high_risk_snps_df)
        upstream_2mb_high_risk_snps = len(upstream_2mb_high_risk_snps_df)
        within_gene_high_risk_snps = len(within_gene_high_risk_snps_df)
        downstream_2mb_high_risk_snps = len(downstream_2mb_high_risk_snps_df)
        total_low_risk_snps = len(total_low_risk_snps_df)
        upstream_2mb_low_risk_snps = len(upstream_2mb_low_risk_snps_df)
        within_gene_low_risk_snps = len(within_gene_low_risk_snps_df)
        downstream_2mb_low_risk_snps = len(downstream_2mb_low_risk_snps_df)

        # For AR cases calculate values for each partner using "snp_inherited_from" info
        if mode_of_inheritance == "autosomal_recessive":

            total_male_high_risk_snps = len(
                filter_snps_by_partner_sex(total_high_risk_snps_df, "male_partner")
            )
            downstream_2mb_male_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    downstream_2mb_high_risk_snps_df, "male_partner"
                )
            )
            within_gene_male_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    within_gene_high_risk_snps_df, "male_partner"
                )
            )
            upstream_2mb_male_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    upstream_2mb_high_risk_snps_df, "male_partner"
                )
            )
            total_female_high_risk_snps = len(
                filter_snps_by_partner_sex(total_high_risk_snps_df, "female_partner")
            )
            downstream_2mb_female_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    downstream_2mb_high_risk_snps_df, "female_partner"
                )
            )
            within_gene_female_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    within_gene_high_risk_snps_df, "female_partner"
                )
            )
            upstream_2mb_female_high_risk_snps = len(
                filter_snps_by_partner_sex(
                    upstream_2mb_high_risk_snps_df, "female_partner"
                )
            )
            downstream_2mb_male_low_risk_snps = len(
                filter_snps_by_partner_sex(
                    downstream_2mb_low_risk_snps_df, "male_partner"
                )
            )
            within_gene_male_low_risk_snps = len(
                filter_snps_by_partner_sex(within_gene_low_risk_snps_df, "male_partner")
            )
            upstream_2mb_male_low_risk_snps = len(
                filter_snps_by_partner_sex(
                    upstream_2mb_low_risk_snps_df, "male_partner"
                )
            )
            total_female_low_risk_snps = len(
                filter_snps_by_partner_sex(total_low_risk_snps_df, "female_partner")
            )
            downstream_2mb_female_low_risk_snps = len(
                filter_snps_by_partner_sex(
                    downstream_2mb_low_risk_snps_df, "female_partner"
                )
            )
            within_gene_female_low_risk_snps = len(
                filter_snps_by_partner_sex(
                    within_gene_low_risk_snps_df, "female_partner"
                )
            )
            upstream_2mb_female_low_risk_snps = len(
                filter_snps_by_partner_sex(
                    upstream_2mb_low_risk_snps_df, "female_partner"
                )
            )

        # Populate summary dataframe
        row = [
            embryo,
            total_high_risk_snps,
            downstream_2mb_high_risk_snps,
            within_gene_high_risk_snps,
            upstream_2mb_high_risk_snps,
            total_low_risk_snps,
            downstream_2mb_low_risk_snps,
            within_gene_low_risk_snps,
            upstream_2mb_low_risk_snps,
        ]

        # For AR cases include values for each partner as calculated above
        if mode_of_inheritance == "autosomal_recessive":
            row = [
                embryo,
                total_high_risk_snps,
                downstream_2mb_high_risk_snps,
                within_gene_high_risk_snps,
                upstream_2mb_high_risk_snps,
                total_low_risk_snps,
                downstream_2mb_low_risk_snps,
                within_gene_low_risk_snps,
                upstream_2mb_low_risk_snps,
                total_male_high_risk_snps,
                downstream_2mb_male_high_risk_snps,
                within_gene_male_high_risk_snps,
                upstream_2mb_male_high_risk_snps,
                total_female_high_risk_snps,
                downstream_2mb_female_high_risk_snps,
                within_gene_female_high_risk_snps,
                upstream_2mb_female_high_risk_snps,
                downstream_2mb_male_low_risk_snps,
                within_gene_male_low_risk_snps,
                upstream_2mb_male_low_risk_snps,
                total_female_low_risk_snps,
                downstream_2mb_female_low_risk_snps,
                within_gene_female_low_risk_snps,
                upstream_2mb_female_low_risk_snps,
            ]

        summary_snps_for_embryos_df = pd.concat(
            [
                summary_snps_for_embryos_df,
                pd.DataFrame([dict(zip(summary_snps_for_embryos_df.columns, row))]),
            ],
            ignore_index=True,
            axis=0,
        )

    return summary_snps_for_embryos_df


def plot_results(
    df,
    summary_df,
    embryo_ids,
    embryo_sex,
    gene_start,
    gene_end,
    mode_of_inheritance,
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

    plots_as_html = []

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
            x=gene_start - 2000000,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene Start",
            annotation_position="left top",
            annotation_textangle=90,
        )
        # add upstream line for 2mb flanking region lines
        fig.add_vline(
            x=gene_end + 2000000,
            line_width=3,
            line_dash="dash",
            line_color="green",
            annotation_text="2Mb from Gene End",
            annotation_position="right top",
            annotation_textangle=90,
        )
        # Set reasonable axis size
        fig.update_xaxes(
            range=[gene_start - 2100000, gene_end + 2100000], exponentformat="none"
        )

        # Functions to add SNP count annotations to plot
        # TODO will refactor to reduce code duplication as calculated in main script

        # Get counts for each risk category per region
        def get_counts(summary_df, embryo, multi_category):
            count = summary_df.loc[
                summary_df["embryo_id"] == embryo,
                multi_category,
            ].item()
            return count

        # Add SNP count annotations to plot - code gets convoluted here as we have to
        # add annotations to three different types of plots (and AR has three faceted plots)
        def add_snp_count_annotation(
            facet_row,
            embryo,
            summary_df,
            annotation_name_high_risk,
            annotation_name_low_risk,
            upstream_high_sum,
            upstream_low_sum,
            within_gene_high_sum,
            within_gene_low_sum,
            downstream_high_sum,
            downstream_low_sum,
        ):
            fig.add_trace(
                go.Scatter(
                    name=annotation_name_high_risk,
                    x=[
                        gene_start - 1900000,
                        (gene_start + gene_end) / 2,
                        gene_end + 1900000,
                    ],
                    y=[
                        3,
                        3,
                        3,
                    ],
                    mode="text",
                    textfont_color="#e60e0e",
                    text=[
                        get_counts(summary_df, embryo, upstream_high_sum),
                        get_counts(summary_df, embryo, within_gene_high_sum),
                        get_counts(summary_df, embryo, downstream_high_sum),
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
                        gene_start - 1900000,
                        (gene_start + gene_end) / 2,
                        gene_end + 1900000,
                    ],
                    y=[-3, -3, -3],
                    mode="text",
                    textfont_color="#0ee60e",
                    text=[
                        get_counts(summary_df, embryo, upstream_low_sum),
                        get_counts(summary_df, embryo, within_gene_low_sum),
                        get_counts(summary_df, embryo, downstream_low_sum),
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
            add_snp_count_annotation(
                1,
                embryo,
                summary_df,
                "High_risk count",
                "Low_risk count",
                "upstream_2mb_high_risk_snps",
                "upstream_2mb_low_risk_snps",
                "within_gene_high_risk_snps",
                "within_gene_low_risk_snps",
                "downstream_2mb_high_risk_snps",
                "downstream_2mb_low_risk_snps",
            )
        # A faceted 3 plot figure is created for AR so that male and female SNPs can be separated out.
        elif mode_of_inheritance == "autosomal_recessive":
            add_snp_count_annotation(  # Male plot
                3,
                embryo,
                summary_df,
                "High_risk count male",
                "Low_risk count male",
                "upstream_2mb_male_high_risk_snps",
                "upstream_2mb_male_low_risk_snps",
                "within_gene_male_high_risk_snps",
                "within_gene_male_low_risk_snps",
                "downstream_2mb_male_high_risk_snps",
                "downstream_2mb_male_low_risk_snps",
            )

            add_snp_count_annotation(  # female plot
                1,
                embryo,
                summary_df,
                "High_risk count female",
                "Low_risk count female",
                "upstream_2mb_female_high_risk_snps",
                "upstream_2mb_female_low_risk_snps",
                "within_gene_female_high_risk_snps",
                "within_gene_female_low_risk_snps",
                "downstream_2mb_female_high_risk_snps",
                "downstream_2mb_female_low_risk_snps",
            )

        fig.update_yaxes(range=[-4, 4], showticklabels=False)
        fig.update_layout(
            height=540,
            width=1800,
            title_text=f"Results for {embryo} (Embryo Sex: {embryo_dict[embryo]})",
        )
        # Convert plot to HTML and add to list of plots for export and insertion in HTML template
        plots_as_html.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))
    return plots_as_html
