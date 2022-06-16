import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


def filter_snps_by_partner_sex(
    df,
    partner_sex="all",
):
    # Check that df is not empty
    if partner_sex != "all":
        filtered_df = df
    else:
        filtered_df = df[df["snp_risk_category"] == partner_sex]

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


# if (
#     lookup_category
#     in (df[df["Position"] >= gene_end][f"{embryo}_risk_category"]).unique()
# ):
#     snp_count = df.loc[
#         (df["Position"] >= gene_end),
#         [
#             f"{embryo}_risk_category",
#         ],
#     ].value_counts()[lookup_category]
# else:
#     snp_count = 0

# return snp_count


def summarise_snps_per_embryo(
    df,
    embryo_ids,
):
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

        total_male_high_risk_snps = len(
            filter_snps_by_partner_sex(total_high_risk_snps_df, "male")
        )
        downstream_2mb_male_high_risk_snps = len(
            filter_snps_by_partner_sex(downstream_2mb_high_risk_snps_df, "male")
        )
        within_gene_male_high_risk_snps = len(
            filter_snps_by_partner_sex(within_gene_high_risk_snps_df, "male")
        )
        upstream_2mb_male_high_risk_snps = len(
            filter_snps_by_partner_sex(upstream_2mb_high_risk_snps_df, "male")
        )
        total_female_high_risk_snps = len(
            filter_snps_by_partner_sex(total_high_risk_snps_df, "female")
        )
        downstream_2mb_female_high_risk_snps = len(
            filter_snps_by_partner_sex(downstream_2mb_high_risk_snps_df, "female")
        )
        within_gene_female_high_risk_snps = len(
            filter_snps_by_partner_sex(within_gene_high_risk_snps_df, "female")
        )
        upstream_2mb_female_high_risk_snps = len(
            filter_snps_by_partner_sex(upstream_2mb_high_risk_snps_df, "female")
        )
        downstream_2mb_male_low_risk_snps = len(
            filter_snps_by_partner_sex(downstream_2mb_low_risk_snps_df, "male")
        )
        within_gene_male_low_risk_snps = len(
            filter_snps_by_partner_sex(within_gene_low_risk_snps_df, "male")
        )
        upstream_2mb_male_low_risk_snps = len(
            filter_snps_by_partner_sex(upstream_2mb_low_risk_snps_df, "male")
        )
        total_female_low_risk_snps = len(
            filter_snps_by_partner_sex(total_low_risk_snps_df, "female")
        )
        downstream_2mb_female_low_risk_snps = len(
            filter_snps_by_partner_sex(downstream_2mb_low_risk_snps_df, "female")
        )
        within_gene_female_low_risk_snps = len(
            filter_snps_by_partner_sex(within_gene_low_risk_snps_df, "female")
        )
        upstream_2mb_female_low_risk_snps = len(
            filter_snps_by_partner_sex(upstream_2mb_low_risk_snps_df, "female")
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
        summary_snps_for_embryos_df = summary_snps_for_embryos_df.append(
            pd.Series(row, summary_snps_for_embryos_df.columns), ignore_index=True
        )

    return summary_snps_for_embryos_df


def plot_results(
    df,
    summary_df,
    embryo_ids,
    gene_start,
    gene_end,
):
    plots_as_html = []

    # preprocess dataframe for input to plotting function
    for embryo in embryo_ids:

        fig = px.scatter(
            df,
            x="Position",
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
            facet_col="snp_inherited_from",
            facet_col_wrap=1,
            color=f"{embryo}_risk_category",
            color_discrete_map={
                "high_risk": "#e60e0e",
                "low_risk": "#0ee60e",
                "NoCall": "#0818a6",
                "miscall": "#f0690a",
                "uninformative": "#52555e",
                "ADO": "#00ccff",
            },
            symbol=f"{embryo}_risk_category",
            symbol_map={
                "high_risk": "line-ns-open",
                "low_risk": "line-ns-open",
                "NoCall": "x",
                "miscall": "line-ns-open",
                "ADO": "line-ns-open",
                "uninformative": "line-ns-open",
            },
            category_orders={
                f"{embryo}_risk_category": [
                    "high_risk",
                    "low_risk",
                    "miscall",
                    "ADO",
                    "NoCall",
                    "uniformative",
                ],
                # TODO Update with future categories
                "snp_inherited_from": [
                    "male_partner",
                    "unassigned",
                    "female_partner",
                ],
            },
            labels={
                "y": "SNP Category",
                "Position": "Genomic coordinates",
            },
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
            hover_data={
                "Position": ":.0f",
                "probeset_id": True,
                "rsID": True,
            },
        )

        # Format facet plot labels
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.for_each_annotation(lambda a: a.update(xshift=-680))

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

        def get_counts(summary_df, embryo, multi_category):
            count = summary_df.loc[
                summary_df["embryo_id"] == embryo,
                multi_category,
            ].item()
            return count

        def add_snp_count_annotation(
            facet_row,
            embryo,
            summary_df,
            annotation_name_high_risk,
            annotation_name_low_risk,
            upstream_high_sum,
            upstream_low_sum,
            downstream_high_sum,
            downstream_low_sum,
        ):
            fig.add_trace(
                go.Scatter(
                    name=annotation_name_high_risk,
                    x=[
                        gene_start - 1900000,
                        gene_end + 1900000,
                    ],
                    y=[
                        2,
                        2,
                    ],
                    mode="text",
                    textfont_color="#e60e0e",
                    text=[
                        get_counts(summary_df, embryo, upstream_high_sum),
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
                        gene_end + 1900000,
                    ],
                    y=[-2, -2],
                    mode="text",
                    textfont_color="#0ee60e",
                    text=[
                        get_counts(summary_df, embryo, upstream_low_sum),
                        get_counts(summary_df, embryo, downstream_low_sum),
                    ],
                    textposition="bottom center",
                ),
                row=facet_row,  # The facet plot to anotate (1=bottom, 2=middle, 3=top)
                col=1,
            )

        add_snp_count_annotation(
            3,
            embryo,
            summary_df,
            "High_risk count male",
            "Low_risk count male",
            "upstream_2mb_male_high_risk_snps",
            "upstream_2mb_male_low_risk_snps",
            "downstream_2mb_male_high_risk_snps",
            "downstream_2mb_male_low_risk_snps",
        )
        # add_snp_count_annotation(
        #     2,
        #     embryo,
        #     "High_risk count unassigned",
        #     "Low_risk count unassigned",
        # )
        add_snp_count_annotation(
            1,
            embryo,
            summary_df,
            "High_risk count female",
            "Low_risk count female",
            "upstream_2mb_female_high_risk_snps",
            "upstream_2mb_female_low_risk_snps",
            "downstream_2mb_female_high_risk_snps",
            "downstream_2mb_female_low_risk_snps",
        )
        fig.update_yaxes(range=[-3, 3], showticklabels=False)
        fig.update_layout(height=540, width=1800, title_text=f"Results for {embryo}")

        plots_as_html.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))
    return plots_as_html
