def lookup_df_result(df, value_column, **lookup_values):
    total = 0
    for lookup_colname, lookup_value in lookup_values.items():
        cell_values = df.query(
            f'{lookup_colname}.str.contains("{lookup_value}")', engine="python"
        )[value_column].sum()
        total = total + cell_values
    return total


def validate_snp_results(
    mode, sample_id, summary_snps_by_region, informative_snps_by_region, all_validation
):
    validation = all_validation[sample_id]
    assert mode == validation["mode"]
    assert sample_id == validation["sample_id"]
    # assert results["num_snps"] == validation["num_snps"] TODO replicate this code

    # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
    # These tests are the same for all modes of inheritance
    assert (
        lookup_df_result(summary_snps_by_region, "snp_count", gene_distance="start")
        == validation["info_snps_upstream_2mb"]
    )
    assert (
        lookup_df_result(
            summary_snps_by_region, "snp_count", gene_distance="within_gene"
        )
        == validation["info_snps_in_gene"]
    )
    assert (
        lookup_df_result(summary_snps_by_region, "snp_count", gene_distance="end")
        == validation["info_snps_downstream_2mb"]
    )
    assert (
        lookup_df_result(
            summary_snps_by_region, "snp_count", gene_distance="total_info_snps"
        )
        == validation["total_info_snps"]
    )

    # Test contents of informative_snps_by_region dataframe which is used to make the informative_snps_table in the report
    # These tests are different depending on the mode of inheritance
    if mode == "autosomal_dominant" or mode == "x_linked":
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="high_risk",
            )
            == validation["high_risk_snps_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="high_risk",
            )
            == validation["high_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="high_risk",
            )
            == validation["high_risk_snps_downstream_2mb"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="low_risk",
            )
            == validation["low_risk_snps_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="low_risk",
            )
            == validation["low_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="low_risk",
            )
            == validation["low_risk_snps_downstream_2mb"]
        )
    elif mode == "autosomal_recessive":
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="high_risk",
                snp_inherited_from="female_partner",
            )
            == validation["high_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="high_risk",
                snp_inherited_from="female_partner",
            )
            == validation["high_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="high_risk",
                snp_inherited_from="female_partner",
            )
            == validation["high_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="low_risk",
                snp_inherited_from="female_partner",
            )
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="low_risk",
                snp_inherited_from="female_partner",
            )
            == validation["low_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="low_risk",
                snp_inherited_from="female_partner",
            )
            == validation["low_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="high_risk",
                snp_inherited_from="male_partner",
            )
            == validation["high_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="high_risk",
                snp_inherited_from="male_partner",
            )
            == validation["high_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="high_risk",
                snp_inherited_from="male_partner",
            )
            == validation["high_risk_snps_downstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="start",
                snp_risk_category="low_risk",
                snp_inherited_from="male_partner",
            )
            == validation["low_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="within_gene",
                snp_risk_category="low_risk",
                snp_inherited_from="male_partner",
            )
            == validation["low_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                informative_snps_by_region,
                "snp_count",
                gene_distance="end",
                snp_risk_category="low_risk",
                snp_inherited_from="male_partner",
            )
            == validation["low_risk_snps_downstream_2mb_from_male"]
        )


def validate_embryo_results(
    mode, sample_id, embryo_id, embryo_count_data_df, all_validation
):
    validation = all_validation[sample_id]
    assert mode == validation["mode"]
    assert sample_id == validation["sample_id"]

    if mode == "autosomal_dominant" or mode == "x_linked":
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                risk_category="high_risk",
            )
            == validation["high_risk_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_gene",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                risk_category="high_risk",
            )
            == validation["high_risk_downstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                risk_category="low_risk",
            )
            == validation["low_risk_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_gene",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                risk_category="low_risk",
            )
            == validation["low_risk_downstream_2mb"]
        )
    elif mode == "autosomal_recessive":
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_gene",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_gene",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_genes",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_downstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="upstream",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="within_gene",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id,
                snp_position="downstream",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_downstream_2mb_from_male"]
        )
