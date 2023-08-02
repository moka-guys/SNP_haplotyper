def lookup_df_result(df, value_column, **lookup_values):
    """
    Filters a given DataFrame based on the criteria provided and sums up the values of a specified column for the filtered rows.  These summarised results can be used
    to compare the output results to theh expected results.
    Parameters:
    df (pandas.DataFrame): The DataFrame to be filtered.
    value_column (str): The name of the column whose values will be summed up after the DataFrame is filtered.
    **lookup_values (dict): Variable length dictionary. Column name(s) and value(s) to be used in the filter.
        Special case for 'snp_inherited_from' column where it filters rows whose value starts with the lookup_value, for other columns, filters rows whose values contain the lookup_value.

    Returns:
    total (float): The sum of the values in the specified column, value_column, for the rows that meet the filter criteria.
    """
    filter_criteria = []
    for lookup_colname, lookup_value in lookup_values.items():
        if lookup_colname == "snp_inherited_from":
            filter_criteria.append(f'{lookup_colname}.str.startswith("{lookup_value}")')
        else:
            filter_criteria.append(f'{lookup_colname}.str.contains("{lookup_value}")')
    filter_string = " and ".join(filter_criteria)
    print(filter_string)
    total = df.query(filter_string, engine="python")[value_column].sum()
    return total


def validate_snp_results(
    mode,
    sample_id,
    number_snps_imported,
    summary_snps_by_region,
    informative_snps_by_region,
    all_validation,
):
    """
    Validates the summary and informative SNP (Single Nucleotide Polymorphisms) results.  Compares the produced results to the expected results as imported from the launch.json.

    Arguments:
        mode (str): The mode of inheritance being tested. Currently supports "autosomal_dominant" and "autosomal_recessive".
        sample_id (str): The unique identifier for the sample.
        number_snps_imported (int): The number of SNPs imported.
        summary_snps_by_region (DataFrame): DataFrame containing summary SNPs information by gene region.
        informative_snps_by_region (DataFrame): DataFrame containing informative SNPs information by gene region.
        all_validation (dict): A dictionary containing validation results.

    Returns:
        None

    Raises:
        AssertionError: If any test fails to validate the results.
    """
    validation = all_validation[sample_id]
    assert mode == validation["mode"]
    assert sample_id == validation["sample_id"]
    assert number_snps_imported == validation["num_snps"]

    if mode == "autosomal_dominant":
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # These tests are the same for all modes of inheritance
        # 'within_gene' is within the genomic region of the gene, 'start' is 2MB or 3Mb upstream of the gene, 'end' is 2MB or 3MB downstream of the gene
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
                summary_snps_by_region, "snp_count", gene_distance="total_snps"
            )
            == validation["total_info_snps"]
        )

        # Test contents of informative_snps_by_region dataframe which is used to make the informative_snps_table in the report
        # These tests are different depending on the mode of inheritance
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
            summary_snps_by_region["snp_count"].sum()
            == validation["info_snps_upstream_2mb"]
            + validation["info_snps_in_gene"]
            + validation["info_snps_downstream_2mb"]
        )

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
    validation = all_validation[sample_id + "_" + embryo_id]
    assert mode == validation["mode"]
    assert sample_id == validation["sample_id"]

    if mode == "autosomal_dominant" or mode == "x_linked":
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                risk_category="high_risk",
            )
            == validation["high_risk_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                risk_category="high_risk",
            )
            == validation["high_risk_downstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                risk_category="low_risk",
            )
            == validation["low_risk_upstream_2mb"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                risk_category="low_risk",
            )
            == validation["low_risk_downstream_2mb"]
        )
    elif mode == "autosomal_recessive":
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                snp_inherited_from="female_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                snp_inherited_from="female_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                snp_inherited_from="male_partner",
                risk_category="high_risk",
            )
            == validation["high_risk_snps_downstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="upstream",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="within_gene",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_within_gene_from_male"]
        )
        assert (
            lookup_df_result(
                embryo_count_data_df,
                embryo_id + ".rhchp",
                snp_position="downstream",
                snp_inherited_from="male_partner",
                risk_category="low_risk",
            )
            == validation["low_risk_snps_downstream_2mb_from_male"]
        )
