import sys


def get_snp_count(df, pattern, category_col, snp_count_col):
    count = int(
        df[
            (df["gene_distance"].str.endswith(pattern))
            & (df[category_col] != "uninformative")
        ][snp_count_col].sum()
    )
    return count


def get_total_snp_count(df, category_col, snp_count_col):
    total = int(df[(df[category_col] != "uninformative")][snp_count_col].sum())
    return total


def get_snp_count_and_filter(
    df,
    pattern,
    category_col,
    snp_count_col,
    risk_filter,
    inherited_from_col="None",  # For Autosomal Recessive samples also filter on the partner SNP inherited from
):
    if inherited_from_col == "None":
        count = int(
            df[
                (df["gene_distance"].str.endswith(pattern))
                & (df[category_col] == risk_filter)
            ][snp_count_col].sum()
        )
    else:
        # For Autosomal Recessive samples also filter on the partner SNP inherited from
        count = int(
            df[
                (df["gene_distance"].str.endswith(pattern))
                & (df[category_col] == risk_filter)
                & (df["snp_inherited_from"] == inherited_from_col)
            ][snp_count_col].sum()
        )
    return count


def stream_autosomal_dominant_output(
    mode_of_inheritance,
    informative_snps_by_region,
    embryo_snps_summary_df,
    number_snps_imported,
    output_prefix,
):
    """Converts autosomal dominant output into JSON compatible format to stream for use in testing

    New column created in the dataframe, df, matching the probes_set IDs to dbSNP rsIDs.

    Args:
        mode_of_inheritance (string): "autosomal_dominant"
        informative_snps_by_region (dataframe): informative_snps_by_region dataframe
        embryo_snps_summary_df (dataframe): embryo_snps_summary_df dataframe
        number_snps_imported (string): Number of SNPs in imported input file
        output_prefix (string): Sample ID
    Returns:
        list: List contains two dictionaries informative_snp_data & embryo_cat_data

    """
    informative_snp_data = {
        "mode": mode_of_inheritance,
        "sample_id": output_prefix,
        "num_snps": number_snps_imported,
        "info_snps_upstream_2mb": get_snp_count(
            informative_snps_by_region, "_start", "snp_risk_category", "snp_count"
        ),
        "info_snps_in_gene": get_snp_count(
            informative_snps_by_region, "within_gene", "snp_risk_category", "snp_count"
        ),
        "info_snps_downstream_2mb": get_snp_count(
            informative_snps_by_region, "_end", "snp_risk_category", "snp_count"
        ),
        "total_info_snps": get_total_snp_count(
            informative_snps_by_region, "snp_risk_category", "snp_count"
        ),
        "high_risk_snps_upstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "high_risk",
        ),
        "high_risk_snps_within_gene": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "high_risk",
        ),
        "high_risk_snps_downstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "high_risk",
        ),
        "low_risk_snps_upstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "low_risk",
        ),
        "low_risk_snps_within_gene": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "low_risk",
        ),
        "low_risk_snps_downstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "low_risk",
        ),
    }
    # Populate stream with additional fields:
    embryo_snps_output_df = embryo_snps_summary_df
    embryo_snps_output_df["mode"] = mode_of_inheritance
    embryo_snps_output_df["sample_id"] = output_prefix
    embryo_cat_data = embryo_snps_output_df.to_dict(orient="record")
    return [informative_snp_data, embryo_cat_data]


def stream_autosomal_recessive_output(
    mode_of_inheritance,
    informative_snps_by_region,
    embryo_snps_summary_df,
    number_snps_imported,
    output_prefix,
):
    """Converts autosomal recessive output into JSON compatible format to stream for use in testing

    New column created in the dataframe, df, matching the probes_set IDs to dbSNP rsIDs.

    Args:
        mode_of_inheritance (string): "autosomal_recessive"
        informative_snps_by_region (dataframe): informative_snps_by_region dataframe
        embryo_snps_summary_df (dataframe): embryo_snps_summary_df dataframe
        number_snps_imported (string): Number of SNPs in imported input file
        output_prefix (string): Sample ID
    Returns:
        list: List contains two dictionaries informative_snp_data & embryo_cat_data

    """
    informative_snp_data = {
        "mode": mode_of_inheritance,
        "sample_id": output_prefix,
        "num_snps": number_snps_imported,
        "info_snps_upstream_2mb": get_snp_count(
            informative_snps_by_region, "_start", "snp_risk_category", "snp_count"
        ),
        "info_snps_in_gene": get_snp_count(
            informative_snps_by_region, "within_gene", "snp_risk_category", "snp_count"
        ),
        "info_snps_downstream_2mb": get_snp_count(
            informative_snps_by_region, "_end", "snp_risk_category", "snp_count"
        ),
        "total_info_snps": get_total_snp_count(
            informative_snps_by_region, "snp_risk_category", "snp_count"
        ),
        "high_risk_snps_upstream_2mb_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "female_partner",
        ),
        "high_risk_within_gene_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "female_partner",
        ),
        "high_risk_snps_downstream_2mb_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "female_partner",
        ),
        "low_risk_snps_upstream_2mb_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "female_partner",
        ),
        "low_risk_within_gene_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "female_partner",
        ),
        "low_risk_snps_downstream_2mb_from_female": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "female_partner",
        ),
        "high_risk_snps_upstream_2mb_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "male_partner",
        ),
        "high_risk_within_gene_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "male_partner",
        ),
        "high_risk_snps_downstream_2mb_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "high_risk",
            "male_partner",
        ),
        "low_risk_snps_upstream_2mb_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "male_partner",
        ),
        "low_risk_within_gene_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "male_partner",
        ),
        "low_risk_snps_downstream_2mb_from_male": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "snp_risk_category",
            "snp_count",
            "low_risk",
            "male_partner",
        ),
    }
    # Populate stream with additional fields:
    embryo_snps_output_df = embryo_snps_summary_df
    embryo_snps_output_df["mode"] = mode_of_inheritance
    embryo_snps_output_df["sample_id"] = output_prefix
    embryo_cat_data = embryo_snps_output_df.to_dict(orient="record")
    return [informative_snp_data, embryo_cat_data]


def stream_x_linked_output(
    mode_of_inheritance,
    informative_snps_by_region,
    embryo_snps_summary_df,
    number_snps_imported,
    output_prefix,
):
    """Converts x-linked output into JSON compatible format to stream for use in testing

    New column created in the dataframe, df, matching the probes_set IDs to dbSNP rsIDs.

    Args:
        mode_of_inheritance (string): "x_linked"
        informative_snps_by_region (dataframe): informative_snps_by_region dataframe
        embryo_snps_summary_df (dataframe): embryo_snps_summary_df dataframe
        number_snps_imported (string): Number of SNPs in imported input file
        output_prefix (string): Sample ID
    Returns:
        list: List cotains two dictionaries informative_snp_data & embryo_cat_data

    """
    informative_snp_data = {
        "mode": mode_of_inheritance,
        "sample_id": output_prefix,
        "num_snps": number_snps_imported,
        "info_snps_upstream_2mb": get_snp_count(
            informative_snps_by_region,
            "_start",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
        ),
        "info_snps_in_gene": get_snp_count(
            informative_snps_by_region,
            "within_gene",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
        ),
        "info_snps_downstream_2mb": get_snp_count(
            informative_snps_by_region,
            "_end",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
        ),
        "total_info_snps": get_total_snp_count(
            informative_snps_by_region,
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
        ),
        "high_risk_snps_upstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "high_risk",
        ),
        "high_risk_snps_within_gene": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "high_risk",
        ),
        "high_risk_snps_downstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "high_risk",
        ),
        "low_risk_snps_upstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_start",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "low_risk",
        ),
        "low_risk_snps_within_gene": get_snp_count_and_filter(
            informative_snps_by_region,
            "within_gene",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "low_risk",
        ),
        "low_risk_snps_downstream_2mb": get_snp_count_and_filter(
            informative_snps_by_region,
            "_end",
            "female_AB_snp_risk_category",
            "female_AB_snp_count",
            "low_risk",
        ),
    }

    # Populate stream with additional fields:
    embryo_snps_output_df = embryo_snps_summary_df
    embryo_snps_output_df["mode"] = mode_of_inheritance
    embryo_snps_output_df["sample_id"] = output_prefix

    embryo_cat_data = embryo_snps_output_df.to_dict(orient="record")
    return [informative_snp_data, embryo_cat_data]
