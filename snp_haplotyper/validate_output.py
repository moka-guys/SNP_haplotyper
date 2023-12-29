import pandas as pd
from EnumDataClasses import InheritanceMode
from pandas.testing import assert_frame_equal


def dict_to_dataframe_AD(data, consanguineous=False):
    """
    Converts an input dictionary of autosomal dominant benchmark data read in from a JSON
    file into a DataFrame in the same format to that produced for autosomal dominant cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]
    risk_categories = ["high_risk", "low_risk"]
    all_categories = [
        "high_risk",
        "high_or_low_risk",
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    organized_data = [
        {
            "snp_position": pos,
            "snp_risk_category_summary": risk,
            "snp_count": data.get(
                f"{risk}_snps_{pos}", 0
            ),  # Using get method to handle missing keys
        }
        for pos in positions
        for risk in risk_categories
    ]

    df = pd.DataFrame(organized_data).set_index(
        ["snp_position", "snp_risk_category_summary"]
    )

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(df.index.names, [positions, all_categories])
        ]
    )
    return df


def dict_to_dataframe_AR(data, consanguineous=False):
    """
    Converts an input dictionary of autosomal recessive benchmark data read in from a JSON
    file into a DataFrame in the same format as that produced for autosomal recessive cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]

    if consanguineous:
        inherited_from_categories = ["male_partner", "both_partners", "female_partner"]
        risk_categories = ["high_risk", "high_or_low_risk", "low_risk"]
    else:
        inherited_from_categories = ["male_partner", "female_partner"]
        risk_categories = ["high_risk", "low_risk"]

    # Initialize organized_data with zero counts
    organized_data = [
        {
            "snp_position": pos,
            "snp_risk_category_summary": risk,
            "snp_inherited_from": inherited_from,
            "snp_count": 0,
        }
        for pos in positions
        for risk in risk_categories
        for inherited_from in inherited_from_categories
    ]

    # Update counts with values from input data if they exist
    for entry in organized_data:
        key = f"{entry['snp_risk_category_summary']}_snps_{entry['snp_position']}_from_{entry['snp_inherited_from']}"
        entry["snp_count"] = data.get(key, entry["snp_count"])

    # Create DataFrame and set multi-index
    df = pd.DataFrame(organized_data).set_index(
        ["snp_inherited_from", "snp_risk_category_summary", "snp_position"]
    )

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(
                df.index.names, [inherited_from_categories, risk_categories, positions]
            )
        ]
    )
    df = df.sort_index()

    if consanguineous:
        # Remove rows where snp_inherited_from is both_partners and snp_risk_category_summary is high_or_low_risk
        df = df.query(
            "not ("
            "(snp_inherited_from == 'male_partner' & snp_risk_category_summary == 'high_or_low_risk') | "
            "(snp_inherited_from == 'female_partner' & snp_risk_category_summary == 'high_or_low_risk') | "
            "(snp_inherited_from == 'both_partners' & (snp_risk_category_summary == 'high_risk' | snp_risk_category_summary == 'low_risk'))"
            ")"
        )

    return df


def dict_to_dataframe_XL(data, consanguineous=False):
    """
    Converts an input dictionary of autosomal dominant benchmark data read in from a JSON
    file into a DataFrame in the same format to that produced for autosomal dominant cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]
    risk_categories = ["high_risk", "low_risk"]
    all_categories = [
        "high_risk",
        "high_or_low_risk",
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    organized_data = [
        {
            "snp_position": pos,
            "snp_risk_category_summary": risk,
            "snp_count": data.get(
                f"{risk}_snps_{pos}", 0
            ),  # Using get method to handle missing keys
        }
        for pos in positions
        for risk in risk_categories
    ]

    df = pd.DataFrame(organized_data).set_index(
        ["snp_position", "snp_risk_category_summary"]
    )

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(df.index.names, [positions, all_categories])
        ]
    )
    return df


def dict_to_dataframe_AD_embryo(data, consanguineous=False):
    """
    Converts an input dictionary of autosomal dominant benchmark data read in from a JSON
    file into a DataFrame in the same format to that produced for autosomal dominant cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]
    risk_categories = ["high_risk", "low_risk"]
    all_categories = [
        "high_risk",
        # "high_or_low_risk", TODO
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    organized_data = [
        {
            "embryo_risk_category": risk,
            "snp_position": pos,
            "snp_count": data.get(
                f"{risk}_snps_{pos}", 0
            ),  # Using get method to handle missing keys
        }
        for risk in risk_categories
        for pos in positions
    ]

    df = pd.DataFrame(organized_data).set_index(
        ["embryo_risk_category", "snp_position"]
    )

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(df.index.names, [all_categories, positions])
        ]
    )
    return df


def dict_to_dataframe_AR_embryo(data, consanguineous=False):
    """
    Converts an input dictionary of autosomal recessive benchmark data read in from a JSON
    file into a DataFrame in the same format as that produced for autosomal recessive cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]

    inherited_from_categories = [
        "male_partner",
        "both_partners",
        "female_partner",
        "unassigned",
    ]
    risk_categories = ["high_risk", "low_risk"]

    # Initialize organized_data with zero counts
    organized_data = [
        {
            "snp_position": pos,
            "snp_risk_category_summary": risk,
            "snp_inherited_from": inherited_from,
            "snp_count": 0,
        }
        for pos in positions
        for risk in risk_categories
        for inherited_from in inherited_from_categories
    ]

    # Update counts with values from input data if they exist
    for entry in organized_data:
        key = f"{entry['snp_risk_category_summary']}_snps_{entry['snp_position']}_from_{entry['snp_inherited_from']}"
        entry["snp_count"] = data.get(key, entry["snp_count"])

    # Create DataFrame and set multi-index
    df = pd.DataFrame(organized_data).set_index(
        ["snp_inherited_from", "snp_risk_category_summary", "snp_position"]
    )

    full_risk_categories = [
        "high_risk",
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(
                df.index.names,
                [inherited_from_categories, full_risk_categories, positions],
            )
        ]
    )
    df = df.sort_index()

    if consanguineous:
        pass
    else:
        # Remove rows where snp_inherited_from is both_partners and snp_risk_category_summary is high_or_low_risk
        # This is necessary summary values of 0 will be provided for each of these categories as they are valid
        # categorical values in the dataframe column, even though they are not viable combinations for this mode of inheritance
        df = df.query(
            "not ("
            "(snp_inherited_from == 'male_partner' & snp_risk_category_summary == 'high_or_low_risk') | "
            "(snp_inherited_from == 'female_partner' & snp_risk_category_summary == 'high_or_low_risk') | "
            "(snp_inherited_from == 'both_partners' & snp_risk_category_summary == 'high_or_low_risk' ) | "
            "(snp_inherited_from == 'unassigned' )"
            ")"
        )
    df.index.names = ["snp_inherited_from", "embryo_risk_category", "snp_position"]

    return df


def dict_to_dataframe_XL_embryo(data, consanguineous=False):
    """
    Converts an input dictionary of x-Linked benchmark data read in from a JSON
    file into a DataFrame in the same format to that produced for autosomal dominant cases.

    Parameters:
    data (dict): Input dictionary with specific keys.

    Returns:
    pd.DataFrame: DataFrame with multi-index.
    """
    positions = ["upstream", "within_gene", "downstream"]
    risk_categories = ["high_risk", "low_risk"]
    all_categories = [
        "high_risk",
        # "high_or_low_risk", TODO
        "low_risk",
        "uninformative",
        "ADO",
        "miscall",
        "NoCall",
        "NoCall_in_trio",
    ]

    organized_data = [
        {
            "embryo_risk_category": risk,
            "snp_position": pos,
            "snp_count": data.get(
                f"{risk}_snps_{pos}", 0
            ),  # Using get method to handle missing keys
        }
        for risk in risk_categories
        for pos in positions
    ]

    df = pd.DataFrame(organized_data).set_index(
        ["embryo_risk_category", "snp_position"]
    )

    # Convert index levels to CategoricalIndex with specified categories and order
    df.index = pd.MultiIndex.from_arrays(
        [
            pd.CategoricalIndex(
                df.index.get_level_values(level), categories=categories, ordered=True
            )
            for level, categories in zip(df.index.names, [all_categories, positions])
        ]
    )
    return df


def validate_snp_results(
    mode,
    sample_id,
    number_snps_imported,
    informative_snps_by_region,
    embryo_snps_by_region,
    all_validation,
    consanguineous=False,
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
    assert mode == InheritanceMode(validation["mode"])
    assert sample_id == validation["sample_id"]
    assert number_snps_imported == validation["num_snps"]

    if mode == InheritanceMode.AUTOSOMAL_DOMINANT:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame

        AD_benchmark_df = dict_to_dataframe_AD(validation)
        assert_frame_equal(informative_snps_by_region, AD_benchmark_df)

    elif mode == InheritanceMode.AUTOSOMAL_RECESSIVE:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame
        AR_benchmark_df = dict_to_dataframe_AR(validation, consanguineous)
        assert_frame_equal(informative_snps_by_region, AR_benchmark_df)

    elif mode == InheritanceMode.X_LINKED:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame
        XL_benchmark_df = dict_to_dataframe_XL(validation)
        XL_benchmark_df = pd.DataFrame(
            XL_benchmark_df.snp_count.tolist(), index=XL_benchmark_df.index
        )
        XL_benchmark_df.columns = [
            "snp_count_female_AB",
            "snp_count_male_AA",
            "snp_count_male_BB",
        ]
        XL_benchmark_df = XL_benchmark_df.sort_index()

        assert_frame_equal(informative_snps_by_region, XL_benchmark_df)


def validate_embryo_results(
    mode: InheritanceMode,
    sample_id: str,
    number_snps_imported: int,
    embryo_id: str,
    embryo_count_data_df: pd.DataFrame,
    all_validation: dict,
    consanguineous=False,
):
    validation = all_validation[sample_id + "_" + embryo_id]
    assert mode == InheritanceMode(validation["mode"])
    assert sample_id == validation["sample_id"]
    # assert number_snps_imported == validation["num_snps"] TODO

    if mode == InheritanceMode.AUTOSOMAL_DOMINANT:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame

        AD_benchmark_df = dict_to_dataframe_AD_embryo(validation)

        test_df = embryo_count_data_df[
            [
                "embryo_risk_category",
                "snp_position",
                f"{embryo_id}.rhchp",
            ]
        ]

        test_df.columns = [
            "embryo_risk_category",
            "snp_position",
            "snp_count",
        ]

        assert_frame_equal(test_df, AD_benchmark_df.reset_index())

    elif mode == InheritanceMode.AUTOSOMAL_RECESSIVE:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame
        AR_benchmark_df = dict_to_dataframe_AR_embryo(validation, consanguineous)

        test_df = embryo_count_data_df[
            [
                "snp_inherited_from",
                "embryo_risk_category",
                "snp_position",
                f"{embryo_id}.rhchp",
            ]
        ]

        test_df.columns = [
            "snp_inherited_from",
            "embryo_risk_category",
            "snp_position",
            "snp_count",
        ]

        AR_benchmark_df = AR_benchmark_df[
            AR_benchmark_df.index.isin(
                ["male_partner", "both_partners", "female_partner"], level=0
            )
        ]

        AR_benchmark_df = AR_benchmark_df.reset_index()

        AR_benchmark_df = AR_benchmark_df.sort_values(
            [
                "snp_inherited_from",
                "embryo_risk_category",
                "snp_position",
            ]
        )

        test_df = test_df.sort_values(
            [
                "snp_inherited_from",
                "embryo_risk_category",
                "snp_position",
            ]
        )

        AR_benchmark_df = AR_benchmark_df.reset_index(drop=True)
        test_df = test_df.reset_index(drop=True)

        assert_frame_equal(test_df, AR_benchmark_df)

    elif mode == InheritanceMode.X_LINKED:
        # Test contents of summary_snps_by_region dataframe which is used to make the summary_snps_table in the report
        # We import the benchmark data as a DataFrame and compare it to the produced DataFrame
        XL_benchmark_df = dict_to_dataframe_XL_embryo(validation)

        test_df = embryo_count_data_df[
            [
                "embryo_risk_category",
                "snp_position",
                f"{embryo_id}.rhchp",
            ]
        ]

        test_df.columns = [
            "embryo_risk_category",
            "snp_position",
            "snp_count",
        ]

        assert_frame_equal(test_df, XL_benchmark_df.reset_index())
