import pytest
import pandas as pd
from pandas import testing as tm
from autosomal_dominant_logic import autosomal_dominant_analysis
from autosomal_recessive_logic import autosomal_recessive_analysis
from snp_haplotype import annotate_distance_from_gene
from exceptions import ArgumentInputError
from merge_array_files import main as merge_array_files_main

# Test Autosomal_dominant logic


# Test Autosomal_recessive logic
def test_ArgumentInputError():
    d = {
        "male_partner": [
            "AB",
        ],
        "female_partner": ["AA"],
        "reference": [
            "AA",
        ],
    }
    test_df = pd.DataFrame(data=d)
    reference_status = "unaffected"
    consanguineous = True

    # Should throw exception if reference_status = "unaffected" and consanguineous = True
    with pytest.raises(ArgumentInputError):
        results_df = autosomal_recessive_analysis(
            test_df,
            "male_partner",
            "female_partner",
            "reference",
            reference_status,
            consanguineous,
        )


@pytest.fixture
def setup_all_combination_of_inputs_AD():
    # Dictionary of all possible combinations of haplotypes for the trio
    d = {
        "affected_partner": [
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
        ],
        "unaffected_partner": [
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
        ],
        "reference": [
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
        ],
    }
    return d


@pytest.fixture
def setup_all_combination_of_inputs_AR():
    # Dictionary of all possible combinations of haplotypes for the trio
    d = {
        "male_partner": [
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
        ],
        "female_partner": [
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
            "AA",
            "BB",
            "AB",
        ],
        "reference": [
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "AA",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "BB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
            "AB",
        ],
    }
    return d


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_affected
@pytest.mark.ref_grandparent
def test_ref_affected_grandparent_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    reference_status = "affected"
    reference_relationship = "grandparent"
    results_df = autosomal_dominant_analysis(
        test_df,
        "affected_partner",
        "unaffected_partner",
        "reference",
        reference_status,
        reference_relationship,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 6)
            + [
                "low_risk",
                "high_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "high_risk",
                "low_risk",
            ]
            + (["uninformative"] * 10),
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_unaffected
@pytest.mark.ref_grandparent
def test_ref_unaffected_grandparent_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    reference_status = "unaffected"
    reference_relationship = "grandparent"
    results_df = autosomal_dominant_analysis(
        test_df,
        "affected_partner",
        "unaffected_partner",
        "reference",
        reference_status,
        reference_relationship,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 6)
            + [
                "high_risk",
                "low_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "low_risk",
                "high_risk",
            ]
            + (["uninformative"] * 10),
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_affected
@pytest.mark.ref_child
def test_ref_affected_child_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    reference_status = "affected"
    reference_relationship = "child"
    results_df = autosomal_dominant_analysis(
        test_df,
        "affected_partner",
        "unaffected_partner",
        "reference",
        reference_status,
        reference_relationship,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 6)
            + [
                "low_risk",
            ]
            + (["uninformative"] * 9)
            + [
                "low_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "high_risk",
                "high_risk",
                "uninformative",
            ],
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.mark.autosomal_dominant_logic
@pytest.mark.ref_unaffected
@pytest.mark.ref_child
def test_ref_unaffected_child_AD(setup_all_combination_of_inputs_AD):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD)
    reference_status = "unaffected"
    reference_relationship = "child"
    results_df = autosomal_dominant_analysis(
        test_df,
        "affected_partner",
        "unaffected_partner",
        "reference",
        reference_status,
        reference_relationship,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 6)
            + [
                "high_risk",
            ]
            + (["uninformative"] * 9)
            + [
                "high_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "low_risk",
                "low_risk",
                "uninformative",
            ],
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.mark.autosomal_recessive_logic
@pytest.mark.ref_affected
def test_ref_affected_AR(setup_all_combination_of_inputs_AR):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AR)
    reference_status = "affected"
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        consanguineous=False,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 2)
            + [
                "low_risk",
            ]
            + (["uninformative"] * 3)
            + [
                "low_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "low_risk",
                "uninformative",
                "low_risk",
            ]
            + (["uninformative"] * 3)
            + [
                "high_risk",
            ]
            + (["uninformative"] * 2)
            + (["high_risk"] * 3)
            + ["uninformative"],
            "snp_inherited_from": (["uninformative"] * 2)
            + [
                "female_partner",
            ]
            + (["uninformative"] * 3)
            + [
                "male_partner",
            ]
            + (["uninformative"] * 7)
            + [
                "female_partner",
                "uninformative",
                "male_partner",
            ]
            + (["uninformative"] * 3)
            + [
                "female_partner",
            ]
            + (
                ["uninformative"] * 2
                + [
                    "female_partner",
                ]
            )
            + (["male_partner"] * 2)
            + ["uninformative"],
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.mark.autosomal_recessive_logic
@pytest.mark.ref_unaffected
def test_ref_unaffected_AR(setup_all_combination_of_inputs_AR):
    test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AR)
    reference_status = "unaffected"
    results_df = autosomal_recessive_analysis(
        test_df,
        "male_partner",
        "female_partner",
        "reference",
        reference_status,
        consanguineous=False,
    )
    expected_results_df = pd.DataFrame(
        data={
            "snp_risk_category": (["uninformative"] * 2)
            + [
                "high_risk",
            ]
            + (["uninformative"] * 3)
            + [
                "high_risk",
            ]
            + (["uninformative"] * 7)
            + [
                "high_risk",
                "uninformative",
                "high_risk",
            ]
            + (["uninformative"] * 3)
            + [
                "low_risk",
            ]
            + (["uninformative"] * 2)
            + (["low_risk"] * 3)
            + ["uninformative"],
            "snp_inherited_from": (["uninformative"] * 2)
            + [
                "female_partner",
            ]
            + (["uninformative"] * 3)
            + [
                "male_partner",
            ]
            + (["uninformative"] * 7)
            + [
                "female_partner",
                "uninformative",
                "male_partner",
            ]
            + (["uninformative"] * 3)
            + [
                "female_partner",
            ]
            + (
                ["uninformative"] * 2
                + [
                    "female_partner",
                ]
            )
            + (["male_partner"] * 2)
            + ["uninformative"],
        }
    )
    tm.assert_series_equal(
        results_df["snp_risk_category"], expected_results_df["snp_risk_category"]
    )


@pytest.fixture
def setup_all_gene_regions():
    # Dictionary of a range of genomic coordinates around a gene (chr1:4000000-4001000), testing edge effects
    d = {
        "probeset_id": list(range(0, 21)),
        "Position": [
            1,
            1000000,
            1999999,
            2000000,
            2000001,
            2999999,
            3000000,
            3000001,
            3999999,
            4000000,
            4000001,
            4000999,
            4001000,
            4001001,
            5000999,
            5001000,
            5001001,
            6000999,
            6001000,
            6001001,
            7000000,
        ],
    }
    return d


@pytest.mark.gene_region
def test_annotate_distance_from_gene(setup_all_gene_regions):
    test_df = pd.DataFrame(data=setup_all_gene_regions)
    chr = 1
    gene_start = 4000000
    gene_end = 4001000
    results_df = annotate_distance_from_gene(test_df, chr, gene_start, gene_end)
    expected_results_df = pd.DataFrame(
        data={
            "gene_distance": (["outside_range"] * 4)
            + (["1-2MB_from_start"] * 3)
            + (["0-1MB_from_start"] * 3)
            + (["within_gene"] * 3)
            + (["0-1MB_from_end"] * 3)
            + (["1-2MB_from_end"] * 3)
            + (["outside_range"] * 2),
        }
    )
    expected_results_df["gene_distance"] = expected_results_df["gene_distance"].astype(
        "category"
    )
    expected_results_df["gene_distance"] = expected_results_df[
        "gene_distance"
    ].cat.set_categories(
        [
            "1-2MB_from_start",
            "0-1MB_from_start",
            "within_gene",
            "0-1MB_from_end",
            "1-2MB_from_end",
        ],
    )
    tm.assert_series_equal(
        results_df["gene_distance"], expected_results_df["gene_distance"]
    )


# def test_criteria_AR1():
#     test_df = pd.DataFrame(data=setup_all_combination_of_inputs_AD())
#     reference_status = "affected"
#     consanguineous = True
#     results_df = autosomal_recessive_analysis(
#         test_df,
#         "male_partner",
#         "female_partner",
#         "reference",
#         reference_status,
#         "ref_status",
#         consanguineous,
#     )
#     assert results_df["snp_risk_category"].unique().tolist() == [
#         "male_partner_low_risk"
#     ]


@pytest.mark.handling_multiple_input_files
def test_merging_input_files():
    truth_df = pd.read_csv("test_data/autosomal_dominant/F4_BRCA2_AD.txt", sep="\t")
    merged_df = merge_array_files_main(
        [
            "test_data/precases/F4_BRCA2_AD_precases.txt",
            "test_data/precases/F4_BRCA2_AD_embryos.txt",
        ]
    )
    assert truth_df.equals(merged_df)
    truth_df = pd.read_csv("test_data/autosomal_recessive/F3_HBB_AR.txt", sep="\t")
    merged_df = merge_array_files_main(
        [
            "test_data/precases/F3_HBB_AR_precases.txt",
            "test_data/precases/F3_HBB_AR_embryos.txt",
        ]
    )
    assert truth_df.equals(merged_df)
    truth_df = pd.read_csv("test_data/x_linked/F10_FMR1_XL.txt", sep="\t")
    merged_df = merge_array_files_main(
        [
            "test_data/precases/F10_FMR1_XL_precases.txt",
            "test_data/precases/F10_FMR1_XL_embryos.txt",
        ]
    )
    assert truth_df.equals(merged_df)


@pytest.mark.handling_multiple_input_files
def test_multifile_input_columns_duplicates():
    # Test that the input files have the same columns, pytest will fail if any exception is raised
    with pytest.raises(Exception):
        merge_array_files_main(
            [
                "test_data/autosomal_dominant/F4_BRCA2_AD.txt",
                "test_data/autosomal_dominant/F4_BRCA2_AD.txt",
            ]
        )


@pytest.mark.handling_multiple_input_files
def test_multifile_input_probesets_are_identical():
    with pytest.raises(Exception):
        merge_array_files_main(
            [
                "test_data/autosomal_dominant/F4_BRCA2_AD.txt",
                "test_data/autosomal_dominant/F5_COL1A1_AD.txt",
            ]
        )
