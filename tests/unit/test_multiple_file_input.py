import pytest
import pandas as pd
from merge_array_files import main as merge_array_files_main


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
