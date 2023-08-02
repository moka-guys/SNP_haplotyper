from argparse import Namespace
from xml.dom import ValidationErr
from snp_haplotype import main as snp_haplotype_main
import os
import pandas as pd
import pytest
import json
from validate_output import validate_snp_results, validate_embryo_results

from excel_parser import main as excel_parser_main
from snp_haplotype import main as snp_haplotype_main
from test_snp_filtering import setup_test_data
import re


def setup_test_data_from_excel(split_by_embryo=False, version=1):
    """
    Creates a dictionary of test data sets based on provided Excel files in specified directories.
    The dictionary keys are run names, and the values are file paths for each run, stored as Namespace objects.

    If 'split_by_embryo' is set to True, the function will split the data sets based on embryo ids from each run.

    The 'version' parameter determines the folder from which the Excel files are retrieved. Currently supports version 1 and 2.

    Parameters:
    split_by_embryo (bool, optional): A flag to determine whether to split the data sets by embryo id or not. Defaults to False.
    version (int, optional): Specifies the version of test data templates to be used. Defaults to 1.

    Raises:
    ValueError: If an unsupported version number is provided.

    Returns:
    dict: A dictionary containing the setup test data. Keys are the run names and values are Namespace objects containing the input spreadsheet file path and optional snp_array_file and run_basher values.
    """
    # We need to get the embryo_ids for each sample from the launch JSON so we can auto generate
    # the test names for each embryo without prereading the excel files
    embryo_data = setup_test_data(False)

    # get list of excel files in folder
    if version == 1:
        excel_folder_path = "test_data/template_test_data_v1"
    elif version == 2:
        excel_folder_path = "test_data/template_test_data_v2"
    else:
        raise ValueError(f"Unsupported version: {version}")

    excel_files = [f for f in os.listdir(excel_folder_path) if f.endswith(".xlsx")]

    # iterate through excel files creating a dictionary with run names as keys and the values are
    # file paths for each run
    run_data_dict = {}
    for excel_file in excel_files:
        run_name = re.sub(
            r"_v\d+\.xlsx$", "", excel_file.replace("excel_test_", "")
        ).replace(".xlsx", "")
        if split_by_embryo == False:
            run_data_dict[run_name] = Namespace(
                input_spreadsheet=os.path.join(
                    excel_folder_path,
                    excel_file,
                ),
                snp_array_file=None,
                run_basher=False,
            )
        elif split_by_embryo == True:
            for embryo_id in embryo_data[run_name].embryo_ids:
                test_name = f"{run_name}_{embryo_id}".replace(".rhchp", "")
                run_data_dict[test_name] = Namespace(
                    input_spreadsheet=os.path.join(
                        excel_folder_path,
                        excel_file,
                    ),
                    snp_array_file=None,
                    run_basher=False,
                )

    return run_data_dict


@pytest.mark.parametrize("name", setup_test_data_from_excel(False, version=1))
def test_informative_snps_excel_v1(name):
    """
    Tests whether the Single Nucleotide Polymorphisms (SNPs) extracted from the provided excel files
    (version 1) are as expected by comparing them with the validation data.

    The function parameters are parametrized using pytest's mark.parametrize decorator. For each test case,
    the function parses the excel file, processes the SNPs, and compares them with the expected output stored in a JSON file.
    If the SNPs extracted do not match the expected SNPs, the test will fail.

    Parameters:
    name (str): The name of the test case. This name is used to fetch the corresponding excel file for testing.

    Raises:
    AssertionError: If the actual output does not match the expected output.
    """
    test_args = setup_test_data_from_excel()
    basher_input_namespace, error_dictionary, input_ok_flag = excel_parser_main(
        test_args[name]
    )

    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = snp_haplotype_main(basher_input_namespace)

    json_file_path = "test_data/informative_snp_validation.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        snp_validation = json.loads(j.read())
        for dict in snp_validation:
            all_validation[dict["sample_id"]] = dict

    validate_snp_results(
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        all_validation,
    )


@pytest.mark.parametrize("name", setup_test_data_from_excel(False, version=2))
def test_informative_snps_excel_v2(name):
    """
    Tests whether the Single Nucleotide Polymorphisms (SNPs) extracted from the provided excel files
    (version 2) are as expected by comparing them with the validation data.

    The function parameters are parametrized using pytest's mark.parametrize decorator. For each test case,
    the function parses the excel file, processes the SNPs, and compares them with the expected output stored in a JSON file.
    If the SNPs extracted do not match the expected SNPs, the test will fail.

    Parameters:
    name (str): The name of the test case. This name is used to fetch the corresponding excel file for testing.

    Raises:
    AssertionError: If the actual output does not match the expected output.
    """
    test_args = setup_test_data_from_excel()
    basher_input_namespace, error_dictionary, input_ok_flag = excel_parser_main(
        test_args[name]
    )

    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = snp_haplotype_main(basher_input_namespace)

    json_file_path = "test_data/informative_snp_validation.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        snp_validation = json.loads(j.read())
        for dict in snp_validation:
            all_validation[dict["sample_id"]] = dict

    validate_snp_results(
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        all_validation,
    )


@pytest.mark.parametrize("name", setup_test_data_from_excel(True, version=1))
def test_embryo_categorization_excel_v1(name):
    """
    Tests the categorization of embryos based on the Single Nucleotide Polymorphisms (SNPs)
    data provided in the excel files (version 1).

    The function parameters are parametrized using pytest's mark.parametrize decorator. For each
    test case, the function parses the excel file, processes the SNPs, and then compares the
    embryo categorization results with the expected results stored in a JSON file.

    Parameters:
    name (str): The name of the test case. The name is used to fetch the corresponding excel file for testing
    and includes the sample_id and embryo_id.

    Raises:
    AssertionError: If the actual output does not match the expected output.
    """
    test_args = setup_test_data_from_excel()
    sample_id, embryo_id = name.rsplit("_", 1)
    (basher_input_namespace, error_dictionary, input_ok_flag) = excel_parser_main(
        test_args[sample_id]
    )

    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = snp_haplotype_main(basher_input_namespace)

    json_file_path = "test_data/embryo_validation_data.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        embryo_validation = json.loads(j.read())
        for dict in embryo_validation:
            all_validation[dict["sample_id"] + "_" + dict["embryo_id"]] = dict

    validate_embryo_results(
        mode_of_inheritance,
        sample_id,
        embryo_id,
        embryo_count_data_df,
        all_validation,
    )


@pytest.mark.parametrize("name", setup_test_data_from_excel(True, version=2))
def test_embryo_categorization_excel_v2(name):
    """
    Tests the categorization of embryos based on the Single Nucleotide Polymorphisms (SNPs)
    data provided in the excel files (version 2).

    The function parameters are parametrized using pytest's mark.parametrize decorator. For each
    test case, the function parses the excel file, processes the SNPs, and then compares the
    embryo categorization results with the expected results stored in a JSON file.

    Parameters:
    name (str): The name of the test case. The name is used to fetch the corresponding excel file for testing
    and includes the sample_id and embryo_id.

    Raises:
    AssertionError: If the actual output does not match the expected output.
    """
    test_args = setup_test_data_from_excel()
    sample_id, embryo_id = name.rsplit("_", 1)
    (basher_input_namespace, error_dictionary, input_ok_flag) = excel_parser_main(
        test_args[sample_id]
    )

    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = snp_haplotype_main(basher_input_namespace)

    json_file_path = "test_data/embryo_validation_data.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        embryo_validation = json.loads(j.read())
        for dict in embryo_validation:
            all_validation[dict["sample_id"] + "_" + dict["embryo_id"]] = dict

    validate_embryo_results(
        mode_of_inheritance,
        sample_id,
        embryo_id,
        embryo_count_data_df,
        all_validation,
    )
