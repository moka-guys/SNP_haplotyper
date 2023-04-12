from argparse import Namespace
from xml.dom import ValidationErr
from snp_haplotype import main as snp_haplotype_main
import os
import pandas as pd
import pytest
import json
from validate_output import validate_snp_results, validate_embryo_results

from excel_parser import main as excel_parser_main
from test_snp_filtering import setup_test_data


def setup_test_data_from_excel(split_by_embryo=False):
    # We need to get the embryo_ids for each sample from the launch JSON so we can auto generate
    # the test names for each embryo without prereading the excel files
    embryo_data = setup_test_data(False)

    # get list of excel files in folder
    excel_folder_path = "test_data/template_test_data"
    excel_files = [f for f in os.listdir(excel_folder_path) if f.endswith(".xlsx")]

    # iterate through excel files creating a dictionary with run names as keys and the values are
    # file paths for each run
    run_data_dict = {}
    for excel_file in excel_files:
        run_name = excel_file.replace("excel_test_", "").replace(".xlsx", "")
        if split_by_embryo == False:
            run_data_dict[run_name] = Namespace(
                input_spreadsheet=os.path.join(
                    excel_folder_path,
                    excel_file,
                ),
                snp_array_file=None,
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
                )

    return run_data_dict


@pytest.mark.parametrize("name", setup_test_data_from_excel(False))
def test_informative_snps_excel(name):
    args = setup_test_data_from_excel()
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = excel_parser_main(args[name])

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


@pytest.mark.parametrize("name", setup_test_data_from_excel(True))
def test_embryo_categorization_excel(name):
    args = setup_test_data_from_excel()
    sample_id, embryo_id = name.rsplit("_", 1)
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = excel_parser_main(args[sample_id])

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
