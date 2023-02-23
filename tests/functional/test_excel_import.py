from argparse import Namespace
from xml.dom import ValidationErr
from snp_haplotype import main, header_to_dict
import os
import pandas as pd
import pytest
import json
from validate_output import validate_snp_results, validate_embryo_results

from excel_parser import parse_excel_input


def setup_test_data_from_excel(split_by_embryo=False):
    # Import launch.json and parse the file for the correct parameters for each of the samples
    excel_folder_path = "test_data/template_test_data"

    # get list of excel files in folder
    excel_files = [f for f in os.listdir(excel_folder_path) if f.endswith(".xlsx")]

    # iterate through excel files calling the parser function on each file, adding the output to a dictionary, with the file name as the key
    run_data_dict = {}
    for excel_file in excel_files:
        run_data_dict[excel_file] = parse_excel_input(
            os.path.join(excel_folder_path, excel_file), False
        )

    # Output dictionaries
    run_data_dictionary = {}
    run_data_per_embryo_dictionary = {}

    # The list returned has elements representing flag arguments (starting --) and their values,
    # which may be a single value or a list of values. We need to amalgamate these values into
    # a single list of values for each flag argument.
    for run_name, arg_dictionary in run_data_dict.items():
        run_name = (
            run_name.replace("excel_test_", "")
            .replace(".xlsx", "")
            .replace("Autosomal_Dominant_", "")
            .replace("Autosomal_Recessive_", "")
            .replace("X_linked_", "")
        )
        embryo_ids_list = arg_dictionary["embryo_data"]["embryo_column_name"].to_list()
        embryo_sex_list = arg_dictionary["embryo_data"]["embryo_sex"].to_list()

        args = Namespace(
            input_file=arg_dictionary["input_file"],
            output_prefix=run_name,
            mode_of_inheritance=arg_dictionary["mode_of_inheritance"],
            male_partner=arg_dictionary["male_partner_col"],
            male_partner_status=arg_dictionary["male_partner_status"],
            female_partner=arg_dictionary["female_partner_col"],
            female_partner_status=arg_dictionary["female_partner_status"],
            reference=arg_dictionary["reference_column_name"],
            reference_status=arg_dictionary["ref_status"],
            reference_relationship=arg_dictionary["ref_relationship"],
            embryo_ids=embryo_ids_list,
            embryo_sex=embryo_sex_list,
            gene_symbol=arg_dictionary["gene"],
            gene_start=int(arg_dictionary["gene_start"]),
            gene_end=int(arg_dictionary["gene_end"]),
            chr=arg_dictionary["chromosome"],
            consanguineous=True
            if arg_dictionary["consanguineous"] == "yes"
            else False,  # arg_dictionary["consanguineous"],
            testing=True,
            trio_only=False,
            header_info=f'PRU={arg_dictionary["pru"]};Hospital No={arg_dictionary["female_partner_hosp_num"]};Biopsy No={arg_dictionary["biopsy_number"]}',
        )

        run_data_dictionary[run_name] = args

        # This dictionary needs to be ammended when passed to the function testing the embryo categorising.
        # We require an entry for each embryo so the mark.parametrize function can create the tests correctly.
        for embryo_id in embryo_ids_list:
            run_data_per_embryo_dictionary[
                run_name + "_" + embryo_id.replace(".rhchp", "")
            ] = args

    if split_by_embryo:
        return run_data_per_embryo_dictionary
    else:
        return run_data_dictionary


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
    ) = main(args[name])

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
    args = setup_test_data_from_excel(True)
    embryo_id = name.rsplit("_", 1)[1]

    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
    ) = main(args[name])

    # Need to select "informative_snp_data" or "embryo_cat_json" as appropriate TODO move into func args
    json_file_path = "test_data/embryo_validation_data.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        embryo_validation = json.loads(j.read())
        for dict in embryo_validation:
            all_validation[dict["sample_id"] + "_" + dict["embryo_id"]] = dict

    validate_embryo_results(
        mode_of_inheritance,
        name.rsplit("_", 1)[0],
        embryo_id,
        embryo_count_data_df,
        all_validation,
    )
