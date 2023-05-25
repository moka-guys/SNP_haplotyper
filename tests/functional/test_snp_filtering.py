from argparse import Namespace
from xml.dom import ValidationErr
from snp_haplotype import main, header_to_dict
import pandas as pd
import pytest
import json
from validate_output import validate_snp_results, validate_embryo_results


def setup_test_data(split_by_embryo=False):
    # Import launch.json and parse the file for the correct parameters for each of the samples
    json_file_path = ".vscode/launch.json"

    with open(json_file_path, "r") as j:
        launch_contents = json.loads(j.read())

    run_data_df = pd.json_normalize(launch_contents["configurations"])

    run_data_df = run_data_df.drop(
        ["pathMappings", "connect.host", "connect.port"], axis=1
    )

    # Run configurations are filtered to exclue docker specific configurations
    run_data_df = run_data_df.drop(
        run_data_df[run_data_df.name == "Python: Remote Attach"].index
    )

    # Run configurations are filtered by the "purpose" field to exclude configurations that are not intended for testing
    run_data_df = run_data_df[
        run_data_df["purpose"].str.contains("debug-test", regex=False)
    ]

    # The name and args columns are zipped into a dictionary for iterating over
    run_data_dict = dict(zip(run_data_df["name"], run_data_df["args"]))

    # Output dictionaries
    run_data_dictionary = {}
    run_data_per_embryo_dictionary = {}

    # The list returned has elements representing flag arguments (starting --) and their values,
    # which may be a single value or a list of values. We need to amalgamate these values into
    # a single list of values for each flag argument.
    for run_name, arg_list in run_data_dict.items():
        arg_dictionary = {}
        dict_values = []
        for i in arg_list:
            if i.startswith("-"):
                dict_values = []
                dict_key = i.strip("-")
            else:
                dict_values.append(i)
            if len(dict_values) == 1:
                arg_dictionary[dict_key] = dict_values[0]  # Don't return a list of 1
            else:
                arg_dictionary[dict_key] = dict_values

        if "trio_only" in arg_dictionary:
            args = Namespace(
                input_file=arg_dictionary["input_file"],
                output_prefix=arg_dictionary["output_prefix"],
                mode_of_inheritance=arg_dictionary["mode_of_inheritance"],
                male_partner=arg_dictionary["male_partner"],
                male_partner_status=arg_dictionary["male_partner_status"],
                female_partner=arg_dictionary["female_partner"],
                female_partner_status=arg_dictionary["female_partner_status"],
                reference=arg_dictionary["reference"],
                reference_status=arg_dictionary["reference_status"],
                reference_relationship=arg_dictionary["reference_relationship"],
                gene_symbol=arg_dictionary["gene_symbol"],
                gene_start=int(arg_dictionary["gene_start"]),
                gene_end=int(arg_dictionary["gene_end"]),
                chr=arg_dictionary["chr"],
                flanking_region_size=arg_dictionary["flanking_region_size"],
                consanguineous=True if "consanguineous" in arg_dictionary else False,
                testing=True,
                trio_only=True,
                header_info=header_to_dict(arg_dictionary["header_info"]),
            )
        else:
            # Ensure that embryo_ids are a list even if only a single embryo is present
            if isinstance(arg_dictionary["embryo_ids"], str):
                embryo_ids_list = [arg_dictionary["embryo_ids"]]
            else:
                embryo_ids_list = arg_dictionary["embryo_ids"]

            if isinstance(arg_dictionary["embryo_sex"], str):
                embryo_sex_list = [arg_dictionary["embryo_sex"]]
            else:
                embryo_sex_list = arg_dictionary["embryo_sex"]

            args = Namespace(
                input_file=arg_dictionary["input_file"],
                output_prefix=arg_dictionary["output_prefix"],
                mode_of_inheritance=arg_dictionary["mode_of_inheritance"],
                male_partner=arg_dictionary["male_partner"],
                male_partner_status=arg_dictionary["male_partner_status"],
                female_partner=arg_dictionary["female_partner"],
                female_partner_status=arg_dictionary["female_partner_status"],
                reference=arg_dictionary["reference"],
                reference_status=arg_dictionary["reference_status"],
                reference_relationship=arg_dictionary["reference_relationship"],
                embryo_ids=embryo_ids_list,
                embryo_sex=embryo_sex_list,
                gene_symbol=arg_dictionary["gene_symbol"],
                gene_start=int(arg_dictionary["gene_start"]),
                gene_end=int(arg_dictionary["gene_end"]),
                chr=arg_dictionary["chr"],
                flanking_region_size=arg_dictionary["flanking_region_size"],
                consanguineous=True if "consanguineous" in arg_dictionary else False,
                testing=True,
                trio_only=False,
                header_info=header_to_dict(arg_dictionary["header_info"]),
            )

        run_data_dictionary[run_name] = args
        # This dictionary needs to be ammended when passed to the function testing the embryo categorising.
        # We require an entry for each embryo so the mark.parametrize function can create the tests correctly.
        if args.trio_only == False:
            for embryo_id in embryo_ids_list:
                run_data_per_embryo_dictionary[
                    run_name + "_" + embryo_id.replace(".rhchp", "")
                ] = args

    if split_by_embryo:
        return run_data_per_embryo_dictionary
    else:
        return run_data_dictionary


@pytest.mark.parametrize("name", setup_test_data(False))
def test_informative_snps(name):
    args = setup_test_data()
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
    ) = main(args[name])
    json_file_path = "test_data/informative_snp_validation.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        snp_validation = json.loads(j.read())
        for dict in snp_validation:
            all_validation[dict["sample_id"]] = dict

    validate_snp_results(
        mode_of_inheritance,
        sample_id.replace("_precase", "").replace(
            "_with_embryo", ""
        ),  # Remove any suffixes from the sample_id which will not be present in the validation data
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        all_validation,
    )


@pytest.mark.parametrize("name", setup_test_data(True))
def test_embryo_categorization(name):
    args = setup_test_data()
    embryo_id = name.rsplit("_", 1)[1]
    name = name.rsplit("_", 1)[0]
    (
        mode_of_inheritance,
        sample_id,
        number_snps_imported,
        summary_snps_by_region,
        informative_snps_by_region,
        embryo_count_data_df,
        html_string,
        pdf_string,
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
        name.split("_", 2)[2],
        embryo_id,
        embryo_count_data_df,
        all_validation,
    )
