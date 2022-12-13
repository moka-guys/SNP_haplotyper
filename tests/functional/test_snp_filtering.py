from argparse import Namespace
from xml.dom import ValidationErr
from snp_haplotype import main
import pandas as pd
import pytest
import json


def setup_test_data(split_by_embryo=False):
    # Import launch.json and parse the file for the correct parameters for each of the samples
    json_file_path = ".vscode/launch.json"

    with open(json_file_path, "r") as j:
        launch_contents = json.loads(j.read())

    run_data_df = pd.json_normalize(launch_contents["configurations"])

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
            consanguineous=True
            if "consanguineous" in arg_dictionary
            else False,  # arg_dictionary["consanguineous"],
            testing=True,
            trio_only=False,
            header_info=json.loads(arg_dictionary["header_info"]),
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


@pytest.mark.parametrize("name", setup_test_data(False))
def test_informative_snps(capsys, name):
    args = setup_test_data()
    main(args[name])
    output = capsys.readouterr().out
    results = json.loads(output)[
        "informative_snp_data"
    ]  # Need to select "informative_snp_data" or "embryo_cat_json" as appropriate TODO move into func args
    json_file_path = "test_data/informative_snp_validation.json"

    all_validation = {}
    with open(json_file_path, "r") as j:
        snp_validation = json.loads(j.read())
        for dict in snp_validation:
            all_validation[dict["sample_id"]] = dict

    if results["mode"] == "autosomal_dominant":
        validation = all_validation[results["sample_id"]]
        assert results["mode"] == validation["mode"]
        assert results["sample_id"] == validation["sample_id"]
        assert results["num_snps"] == validation["num_snps"]
        assert results["info_snps_upstream_2mb"] == validation["info_snps_upstream_2mb"]
        assert results["info_snps_in_gene"] == validation["info_snps_in_gene"]
        assert (
            results["info_snps_downstream_2mb"]
            == validation["info_snps_downstream_2mb"]
        )
        assert results["total_info_snps"] == validation["total_info_snps"]
        assert (
            results["high_risk_snps_upstream_2mb"]
            == validation["high_risk_snps_upstream_2mb"]
        )
        assert (
            results["high_risk_snps_within_gene"] == validation["high_risk_within_gene"]
        )
        assert (
            results["high_risk_snps_downstream_2mb"]
            == validation["high_risk_snps_downstream_2mb"]
        )
        assert (
            results["low_risk_snps_upstream_2mb"]
            == validation["low_risk_snps_upstream_2mb"]
        )
        assert (
            results["low_risk_snps_within_gene"] == validation["low_risk_within_gene"]
        )
        assert (
            results["low_risk_snps_downstream_2mb"]
            == validation["low_risk_snps_downstream_2mb"]
        )
    elif results["mode"] == "x_linked":
        # TODO add check for male embryo informative SNPs
        validation = all_validation[results["sample_id"]]
        assert results["mode"] == validation["mode"]
        assert results["sample_id"] == validation["sample_id"]
        assert results["num_snps"] == validation["num_snps"]
        assert results["info_snps_upstream_2mb"] == validation["info_snps_upstream_2mb"]
        assert results["info_snps_in_gene"] == validation["info_snps_in_gene"]
        assert (
            results["info_snps_downstream_2mb"]
            == validation["info_snps_downstream_2mb"]
        )
        assert results["total_info_snps"] == validation["total_info_snps"]
        assert (
            results["high_risk_snps_upstream_2mb"]
            == validation["high_risk_snps_upstream_2mb"]
        )
        assert (
            results["high_risk_snps_within_gene"] == validation["high_risk_within_gene"]
        )
        assert (
            results["high_risk_snps_downstream_2mb"]
            == validation["high_risk_snps_downstream_2mb"]
        )
        assert (
            results["low_risk_snps_upstream_2mb"]
            == validation["low_risk_snps_upstream_2mb"]
        )
        assert (
            results["low_risk_snps_within_gene"] == validation["low_risk_within_gene"]
        )
        assert (
            results["low_risk_snps_downstream_2mb"]
            == validation["low_risk_snps_downstream_2mb"]
        )
    elif results["mode"] == "autosomal_recessive":
        validation = all_validation[results["sample_id"]]
        assert results["mode"] == validation["mode"]
        assert results["sample_id"] == validation["sample_id"]
        assert results["num_snps"] == validation["num_snps"]
        assert results["total_info_snps"] == validation["total_info_snps"]
        assert results["info_snps_upstream_2mb"] == validation["info_snps_upstream_2mb"]
        assert results["info_snps_in_gene"] == validation["info_snps_in_gene"]
        assert (
            results["info_snps_downstream_2mb"]
            == validation["info_snps_downstream_2mb"]
        )
        assert (
            results["high_risk_snps_upstream_2mb_from_female"]
            == validation["high_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            results["high_risk_within_gene_from_female"]
            == validation["high_risk_within_gene_from_female"]
        )
        assert (
            results["high_risk_snps_downstream_2mb_from_female"]
            == validation["high_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            results["low_risk_snps_upstream_2mb_from_female"]
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            results["low_risk_within_gene_from_female"]
            == validation["low_risk_within_gene_from_female"]
        )
        assert (
            results["low_risk_snps_downstream_2mb_from_female"]
            == validation["low_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            results["high_risk_snps_upstream_2mb_from_male"]
            == validation["high_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            results["high_risk_within_gene_from_male"]
            == validation["high_risk_within_gene_from_male"]
        )
        assert (
            results["high_risk_snps_downstream_2mb_from_male"]
            == validation["high_risk_snps_downstream_2mb_from_male"]
        )
        assert (
            results["low_risk_snps_upstream_2mb_from_male"]
            == validation["low_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            results["low_risk_within_gene_from_male"]
            == validation["low_risk_within_gene_from_male"]
        )
        assert (
            results["low_risk_snps_downstream_2mb_from_male"]
            == validation["low_risk_snps_downstream_2mb_from_male"]
        )


@pytest.mark.parametrize("name", setup_test_data(True))
def test_embryo_categorization(capsys, name):
    args = setup_test_data()
    embryo_id = name.rsplit("_", 1)[1]
    name = name.rsplit("_", 1)[0]

    main(args[name])
    output = capsys.readouterr().out
    results = json.loads(output)[
        "embryo_cat_json"
    ]  # Need to select "informative_snp_data" or "embryo_cat_json" as appropriate TODO move into func args
    json_file_path = "test_data/embryo_validation_data.json"

    all_validation = {}
    sample_validation = {}
    with open(json_file_path, "r") as j:
        embryo_validation = json.loads(j.read())
        for dict in embryo_validation:
            all_validation[dict["sample_id"] + "_" + dict["embryo_id"]] = dict

    results_dict = {}
    # TODO rename dict variable to something more meaningful
    # TODO refactor all tests to use x_linked format
    try:
        # Required for legacy method testing AR and AD
        for result_dictionary in results:
            results_dict[
                result_dictionary["embryo_id"].replace(".rhchp", "")
            ] = result_dictionary
    except TypeError:
        # Required for new method testing XL
        pass

    try:
        # Required for legacy method testing AR and AD
        result = results_dict[embryo_id]
    except KeyError:
        result = results[embryo_id]
        result["mode"] = results["mode"]
        result["sample_id"] = results["sample_id"]

    validation = all_validation[name.split("_", 2)[2] + "_" + embryo_id]

    if result["mode"] == "autosomal_dominant":
        assert result["mode"] == validation["mode"]
        assert result["sample_id"] == validation["sample_id"]
        assert (
            result["upstream_2mb_high_risk_snps"]
            == validation["high_risk_upstream_2mb"]
        )
        assert (
            result["within_gene_high_risk_snps"] == validation["high_risk_within_gene"]
        )
        assert (
            result["downstream_2mb_high_risk_snps"]
            == validation["high_risk_downstream_2mb"]
        )
        assert (
            result["upstream_2mb_low_risk_snps"] == validation["low_risk_upstream_2mb"]
        )
        assert result["within_gene_low_risk_snps"] == validation["low_risk_within_gene"]
        assert (
            result["downstream_2mb_low_risk_snps"]
            == validation["low_risk_downstream_2mb"]
        )
    elif result["mode"] == "x_linked":
        assert result["mode"] == validation["mode"]
        assert result["sample_id"] == validation["sample_id"]
        assert result["high_risk,upstream"] == validation["high_risk_upstream_2mb"]
        assert result["high_risk,within_gene"] == validation["high_risk_within_gene"]
        assert result["high_risk,downstream"] == validation["high_risk_downstream_2mb"]
        assert result["low_risk,upstream"] == validation["low_risk_upstream_2mb"]
        assert result["low_risk,within_gene"] == validation["low_risk_within_gene"]
        assert result["low_risk,downstream"] == validation["low_risk_downstream_2mb"]
    elif result["mode"] == "autosomal_recessive":
        assert result["mode"] == validation["mode"]
        assert result["sample_id"] == validation["sample_id"]
        assert (
            result["upstream_2mb_female_high_risk_snps"]
            == validation["high_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            result["within_gene_female_high_risk_snps"]
            == validation["high_risk_within_gene_from_female"]
        )
        assert (
            result["downstream_2mb_female_high_risk_snps"]
            == validation["high_risk_snps_downstream_2mb_from_female"]
        )
        assert (
            result["upstream_2mb_female_low_risk_snps"]
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            result["within_gene_female_low_risk_snps"]
            == validation["low_risk_within_gene_from_female"]
        )
        assert (
            result["upstream_2mb_female_low_risk_snps"]
            == validation["low_risk_snps_upstream_2mb_from_female"]
        )
        assert (
            result["upstream_2mb_male_high_risk_snps"]
            == validation["high_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            result["within_gene_male_high_risk_snps"]
            == validation["high_risk_within_gene_from_male"]
        )
        assert (
            result["downstream_2mb_male_high_risk_snps"]
            == validation["high_risk_snps_downstream_2mb_from_male"]
        )
        assert (
            result["upstream_2mb_male_low_risk_snps"]
            == validation["low_risk_snps_upstream_2mb_from_male"]
        )
        assert (
            result["within_gene_male_low_risk_snps"]
            == validation["low_risk_within_gene_from_male"]
        )
        assert (
            result["downstream_2mb_male_low_risk_snps"]
            == validation["low_risk_snps_downstream_2mb_from_male"]
        )
