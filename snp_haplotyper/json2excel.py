import argparse
from openpyxl import load_workbook
from openpyxl.utils import get_column_interval
import pandas as pd
import json

# Import command line arguments
parser = argparse.ArgumentParser(
    description="Populate excel spreadsheet template with test data from JSON file"
)

# Template file for creating excel file
parser.add_argument(
    "-t",
    "--excel_template",
    type=str,
    help="Excel template file for entering SNP Array meta data",
)

# JSON file containing test data
parser.add_argument(
    "-j",
    "--json_file",
    type=str,
    help="JSON file containing SNP Array parameters for test data",
)

args = parser.parse_args()


def defined_name_to_cell_location(workbook):
    """
    Convert defined name to cell location

    Defined names are in the format data_entry!$B$31 or data_entry!$F$22:$L$22 (merged cells)

    args:   defined_name (str) - Defined name in the format data_entry!$B$31 or data_entry!$F$22:$L$22 (merged cells)
            workbook (openpyxl.workbook.workbook.Workbook) - Excel workbook

    return: cell_location (str) - Cell location in the format B31 or F22:L22
    """
    wb = workbook
    location_lookup = {}
    for dn in wb.defined_names.definedName:
        cell_location = dn.attr_text.split("!")[1].replace("$", "")
        location_lookup[dn.name] = cell_location

    return location_lookup


def read_launch_json(json_file_path):

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

            args = {
                "input_file": arg_dictionary["input_file"],
                "output_prefix": arg_dictionary["output_prefix"],
                "mode_of_inheritance": arg_dictionary["mode_of_inheritance"],
                "male_partner": arg_dictionary["male_partner"],
                "male_partner_status": arg_dictionary["male_partner_status"],
                "female_partner": arg_dictionary["female_partner"],
                "female_partner_status": arg_dictionary["female_partner_status"],
                "reference": arg_dictionary["reference"],
                "reference_status": arg_dictionary["reference_status"],
                "reference_relationship": arg_dictionary["reference_relationship"],
                "embryo_ids": embryo_ids_list,
                "embryo_sex": embryo_sex_list,
                "gene_symbol": arg_dictionary["gene_symbol"],
                "gene_start": int(arg_dictionary["gene_start"]),
                "gene_end": int(arg_dictionary["gene_end"]),
                "chr": arg_dictionary["chr"],
                "consanguineous": "No",  # arg_dictionary["consanguineous"],
            }
        run_data_dictionary[run_name] = args

    return run_data_dictionary


def write_data_to_excel(run_data_dict, template_file_path):
    """
    Writes the following defined cells/ranges from the provided excel file:
    biopsy_number
    chromosome
    consanguineous
    embryo_data - range of cells containing embryo data:
        biopsy_no
        embryo_id
        embryo_sex
        embryo_column_name
    gene
    gene_end
    gene_start
    input_file
    mode_of_inheritance
    partner1_details - range of cells containing partner1 details:
        partner1_type
        partner1_sex
        partner1_column_name
    partner2_details - range of cells containing partner2 details:
        partner2_type
        partner2_sex
        partner2_column_name
    paste_gene
    reference
    ref_relationship
    ref_relationship_to_couple
    ref_seq
    ref_status
    """

    wb = load_workbook(
        filename=template_file_path,
        keep_vba=True,
        data_only=True,
        keep_links=True,
    )

    data_entry_sheet = wb["data_entry"]  #

    cell_lookup = defined_name_to_cell_location(wb)

    # Iterate over the run_data_dict and write the data to the excel file
    for run_dict in run_data_dict:
        argument_dict = run_data_dict[run_dict]
        data_entry_sheet[
            "B9"
        ] = "12345"  # data_entry_sheet["biopsy_number"] = "12345" - not recognizing key for some reason
        data_entry_sheet[cell_lookup["chromosome"]].value = argument_dict["chr"]
        data_entry_sheet[cell_lookup["consanguineous"]].value = argument_dict[
            "consanguineous"
        ]
        data_entry_sheet[cell_lookup["female_partner_hosp_num"]].value = "12345"
        data_entry_sheet[cell_lookup["gene"]].value = argument_dict["gene_symbol"]
        data_entry_sheet[cell_lookup["gene_end"]].value = argument_dict["gene_end"]
        data_entry_sheet[cell_lookup["gene_start"]].value = argument_dict["gene_start"]
        data_entry_sheet[cell_lookup["input_file"]].value = argument_dict["input_file"]
        data_entry_sheet[cell_lookup["mode_of_inheritance"]].value = argument_dict[
            "mode_of_inheritance"
        ]
        data_entry_sheet[
            cell_lookup["paste_gene"]
        ].value = f'{argument_dict["chr"]}:{argument_dict["gene_start"]}-{argument_dict["gene_end"]}'  # genomic range in format chr3:100000-200000 used to populate other fields
        data_entry_sheet[cell_lookup["pgd_worksheet"]].value = "12345"
        data_entry_sheet[cell_lookup["pru"]].value = "123456"
        data_entry_sheet["M8"].value = argument_dict["reference"]
        data_entry_sheet["H12"].value = argument_dict["reference_relationship"]
        data_entry_sheet["H11"].value = argument_dict["reference_status"]

        if argument_dict["mode_of_inheritance"] == "autosomal_dominant":
            if (
                argument_dict["male_partner_status"] == "affected"
                and argument_dict["female_partner_status"] == "unaffected"
            ):
                data_entry_sheet["G4"].value = "affected_partner"
                data_entry_sheet["M4"].value = argument_dict["male_partner"]
                data_entry_sheet["H4"].value = "Male"
                data_entry_sheet["G6"].value = "unaffected_partner"
                data_entry_sheet["M6"].value = argument_dict["female_partner"]
                data_entry_sheet["H6"].value = "Female"
            elif (
                argument_dict["male_partner_status"] == "unaffected"
                and argument_dict["female_partner_status"] == "affected"
            ):
                data_entry_sheet["G4"].value = "affected_partner"
                data_entry_sheet["M4"].value = argument_dict["female_partner"]
                data_entry_sheet["H4"].value = "Female"
                data_entry_sheet["G6"].value = "unaffected_partner"
                data_entry_sheet["M6"].value = argument_dict["male_partner"]
                data_entry_sheet["H6"].value = "Male"
        elif argument_dict["mode_of_inheritance"] == "autosomal_recessive":
            data_entry_sheet["G4"].value = "carrier_partner"
            data_entry_sheet["M4"].value = argument_dict["female_partner"]
            data_entry_sheet["H4"].value = "Female"
            data_entry_sheet["G6"].value = "carrier_partner"
            data_entry_sheet["M6"].value = argument_dict["male_partner"]
            data_entry_sheet["H6"].value = "Male"
        elif argument_dict["mode_of_inheritance"] == "x_linked":
            data_entry_sheet["G4"].value = "carrier_female_partner"
            data_entry_sheet["M4"].value = argument_dict["female_partner"]
            data_entry_sheet["H4"].value = "Female"
            data_entry_sheet["G6"].value = "unaffected_male_partner"
            data_entry_sheet["M6"].value = argument_dict["male_partner"]
            data_entry_sheet["H6"].value = "Male"
        else:
            pass

        # Clear embryo data from previous run, necessary because the number of embryos is variable
        # and is not guraranteed to be overwritten in the next step
        for count in range(4, 28):
            data_entry_sheet[f"R{count}"].value = ""
            data_entry_sheet[f"S{count}"].value = ""
            data_entry_sheet[f"T{count}"].value = ""
            data_entry_sheet[f"Q{count}"].value = ""

        # Populate embryo data
        for count, embryo_id in enumerate(argument_dict["embryo_ids"], start=4):
            data_entry_sheet[f"R{count}"].value = embryo_id
            data_entry_sheet[f"T{count}"].value = embryo_id
            data_entry_sheet[f"Q{count}"].value = "12345"

        for count, embryo_sex in enumerate(argument_dict["embryo_sex"], start=4):
            data_entry_sheet[f"S{count}"].value = embryo_sex

        wb.save(f"test_data/template_test_data/excel_test_{run_dict}.xlsx")


def main(args=None):
    if args is None:
        args = parser.parse_args()

    test_data_parameters = read_launch_json(args.json_file)
    write_data_to_excel(test_data_parameters, args.excel_template)


if __name__ == "__main__":
    main()
