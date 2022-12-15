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
    "--template_file",
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

    return run_data_dict


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

    # Get worksheet "data_entry" from template file
    data_entry_sheet = wb["data_entry"]

    # Iterate over the run_data_dict and write the data to the excel file
    for argument_dict in run_data_dict:
        data_entry_sheet["biopsy_number"] = "12345"
        data_entry_sheet["chromosome"] = argument_dict["chromosome"]
        data_entry_sheet["consanguineous"] = argument_dict["consanguineous"]
        data_entry_sheet["female_partner_hosp_num"] = "12345"
        data_entry_sheet["gene_symbol"] = argument_dict["gene_symbol"]
        data_entry_sheet["gene_end"] = argument_dict["gene_end"]
        data_entry_sheet["gene_omim"] = argument_dict["gene_omim"]
        data_entry_sheet["gene_start"] = argument_dict["gene_start"]
        data_entry_sheet["input_file"] = argument_dict["input_file"]
        data_entry_sheet["mode_of_inheritance"] = argument_dict["mode_of_inheritance"]
        data_entry_sheet[
            "paste_gene"
        ] = f'{argument_dict["chromosome"]}:{argument_dict["gene_start"]}-{argument_dict["gene_end"]}'  # genomic range in format chr3:100000-200000 used to populate other fields
        data_entry_sheet["pgd_worksheet"] = "12345"
        data_entry_sheet["pgd_worksheet_denovo"] = argument_dict["pgd_worksheet_denovo"]
        data_entry_sheet["pru"] = "123456"
        data_entry_sheet["reference"] = argument_dict["reference"]
        data_entry_sheet["ref_relationship"] = argument_dict["ref_relationship"]
        data_entry_sheet["ref_relationship_to_couple"] = argument_dict[
            "ref_relationship_to_couple"
        ]
        data_entry_sheet["ref_seq"] = argument_dict["ref_seq"]
        data_entry_sheet["ref_status"] = argument_dict["ref_status"]
        data_entry_sheet["template_version"] = argument_dict["template_version"]
        data_entry_sheet["embryo_data_df"] = argument_dict["embryo_data"]


def main(args=None):  # default argument allows pytest to override argparse for testing
    if args is None:
        args = parser.parse_args()

    test_data_parameters = read_launch_json(args.json_file)
    write_data_to_excel(test_data_parameters, args.template_file_path)


if __name__ == "__main__":
    main()
