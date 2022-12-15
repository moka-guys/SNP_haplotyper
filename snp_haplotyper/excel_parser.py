import argparse
from openpyxl import load_workbook
from openpyxl.utils import get_column_interval
import pandas as pd
import re
import subprocess
import config as config

# Import command line arguments (these can be automatically generated from the sample sheet using sample_sheet_reader.py)
parser = argparse.ArgumentParser(
    description="Parses meta data from provided excel file and runs SNP haplotyper"
)

# File input/output data
parser.add_argument(
    "-i",
    "--input_spreadsheet",
    type=str,
    help="Excel file containing SNP Array meta data",
)

args = parser.parse_args()


def load_workbook_range(range_string, worksheet):
    """
    Import rows in the range with text entered
    """
    col_start, col_end = re.findall("[A-Z]+", range_string)

    data_rows = []
    for row in worksheet[range_string]:
        data_rows.append([cell.value for cell in row])

    return pd.DataFrame(data_rows, columns=get_column_interval(col_start, col_end))


def parse_excel_input(input_file):
    """
    Imports the following defined cells/ranges from the provided excel file:
        biopsy_number
        chromosome
        consanguineous
        de_novo
        disease
        disease_omim
        embryo_data - range of cells containing embryo data:
            biopsy_no
            embryo_id
            embryo_sex
            embryo_column_name
        exclusion
        female_partner_hosp_num
        flanking_region_size
        gene
        gene_end
        gene_omim
        gene_start
        input_file
        maternal_mutation
        mode_of_inheritance
        multi_analysis
        partner1_details - range of cells containing partner1 details:
            partner1_type
            partner1_sex
            partner1_forename
            partner1_surname
            partner1_dob
            partner1_DNA_number
            partner1_column_name
        partner2_details - range of cells containing partner2 details:
            partner2_type
            partner2_sex
            partner2_forename
            partner2_surname
            partner2_dob
            partner2_DNA_number
            partner2_column_name
        paste_gene
        paternal_mutation
        pgd_worksheet
        pgd_worksheet_denovo
        pru
        reference
        ref_relationship
        ref_relationship_to_couple
        ref_seq
        ref_status
        template_version
    """
    wb = load_workbook(
        filename=input_spreadsheet,
        keep_vba=False,
        data_only=True,
        keep_links=True,
    )

    argument_dict = {}

    # Get list of defined ranges from provided excel sheet
    defined_ranges = wb.defined_names
    # TODO Add to logger when implemented
    #  print(defined_ranges)
    data_entry_sheet = wb["data_entry"]

    for dn in wb.defined_names.definedName:
        # TODO Add to logger when implemented
        # print(dn.name)
        # print(dn.attr_text)
        input_name = dn.name
        if input_name in ["embryo_data", "partner1_details", "partner2_details"]:
            # Import excel ranges
            df = load_workbook_range(
                dn.attr_text.split("!")[1].replace("$", ""), data_entry_sheet
            )
            argument_dict[input_name] = df.dropna(how="all")  # Remove empty rows
        else:
            # Process cell locations in the format data_entry!$B$31 or data_entry!$F$22:$L$22 (merged cells)
            cell_location = dn.attr_text.split(":")[0].split("!")[1].replace("$", "")
            cell_value = data_entry_sheet[cell_location].value
            argument_dict[input_name] = cell_value

    biopsy_number = argument_dict["biopsy_number"]
    chromosome = argument_dict["chromosome"]
    consanguineous = argument_dict["consanguineous"]
    de_novo = argument_dict["de_novo"]
    disease = argument_dict["disease"]
    disease_omim = argument_dict["disease_omim"]
    exclusion = argument_dict["exclusion"]
    female_partner_hosp_num = argument_dict["female_partner_hosp_num"]
    flanking_region_size = argument_dict["flanking_region_size"]
    gene_symbol = argument_dict["gene"]
    gene_end = argument_dict["gene_end"]
    gene_omim = argument_dict["gene_omim"]
    gene_start = argument_dict["gene_start"]
    input_file = argument_dict["input_file"]
    maternal_mutation = argument_dict["maternal_mutation"]
    mode_of_inheritance = argument_dict["mode_of_inheritance"]
    multi_analysis = argument_dict["multi_analysis"]
    paste_gene = argument_dict[
        "paste_gene"
    ]  # genomic range in format chr3:100000-200000 used to populate other fields
    paternal_mutation = argument_dict["paternal_mutation"]
    pgd_worksheet = argument_dict["pgd_worksheet"]
    pgd_worksheet_denovo = argument_dict["pgd_worksheet_denovo"]
    pru = argument_dict["pru"]
    reference = argument_dict["reference"]
    ref_relationship = argument_dict["ref_relationship"]
    ref_relationship_to_couple = argument_dict["ref_relationship_to_couple"]
    ref_seq = argument_dict["ref_seq"]
    ref_status = argument_dict["ref_status"]
    template_version = argument_dict["template_version"]

    embryo_data_df = argument_dict["embryo_data"]
    embryo_data_df.columns = [
        "biopsy_no",
        "embryo_id",
        "embryo_sex",
        "embryo_column_name",
    ]

    # Filter embryo data to only include embryos from the current biopsy
    filtered_embryo_data_df = embryo_data_df[
        embryo_data_df["biopsy_no"] == biopsy_number
    ]

    # Check if there are any embryos in the current biopsy
    if filtered_embryo_data_df.empty:
        trio_only = True
    else:
        trio_only = False

    (
        partner1_type,
        partner1_sex,
        partner1_forename,
        partner1_surname,
        partner1_dob,
        partner1_DNA_number,
        partner1_column_name,
    ) = argument_dict["partner1_details"].values.tolist()[0]

    (
        partner2_type,
        partner2_sex,
        partner2_forename,
        partner2_surname,
        partner2_dob,
        partner2_DNA_number,
        partner2_column_name,
    ) = argument_dict["partner2_details"].values.tolist()[0]

    # Clean imported data
    partner1_sex = partner1_sex.lower()
    partner2_sex = partner2_sex.lower()

    if partner1_sex == "male" and partner2_sex == "female":
        male_partner_status = partner1_type
        male_partner_col = partner1_column_name
        female_partner_status = partner2_type
        female_partner_col = partner2_column_name
    elif partner1_sex == "female" and partner2_sex == "male":
        male_partner_status = partner2_type
        male_partner_col = partner2_column_name
        female_partner_status = partner1_type
        female_partner_col = partner1_column_name

    output_prefix = input_file

    if mode_of_inheritance == "autosomal_dominant":
        female_partner_status = female_partner_status.split("_")[0]
        male_partner_status = male_partner_status.split("_")[0]
    elif mode_of_inheritance == "autosomal_recessive":
        female_partner_status = (
            "carrier"
            if female_partner_status == "carrier_partner"
            else female_partner_status
        )
        male_partner_status = (
            "carrier"
            if male_partner_status == "carrier_partner"
            else male_partner_status
        )
    elif mode_of_inheritance == "x_linked":
        female_partner_status = (
            "carrier"
            if female_partner_status == "carrier_female_partner"
            else female_partner_status
        )
        male_partner_status = (
            "unaffected"
            if male_partner_status == "unaffected_male_partner"
            else male_partner_status
        )

    cmd = (
        f" {config.python_location} -i {config.snp_haplotype_script}"
        f" --input_file {input_file} --output_folder {config.output_folder}"
        f" --output_prefix {output_prefix} --mode_of_inheritance {mode_of_inheritance}"
        f" --male_partner {male_partner_col} --male_partner_status {male_partner_status}"
        f" --female_partner {female_partner_col} --female_partner_status {female_partner_status}"
        f" --reference {reference} --reference_status {ref_status}"
        f" --reference_relationship {ref_relationship}"
        f" --gene_symbol {gene_symbol} --gene_start {gene_start}"
        f" --gene_end {gene_end} --chr {chromosome}"
    )

    if trio_only == False:
        cmd = cmd + (
            f" --embryo_ids {' '.join(filtered_embryo_data_df.embryo_column_name.to_list())}"
            f" --embryo_sex {' '.join(filtered_embryo_data_df.embryo_sex.to_list())}"
        )

    else:
        cmd = cmd + f" --trio_only"

    # Add header info to cmd string
    cmd = (
        cmd
        + f" --header 'PRU:{pru},Hospital No:{female_partner_hosp_num},Biopsy No:{biopsy_number}'"
    )

    # Run SNP haplotyping script. (if statement to allow for future development of mode dependent arguments, currently no difference)
    # TODO if not required remove if statement
    if mode_of_inheritance == "autosomal_dominant":
        subprocess.run(cmd, shell=True)
    elif mode_of_inheritance == "autosomal_recessive":
        subprocess.run(cmd, shell=True)
    elif mode_of_inheritance == "x_linked":
        subprocess.run(cmd, shell=True)


def main():
    parse_excel_input(args.input_spreadsheet)


if __name__ == "__main__":
    main()
