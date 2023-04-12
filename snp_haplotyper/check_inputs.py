import logging
import re

logger = logging.getLogger("BASHer_logger")


def check_input(input_namespace, input_file):
    """Check the input arguments for the script.

    Args:
        input_namespace (argparse.Namespace): Input arguments.
        input_files (list): List of input files.
        column_names (list): List of column names from input files.

    Returns:
        bool: True if input arguments are valid, False otherwise.

    Raises:
        ValueError: If input arguments are invalid.
    """

    return_value = True

    # read in input_file to get column names (used to check columns listed in excel sheet match those in files)
    with open(input_file, "r") as f:
        column_names = f.readline().strip().split("\t")

    if input_namespace.mode_of_inheritance not in [
        "autosomal_dominant",
        "autosomal_recessive",
        "x_linked",
    ]:
        logging.error(
            f"Invalid Mode of Inheritance {input_namespace.mode_of_inheritance} entered as argument, should be 'Autosomal Dominant', 'Autosomal Recessive', or 'X-Linked'"
        )
        return_value = False

    if input_namespace.chr not in [
        "x",
        "y",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
    ]:
        logging.error(
            f"Invalid Chromosome '{input_namespace.chr}' entered as an argument - must be one of 1-22, X, or Y."
        )
        return_value = False

    if input_namespace.mode_of_inheritance == "x-Linked":
        if input_namespace.chr != "x":
            logging.error(
                f"Invalid Chromosome '{input_namespace.chr}' entered as an argument - must be x for X-Linked inheritance."
            )
            return_value = False

    # TODO a list
    # if input_namespace.input_file.endswith(".txt"):
    #     logging.error(
    #         f"Invalid Chromosome '{input_namespace.chr}' entered as an argument - must be one of 1-22, X, or Y."
    #     )
    #     return_value = False

    # Check if provided column names match those in the input file.
    # for key, column in {
    #     "partner1_column_name": input_namespace.partner1_column_name,
    #     "partner2_column_name": input_namespace.partner2_column_name,
    #     "reference_column_name": input_namespace.reference_column_name,
    # }:
    #     if column not in column_names:
    #         logging.error(
    #             f"{key} '{column}' not found in input file columns {column_names}."
    #         )
    #         return_value = False

    # if input_namespace.mode_of_inheritance == "Autosomal Dominant":
    #     for key, value in {
    #         "partner1_status": input_namespace.partner1_status,
    #         "partner2_status": input_namespace.partner2_status,
    #         "ref_status": input_namespace.ref_status,
    #     }:
    #         if value not in ["Affected", "Unaffected"]:
    #             logging.error(
    #                 f"{key} '{value}' not valid, must be either 'Affected' or 'Unaffected'."
    #             )
    #             return_value = False
    #     if input_namespace.partner1_status == input_namespace.partner2_status:
    #         logging.error(
    #             f"partner1_status '{input_namespace.partner1_status}' and partner2_status '{input_namespace.partner2_status}' cannot be the same as if they are this negates the need for testing"
    #         )
    #         return_value = False

    # elif input_namespace.mode_of_inheritance == "Autosomal Recessive":
    #     input_namespace.ref_status
    # elif input_namespace.mode_of_inheritance == "X-Linked":
    #     input_namespace.ref_status

    # # Check partners sex are valid.
    # for key, value in {
    #     "partner1_sex": input_namespace.partner1_sex,
    #     "partner2_sex": input_namespace.partner2_sex,
    # }:
    #     if value.title() not in ["Male", "Female"]:
    #         logging.error(f"{key} '{value}' must be either 'Male' or 'Female'.")
    #         return_value = False

    # Check reference sex is valid

    # "reference_sex" : input_namespace.reference_sex

    # Sanity check on genomic co-ordinates
    if input_namespace.gene_start > input_namespace.gene_end:
        logging.error(
            f"Gene: {input_namespace.gene} gene_start {input_namespace.gene_start} is greater than gene_end '{input_namespace.gene_end}'"
        )
        return_value = False

    # version_pattern = re.compile("^(v|V)[0-9]+\.[0-9]+\.[0-9]+$")

    # check if input_namespace contains the argument input_namespace:
    # if hasattr(input_namespace, "template_version"):
    #     # Check semantic version has been provided for the input spreadsheet
    #     if version_pattern.match(input_namespace.template_version) == False:
    #         logging.error(
    #             f"Invalid template version '{input_namespace.template_version}' entered as an argument - must be in the format v1.0.0"
    #         )
    #         return_value = False
    # else:
    #     logging.error(
    #         f"template_version not provided as an argument in excel template - must be in the format v1.0.0"
    #     )
    #     return_value = False

    # # check flanking_region_size is '2mb' or '3mb'
    # if input_namespace.flanking_region_size not in ["2mb", "3mb"]:
    #     logging.error(
    #         f"Invalid flanking_region_size '{input_namespace.flanking_region_size}' entered as an argument - must be either '2mb' or '3mb'"
    #     )

    # ref_relationship
    # ref_relationship_to_couple

    # embryo_id
    # embryo_sex
    # embryo_column_name

    # Check if any header information out of Pru, female_partner_hosp_num, and biopsy_no are missing:

    # consanguineous
    # de_novo
    # template_version

    # exclusion
    # multi_analysis
    # partner1_type
    # partner2_type

    # print other metadata which is not used in the analysis
    # for arg in vars(input_namespace):
    #     print(arg, getattr(input_namespace, arg))

    return return_value
