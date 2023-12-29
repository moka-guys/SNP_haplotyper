import logging
import re


# Custom error handler which saves errors to a dictionary for feedback to user
class DictErrorHandler(logging.Handler):
    def __init__(self, error_dict):
        super().__init__()
        self.error_dict = error_dict

    def emit(self, record):
        if record.levelno == logging.ERROR:
            error_msg = self.format(record)
            if error_msg not in self.error_dict:
                self.error_dict[error_msg] = 0
            self.error_dict[error_msg] += 1


# TODO I originally had this script called using os - I now have it imported as a module  - so this code is not used to set the error_dict which is now done in app.py
# I'll keep this code for present as I need to check what the commandline version uses.  Remove as appropriate.

# Create an error dictionary
error_dict = {}

# Initialize the custom error handler
dict_error_handler = DictErrorHandler(error_dict)

# Configure the logger to use the custom handler
logger = logging.getLogger("BASHer_logger")
logger.setLevel(logging.ERROR)
logger.addHandler(dict_error_handler)


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

    input_ok_flag = True

    # read in input_file to get column names (used to check columns listed in excel sheet match those in files)
    # with open(input_file, "r") as f:
    #    column_names = f.readline().strip().split("\t")

    if input_namespace.mode_of_inheritance not in [
        "autosomal_dominant",
        "autosomal_recessive",
        "x_linked",
    ]:
        logger.error(
            f"Invalid Mode of Inheritance {input_namespace.mode_of_inheritance} entered as argument, should be 'Autosomal Dominant', 'Autosomal Recessive', or 'X_Linked'"
        )
        input_ok_flag = False

    if input_namespace.chr not in [
        "X",
        "Y",
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
        logger.error(
            f"Invalid Chromosome '{input_namespace.chr}' entered as an argument - must be one of 1-22, X, or Y."
        )
        input_ok_flag = False

    if input_namespace.mode_of_inheritance == "x_linked":
        if input_namespace.chr != "x":
            logger.error(
                f"Invalid Chromosome '{input_namespace.chr}' entered as an argument - must be x for X-Linked inheritance."
            )
            input_ok_flag = False

    # Check if consanguineous is boolean
    if not isinstance(input_namespace.consanguineous, bool):
        logger.error("Invalid consanguineous: must be True or False.")
        input_ok_flag = False

    # Check if trio_only is boolean
    if not isinstance(input_namespace.trio_only, bool):
        logger.error("trio_only: must be True or False.")
        input_ok_flag = False

    # Check if flanking_region_size is "2mb" or "3mb"
    if input_namespace.flanking_region_size.value not in [2, 3]:
        logger.error("Invalid flanking_region_size: must be either '2mb' or '3mb'.")
        input_ok_flag = False

    # Check if gene_symbol is a non-empty string
    if not input_namespace.gene_symbol or input_namespace.gene_symbol.strip() == "":
        logger.error("Invalid gene_symbol: must be a non-empty string.")
        input_ok_flag = False

    # Check if input_file is a non-empty string and ends with csv or txt
    if not input_namespace.input_file or input_namespace.input_file.strip() == "":
        logger.error("Invalid input_file: must be a non-empty string.")
        input_ok_flag = False
    elif not (
        input_namespace.input_file.endswith(".csv")
        or input_namespace.input_file.endswith(".txt")
    ):
        logger.error(
            "Invalid input_file: must have a file extension of '.csv' or '.txt'."
        )
        input_ok_flag = False

    # # Check if ref_relationship is a non-empty string
    # if (
    #     not input_namespace.reference_relationship
    #     or input_namespace.reference_relationship.strip() == ""
    # ):
    #     logger.error("Invalid reference_relationship: must be a non-empty string.")
    #     input_ok_flag = False

    # Check if reference_status is a non-empty string
    # if (
    #     not input_namespace.reference_status
    #     or input_namespace.reference_status.strip() == ""
    # ):
    #     logger.error("Invalid reference_status: must be a non-empty string.")
    #     input_ok_flag = False

    # Sanity check on genomic coordinates
    if (
        not isinstance(input_namespace.gene_start, int)
        or input_namespace.gene_start < 0
    ):
        logger.error("Gene_start must be a non-negative integer.")
        input_ok_flag = False

    if not isinstance(input_namespace.gene_end, int) or input_namespace.gene_end < 0:
        logger.error("Gene_end must be a non-negative integer.")
        input_ok_flag = False

    if input_namespace.gene_start > input_namespace.gene_end:
        logger.error(
            f"Gene: {input_namespace.gene} gene_start {input_namespace.gene_start} is greater than gene_end '{input_namespace.gene_end}'"
        )
        input_ok_flag = False

    # Check if the length of embryo_ids and embryo_sex lists are the same
    if len(input_namespace.embryo_ids) != len(input_namespace.embryo_sex):
        logger.error("The lengths of embryo_ids and embryo_sex lists must be the same.")
        input_ok_flag = False

    # Check if embryo_sex list contains only allowable values
    allowed_embryo_sex_values = ["male", "female", "unknown"]
    for sex in input_namespace.embryo_sex:
        if sex.lower() not in allowed_embryo_sex_values:
            logger.error(
                f"Invalid value '{sex}' in embryo_sex list. Allowed values are 'male', 'female', and 'unknown'."
            )
            input_ok_flag = False

    allowable_values = {
        "x_linked": {
            "reference_sex": {"female", "male"},
            "reference_status": {"carrier", "affected", "unaffected"},
            "reference_relationship": {
                "father",
                "mother",
                "son",
                "daughter",
                "child",
                "embryo",
                "prenatal",
            },
        },
        "autosomal_dominant": {
            "reference_sex": {"male", "female", "unknown"},
            "reference_status": {"affected", "unaffected"},
            "reference_relationship": {
                "grandparent",
                "father",
                "mother",
                "son",
                "daughter",
                "child",
                "embryo",
                "prenatal",
                "both",
            },
        },
        "autosomal_recessive": {
            "reference_sex": {"male", "female", "unknown"},
            "reference_status": {"affected", "unaffected"},
            "reference_relationship": {
                "son",
                "daughter",
                "child",
                "embryo",
                "prenatal",
            },
        },
    }

    # # Check if reference_status is an allowable value based on the mode_of_inheritance
    # if (
    #     input_namespace.reference_status.lower()
    #     not in allowable_values[input_namespace.mode_of_inheritance]["reference_status"]
    # ):
    #     logger.error(
    #         f"Invalid reference_status '{input_namespace.reference_status}' for mode_of_inheritance '{input_namespace.mode_of_inheritance}'"
    #     )
    #     input_ok_flag = False

    # # Check if reference_relationship is an allowable value based on the mode_of_inheritance
    # if (
    #     input_namespace.reference_relationship.lower()
    #     not in allowable_values[input_namespace.mode_of_inheritance][
    #         "reference_relationship"
    #     ]
    # ):
    #     logger.error(
    #         f"Invalid reference_relationship '{input_namespace.reference_relationship}' for mode_of_inheritance '{input_namespace.mode_of_inheritance}'"
    #     )
    #     input_ok_flag = False

    # Define a dictionary with the allowable partner statuses for each mode of inheritance
    allowable_partner_statuses = {
        "x_linked": {
            "male_partner_status": ["unaffected"],
            "female_partner_status": ["carrier"],
        },
        "autosomal_dominant": {
            "male_partner_status": ["affected", "unaffected"],
            "female_partner_status": ["affected", "unaffected"],
        },
        "autosomal_recessive": {
            "male_partner_status": ["carrier"],
            "female_partner_status": ["carrier"],
        },
    }

    # Get the allowable statuses for the current mode of inheritance
    allowable_statuses = allowable_partner_statuses.get(
        input_namespace.mode_of_inheritance, {}
    )

    # # Check if the male_partner_status is allowable for the current mode of inheritance
    # if input_namespace.male_partner_status not in allowable_statuses.get(
    #     "male_partner_status", []
    # ):
    #     raise ValueError(
    #         f"Invalid male_partner_status '{input_namespace.male_partner_status}' for mode_of_inheritance '{input_namespace.mode_of_inheritance}'."
    #     )

    # # Check if the female_partner_status is allowable for the current mode of inheritance
    # if input_namespace.female_partner_status not in allowable_statuses.get(
    #     "female_partner_status", []
    # ):
    #     raise ValueError(
    #         f"Invalid female_partner_status '{input_namespace.female_partner_status}' for mode_of_inheritance '{input_namespace.mode_of_inheritance}'."
    #     )

    # # Check if both partners are "unaffected" in autosomal_dominant cases
    # if (
    #     input_namespace.mode_of_inheritance == "autosomal_dominant"
    #     and input_namespace.male_partner_status == "unaffected_partner"
    #     and input_namespace.female_partner_status == "unaffected_partner"
    # ):
    #     raise ValueError(
    #         "In autosomal_dominant cases, both partners cannot be 'unaffected'."
    #     )

    return error_dict, input_ok_flag
