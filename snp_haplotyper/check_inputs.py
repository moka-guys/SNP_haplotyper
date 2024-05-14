import logging

from EnumDataClasses import InheritanceMode, Sex, Status


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
        InheritanceMode.AUTOSOMAL_DOMINANT,
        InheritanceMode.AUTOSOMAL_RECESSIVE,
        InheritanceMode.X_LINKED,
    ]:
        error_msg = (
            "Invalid Mode of Inheritance %s entered as argument, should be 'Autosomal Dominant', "
            "'Autosomal Recessive', or 'X-Linked'" % input_namespace.mode_of_inheritance
        )
        logger.error(error_msg)
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
            "Invalid Chromosome '%s' entered as an argument - must be one of 1-22, X, or Y.", input_namespace.chr
        )
        input_ok_flag = False

    if input_namespace.mode_of_inheritance == "x_linked":
        if input_namespace.chr != "x":
            error_msg = (
                "Invalid Chromosome '%s' entered as an argument - must be x for X-Linked inheritance."
                % input_namespace.chr
            )
            logger.error(error_msg)
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
        logger.error("Invalid flanking_region_size in mb: must be either  2 or 3.")
        input_ok_flag = False

    # Check if gene_symbol is a non-empty string
    if not input_namespace.gene_symbol or input_namespace.gene_symbol.strip() == "":
        logger.error("Invalid gene_symbol: must be a non-empty string.")
        input_ok_flag = False

    # Check if input_file is a non-empty string and ends with csv or txt
    if not input_namespace.input_file or input_namespace.input_file.strip() == "":
        logger.error("Invalid input_file: must be a non-empty string.")
        input_ok_flag = False
    elif not (input_namespace.input_file.endswith(".csv") or input_namespace.input_file.endswith(".txt")):
        logger.error("Invalid input_file: must have a file extension of '.csv' or '.txt'.")
        input_ok_flag = False

    # Sanity check on genomic coordinates
    if not isinstance(input_namespace.gene_start, int) or input_namespace.gene_start < 0:
        logger.error("Gene_start must be a non-negative integer.")
        input_ok_flag = False

    if not isinstance(input_namespace.gene_end, int) or input_namespace.gene_end < 0:
        logger.error("Gene_end must be a non-negative integer.")
        input_ok_flag = False

    if input_namespace.gene_start > input_namespace.gene_end:
        logger.error(
            "Gene: %s gene_start %s is greater than gene_end '%s'",
            input_namespace.gene,
            input_namespace.gene_start,
            input_namespace.gene_end,
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
                "Invalid value '%s' in embryo_sex list. Allowed values are 'male', 'female', and 'unknown'.", sex
            )
            input_ok_flag = False

    allowable_values = {
        InheritanceMode.X_LINKED: {
            "reference_sex": {Sex.FEMALE, Sex.MALE},
            "reference_status": {Status.CARRIER, Status.AFFECTED, Status.UNAFFECTED},
        },
        InheritanceMode.AUTOSOMAL_DOMINANT: {
            "reference_sex": {Sex.FEMALE, Sex.MALE, Sex.UNKNOWN},
            "reference_status": {Status.AFFECTED, Status.UNAFFECTED},
        },
        InheritanceMode.AUTOSOMAL_RECESSIVE: {
            "reference_sex": {Sex.FEMALE, Sex.MALE, Sex.UNKNOWN},
            "reference_status": {Status.AFFECTED, Status.UNAFFECTED},
        },
    }

    # Check if reference_status is an allowable value based on the mode_of_inheritance
    if (
        input_namespace.reference_status
        not in allowable_values[input_namespace.mode_of_inheritance]["reference_status"]
    ):
        logger.error(
            f"Invalid reference_status '{input_namespace.reference_status}' for mode_of_inheritance '{input_namespace.mode_of_inheritance}'"
        )
        input_ok_flag = False

    # Define a dictionary with the allowable partner statuses for each mode of inheritance
    allowable_partner_statuses = {
        InheritanceMode.X_LINKED: {
            "male_partner_status": [Status.UNAFFECTED],
            "female_partner_status": [Status.CARRIER],
        },
        InheritanceMode.AUTOSOMAL_DOMINANT: {
            "male_partner_status": [Status.AFFECTED, Status.UNAFFECTED],
            "female_partner_status": [Status.AFFECTED, Status.UNAFFECTED],
        },
        InheritanceMode.AUTOSOMAL_RECESSIVE: {
            "male_partner_status": [Status.CARRIER],
            "female_partner_status": [Status.CARRIER],
        },
    }

    # Get the allowable statuses for the current mode of inheritance
    allowable_statuses = allowable_partner_statuses.get(input_namespace.mode_of_inheritance, {})

    # Check if the male_partner_status is allowable for the current mode of inheritance
    if input_namespace.male_partner_status not in allowable_statuses.get("male_partner_status", []):
        raise ValueError(
            f"Invalid male_partner_status '{input_namespace.male_partner_status}' "
            f"for mode_of_inheritance '{input_namespace.mode_of_inheritance}'."
        )

    # Check if the female_partner_status is allowable for the current mode of inheritance
    if input_namespace.female_partner_status not in allowable_statuses.get("female_partner_status", []):
        raise ValueError(
            f"Invalid female_partner_status '{input_namespace.female_partner_status}' "
            f"for mode_of_inheritance '{input_namespace.mode_of_inheritance}'."
        )

    # Check if both partners are "unaffected" in autosomal_dominant cases
    if (
        input_namespace.mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT
        and input_namespace.male_partner_status == Status.UNAFFECTED
        and input_namespace.female_partner_status == Status.UNAFFECTED
    ):
        raise ValueError("In autosomal_dominant cases, both partners cannot be 'unaffected'.")

    return error_dict, input_ok_flag
