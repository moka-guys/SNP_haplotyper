"""
Module contains configuration settings for the SNP haplotyper program.
"""

# NOTE: CHECK DEPLOYMENT DOCS FOR MORE INFO ON CONFIGURATION ON TRUST

# The genome build used by the SNP array
GENOME_BUILD = "GRCh38"

# Current version of BASHer used to generate the report
BASHER_VERSION = "2.0.0"

# These flags can be used to prevent the program from running certain parts of the analysis
# for example, if you have validated specific modes of inheritance, you can set the flags
# to skip the analysis of other modes of inheritance
ALLOW_AUTOSOMAL_DOMINANT_CASES = True
ALLOW_AUTOSOMAL_RECESSIVE_CASES = True
ALLOW_AUTOSOMAL_X_LINKED_CASES = True
ALLOW_CONSANGUINEOUS_CASES = True
ALLOW_TRIO_ONLY_ANALYSIS = True

# The following flag adds a warning to the report if the version of BASHer used to generate it
# is still in development
RELEASED_TO_PRODUCTION = False


# Filepaths used by the excel_parser.py script
OUTPUT_FOLDER = "/home/graeme/Desktop/SNP_haplotyper/output"

INPUT_FOLDER = "/home/graeme/Desktop/SNP_haplotyper/"

# File containing details of the probesets on the array (currently Thermo Fisher Scientific HT-CMA_96.r3 SNP array)
PROBESETS_MAPPING_FILE = "../test_data/AffyID2rsid.txt"
FIGURE_NUM = 5
