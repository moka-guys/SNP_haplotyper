# NOTE: CHECK DEPLOYMENT DOCS FOR MORE INFO ON CONFIGURATION ON TRUST

# The genome build used by the SNP array
genome_build = "GRCh38"

# Current version of BASHer used to generate the report
basher_version = "1.0.0"

# These flags can be used to prevent the program from running certain parts of the analysis
# for example, if you have validated specific modes of inheritance, you can set the flags
# to skip the analysis of other modes of inheritance
allow_autosomal_dominant_cases = True
allow_autosomal_recessive_cases = True
allow_x_linked_cases = True
allow_consanguineous_cases = False
allow_trio_only_analysis = True

# The following flag adds a warning to the report if the version of BASHer used to generate it
# is still in development
released_to_production = False


# Filepaths used by the excel_parser.py script
output_folder = "/home/graeme/Desktop/SNP_haplotyper/output"

input_folder = "/home/graeme/Desktop/SNP_haplotyper/"
