# The genome build used by the SNP array
genome_build = "GRCh38"

# Current version of BASHer used to generate the report
basher_version = "1.0.0"

# If this flag is set to True, the program will stream the results of the analysis in JSON
# format to the standard output intead of creating a HTML report. This is useful for
# debugging issues with the test suite as pytest captures this output during testing and
# compares it to the expected output

stream_results = False

# These flags can be used to prevent the program from running certain parts of the analysis
# for example, if you have validated specific modes of inheritance, you can set the flags
# to skip the analysis of other modes of inheritance
allow_autosomal_dominant_cases = True
allow_autosomal_recessive_cases = True
allow_x_linked_cases = True
allow_consanguineous_cases = True
allow_trio_only_analysis = False

# The following flag adds a warning to the report if the version of BASHer used to generate it
# is still in development
released_to_production = False
