import argparse

# Import csv with sample meta data
parser = argparse.ArgumentParser(
    description="Import SNP Haplotying meta data from Excel spreadsheet"
)

parser.add_argument(
    "-s",
    "--sample_sheet",
    help="Excel input file containing sample meta data",
    type=argparse.FileType("r"),
)

# Manually assign samples to categories
reference = "E08_P2C61.rhchp"
unaffected_partner = "C08_P2C59.rhchp"
affected_partner = "D08_P2C60.rhchp"
reference_type = "ref_affected"
embryo_samples = [
    "F08_P2C62.rhchp",
    "G08_P2C63.rhchp",
    "H08_P2C64.rhchp",
]
mode_of_inheritance = "placeholder"


def parse_experiment_design(metadata_csv):
    """
    Parses CSV with the experimental data for the run
    """
    metadata_df = pd.read_excel(open(args.sample_sheet, "rb"), sheet_name="data_entry")

    # Read in variables from samplesheet
    mode_of_inheritance = ""
    reference = "E08_P2C61.rhchp"
    unaffected_partner = "C08_P2C59.rhchp"
    affected_partner = "D08_P2C60.rhchp"
    reference_type = "ref_affected"
    embryo_samples = [
        "F08_P2C62.rhchp",
        "G08_P2C63.rhchp",
        "H08_P2C64.rhchp",
    ]

    pass
