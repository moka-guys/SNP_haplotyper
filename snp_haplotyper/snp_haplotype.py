import argparse
import json
import logging
import os
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum, auto
from io import IOBase
from pathlib import Path
from typing import Any, Dict, List, Optional

import config as config  # TODO add code to use this dependency
import numpy as np
import pandas as pd
import pdfkit
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from helper_functions import custom_order_generator, import_haplotype_data
from jinja2 import Environment, PackageLoader
from pydantic import (
    BaseModel,
    ValidationError,
    conlist,
    constr,
    fields,
    root_validator,
    validator,
)
from pydantic.dataclasses import dataclass

# Import mode of inheritance specific code

logger = logging.getLogger("BASHer_logger")

# Add the directory containing this script to the PYTHOPATH
sys.path.append(os.path.dirname(__file__))
mod_path = Path(__file__).parent

import config as config
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from exceptions import ArgumentInputError, InvalidParameterSelectedError
from FamilyDataClass import FamilyData
from ReportDataClass import ReportData, ReportGenerator
from SNPAnalysisClass import SNPAnalysis

# TODO Copy rsID from hover tap
# TODO Check telomeric/centromeric genes work with 2mb window (FHSD1 - D4Z4 repeat, PKD1)
# TODO Add support for no embryos (just TRIOs being run to check if enough informative SNPs)
# TODO Add ADO % to table


def header_to_dict(header_str: str) -> Dict[str, str]:
    """
    Converts a string of header_info into a dictionary
    Args:
        header_str (str): A string in the key=value pairs like "PRU=1234;Hospital No=1234;Biopsy No:111", where the keys will be the titles of the fields in the header
    Returns:
        dict: A dictionary of the header info with field titles as keys and values as values
    """
    if header_str is None:
        return {}
    else:
        d = dict(x.split("=") for x in header_str.split(";"))
        return d


class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums
    """

    def __init__(self, **kwargs):
        # Pop off the type value
        enum_type = kwargs.pop("type", None)

        # Ensure an Enum subclass is provided
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, Enum):
            raise TypeError("type must be an Enum when using EnumAction")

        # Generate choices from the Enum
        kwargs.setdefault("choices", tuple(e.value for e in enum_type))

        super(EnumAction, self).__init__(**kwargs)

        self._enum = enum_type

    def __call__(self, parser, namespace, values, option_string=None):
        # Convert value back into an Enum
        value = self._enum(values)
        setattr(namespace, self.dest, value)


# Import environment variables set by docker-compose
UPLOAD_FOLDER = os.getenv("UPLOAD_FOLDER")

# Import command line arguments (these can be automatically generated from the sample sheet using sample_sheet_reader.py)
parser = argparse.ArgumentParser(description="SNP Haplotying from SNP Array data")


def flanking_region_size_type(value):
    try:
        return FlankingRegions(int(value))
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid flanking region size")


# File input/output data
parser.add_argument(
    "-i",
    "--input_file",
    type=str,
    help="Input txt file containing SNP Array output",
)

parser.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    help="Output filename prefix",
)

parser.add_argument(
    "-f",
    "--output_folder",
    type=str,
    help="Output folder path",
)

# Patient data
parser.add_argument(
    "-m",
    "--mode_of_inheritance",
    type=InheritanceMode,
    action=EnumAction,
    help="The mode of inheritance",
)

parser.add_argument(
    "-mp",
    "--male_partner",
    type=str,
    help="ID in input table for male_partner",
)

parser.add_argument(
    "-mps",
    "--male_partner_status",
    type=Status,
    action=EnumAction,
    help="Status of male_partner",
)

parser.add_argument(
    "-fp",
    "--female_partner",
    type=str,
    help="ID in input table for male_partner",
)

parser.add_argument(
    "-fps",
    "--female_partner_status",
    type=Status,
    action=EnumAction,
    help="ID in input table for female_partner",
)

parser.add_argument(
    "-consang",
    "--consanguineous",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to indicate that partners are consanguineous",
)

parser.add_argument(
    "-r",
    "--reference",
    type=str,
    help="ID in input table for reference sample",
)

parser.add_argument(
    "-rs",
    "--reference_status",
    type=Status,
    action=EnumAction,
    help="Status of Reference",
)

parser.add_argument(
    "-rr",
    "--reference_relationship",
    type=Relationship,
    action=EnumAction,
    help="Reference relationship to pro-band",
)

parser.add_argument(
    "-rx",
    "--reference_sex",
    type=Sex,
    action=EnumAction,
    help="Reference sex",
)

parser.add_argument(
    "-e",
    "--embryo_ids",
    nargs="+",
    type=str,
    help="IDs of embryos in the input table",
)

parser.add_argument(
    "-es",
    "--embryo_sex",
    nargs="+",
    type=str,
    choices=[sex.value for sex in Sex],
    help="Embryo sex - must be in same order as embryo_ids (cannot be unknown for X-linked diseases)",
)

# Gene/ROI data
parser.add_argument(
    "-g",
    "--gene_symbol",
    type=str,
    help="Gene Symbol",
)

parser.add_argument(
    "-gs",
    "--gene_start",
    type=int,
    help="Gene Start genomic co-ordinate (1-based referencing)",
)

parser.add_argument(
    "-ge",
    "--gene_end",
    type=int,
    help="Gene End genomic co-ordinate (1-based referencing)",
)

parser.add_argument(
    "-c",
    "--chr",
    type=str,
    help="Chromosome of ROI/gene",
)

parser.add_argument(
    "--flanking_region_size",
    type=flanking_region_size_type,
    nargs="?",
    choices=list(FlankingRegions),
    default=FlankingRegions.FLANK_2MB,
    help="Size of the flanking region either side of the gene in Mb",
)

parser.add_argument(
    "--trio_only",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Flag to produce a preliminary report without looking at embryos, must be used if not embryo data is provided.",
)

parser.add_argument(
    "--header_info",
    type=str,
    required=True,
    help="Pass a string to populate the report header. A field will be created for each entry field_title=field_value separated by ';', for example 'PRU=1234;Hospital No=1234;Biopsy No=111' will produce 3 fields in the header with the titles PRU, Hospital No, and Biopsy No.",
)


def main(args):
    number_snps_to_import = import_haplotype_data(args.input_file).shape[0]

    logger.info(
        f"Number of SNPs to import from SNP Array File = {number_snps_to_import}."
    )

    def create_family_data_from_args(args):
        return FamilyData(
            mode_of_inheritance=args.mode_of_inheritance,
            male_partner=args.male_partner,
            male_partner_status=args.male_partner_status,
            female_partner=args.female_partner,
            female_partner_status=args.female_partner_status,
            consanguineous=args.consanguineous,
            reference=args.reference,
            reference_status=args.reference_status,
            reference_relationship=args.reference_relationship,
            reference_sex=args.reference_sex,
            embryo_ids=args.embryo_ids,
            embryo_sex=args.embryo_sex,
            gene_symbol=args.gene_symbol,
            gene_start=args.gene_start,
            gene_end=args.gene_end,
            chr=args.chr,
            flanking_region_size=args.flanking_region_size,
            trio_only=args.trio_only,
            CHaS_input_fields=list(import_haplotype_data(args.input_file).columns),
            report_header_info=header_to_dict(args.header_info),
            input_ChAS_filepath=args.input_file,
        )

    # Instantiate SNPAnalysis object
    snp_pipeline = SNPAnalysis(
        create_family_data_from_args(args), import_haplotype_data(args.input_file)
    )

    logger.info(f"SNP Analysis complete - getting ready to prepare report.")

    report_gen = ReportGenerator(snp_pipeline.report_data)

    html_string = report_gen.render("html")
    pdf_string = report_gen.render("pdf")

    # Save HTML report to file in output folder, including timestamp in filename

    timestr = datetime.now().strftime("%Y%m%d-%H%M%S")
    output_dir = args.output_folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(
        os.path.join(output_dir, args.output_prefix + "_" + timestr + ".html"),
        "w",
    ) as f:
        f.write(html_string)

    # Convert HTML report to PDF
    pdfkit.from_string(
        pdf_string,
        os.path.join(args.output_folder, args.output_prefix + "_" + timestr + ".pdf"),
    )

    return (
        args.mode_of_inheritance,
        args.output_prefix,
        snp_pipeline.snp_data.number_snps_imported,
        snp_pipeline.snp_data.pytest_format_snp_df,
        pd.DataFrame()
        if args.trio_only
        else snp_pipeline.pytest_format_embryo_df,  # TODO Change to pytest specific df
        html_string,
        pdf_string,
    )


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
    # TODO move report saving code to here and capture output from main function
