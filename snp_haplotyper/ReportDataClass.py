from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from datetime import datetime
from typing import IO, Any, Union

import pandas as pd
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from helper_functions import (
    custom_order_generator,
    dict2html,
    format_plot_html_str,
    generate_html_plot,
    generate_pdf_plot,
    generate_plots,
    produce_html_table,
)
from jinja2 import Environment, PackageLoader


@dataclass
class ReportData:
    html_text_for_plots: str
    pdf_text_for_plots: str
    text_for_plots: str
    warning: str
    header_html: str
    mode_of_inheritance: str
    gene_symbol: str
    chromosome: str
    gene_start: str
    gene_end: str
    genome_build: str
    basher_version: str
    input_file: str
    male_partner: Any
    male_partner_status: Any
    female_partner: Any
    female_partner_status: Any
    reference: Any
    reference_status: Any
    reference_relationship: Any
    results_table_1: str
    qc_table: str
    summary_snps_table: str
    summary_embryo_table: str
    summary_embryo_by_region_table: str
    consanguinity_flag: str
    report_date: str = datetime.today().strftime("%Y-%m-%d %H:%M:%S")

    def __post_init__(self):
        # Format gene_start and gene_end with 1000s comma separator
        # self.gene_start = f"{int(self.gene_start):,}" # TODO check if already comma separated
        # self.gene_end = f"{int(self.gene_end):,}" # TODO check if already comma separated
        # Check if input_file is a file object or a string
        if isinstance(self.input_file, IO):
            self.input_file = self.input_file.name


class ReportGenerator:
    def __init__(self, report_data: ReportData):
        # Initialize with the provided ReportData
        self.report_data = report_data
        # Load the Jinja2 template
        env = Environment(loader=PackageLoader("snp_haplotype", "templates"))
        self.template = env.get_template("report_template.html")

    def render(self, file_type: str) -> str:
        if file_type == "html":
            # HTML reports have dynamic plots
            self.report_data.text_for_plots = self.report_data.html_text_for_plots
            return self.template.render(asdict(self.report_data))
        elif file_type == "pdf":
            # PDF reports have static plots
            self.report_data.text_for_plots = self.report_data.pdf_text_for_plots
            return self.template.render(asdict(self.report_data))
        else:
            raise ValueError(f"Unsupported file_type: {file_type}")
