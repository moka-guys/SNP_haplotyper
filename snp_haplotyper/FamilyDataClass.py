"""
FamilyData dataclass
=================

This class represents the family metada provided for the SNP haplotyping and contains relevant data validation of
the input.

"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from pydantic import BaseModel, root_validator

logger = logging.getLogger("BASHer_logger")


class FamilyData(BaseModel):
    """
    Represents detailed genetic and familial data necessary for SNP haplotyping and genetic risk assessment. This class
    encapsulates information about the mode of inheritance, partner statuses, consanguinity, reference data, and
    specific gene information. It includes extensive validation logic to ensure the integrity and consistency of the
    data provided for genetic analysis.

    Attributes:
        mode_of_inheritance (InheritanceMode): Specifies the genetic mode of inheritance (e.g., autosomal dominant,
                                               autosomal recessive, X-linked).
        male_partner (str): Identifier for the male partner, expected to end with ".rhchp".
        male_partner_status (Status): Genetic status of the male partner (e.g., affected, unaffected, carrier).
        female_partner (str): Identifier for the female partner, expected to end with ".rhchp".
        female_partner_status (Status): Genetic status of the female partner.
        consanguineous (bool): Indicates whether the family has a consanguineous background.
        reference (str): Reference identifier used in genetic analysis.
        reference_status (Status): Genetic status of the reference individual.
        reference_relationship (Relationship): Relationship of the reference individual to the family.
        reference_sex (Sex): Sex of the reference individual.
        embryo_ids (List[str]): List of embryo identifiers, each expected to end with ".rhchp".
        embryo_sex (List[Sex]): List of sexes corresponding to each embryo ID.
        gene_symbol (str): Symbol of the gene of interest.
        gene_start (int): Start position of the gene on the chromosome.
        gene_end (int): End position of the gene on the chromosome.
        chr (Chromosome): Chromosome on which the gene is located.
        flanking_region_size (FlankingRegions): Size of the flanking regions around the gene of interest.
        trio_only (bool): Indicates if the analysis is restricted to the trio (mother, father, child) only.
        CHaS_input_fields (List[str]): Specific fields required for CHaS input.
        report_header_info (Dict[str, str]): Information to be included in the report header.
        input_ChAS_filepath (Path): Filepath to the input CHaS file.
        expected_CHaS_file_name (Optional[str]): Expected filename for the CHaS input file, used for validation.

    The class includes several `@root_validator` methods to enforce constraints such as gene start and end positions,
    matching lengths of embryo IDs and sexes, unique reference identification, and proper file extensions for genetic
    data files. It also validates the appropriateness of reference statuses and relationships based on the mode of
    inheritance and ensures that partner statuses are consistent with the inheritance mode.

    Properties:
        unaffected_partner: Returns the identifier of the unaffected partner in autosomal dominant cases.
        affected_partner: Returns the identifier of the affected partner in autosomal dominant cases.
        get_partner1: Logic to determine the first partner based on the mode of inheritance.
        get_partner2: Logic to determine the second partner based on the mode of inheritance.
    """

    mode_of_inheritance: InheritanceMode
    male_partner: str
    male_partner_status: Status
    female_partner: str
    female_partner_status: Status
    consanguineous: bool
    reference: str
    reference_status: Status
    reference_relationship: Relationship
    reference_sex: Sex
    embryo_ids: List[str] = []
    embryo_sex: List[Sex] = []
    gene_symbol: str
    gene_start: int
    gene_end: int
    chr: Chromosome
    flanking_region_size: FlankingRegions
    trio_only: bool
    CHaS_input_fields: List[str] = []
    report_header_info: Dict[str, str] = {}
    input_ChAS_filepath: Path
    expected_CHaS_file_name: Optional[str] = None

    class Config:
        str_strip_whitespace = True  # Auto strip leading/trailing white spaces

    @root_validator(pre=True, allow_reuse=True)
    def check_gene_end_after_gene_start(cls, values):
        """
        Check that gene_end is after gene_start"""
        gene_start = values.get("gene_start")
        gene_end = values.get("gene_end")
        if gene_start > gene_end:
            raise ValueError("gene_end should not be before gene_start")
        return values

    @root_validator
    def check_embryo_ids_and_sex_length(cls, values):
        """
        Check that every embryo_id has the embryo sex provided,
        either 'UNKNOWN', 'MALE' or 'FEMALE'"""
        embryo_ids = values.get("embryo_ids")
        embryo_sex = values.get("embryo_sex")
        if len(embryo_ids) != len(embryo_sex):
            raise ValueError("Number of embryo_ids should match number of embryo_sex entries")
        return values

    @root_validator
    def check_reference_not_partner(cls, values):
        """
        Check that reference is has not been duplicated as a partner"""
        reference = values.get("reference")
        if reference == values.get("male_partner") or reference == values.get("female_partner"):
            raise ValueError("Reference should not be the same as either male_partner or female_partner")
        return values

    @root_validator
    def check_x_linked_inheritance_and_chromosome(cls, values):
        """Ensure that for X-linked inheritance, the gene is on chromosome X"""
        mode_of_inheritance = values.get("mode_of_inheritance").value
        chromosome = values.get("chr").value
        if mode_of_inheritance == InheritanceMode.X_LINKED and chromosome != Chromosome.CHR_X:
            raise ValueError(f"If inheritance mode is X_LINKED, chromosome should be X, not {chromosome}")
        return values

    @root_validator
    def check_rhchp_extension(scls, values):
        """
        Check that all the identifiers end with '.rhchp'"""
        male_partner = values.get("male_partner")
        female_partner = values.get("female_partner")
        embryo_ids = values.get("embryo_ids")

        if not male_partner.endswith(".rhchp"):
            raise ValueError(f"male_partner {male_partner} should end with '.rhchp'")
        if not female_partner.endswith(".rhchp"):
            raise ValueError(f"female_partner {female_partner} should end with '.rhchp'")
        for embryo_id in embryo_ids:
            if not embryo_id.endswith(".rhchp"):
                raise ValueError(f"embryo_id {embryo_id} should end with '.rhchp'")

        return values

    @root_validator
    def check__mode_of_inheritance(cls, values):
        """
        Check that the mode of inheritance is valid.
        """
        mode_of_inheritance = values.get("mode_of_inheritance").value

        if mode_of_inheritance not in [
            "autosomal_dominant",
            "autosomal_recessive",
            "x_linked",
        ]:
            raise ValueError(
                f"Invalid Mode of Inheritance {mode_of_inheritance } entered as argument,  ",
                "should be 'Autosomal Dominant', 'Autosomal Recessive', or 'X-Linked'",
            )

        return values

    @root_validator
    def check_gene_symbol(cls, values):
        """
        Check that gene_symbol is a non-empty string.
        """
        gene_symbol = values.get("gene_symbol")

        if not gene_symbol or gene_symbol.strip() == "":
            raise ValueError("Invalid gene_symbol: must be a non-empty string.")

        return values

    @root_validator
    def check_reference_values_for_mode_of_inheritance(cls, values):
        """
        Check that reference values are valid for the mode of inheritance.
        """
        mode_of_inheritance = values.get("mode_of_inheritance").value
        # reference_sex = values.get("reference_sex").value
        reference_status = values.get("reference_status").value
        reference_sex = values.get("reference_sex").value

        # Define allowable values based on mode_of_inheritance
        allowable_values = {
            "x_linked": {
                "reference_sex": {"female", "male"},
                "reference_status": {"carrier", "affected", "unaffected"},
            },
            "autosomal_dominant": {
                "reference_sex": {"male", "female", "unknown"},
                "reference_status": {"affected", "unaffected"},
            },
            "autosomal_recessive": {
                "reference_sex": {"male", "female", "unknown"},
                "reference_status": {"affected", "unaffected"},
            },
        }

        # Validate reference_sex for mode of inheritance
        if reference_sex not in allowable_values[mode_of_inheritance]["reference_sex"]:
            raise ValueError(f"Invalid reference_sex '{reference_sex}' for mode_of_inheritance '{mode_of_inheritance}'")

        # Validate reference_status for mode of inheritance
        if reference_status not in allowable_values[mode_of_inheritance]["reference_status"]:
            raise ValueError(
                f"Invalid reference_status '{reference_status}' for mode_of_inheritance '{mode_of_inheritance}'"
            )

        return values

    @root_validator
    def check_fields_in_CHaS_input_fields(cls, values):
        """
        Check that all the fields in CHaS_input_fields are present in the input fields.
        """
        male_partner = values.get("male_partner")
        female_partner = values.get("female_partner")
        embryo_ids = values.get("embryo_ids")
        CHaS_input_fields = values.get("CHaS_input_fields")

        if male_partner not in CHaS_input_fields:
            raise ValueError(f"male_partner '{male_partner}' not found in CHaS_input_fields")
        if female_partner not in CHaS_input_fields:
            raise ValueError(f"female_partner '{female_partner}' not found in CHaS_input_fields")
        for embryo_id in embryo_ids:
            if embryo_id not in CHaS_input_fields:
                raise ValueError(f"embryo_id '{embryo_id}' not found in CHaS_input_fields")

        return values

    @root_validator
    def check_partners_in_autosomal_dominant(cls, values):
        """
        Check that in autosomal dominant cases, both partners are not unaffected.
        """
        mode_of_inheritance = values.get("mode_of_inheritance").value
        male_partner_status = values.get("male_partner_status").value
        female_partner_status = values.get("female_partner_status").value

        if (
            mode_of_inheritance == "autosomal_dominant"
            and male_partner_status == "unaffected"
            and female_partner_status == "unaffected"
        ):
            raise ValueError("In autosomal_dominant cases, both partners cannot be 'unaffected'.")

        return values

    @root_validator
    def check_filename(cls, values):
        """
        Check that the filename matches the expected filename.
        """
        filepath = values.get("input_ChAS_filepath")
        expected_filename = values.get("expected_CHaS_file_name")
        if expected_filename and filepath.name != expected_filename:
            raise ValueError(
                f"The filename '{filepath.name}' does not match the expected filename '{expected_filename}'!",
                "As specified in provided excel spreadsheet",
            )
        return values

    @property
    def unaffected_partner(self) -> Optional[str]:
        """
        Returns the identifier of the unaffected partner in autosomal dominant cases."""
        if self.mode_of_inheritance != InheritanceMode.AUTOSOMAL_DOMINANT:
            return None

        if self.male_partner_status == Status.UNAFFECTED:
            return self.male_partner
        elif self.female_partner_status == Status.UNAFFECTED:
            return self.female_partner
        else:
            raise ValueError("Both partners cannot be affected in autosomal dominant inheritance.")

    @property
    def affected_partner(self) -> Optional[str]:
        """
        Returns the identifier of the affected partner in autosomal dominant cases.
        """
        if self.mode_of_inheritance != InheritanceMode.AUTOSOMAL_DOMINANT:
            return None

        if self.male_partner_status == Status.AFFECTED:
            return self.male_partner
        elif self.female_partner_status == Status.AFFECTED:
            return self.female_partner
        else:
            raise ValueError("Both partners cannot be unaffected in autosomal dominant inheritance.")

    def get_partner1(self):
        """
        Logic to determine the first partner based on the mode of inheritance.
        """
        if self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT:
            return self.affected_partner
        elif (
            self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE
            or self.mode_of_inheritance == InheritanceMode.X_LINKED
        ):
            return self.female_partner
        else:
            raise ValueError("Invalid mode of inheritance")

    def get_partner2(self):
        """
        Logic to determine the second partner based on the mode of inheritance.
        """
        if self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT:
            return self.unaffected_partner
        elif (
            self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE
            or self.mode_of_inheritance == InheritanceMode.X_LINKED
        ):
            return self.male_partner
        else:
            raise ValueError("Invalid mode of inheritance")
