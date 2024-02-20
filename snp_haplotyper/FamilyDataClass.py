import logging
import pandas as pd
from pydantic import BaseModel, root_validator
from typing import Dict, List, Optional
from pathlib import Path

from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)

logger = logging.getLogger("BASHer_logger")


class FamilyData(BaseModel):
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
    expected_CHaS_file_name: Optional[str] = None  # TODO Link to front end import
    """
    FamilyData dataclass
    =================

    This class represents the family metada provided for the SNP haplotyping and contains relevant data validation of the input.

    """

    @root_validator(pre=True)
    def check_gene_end_after_gene_start(cls, values):
        gene_start = values.get("gene_start")
        gene_end = values.get("gene_end")
        assert gene_start < gene_end, "gene_end should not be before gene_start"
        return values

    @root_validator
    def check_embryo_ids_and_sex_length(cls, values):
        embryo_ids = values.get("embryo_ids")
        embryo_sex = values.get("embryo_sex")
        if len(embryo_ids) != len(embryo_sex):
            raise ValueError(
                "Number of embryo_ids should match number of embryo_sex entries"
            )
        return values

    @root_validator
    def check_reference_not_partner(cls, values):
        reference = values.get("reference")
        if reference == values.get("male_partner") or reference == values.get(
            "female_partner"
        ):
            raise ValueError(
                "Reference should not be the same as either male_partner or female_partner"
            )
        return values

    @root_validator
    def check_x_linked_inheritance_and_chromosome(cls, values):
        mode_of_inheritance = values.get("mode_of_inheritance").value
        chr = values.get("chr").value
        if mode_of_inheritance == InheritanceMode.X_LINKED and chr != Chromosome.CHR_X:
            raise ValueError("If inheritance mode is X_LINKED, chromosome should be X")
        return values

    @root_validator
    def check_rhchp_extension(cls, values):
        male_partner = values.get("male_partner")
        female_partner = values.get("female_partner")
        embryo_ids = values.get("embryo_ids")

        if not male_partner.endswith(".rhchp"):
            raise ValueError(f"male_partner {male_partner} should end with '.rhchp'")
        if not female_partner.endswith(".rhchp"):
            raise ValueError(
                f"female_partner {female_partner} should end with '.rhchp'"
            )
        for embryo_id in embryo_ids:
            if not embryo_id.endswith(".rhchp"):
                raise ValueError(f"embryo_id {embryo_id} should end with '.rhchp'")

        return values

    @root_validator
    def check_gene_symbol(cls, values):
        gene_symbol = values.get("gene_symbol")

        if not gene_symbol or gene_symbol.strip() == "":
            raise ValueError("Invalid gene_symbol: must be a non-empty string.")

        return values

    @root_validator
    def check_reference_values_for_mode_of_inheritance(cls, values):
        mode_of_inheritance = values.get("mode_of_inheritance").value
        # reference_sex = values.get("reference_sex").value
        reference_status = values.get("reference_status").value

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
        # TODO
        # Validate reference_sex for mode of inheritance
        # if reference_sex not in allowable_values[mode_of_inheritance]["reference_sex"]:
        # raise ValueError(
        #     f"Invalid reference_sex '{reference_sex}' for mode_of_inheritance '{mode_of_inheritance}'"
        # )

        # Validate reference_status for mode of inheritance
        if (
            reference_status
            not in allowable_values[mode_of_inheritance]["reference_status"]
        ):
            raise ValueError(
                f"Invalid reference_status '{reference_status}' for mode_of_inheritance '{mode_of_inheritance}'"
            )

        return values

    @root_validator
    def check_fields_in_CHaS_input_fields(cls, values):
        male_partner = values.get("male_partner")
        female_partner = values.get("female_partner")
        embryo_ids = values.get("embryo_ids")
        CHaS_input_fields = values.get("CHaS_input_fields")

        if male_partner not in CHaS_input_fields:
            raise ValueError(
                f"male_partner '{male_partner}' not found in CHaS_input_fields"
            )
        if female_partner not in CHaS_input_fields:
            raise ValueError(
                f"female_partner '{female_partner}' not found in CHaS_input_fields"
            )
        for embryo_id in embryo_ids:
            if embryo_id not in CHaS_input_fields:
                raise ValueError(
                    f"embryo_id '{embryo_id}' not found in CHaS_input_fields"
                )

        return values

    @root_validator
    def check_partners_in_autosomal_dominant(cls, values):
        mode_of_inheritance = values.get("mode_of_inheritance").value
        male_partner_status = values.get("male_partner_status").value
        female_partner_status = values.get("female_partner_status").value

        if (
            mode_of_inheritance == "autosomal_dominant"
            and male_partner_status == "unaffected"
            and female_partner_status == "unaffected"
        ):
            raise ValueError(
                "In autosomal_dominant cases, both partners cannot be 'unaffected'."
            )

        return values

    @root_validator
    def check_filename(cls, values):
        filepath = values.get("input_ChAS_filepath")
        expected_filename = values.get("expected_CHaS_file_name")
        if expected_filename and filepath.name != expected_filename:
            raise ValueError(
                f"The filename '{filepath.name}' does not match the expected filename '{expected_filename}'!"
            )
        return values

    @property
    def unaffected_partner(self) -> Optional[str]:
        if self.mode_of_inheritance != InheritanceMode.AUTOSOMAL_DOMINANT:
            return None

        if self.male_partner_status == Status.UNAFFECTED:
            return self.male_partner
        elif self.female_partner_status == Status.UNAFFECTED:
            return self.female_partner
        else:
            raise ValueError(
                "Both partners cannot be affected in autosomal dominant inheritance."
            )

    @property
    def affected_partner(self) -> Optional[str]:
        if self.mode_of_inheritance != InheritanceMode.AUTOSOMAL_DOMINANT:
            return None

        if self.male_partner_status == Status.AFFECTED:
            return self.male_partner
        elif self.female_partner_status == Status.AFFECTED:
            return self.female_partner
        else:
            raise ValueError(
                "Both partners cannot be unaffected in autosomal dominant inheritance."
            )

    def get_partner1(self):
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
        if self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_DOMINANT:
            return self.unaffected_partner
        elif (
            self.mode_of_inheritance == InheritanceMode.AUTOSOMAL_RECESSIVE
            or self.mode_of_inheritance == InheritanceMode.X_LINKED
        ):
            return self.male_partner
        else:
            raise ValueError("Invalid mode of inheritance")
