import pydantic
import pytest
from EnumDataClasses import (
    Chromosome,
    FlankingRegions,
    InheritanceMode,
    Relationship,
    Sex,
    Status,
)
from FamilyDataClass import FamilyData


def valid_family_data():
    # A function that provides valid family data for testing purposes.
    return {
        "mode_of_inheritance": InheritanceMode.X_LINKED,
        "male_partner": "male_sample.rhchp",
        "male_partner_status": Status.CARRIER,
        "female_partner": "female_sample.rhchp",
        "female_partner_status": Status.AFFECTED,
        "consanguineous": False,
        "reference": "reference_sample",
        "reference_status": Status.UNAFFECTED,
        "reference_relationship": Relationship.GRANDPARENT,
        "embryo_ids": ["embryo1.rhchp", "embryo2.rhchp"],
        "embryo_sex": [Sex.MALE, Sex.FEMALE],
        "gene_symbol": "ABCD",
        "gene_start": 100,
        "gene_end": 200,
        "chr": Chromosome.CHR_X,
        "flanking_region_size": FlankingRegions.FLANK_2MB,
        "trio_only": False,
        "CHaS_input_fields": [
            "male_sample.rhchp",
            "female_sample.rhchp",
            "embryo1.rhchp",
            "embryo2.rhchp",
        ],
    }


def test_check_gene_end_after_gene_start():
    with pytest.raises(
        pydantic.error_wrappers.ValidationError,
    ):
        FamilyData(**{**valid_family_data(), "gene_start": 300, "gene_end": 200})


def test_check_embryo_ids_and_sex_length():
    with pytest.raises(
        pydantic.error_wrappers.ValidationError,
        match="Number of embryo_ids should match number of embryo_sex entries",
    ):
        FamilyData(**{**valid_family_data(), "embryo_ids": ["embryo1.rhchp"]})

    with pytest.raises(
        pydantic.error_wrappers.ValidationError,
        match="Number of embryo_ids should match number of embryo_sex entries",
    ):
        FamilyData(**{**valid_family_data(), "embryo_sex": [Sex.MALE]})


def test_check_reference_not_partner():
    with pytest.raises(
        ValueError,
        match="Reference should not be the same as either male_partner or female_partner",
    ):
        FamilyData(**{**valid_family_data(), "reference": "male_sample.rhchp"})

    with pytest.raises(
        ValueError,
        match="Reference should not be the same as either male_partner or female_partner",
    ):
        FamilyData(**{**valid_family_data(), "reference": "female_sample.rhchp"})
