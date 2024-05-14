"""
This module contains Enum classes for the SNP haplotyper package.
The Enum classes are used to define constants for the package."""

from enum import Enum

# Define constants


class FlankingRegions(int, Enum):
    """Enum class for flanking regions. The values are the number of megabases to flank the SNP with."""

    FLANK_2MB = 2
    FLANK_3MB = 3


class InheritanceMode(Enum):
    """Enum class for inheritance modes."""

    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    X_LINKED = "x_linked"


class Status(Enum):
    """Status enum class for the status of a family member."""

    AFFECTED = "affected"
    UNAFFECTED = "unaffected"
    CARRIER = "carrier"


class Relationship(Enum):
    """Enum class for the relationship between family members."""

    GRANDPARENT = "grandparent"
    CHILD = "child"


class Sex(Enum):
    """Enum class for the sex of the family member."""

    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"


class Chromosome(str, Enum):
    """Enum class for the chromosome. The values are the standard chromosome names."""

    CHR_1 = "1"
    CHR_2 = "2"
    CHR_3 = "3"
    CHR_4 = "4"
    CHR_5 = "5"
    CHR_6 = "6"
    CHR_7 = "7"
    CHR_8 = "8"
    CHR_9 = "9"
    CHR_10 = "10"
    CHR_11 = "11"
    CHR_12 = "12"
    CHR_13 = "13"
    CHR_14 = "14"
    CHR_15 = "15"
    CHR_16 = "16"
    CHR_17 = "17"
    CHR_18 = "18"
    CHR_19 = "19"
    CHR_20 = "20"
    CHR_21 = "21"
    CHR_22 = "22"
    CHR_X = "X"  # Capitalized X and Y for standard naming conventions
    CHR_Y = "Y"
