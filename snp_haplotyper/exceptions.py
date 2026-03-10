"""
This module contains custom exceptions for the snp_haplotyper package.
"""


class Error(Exception):
    """Base class for other exceptions"""


class ArgumentInputError(Error):
    """Raised when the input parameters do not make sense considering the underlying biology"""


class InvalidParameterSelectedError(Error):
    """The config.py has a flag prohibiting the script from running with the selected parameters.  This is to prevent
    the user from running the script for options which have not yet been validated."""
