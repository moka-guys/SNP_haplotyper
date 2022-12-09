# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""

    pass


class ArgumentInputError(Error):
    """Raised when the input parameters do not make sense considering the underlying biology"""

    pass


class InvalidParameterSelectedError(Error):
    """The config.py has a flag prohibiting the script from running with the selected parameters.  This is to prevent the user from running the script for options which have not yet been validated."""

    pass
