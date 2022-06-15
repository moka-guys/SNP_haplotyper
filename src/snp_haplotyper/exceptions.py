# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""

    pass


class ArgumentInputError(Error):
    """Raised when the input parameters do not make sense considering the underlying biology"""

    pass
