"""
Where the package looks for coefficients, and how to tell it otherwise.

Configuration used to live in profiles/conf.py, which meant changing it
required editing a file inside the installed package - so every machine's
install drifted and nothing was reproducible. Resolution order is now:

1. an explicit path passed to the object that needs it
2. the WXUAS_DIR environment variable
3. ~/.wxuas

profiles.conf still exists and still works, with a DeprecationWarning.
"""
import os
from pathlib import Path

#: Environment variable naming the coefficient directory.
ENV_VAR = 'WXUAS_DIR'

#: Where to look when nothing else says otherwise.
DEFAULT_DIR = Path.home() / '.wxuas'

#: Files expected in that directory.
COEF_FILE = 'MasterCoefList.csv'
COPTER_ID_FILE = 'copterID.csv'
SCOOPS_FILE = 'scoops.csv'


def coefficient_dir(explicit=None):
    """ The directory holding the coefficient tables.

    :param explicit: a path that overrides everything else, or None
    :rtype: pathlib.Path
    """
    if explicit is not None:
        return Path(explicit)

    # profiles.conf is the deprecated route, but while it exists a value
    # set there has to keep winning over the default - test suites and
    # existing scripts repoint it at runtime.
    from profiles.conf import coef_info
    if getattr(coef_info, 'FILE_PATH', None):
        return Path(coef_info.FILE_PATH)

    from_env = os.environ.get(ENV_VAR)
    if from_env:
        return Path(from_env)

    return DEFAULT_DIR


def describe_lookup(explicit=None):
    """ Human-readable account of where the directory came from.

    Used in error messages, so that "no coefficients found" says which
    directory was searched and why.

    :rtype: str
    """
    if explicit is not None:
        return f'{explicit} (passed explicitly)'

    from profiles.conf import coef_info
    if getattr(coef_info, 'FILE_PATH', None):
        return f'{coef_info.FILE_PATH} (from profiles.conf.coef_info, deprecated)'

    from_env = os.environ.get(ENV_VAR)
    if from_env:
        return f'{from_env} (from ${ENV_VAR})'

    return f'{DEFAULT_DIR} (default)'
