"""
Deprecated. Use profiles.config, or set $WXUAS_DIR.

Configuration used to live here, inside the package, so changing it meant
editing installed source and every machine drifted. This module remains
because scripts and test suites repoint coef_info.FILE_PATH at runtime, and
profiles.config still honours that when it is set.
"""
import os
import warnings
from types import SimpleNamespace

warnings.warn(
    'profiles.conf is deprecated; set the WXUAS_DIR environment variable or '
    'pass a directory to TableCalibration(). See profiles.config.',
    DeprecationWarning, stacklevel=2)

wxuas_dir = os.path.join(os.path.expanduser("~"), ".wxuas")

#: Legacy configuration namespace. FILE_PATH is still read by
#: profiles.config.coefficient_dir and still takes precedence when set.
#: USE_AZURE and AZURE_CONNECTION_STRING are ignored - the Azure backend
#: was entirely commented out and has been removed.
coef_info = SimpleNamespace(USE_AZURE="NO", AZURE_CONNECTION_STRING=None,
                            FILE_PATH=None)
