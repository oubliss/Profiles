"""
Deprecated import location. Use ``profiles.flight`` instead.

The class is now called FlightLog, because it holds a whole flight at native
rate rather than a profile; the old name made `Profile` and `Raw_Profile`
read as variations on one thing when they are different objects.
"""
import warnings

from profiles.flight import FlightLog, Raw_Profile  # noqa: F401

warnings.warn(
    'profiles.Raw_Profile is deprecated; import FlightLog from '
    'profiles.flight instead.', DeprecationWarning, stacklevel=2)
