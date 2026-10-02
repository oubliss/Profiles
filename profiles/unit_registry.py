"""
The shared pint registry, with this package's extra units.

These definitions used to run as an import side effect of Raw_Profile and
again inside Profile_Set.read_netCDF, which meant anything wanting `gPerKg`
had to import a data class first to get it. Import this module instead.

Re-defining a unit pint already knows raises, and modern pint and MetPy
ship `percent`, so each definition is guarded.
"""
from metpy.units import units

#: Units this package adds on top of MetPy's registry.
EXTRA_UNITS = (
    'percent = 0.01*count = %',
    'gPerKg = 0.001*count = g/Kg',
)

for _definition in EXTRA_UNITS:
    _name = _definition.split(' ', 1)[0]
    if _name not in units:
        units.define(_definition)

__all__ = ['units', 'EXTRA_UNITS']
