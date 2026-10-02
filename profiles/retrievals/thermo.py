"""
Thermodynamic quantities derived from temperature, humidity and pressure.

Thin wrappers over MetPy, gathered here so that the derivation is one
testable function rather than a block in the middle of a constructor.
"""
from metpy import calc
from profiles.unit_registry import units


def derive(pressure, temperature, relative_humidity):
    """ Mixing ratio, potential temperature, dewpoint and specific humidity.

    :param pressure: a pint Quantity in pressure units
    :param temperature: a pint Quantity in temperature units
    :param relative_humidity: a pint Quantity in percent
    :rtype: dict
    :return: {'mixing_ratio', 'theta', 'T_d', 'q'} as pint Quantities
    """
    mixing_ratio = calc.mixing_ratio_from_relative_humidity(
        pressure, temperature, relative_humidity.magnitude / 100)

    return {
        'mixing_ratio': mixing_ratio,
        'theta': calc.potential_temperature(pressure, temperature),
        'T_d': calc.dewpoint_from_relative_humidity(temperature,
                                                    relative_humidity),
        # NOTE: MetPy returns kg/kg and this attaches the g/kg label
        # without rescaling, so the magnitude is kg/kg (~0.016) while the
        # units read gPerKg. Profile.save_netcdf compensates with an
        # explicit * 1e3; Thermo_Profile._save_netCDF does not, so the
        # thermo_* intermediate is mislabelled by a factor of 1000.
        # Preserved here rather than silently corrected - deciding whether
        # to rescale the value or fix the label is Stage 5's writer work.
        'q': calc.specific_humidity_from_mixing_ratio(mixing_ratio)
             * units.gPerKg,
    }
