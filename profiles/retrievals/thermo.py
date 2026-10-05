"""
Thermodynamic quantities derived from temperature, humidity and pressure.

Thin wrappers over MetPy, gathered here so that the derivation is one
testable function rather than a block in the middle of a constructor.
"""
from metpy import calc


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
        # MetPy returns a dimensionless fraction (kg/kg). It is kept in
        # those units, so q.to('g/kg') is the right conversion for display
        # and the writers do it explicitly. This used to attach the gPerKg
        # label without rescaling, which left the magnitude 1000x too small
        # for its units.
        'q': calc.specific_humidity_from_mixing_ratio(mixing_ratio)
             .to('kg/kg'),
    }
