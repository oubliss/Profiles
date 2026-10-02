"""
Raw sensor readings to calibrated geophysical values.

This block existed twice, byte-identical apart from whitespace, in
FlightLog.apply_thermo_coeffs and Thermo_Profile._init2 - one needed
calibrated values for the a0 file and the other for the gridded profile, so
each grew its own copy.

The original selected sensors by sniffing substrings out of the thermo_data
dictionary's keys ("resi" in key, "temp" in key and "_" not in key, ...),
which made the result depend on dictionary insertion order. The selection is
explicit here; the outcome is unchanged because every CopterSonde log
carries resistances.
"""
import numpy as np

import profiles.utils as utils
from profiles import schema


def _sensor_series(thermo_data, prefix):
    """Per-sensor arrays for a prefix, in sensor order, skipping absent ones."""
    series = []
    for number in range(1, schema.N_SENSORS + 1):
        key = f'{prefix}{number}'
        if key in thermo_data:
            series.append(thermo_data[key].magnitude)
    return series


def calibrate_temperature(thermo_data, serial_numbers):
    """ Per-sensor temperature in K.

    Resistance is preferred where the log carries it, because the
    Steinhart-Hart conversion is per-sensor and non-linear - it has to be
    applied to each thermistor before any averaging. Logs without
    resistances fall back to the temperature the autopilot recorded.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :rtype: list[np.ndarray]
    """
    resistances = _sensor_series(thermo_data, 'resi')

    if resistances:
        return [utils.temp_calib(values, serial_numbers[f'imet{i + 1}'])
                for i, values in enumerate(resistances)]

    return _sensor_series(thermo_data, 'temp')


def calibrate_humidity(thermo_data, serial_numbers):
    """ Per-sensor relative humidity in percent.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :rtype: list[np.ndarray]
    """
    return [utils.rh_calib(values, serial_numbers[f'rh{i + 1}'])
            for i, values in enumerate(_sensor_series(thermo_data, 'rh'))]
