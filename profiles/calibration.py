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


def calibrate_temperature(thermo_data, serial_numbers, record=None):
    """ Per-sensor temperature in K.

    Resistance is preferred where the log carries it, because the
    Steinhart-Hart conversion is per-sensor and non-linear - it has to be
    applied to each thermistor before any averaging. Logs without
    resistances fall back to the temperature the autopilot recorded.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :param dict record: if given, the coefficient row used for each sensor
       is stored here under 'imet<n>', for output provenance
    :rtype: list[np.ndarray]
    """
    resistances = _sensor_series(thermo_data, 'resi')

    if not resistances:
        if record is not None:
            record['temperature_source'] = 'logged (no resistances present)'
        return _sensor_series(thermo_data, 'temp')

    if record is not None:
        record['temperature_source'] = 'Steinhart-Hart from logged resistance'

    calibrated = []
    for i, values in enumerate(resistances):
        serial = serial_numbers[f'imet{i + 1}']
        calibrated.append(utils.temp_calib(values, serial))
        if record is not None:
            try:
                record[f'imet{i + 1}'] = utils.get_coef_manager().get_coefs(
                    'Imet', serial)
            except Exception as exc:          # provenance must never break
                record[f'imet{i + 1}'] = {'error': str(exc)}
    return calibrated


def calibrate_humidity(thermo_data, serial_numbers):
    """ Per-sensor relative humidity in percent.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :rtype: list[np.ndarray]
    """
    return [utils.rh_calib(values, serial_numbers[f'rh{i + 1}'])
            for i, values in enumerate(_sensor_series(thermo_data, 'rh'))]
