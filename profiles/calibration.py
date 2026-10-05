"""
Raw sensor readings to calibrated geophysical values.

This block existed twice, byte-identical apart from whitespace, in
FlightLog.apply_thermo_coeffs and Thermo_Profile._init2 - one needed
calibrated values for the a0 file and the other for the gridded profile, so
each grew its own copy.

The original selected sensors by sniffing substrings out of the thermo_data
dictionary's keys ("resi" in key, "temp" in key and "_" not in key, ...),
which made the result depend on dictionary insertion order. The selection is
explicit here.

Both copies also chose the temperature path by looking for resistances in
the log. That choice now belongs to the flight's CalibrationSource - see
calibrate_temperature.
"""
import numpy as np

import profiles.utils as utils
from profiles import Coef_Manager, schema


def _sensor_series(thermo_data, prefix):
    """ Per-sensor arrays for a prefix, one per sensor slot.

    The list index is the sensor number minus one, always: a sensor missing
    from thermo_data comes back as an all-NaN array of the same length as
    the others, not as a gap that shifts every later sensor down a slot.
    Callers pair these with serial numbers and output variables by slot, and
    the ensemble QC treats an all-NaN series as an empty position.

    :rtype: list[np.ndarray]
    :return: N_SENSORS arrays, or [] if no sensor with this prefix is present
    """
    present = {number: np.asarray(thermo_data[f'{prefix}{number}'].magnitude,
                                  dtype=float)
               for number in range(1, schema.N_SENSORS + 1)
               if f'{prefix}{number}' in thermo_data}
    if not present:
        return []
    length = len(next(iter(present.values())))
    return [present.get(number, np.full(length, np.nan))
            for number in range(1, schema.N_SENSORS + 1)]


def calibrate_temperature(thermo_data, serial_numbers, record=None,
                          source=None, when=None):
    """ Per-sensor temperature in K.

    Which path is taken is decided by the calibration source, not by what
    the log happens to contain. Under table calibration, resistance is
    preferred because the Steinhart-Hart conversion is per-sensor and
    non-linear - it has to be applied to each thermistor before any
    averaging. Under onboard calibration the autopilot has already done
    exactly that, per sensor, with the coefficients for the thermistor
    actually fitted, so the logged temperature is taken as it stands.

    Deciding by availability was wrong once firmware began logging both:
    current logs carry resistances and calibrated temperatures, so the
    resistance branch won, and with no serial numbers to look up it applied
    the catch-all `Imet,0` coefficients to every sensor.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :param dict record: if given, the coefficient row used for each sensor
       is stored here under 'imet<n>', for output provenance
    :param source: a CalibrationSource. Resolved from the log's serial
       numbers when omitted. Every coefficient lookup goes through it, so
       the directory it was built with is the one that decides the numbers.
    :param when: the flight's start time, to select among dated coefficient
       rows for a recalibrated sensor
    :rtype: list[np.ndarray]
    """
    if source is None:
        source = Coef_Manager.source_for_flight(serial_numbers)

    if source.temperature_from == 'logged':
        if record is not None:
            record['temperature_source'] = (
                'calibrated onboard; IMET.T used as logged')
        return _sensor_series(thermo_data, 'temp')

    resistances = _sensor_series(thermo_data, 'resi')

    if not resistances:
        if record is not None:
            record['temperature_source'] = 'logged (no resistances present)'
        return _sensor_series(thermo_data, 'temp')

    if record is not None:
        record['temperature_source'] = 'Steinhart-Hart from logged resistance'

    calibrated = []
    for i, values in enumerate(resistances):
        if f'resi{i + 1}' not in thermo_data:
            # Absent slot: stays NaN so the list index remains the sensor.
            calibrated.append(values)
            continue
        serial = serial_numbers[f'imet{i + 1}']
        # One lookup serves both the arithmetic and the provenance record,
        # so what is recorded is what was applied.
        coefs = source.get_coefs('Imet', serial, when=when)
        calibrated.append(utils.steinhart_hart(values, coefs))
        if record is not None:
            record[f'imet{i + 1}'] = coefs
    return calibrated


def calibrate_humidity(thermo_data, serial_numbers):
    """ Per-sensor relative humidity in percent.

    :param dict thermo_data: as returned by FlightLog.thermo_data()
    :param dict serial_numbers: sensor serials, 0 where unknown
    :rtype: list[np.ndarray]
    """
    return [utils.rh_calib(values, serial_numbers[f'rh{i + 1}'])
            if f'rh{i + 1}' in thermo_data else values
            for i, values in enumerate(_sensor_series(thermo_data, 'rh'))]
