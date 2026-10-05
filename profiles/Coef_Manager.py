"""
Access to the sensor coefficient tables.

The public surface is unchanged - get_tail_n, get_sensors, get_coefs - but
the implementation is now profiles.coef_table, which indexes the tables
once instead of copying a DataFrame per lookup, and can select a row by the
flight's date.

The Azure backend is gone. Every line of it had been commented out for
years, and the branch that would have selected it fell through to `pass`,
leaving sub_manager as None and failing with AttributeError on first use.
"""
import os
from abc import ABC, abstractmethod

import pandas as pd

from profiles import config
from profiles.coef_table import (AmbiguousCoefficients, CoefTable,
                                 CopterRegistry, MissingCoefficients)

__all__ = ['Coef_Manager', 'CalibrationSource', 'TableCalibration',
           'OnboardCalibration', 'AmbiguousCoefficients',
           'MissingCoefficients']


class CalibrationSource(ABC):
    """ Where calibrated values come from for a given flight.

    Two implementations. Which applies depends on the firmware the flight
    was flown with, not on a preference.

    :cvar str temperature_from: 'resistance' to apply Steinhart-Hart to the
       logged resistances, 'logged' to take IMET.T as it stands. This is a
       property of the flight, not of the data present in the log: current
       logs carry resistances *and* calibrated temperatures, so picking by
       what is available silently recomputes already-calibrated values.
    """

    temperature_from = 'resistance'

    @abstractmethod
    def get_tail_n(self, copterID, when=None):
        """ Tail number for a short vehicle ID.

        :param copterID: the short ID logged as SYSID_THISMAV
        :param when: flight time
        :rtype: str
        """

    @abstractmethod
    def get_coefs(self, type, serial_number, equation=None, when=None):
        """ Coefficients for one sensor.

        :param str type: 'Imet', 'RH' or 'Wind'
        :param serial_number: the sensor serial or airframe tail number
        :param str equation: disambiguates when several rows match
        :param when: flight time
        :rtype: dict
        """

    @abstractmethod
    def get_sensors(self, scoopID):
        """ Sensor serial numbers fitted to a scoop.

        :param str scoopID: the scoop's identifier
        :rtype: dict
        """


class TableCalibration(CalibrationSource):
    """ Coefficients from the CSV tables.

    Correct for flights whose firmware logged raw resistances and sensor
    serial numbers, i.e. everything up to the move to onboard calibration.
    """

    def __init__(self, directory=None):
        """
        :param directory: where the tables live; resolved through
           profiles.config when omitted
        """
        self.directory = config.coefficient_dir(directory)

        if not os.path.isdir(self.directory):
            raise MissingCoefficients(
                f'coefficient directory not found: '
                f'{config.describe_lookup(directory)}')

        self._coefs = CoefTable(self.directory / config.COEF_FILE)
        self._copters = CopterRegistry(self.directory / config.COPTER_ID_FILE)

    def get_tail_n(self, copterID, when=None):
        return self._copters.tail_number(copterID, when=when)

    def get_coefs(self, type, serial_number, equation=None, when=None):
        return self._coefs.lookup(type, serial_number, equation=equation,
                                  when=when)

    def get_sensors(self, scoopID):
        scoop_file = self.directory / config.SCOOPS_FILE
        if not os.path.exists(scoop_file):
            raise MissingCoefficients(f'{scoop_file} does not exist')

        table = pd.read_csv(scoop_file)
        matched = table[table.name == scoopID]
        if matched.empty:
            raise MissingCoefficients(
                f'scoop {scoopID!r} is not listed in {scoop_file}')

        return {'imet1': str(matched.imet1.values[0]),
                'imet2': str(matched.imet2.values[0]),
                'imet3': str(matched.imet3.values[0]),
                'imet4': None,
                'rh1': str(matched.rh1.values[0]),
                'rh2': str(matched.rh2.values[0]),
                'rh3': str(matched.rh3.values[0]),
                'rh4': None}


class OnboardCalibration(CalibrationSource):
    """ Trust the values the CopterSonde already calibrated in flight.

    Current firmware applies the thermistor calibration itself and logs the
    result as IMET.T1..T4, and no longer logs USER_SENSORS parameters - so
    the serial numbers a table lookup needs are not there. Asking this
    source for thermodynamic coefficients is a mistake and it says so,
    rather than quietly handing back generic ones.

    Wind is still a table lookup: the airframe calibration is per tail
    number and is not applied onboard.

    This source does no arithmetic on temperature, which is the point: the
    calibration has already happened. Its job is to stop the pipeline
    redoing it. Before this existed, a current-generation flight still took
    the resistance path, and because those logs report no serial numbers
    every lookup fell to the catch-all `Imet,0` row - one shared
    Steinhart-Hart curve applied to all four thermistors. On
    KAEFS_flight2859 that shifted imet2 by +0.100 K and imet3 by -0.109 K
    and, worse, compressed the inter-sensor spread from 0.273 K to 0.064 K,
    manufacturing agreement between sensors and feeding the ensemble QC a
    flattened signal.
    """

    temperature_from = 'logged'

    def __init__(self, directory=None):
        self._directory = directory
        self._table_source = None

    @property
    def _table(self):
        """ The wind/tail-number tables, read only if something asks.

        Built lazily so that processing temperature from an onboard-
        calibrated flight needs no coefficient directory at all - which is
        the whole claim this class makes.
        """
        if self._table_source is None:
            self._table_source = TableCalibration(self._directory)
        return self._table_source

    @property
    def directory(self):
        return self._table.directory

    def get_tail_n(self, copterID, when=None):
        return self._table.get_tail_n(copterID, when=when)

    def get_coefs(self, type, serial_number, equation=None, when=None):
        if type in ('Imet', 'RH'):
            raise MissingCoefficients(
                f'{type} values from this flight were calibrated onboard; '
                f'there are no table coefficients to apply. Use the logged '
                f'temperatures rather than recomputing from resistance.')
        return self._table.get_coefs(type, serial_number, equation=equation,
                                     when=when)

    def get_sensors(self, scoopID):
        return self._table.get_sensors(scoopID)


CALIBRATION_MODES = ('auto', 'table', 'onboard')


def logs_sensor_serials(serial_numbers):
    """ Whether the log names the thermodynamic sensors fitted to it.

    :param dict serial_numbers: as parsed from the log's PARM records
    :rtype: bool
    """
    return any(serial_numbers.get(f'{kind}{n}')
               for kind in ('imet', 'rh')
               for n in range(1, 5))


def source_for_flight(serial_numbers, directory=None, mode='auto'):
    """ Pick the calibration source a flight's log implies.

    A log that reports sensor serial numbers was flown with firmware that
    expected table calibration. One that reports none calibrated onboard.

    :param dict serial_numbers: as parsed from the log's PARM records
    :param directory: coefficient directory, or None to resolve it
    :param str mode: 'auto' to decide from the log, or 'table'/'onboard'
       to force one. Forcing is for logs that are wrong about themselves -
       a reprocess of an early onboard-calibrated flight whose firmware
       still wrote the parameters, say.
    :rtype: CalibrationSource
    """
    if mode not in CALIBRATION_MODES:
        raise ValueError(
            f'unknown calibration mode {mode!r}; '
            f'available: {list(CALIBRATION_MODES)}')

    if mode == 'auto':
        mode = 'table' if logs_sensor_serials(serial_numbers) else 'onboard'

    return (TableCalibration(directory) if mode == 'table'
            else OnboardCalibration(directory))


class Coef_Manager(TableCalibration):
    """ Backwards-compatible entry point.

    Equivalent to TableCalibration. Kept because utils.coef_manager and a
    good deal of calling code still name it.
    """
