"""
The flight's calibration source decides the coefficients, not just the path.

Before this, every lookup went through a process-global manager built from
whatever directory was configured first, so a source built for a particular
directory changed nothing, dated rows were never selected (no caller passed
a date), and an explicit tail number was overwritten by the registry.
"""
import shutil

import netCDF4
import numpy as np
import pandas as pd
import pytest

from profiles import calibration
from profiles.Coef_Manager import OnboardCalibration, TableCalibration
from profiles.coef_table import AmbiguousCoefficients, MissingCoefficients
from profiles.flight import FlightLog
from profiles.processing import ProcessingConfig, all_profiles, process_flights
from test import COEF_PATH
from test.harness import (PROFILE_START_HEIGHT, RESOLUTION, RES_UNITS,
                          staged_bin)

#: imet1 of the reference flight, and the test table's own value for it.
SENSOR = '62275'
ORIGINAL_A = '9.89E-04'


class _Magnitude:
    """Stands in for a pint Quantity: calibration reads only .magnitude."""

    def __init__(self, values):
        self.magnitude = values


def thermo_data():
    return {'resi1': _Magnitude(np.array([10000.0]))}


def config(tmp, **extra):
    return ProcessingConfig(
        resolution=RESOLUTION, res_units=RES_UNITS, dev=True, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT, min_levels=50,
        coefficient_dir=str(tmp / 'coefs'), **extra)


def coefs_copy(tmp):
    shutil.copytree(COEF_PATH, tmp / 'coefs')
    return tmp / 'coefs'


def edit_table(directory, edit):
    path = directory / 'MasterCoefList.csv'
    table = pd.read_csv(path, dtype=str).fillna('')
    table = edit(table)
    table.to_csv(path, index=False)


def split_sensor(directory, pre_a, post_a, split='2022-01-01'):
    """Replace one Imet sensor's row with a pre- and a post-`split` pair."""
    def edit(table):
        table['ValidFrom'] = ''
        table['ValidTo'] = ''
        mask = (table.SensorType == 'Imet') & (table.SerialNumber == SENSOR)
        original = table[mask].iloc[0].copy()
        table = table[~mask]

        before, after = original.copy(), original.copy()
        before['A'], before['ValidTo'] = pre_a, split
        after['A'], after['ValidFrom'] = post_a, split
        return pd.concat([table, before.to_frame().T, after.to_frame().T])

    edit_table(directory, edit)


def run(tmp, **extra):
    results = process_flights([str(staged_bin(tmp))], config(tmp, **extra))
    assert all(result.ok for result in results), \
        [result.error for result in results]
    return all_profiles(results)[0]


class TestDatedCoefficients:
    def test_flight_date_selects_the_row_end_to_end(self, tmp_path):
        directory = coefs_copy(tmp_path)
        split_sensor(directory, pre_a=ORIGINAL_A, post_a='5.00E-04')

        profile = run(tmp_path)   # flown 2021-06-09: before the split

        applied = profile.calibration_record['imet1']
        assert applied['A'] == ORIGINAL_A
        assert applied['ValidTo'] == '2022-01-01'

        c1 = tmp_path / 'c1.cdf'
        profile.save_netcdf(str(c1))
        with netCDF4.Dataset(c1) as handle:
            assert handle.getncattr('coef_imet1_a') == ORIGINAL_A
            assert handle.getncattr('coef_imet1_validto') == '2022-01-01'
            assert handle.getncattr('coefficient_lookup_time').startswith(
                '2021-06-09')

    def test_a_later_flight_gets_the_other_row(self, tmp_path):
        directory = coefs_copy(tmp_path)
        split_sensor(directory, pre_a=ORIGINAL_A, post_a='5.00E-04')
        source = TableCalibration(directory)
        serials = {'imet1': SENSOR}

        before, after = {}, {}
        calibration.calibrate_temperature(
            thermo_data(), serials, record=before, source=source,
            when='2021-06-09')
        calibration.calibrate_temperature(
            thermo_data(), serials, record=after, source=source,
            when='2023-01-01')
        assert before['imet1']['A'] == ORIGINAL_A
        assert after['imet1']['A'] == '5.00E-04'

    def test_undated_lookup_of_a_split_sensor_is_still_ambiguous(
            self, tmp_path):
        """The failure the flight date fixes - nobody used to pass one."""
        directory = coefs_copy(tmp_path)
        split_sensor(directory, pre_a=ORIGINAL_A, post_a='5.00E-04')
        with pytest.raises(AmbiguousCoefficients):
            TableCalibration(directory).get_coefs('Imet', SENSOR)

    def test_start_time_is_the_first_logged_timestamp(self, tmp_path):
        flight = FlightLog(str(staged_bin(tmp_path)), dev=True,
                           nc_level=None, coefficient_dir=str(COEF_PATH))
        assert flight.start_time.date().isoformat() == '2021-06-09'


class TestSourceDecidesTheNumbers:
    def test_directory_changes_the_temperature(self, tmp_path):
        directory = coefs_copy(tmp_path)

        def edit(table):
            mask = ((table.SensorType == 'Imet')
                    & (table.SerialNumber == SENSOR))
            table.loc[mask, 'A'] = '1.50E-03'
            return table
        edit_table(directory, edit)

        serials = {'imet1': SENSOR}
        default = calibration.calibrate_temperature(
            thermo_data(), serials, source=TableCalibration(COEF_PATH))
        edited = calibration.calibrate_temperature(
            thermo_data(), serials, source=TableCalibration(directory))
        assert not np.allclose(default[0], edited[0])

    def test_onboard_source_refuses_resistance_lookups(self):
        """Its refusal used to be unreachable: lookups bypassed the source."""
        source = OnboardCalibration(COEF_PATH)
        source.temperature_from = 'resistance'     # force the lookup path
        with pytest.raises(MissingCoefficients, match='calibrated onboard'):
            calibration.calibrate_temperature(
                thermo_data(), {'imet1': SENSOR}, source=source)

    def test_provenance_reports_the_sources_directory(self, tmp_path):
        directory = coefs_copy(tmp_path)
        profile = run(tmp_path)
        c1 = tmp_path / 'c1.cdf'
        profile.save_netcdf(str(c1))
        with netCDF4.Dataset(c1) as handle:
            assert handle.getncattr('coefficient_directory') == str(directory)
            assert handle.getncattr('calibration_source') == 'TableCalibration'

    def test_flights_source_is_used_for_wind(self, tmp_path):
        directory = coefs_copy(tmp_path)

        def edit(table):
            mask = ((table.SensorType == 'Wind')
                    & (table.SerialNumber == 'N934UA'))
            table.loc[mask, 'A'] = '2.00E+01'
            return table
        edit_table(directory, edit)

        profile = run(tmp_path)
        assert profile.calibration_record['wind']['A'] == '2.00E+01'


class TestExplicitTailNumber:
    def test_explicit_tail_number_beats_the_registry(self, tmp_path):
        directory = coefs_copy(tmp_path)
        with open(directory / 'MasterCoefList.csv', 'a') as handle:
            handle.write('Wind,N944UA,na,E1,2.00E+01,-1.00E+00,na,na,na,'
                         'Active\n')
        # The log's copterID (6) maps to N934UA in copterID.csv.
        profile = run(tmp_path, tail_number='N944UA')

        assert profile.tail_number == 'N944UA'
        assert profile.calibration_record['wind']['SerialNumber'] == 'N944UA'
        assert profile.calibration_record['wind']['A'] == '2.00E+01'

    def test_registry_is_used_when_none_is_given(self, tmp_path):
        coefs_copy(tmp_path)
        assert run(tmp_path).tail_number == 'N934UA'
