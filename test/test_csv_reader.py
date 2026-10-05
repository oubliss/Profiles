"""
The legacy CIWRO-plotter CSV reader (profiles-l9l.25), on a synthetic file.

No real file in this format exists in the repository or alongside it, so
these tests pin only what is unambiguous: columns land in the slots the
other readers use, every IMU channel is read (the block used to sit outside
its loop and populated only the last), VD comes from the vz column, and the
read works without the pandas option removed in pandas 3.
"""
import warnings

import numpy as np
import pytest

from profiles.flight import FlightLog
from profiles.schema import PRESSURE

HEADER = ['date', 'lat', 'lon', 'alt', 'pressure', 'roll', 'pitch', 'yaw',
          'gyry', 'gyrx', 'gyrz', 'vx', 'vy', 'vz', 'accx', 'accy', 'accz',
          'temp1', 'temp2', 'temp3', 'temp4', 'temp5',
          'rh1', 'rh2', 'rh3', 'rh4', 'rh5',
          'gpsboottime', 'pressureboottime', 'attitudeboottime',
          'imetboottime',
          'temp_r1', 'temp_r2', 'temp_r3', 'temp_r4', 'temp_r5']

ROWS = 4


def row(i):
    values = {name: float(index + 1) for index, name in enumerate(HEADER[1:])}
    values.update({'gyrx': 10.0 + i, 'gyry': 20.0 + i, 'gyrz': 30.0 + i,
                   'accx': 40.0 + i, 'accy': 50.0 + i, 'accz': 60.0 + i,
                   'vx': 1.0, 'vy': 2.0, 'vz': 3.0 + i,
                   'temp1': 290.0 + i})
    # Seven fractional digits, as the plotter writes; the reader trims two
    # characters to get the six datetime can parse.
    stamp = f'2025-10-15T13:43:{34 + i:02d}.1234567Z'
    return stamp + ',' + ','.join(repr(values[n]) for n in HEADER[1:])


@pytest.fixture(scope='module')
def flight(tmp_path_factory):
    path = tmp_path_factory.mktemp('csv') / 'plotter.csv'
    path.write_text('\n'.join(row(i) for i in range(ROWS)) + '\n')
    with warnings.catch_warnings():
        warnings.simplefilter('error')  # a deprecated read_csv option fails
        return FlightLog(str(path), dev=True, nc_level=None)


def test_every_imu_channel_is_read(flight):
    gyr_x, gyr_y, gyr_z, acc_x, acc_y, acc_z, times = flight.imu
    assert list(gyr_x) == [10.0, 11.0, 12.0, 13.0]
    assert list(gyr_y) == [20.0, 21.0, 22.0, 23.0]
    assert list(gyr_z) == [30.0, 31.0, 32.0, 33.0]
    assert list(acc_x) == [40.0, 41.0, 42.0, 43.0]
    assert list(acc_y) == [50.0, 51.0, 52.0, 53.0]
    assert list(acc_z) == [60.0, 61.0, 62.0, 63.0]
    assert len(times) == ROWS


def test_down_velocity_comes_from_vz(flight):
    speed_down = flight.rotation[2]
    assert list(speed_down.magnitude) == [3.0, 4.0, 5.0, 6.0]
    assert str(speed_down.units) == 'meter / second'


def test_other_columns_land_in_their_slots(flight):
    assert list(flight.temp[0].magnitude) == [290.0, 291.0, 292.0, 293.0]
    assert len(flight.temp) == 10 and len(flight.rh) == 9 + 1
    assert len(flight.pos[-1]) == ROWS
    assert flight.pos[-1][0].year == 2025


def test_barometer_temperature_is_celsius(flight):
    assert str(flight.pres[1].units) == 'degree_Celsius'


def test_csv_flight_has_no_parsed_groups_and_needs_a_tail_number(flight):
    assert flight.data == {}
    assert 'copterID' not in flight.serial_numbers
    with pytest.raises(KeyError):
        flight.resolve_tail_number()


def test_baro_fields_declare_celsius():
    """ArduPilot's BARO message logs Temp and GndTemp in degC."""
    declared = {f.source: f.units for f in PRESSURE.fields}
    assert declared['Temp'] == 'degC'
    assert declared['GndTemp'] == 'degC'
