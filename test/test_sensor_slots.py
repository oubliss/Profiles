"""
Per-sensor series are slot-aligned (profiles-l9l.31).

calibration._sensor_series used to skip absent keys, so with sensor 2
missing the list was [s1, s3, s4] and everything that used the list index as
the sensor number was off by one: calibrate_temperature paired resistance
index 1 (sensor 3) with serial imet2, and the a0 writer stored sensor 3's
calibrated temperature as calib_temp2.
"""
import numpy as np
import pytest
from metpy.units import units as u

from profiles import calibration, qc
from profiles.flight import FlightLog
from test.test_unset_clock import stream


class RecordingSource:
    """A calibration source that records which serial it was asked for."""

    def __init__(self, temperature_from='table'):
        self.temperature_from = temperature_from
        self.serials = []

    def get_coefs(self, kind, serial, when=None):
        self.serials.append(serial)
        # Steinhart-Hart A, B, C; the serial is irrelevant to the arithmetic.
        return {'A': 1.0e-3, 'B': 2.5e-4, 'C': 1.0e-7}


def thermo(missing=(), n=6):
    data = {}
    for number in range(1, 5):
        if number in missing:
            continue
        data[f'temp{number}'] = np.full(n, 280.0 + number) * u.kelvin
        data[f'resi{number}'] = np.full(n, 1.0e4 * number) * u.ohm
        data[f'rh{number}'] = np.full(n, 40.0 + number) * u.percent
    return data


SERIALS = {'imet1': 101, 'imet2': 102, 'imet3': 103, 'imet4': 104,
           'rh1': 1, 'rh2': 2, 'rh3': 3, 'rh4': 4}


def test_absent_middle_sensor_keeps_its_slot():
    series = calibration._sensor_series(thermo(missing=(2,)), 'temp')
    assert len(series) == 4
    assert series[0][0] == 281.0
    assert np.isnan(series[1]).all() and len(series[1]) == 6
    assert series[2][0] == 283.0
    assert series[3][0] == 284.0


def test_no_sensor_of_a_kind_gives_empty_list():
    assert calibration._sensor_series({}, 'temp') == []


def test_resistances_pair_with_the_right_serial():
    source = RecordingSource()
    record = {}
    result = calibration.calibrate_temperature(
        thermo(missing=(2,)), SERIALS, record=record, source=source)
    assert source.serials == [101, 103, 104]
    assert len(result) == 4
    assert np.isnan(result[1]).all()
    assert 'imet2' not in record
    assert {'imet1', 'imet3', 'imet4'} <= set(record)


def test_logged_temperature_slots_align():
    source = RecordingSource(temperature_from='logged')
    result = calibration.calibrate_temperature(
        thermo(missing=(2,)), SERIALS, source=source)
    assert [float(s[0]) for s in result[::2]] == [281.0, 283.0]
    assert np.isnan(result[1]).all()


def test_humidity_slots_align():
    result = calibration.calibrate_humidity(thermo(missing=(1,)), SERIALS)
    assert np.isnan(result[0]).all()
    assert [float(s[0]) for s in result[1:]] == [42.0, 43.0, 44.0]


def test_ensemble_qc_treats_absent_slot_as_empty():
    series = calibration._sensor_series(thermo(missing=(2,)), 'temp')
    flags = qc.qc(series, 5.0, 5.0)
    assert flags[1] == qc.EMPTY
    assert flags[0] == flags[2] == flags[3] == qc.GOOD


@pytest.fixture
def flight(tmp_path):
    import json
    path = tmp_path / 'slots.json'
    path.write_text('\n'.join(json.dumps(m) for m in stream(unset=0)))
    return FlightLog(str(path), dev=True, nc_level=None)


def test_a0_writes_calibrated_values_to_the_right_slot(flight, tmp_path):
    import netCDF4
    n = len(flight.temp[-1])
    flight.calib_temp = [np.full(n, 281.0), np.full(n, np.nan),
                         np.full(n, 283.0), np.full(n, np.nan)]
    flight.calib_rh = [np.full(n, 51.0), np.full(n, 52.0),
                       np.full(n, np.nan), np.full(n, np.nan)]
    path = tmp_path / 'slots_a0.nc'
    flight._save_netCDF(str(path))

    with netCDF4.Dataset(str(path)) as dataset:
        temp = dataset['temp'].variables
        assert sorted(v for v in temp if v.startswith('calib_temp')) \
            == ['calib_temp1', 'calib_temp3']
        assert float(temp['calib_temp3'][0]) == 283.0
        rh = dataset['rh'].variables
        assert sorted(v for v in rh if v.startswith('calib_rh')) \
            == ['calib_rh1', 'calib_rh2']

    reloaded = FlightLog(str(path), dev=True)
    assert len(reloaded.calib_temp) == 4
    assert reloaded.calib_temp[0][0] == 281.0
    assert np.isnan(reloaded.calib_temp[1]).all()
    assert reloaded.calib_temp[2][0] == 283.0
    assert np.isnan(reloaded.calib_temp[3]).all()
    assert reloaded.calib_rh[1][0] == 52.0
    assert np.isnan(reloaded.calib_rh[2]).all()


def test_a0_without_calibrated_values_reads_none(flight, tmp_path):
    flight.calib_temp = flight.calib_rh = None
    path = tmp_path / 'plain_a0.nc'
    flight._save_netCDF(str(path))
    reloaded = FlightLog(str(path), dev=True)
    assert reloaded.calib_temp is None and reloaded.calib_rh is None
