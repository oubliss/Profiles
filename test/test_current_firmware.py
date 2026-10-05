"""
A trimmed current-firmware (2026) log, flight2862, through the whole pipeline.

Current firmware differs from the 2021 reference flight in every way the
parser cares about: barometers are BARO tagged with I, EKF cores are XKF1
tagged with C, ESC instances are output channels rather than 0-3, and there
are no USER_SENSORS parameters (the thermistors are calibrated onboard).

test/data/flight2862_ascent_trim.BIN is a byte-for-byte excerpt of
KAEFS_flight2862_20260519_141710.BIN: every FMT/FMTU/UNIT/MULT/PARM/MSG/EV
record, plus the records profiles reads for t = 322..492 s of the log (the
whole climb 0 -> ~520 m), with ESC and IMU thinned to one record in four.
It exists so these checks do not depend on the 22 MB original, which lives
outside the repository.
"""
import os
import shutil

import numpy as np
import pytest

from profiles.Coef_Manager import OnboardCalibration
from profiles.flight import FlightLog
from profiles.parsing import parse
from profiles.processing import ProcessingConfig, all_profiles, process_flights
from profiles.readers import iter_bin
from test import BASE_TEST_PATH, COEF_PATH

FIXTURE = BASE_TEST_PATH / 'flight2862_ascent_trim.BIN'

#: How many samples each series has in the trimmed window (BARO and XKF1
#: log every instance, so the raw file holds two and three times as many).
N_SAMPLES = 1700

#: A tail number the test coefficient table knows a wind calibration for.
#: Passed explicitly because the copterID registry is not what is under test.
TAIL = 'FA3TANE3MF'


@pytest.fixture(scope='module')
def parsed():
    return parse(iter_bin(str(FIXTURE)))


@pytest.fixture(scope='module')
def staged(tmp_path_factory):
    path = tmp_path_factory.mktemp('fw2026') / FIXTURE.name
    shutil.copy(FIXTURE, path)
    return str(path)


def test_fixture_is_small():
    assert os.path.getsize(FIXTURE) < 3 * 1024 * 1024


def test_pymavlink_reads_the_fixture():
    from pymavlink import mavutil
    log = mavutil.mavlink_connection(str(FIXTURE))
    seen = set()
    while True:
        message = log.recv_match(blocking=False)
        if message is None:
            break
        seen.add(message.get_type())

    assert {'BARO', 'XKF1', 'ESC', 'IMU', 'IMET', 'RHUM', 'POS',
            'WIND', 'PARM'} <= seen
    # The FMT table is whole even though the bulk records were dropped.
    assert 'ISBD' in log.name_to_id and 'XKY0' in log.name_to_id
    assert 'ISBD' not in seen


def test_parse_succeeds_with_every_group(parsed):
    assert set(parsed['groups']) == {'temp', 'rh', 'pos', 'pres', 'rotation',
                                     'wind', 'imu'}
    for name, dataset in parsed['groups'].items():
        times = dataset[f'{name}_time'].values
        assert not np.isnat(times).any(), name
        assert np.all(np.diff(times) > np.timedelta64(0, 'ns')), name


def test_one_pressure_series_not_one_per_barometer(parsed):
    pres = parsed['groups']['pres']
    assert pres.sizes['pres_time'] == N_SAMPLES    # raw file has 2x that
    assert pres.attrs['source_message_type'] == 'BARO'
    assert pres.attrs['source_instance'] == 1
    assert np.all(np.isfinite(pres['pres'].values))


def test_one_ekf_core_is_kept(parsed):
    rotation = parsed['groups']['rotation']
    assert rotation.sizes['rotation_time'] == N_SAMPLES   # 3 cores logged
    assert rotation.attrs['source_instance'] == 0
    assert parsed['groups']['imu'].attrs['source_instance'] == 0


def test_four_rpm_motors_on_remapped_esc_channels(parsed):
    rpm = parsed['rpm']
    assert len(rpm) == 5                   # four motors plus the time row
    lengths = {len(row) for row in rpm}
    assert len(lengths) == 1 and lengths.pop() > 100
    for motor in rpm[:4]:
        assert np.nanmax(motor) > 1000     # every motor actually spun


def test_no_user_sensors_serials(parsed):
    assert set(parsed['serial_numbers']) == {'copterID'}


def test_flight_selects_onboard_calibration(staged):
    flight = FlightLog(staged, dev=True, nc_level=None,
                       coefficient_dir=str(COEF_PATH))
    assert isinstance(flight.calibration_source, OnboardCalibration)
    assert flight.start_time.date().isoformat() == '2026-05-19'


def test_process_flights_makes_a_profile(staged):
    config = ProcessingConfig(dev=True, nc_level=None, tail_number=TAIL,
                              coefficient_dir=str(COEF_PATH), min_levels=20)
    results = process_flights([staged], config)
    assert all(result.ok for result in results), \
        [result.error for result in results]

    profiles = all_profiles(results)
    assert len(profiles) == 1
    profile = profiles[0]

    assert profile.tail_number == TAIL
    assert 'calibrated onboard' in profile.calibration_record[
        'temperature_source']
    assert profile.calibration_record['wind']['SerialNumber'] == TAIL

    # An ascent through ~500 m of a May afternoon: plausible and finite.
    temp = profile.temp.magnitude
    assert profile.n_populated_levels >= 20
    assert np.isfinite(temp).sum() >= 20
    assert 260 < np.nanmin(temp) and np.nanmax(temp) < 320
    alt = profile.alt.magnitude
    assert np.all(np.diff(alt) > 0)
    assert alt[-1] - alt[0] > 300
    assert np.isfinite(profile.theta.magnitude).any()
    assert np.isfinite(profile.u.magnitude).any()
