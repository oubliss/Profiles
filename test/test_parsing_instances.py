"""
Instance selection and ESC motor mapping in profiles.parsing.

Current ArduPilot logs every barometer as BARO and every EKF core as XKF1,
tagged with an instance field (I, C); older logs have one of each, or name
the extras BAR2. These tests use synthetic message streams so the behaviour
is pinned without a large log, plus an optional check against a real
current-firmware flight that is skipped when the file is absent.
"""
import os

import numpy as np
import pytest

from profiles.parsing import parse
from profiles.readers import iter_bin

T0 = 1_780_000_000.0


def msg(kind, t, **data):
    return {'meta': {'type': kind, 'timestamp': T0 + t}, 'data': data}


def baro(kind, t, press, **extra):
    return msg(kind, t, Press=press, Temp=20.0, GndTemp=20.0, Alt=100.0, **extra)


def pos(t):
    return msg('POS', t, Lat=35.0, Lng=-97.0, Alt=350.0,
               RelHomeAlt=0.0, RelOriginAlt=0.0)


def ekf(kind, t, **extra):
    return msg(kind, t, VE=1.0, VN=2.0, VD=3.0, Roll=0.0, Pitch=0.0,
               Yaw=0.0, PN=0.0, PE=0.0, PD=0.0, **extra)


def test_baro_instances_keep_one_series():
    stream = []
    for k in range(5):
        stream.append(baro('BARO', k, 95000.0, I=0))
        stream.append(baro('BARO', k, 94900.0, I=1))
    pres = parse(stream)['groups']['pres']

    assert pres.sizes['pres_time'] == 5
    assert np.all(pres['pres'].values == 94900.0)       # default is instance 1
    assert pres.attrs['source_message_type'] == 'BARO'
    assert pres.attrs['source_instance'] == 1


def test_baro_instance_is_configurable():
    stream = []
    for k in range(3):
        stream.append(baro('BARO', k, 95000.0, I=0))
        stream.append(baro('BARO', k, 94900.0, I=1))
    pres = parse(stream, instances={'pres': 0})['groups']['pres']

    assert np.all(pres['pres'].values == 95000.0)
    assert pres.attrs['source_instance'] == 0


def test_legacy_bar2_still_beats_baro_without_instance_field():
    stream = []
    for k in range(4):
        stream.append(baro('BARO', k, 90000.0))
        stream.append(baro('BAR2', k, 90010.0))
    pres = parse(stream)['groups']['pres']

    assert pres.sizes['pres_time'] == 4
    assert np.all(pres['pres'].values == 90010.0)
    assert pres.attrs['source_message_type'] == 'BAR2'
    assert 'source_instance' not in pres.attrs


def test_legacy_baro_alone_has_no_instance_and_is_kept():
    pres = parse([baro('BARO', k, 90000.0) for k in range(3)])['groups']['pres']
    assert pres.sizes['pres_time'] == 3


def test_xkf1_keeps_primary_core_only():
    stream = [ekf('XKF1', k, C=core) for k in range(4) for core in range(3)]
    rotation = parse(stream)['groups']['rotation']

    assert rotation.sizes['rotation_time'] == 4
    assert rotation.attrs['source_message_type'] == 'XKF1'
    assert rotation.attrs['source_instance'] == 0
    times = rotation['rotation_time'].values
    assert np.all(np.diff(times) > np.timedelta64(0, 'ns'))


def test_xkf1_core_is_configurable():
    stream = [ekf('XKF1', k, C=core) for k in range(4) for core in range(3)]
    rotation = parse(stream, instances={'rotation': 2})['groups']['rotation']
    assert rotation.sizes['rotation_time'] == 4
    assert rotation.attrs['source_instance'] == 2


def test_nkf1_without_core_field_is_unfiltered():
    rotation = parse([ekf('NKF1', k) for k in range(6)])['groups']['rotation']
    assert rotation.sizes['rotation_time'] == 6
    assert rotation.attrs['source_message_type'] == 'NKF1'
    assert 'source_instance' not in rotation.attrs


def test_imu_first_instance_only_is_still_enforced():
    def imu(t, i):
        return msg('IMU', t, I=i, GyrX=0.0, GyrY=0.0, GyrZ=0.0,
                   AccX=0.0, AccY=0.0, AccZ=0.0)

    stream = [imu(k, i) for k in range(3) for i in range(3)]
    assert parse(stream)['groups']['imu'].sizes['imu_time'] == 3


def esc(t, instance, rpm):
    return msg('ESC', t, Instance=instance, RPM=rpm)


def test_esc_instances_map_to_motors_in_sorted_order():
    instances = (8, 9, 11, 12)
    stream = []
    for k in range(5):
        for n, inst in enumerate(instances):
            stream.append(esc(k + 0.001 * n, inst, 1000.0 * (n + 1) + k))
    rpm = parse(stream)['rpm']

    assert len(rpm) == 5                                  # 4 motors + times
    for n in range(4):
        assert rpm[n] == [1000.0 * (n + 1) + k for k in range(5)]
    assert len(rpm[-1]) == 5


def test_esc_legacy_instances_zero_to_three_unchanged():
    stream = [esc(k, i, 100.0 * i + k) for k in range(3) for i in range(4)]
    rpm = parse(stream)['rpm']
    assert [rpm[i][1] for i in range(4)] == [1.0, 101.0, 201.0, 301.0]


def test_esc_dropped_message_does_not_shift_later_motors():
    stream = [esc(0, 8, 1.0), esc(0, 9, 2.0), esc(0, 11, 3.0), esc(0, 12, 4.0),
              esc(1, 8, 5.0), esc(1, 9, 6.0), esc(1, 12, 8.0),      # 11 lost
              esc(2, 8, 9.0), esc(2, 9, 10.0), esc(2, 11, 11.0), esc(2, 12, 12.0)]
    rpm = parse(stream)['rpm']

    assert rpm[2][0] == 3.0 and np.isnan(rpm[2][1]) and rpm[2][2] == 11.0
    assert rpm[3] == [4.0, 8.0, 12.0]


def test_esc_absent_gives_none():
    assert parse([pos(0)])['rpm'] is None


OK3DM = os.path.expanduser('~/Data/OK3DM/KAEFS_flight2862_20260519_141710.BIN')


@pytest.mark.skipif(not os.path.exists(OK3DM), reason='OK3DM flight2862 not present')
def test_current_firmware_flight_has_strictly_increasing_times():
    parsed = parse(iter_bin(OK3DM))
    for name, dataset in parsed['groups'].items():
        times = dataset[f'{name}_time'].values
        times = times[~np.isnat(times)]
        assert np.all(np.diff(times) > np.timedelta64(0, 'ns')), name

    assert parsed['groups']['pres'].attrs['source_instance'] == 1
    assert parsed['groups']['rotation'].attrs['source_instance'] == 0
    assert len(parsed['rpm']) == 5
