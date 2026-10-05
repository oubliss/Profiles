"""
Barometer/EKF instance choice reaching FlightLog, a0 and c1 files
(profiles-l9l.29), plus the ESC unset-clock and is_equal fixes
(profiles-l9l.36). Synthetic streams only.
"""
import json
from types import SimpleNamespace

import pytest

from profiles.flight import FlightLog
from profiles.io.provenance import provenance_attributes
from profiles.parsing import parse
from profiles.processing import ProcessingConfig

from test.test_unset_clock import CLOCK_UNSET, T0, msg, stream


def two_baro_stream(instances=(0, 1), n=6):
    out = [m for m in stream(unset=0, valid=n) if m['meta']['type'] != 'BARO']
    for i in range(n):
        for inst in instances:
            out.append(msg('BARO', T0 + i, I=inst,
                           Press=90000.0 + 1000 * inst - i,
                           Temp=20.0, GndTemp=20.0, Alt=float(i)))
    return out


def write(tmp_path, messages, name='log.json'):
    path = tmp_path / name
    path.write_text('\n'.join(json.dumps(m) for m in messages))
    return str(path)


def test_config_defaults():
    config = ProcessingConfig()
    assert (config.baro_instance, config.ekf_core) == (1, 0)


@pytest.mark.parametrize('instance', [0, 1])
def test_flightlog_uses_requested_barometer(tmp_path, instance):
    path = write(tmp_path, two_baro_stream())
    flight = FlightLog(path, dev=True, nc_level=None, baro_instance=instance)
    assert flight.baro_instance == instance
    assert flight.pres[0].magnitude[0] == pytest.approx(
        90000.0 + 1000 * instance)


def test_flightlog_default_is_instance_one(tmp_path):
    flight = FlightLog(write(tmp_path, two_baro_stream()), dev=True,
                       nc_level=None)
    assert flight.baro_instance == 1


def test_missing_instance_is_an_error_naming_those_present(tmp_path):
    path = write(tmp_path, two_baro_stream(instances=(0, 2)))
    with pytest.raises(ValueError, match=r"instance 1.*present: \[0, 2\]"):
        FlightLog(path, dev=True, nc_level=None)
    with pytest.raises(ValueError, match=r"instance 5.*present: \[0, 2\]"):
        FlightLog(path, dev=True, nc_level=None, baro_instance=5)


def test_missing_ekf_core_is_an_error():
    messages = [msg('XKF1', T0 + i, C=0, VE=0.0, VN=0.0, VD=0.0, Roll=0.0,
                    Pitch=0.0, Yaw=0.0, PN=0.0, PE=0.0, PD=0.0)
                for i in range(3)]
    with pytest.raises(ValueError, match=r"rotation instance 3.*\[0\]"):
        parse(messages, instances={'rotation': 3})


def test_legacy_log_ignores_the_setting(tmp_path):
    legacy = [m for m in stream(unset=0, valid=5)
              if m['meta']['type'] != 'BARO']
    for i in range(5):
        for kind, press in (('BARO', 90000.0), ('BAR2', 90010.0)):
            legacy.append(msg(kind, T0 + i, Press=press, Temp=20.0,
                              GndTemp=20.0, Alt=0.0))
    path = write(tmp_path, legacy)
    for instance in (0, 1, 7):
        flight = FlightLog(path, dev=True, nc_level=None,
                           baro_instance=instance)
        assert flight.baro == 'BAR2'
        assert flight.baro_instance is None
        assert flight.pres[0].magnitude[0] == 90010.0


def test_instance_in_a0_and_c1(tmp_path):
    path = write(tmp_path, two_baro_stream())
    flight = FlightLog(path, dev=True, nc_level=None, baro_instance=0)
    a0 = tmp_path / 'a0.nc'
    flight._save_netCDF(str(a0))
    assert FlightLog(str(a0), dev=True).baro_instance == 0

    profile = SimpleNamespace(_raw_profile=flight, resolution=10, ascent=True)
    attrs = provenance_attributes(profile)
    assert attrs['baro_instance'] == 0
    assert attrs['baro_message_type'] == 'BARO'

    flight.baro_instance = None
    assert 'baro_instance' not in provenance_attributes(profile)


def test_esc_before_gps_clock_is_dropped():
    def esc(t, inst, rpm):
        return msg('ESC', t, Instance=inst, RPM=rpm)

    messages = [esc(CLOCK_UNSET, 0, 1.0), esc(CLOCK_UNSET, 1, 2.0),
                esc(T0, 0, 3.0), esc(T0 + 0.01, 1, 4.0),
                esc(T0 + 1, 0, 5.0), esc(T0 + 1.01, 1, 6.0)]
    rpm = parse(messages)['rpm']
    assert rpm[0] == [3.0, 5.0] and rpm[1] == [4.0, 6.0]
    assert len(rpm[-1]) == 2
    assert all(t.astype('datetime64[Y]').astype(int) + 1970 > 2000
               for t in rpm[-1])
    assert parse(messages[:2])['rpm'] is None


def test_is_equal_works_with_quantities(tmp_path):
    path = write(tmp_path, two_baro_stream())
    first = FlightLog(path, dev=True, nc_level=None)
    second = FlightLog(path, dev=True, nc_level=None)
    assert first.is_equal(second)

    other = FlightLog(path, dev=True, nc_level=None, baro_instance=0)
    assert not first.is_equal(other)
