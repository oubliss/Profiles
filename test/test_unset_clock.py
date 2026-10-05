"""
Samples logged before the GPS clock was set (profiles-l9l.30).

Their envelope timestamp is 1970-ish, which the parser maps to NaT. Carrying
that as NaN in the time arrays broke the a0 writer (date2num), Profile._trim
and regrid_data (datetime compared with float). The rows are now dropped at
the source, whole, so every group stays aligned with its own time.
"""
import json

import numpy as np
import pytest

from profiles.flight import FlightLog
from profiles.parsing import parse

T0 = 1_780_000_000.0
CLOCK_UNSET = 12.0  # seconds since boot, as the envelope carries it


def msg(kind, t, **data):
    return {'meta': {'type': kind, 'timestamp': t}, 'data': data}


def stream(unset=3, valid=5):
    """Messages of every required type; the first `unset` of each are early."""
    out = []
    for i in range(unset + valid):
        t = CLOCK_UNSET + i if i < unset else T0 + i
        out.append(msg('IMET', t, T1=290.0 + i, T2=291.0, T3=292.0, T4=0.0,
                       R1=10e3, R2=10e3, R3=10e3, R4=0.0, Fan=1.0))
        out.append(msg('RHUM', t, H1=50.0 + i, H2=51.0, H3=52.0, H4=0.0,
                       T1=290.0, T2=290.0, T3=290.0, T4=0.0))
        out.append(msg('POS', t, Lat=35.0, Lng=-97.0, Alt=300.0 + i,
                       RelHomeAlt=float(i), RelOriginAlt=float(i)))
        out.append(msg('BARO', t, I=1, Press=97000.0 - i, Temp=21.5,
                       GndTemp=20.0, Alt=float(i)))
        out.append(msg('NKF1', t, VE=1.0, VN=2.0, VD=3.0, Roll=0.0,
                       Pitch=0.0, Yaw=0.0, PN=0.0, PE=0.0, PD=0.0))
    out.append(msg('EV', CLOCK_UNSET, Id=10))
    out.append(msg('EV', T0 + 2, Id=11))
    out.append(msg('MSG', CLOCK_UNSET, Message='early'))
    out.append(msg('MSG', T0 + 2, Message='late'))
    return out


@pytest.fixture(scope='module')
def flight(tmp_path_factory):
    path = tmp_path_factory.mktemp('unset') / 'early.json'
    path.write_text('\n'.join(json.dumps(m) for m in stream()))
    return FlightLog(str(path), dev=True, nc_level=None)


def test_unset_clock_rows_are_dropped_whole():
    groups = parse(stream(unset=3, valid=5))['groups']
    for name, dataset in groups.items():
        times = dataset[f'{name}_time'].values
        assert not np.isnat(times).any(), name
        assert all(dataset.sizes[d] == 5 for d in dataset.dims), name
    # the surviving rows are the later ones, not a mixed-up selection
    assert list(groups['temp']['temp1'].values) == [293.0, 294.0, 295.0,
                                                    296.0, 297.0]
    assert list(groups['pos']['alt_MSL'].values) == [303.0, 304.0, 305.0,
                                                     306.0, 307.0]


def test_all_unset_group_is_absent():
    assert parse(stream(unset=4, valid=0))['groups'] == {}


def test_events_and_messages_drop_unset_entries():
    parsed = parse(stream())
    assert parsed['events'][0] == [11]
    assert parsed['messages'][0] == ['late']
    assert not np.isnat(np.array(parsed['events'][1])).any()


def test_time_arrays_hold_only_datetimes(flight):
    for series in (flight.temp, flight.rh, flight.pos, flight.pres,
                   flight.rotation):
        assert len(series[-1]) == 5
        assert all(hasattr(stamp, 'year') for stamp in series[-1])
    assert len(flight.events[1]) == 1
    assert all(hasattr(stamp, 'year') for stamp in flight.events[1])


def test_a0_writer_and_reader_survive_early_samples(flight, tmp_path):
    nc_path = tmp_path / 'early_a0.nc'
    flight._save_netCDF(str(nc_path))
    reloaded = FlightLog(str(nc_path), dev=True)
    for original, back in ((flight.temp, reloaded.temp),
                           (flight.pos, reloaded.pos),
                           (flight.pres, reloaded.pres)):
        assert list(back[-1]) == list(original[-1])
        for a, b in zip(original[:-1], back[:-1]):
            np.testing.assert_allclose(getattr(a, 'magnitude', a),
                                       getattr(b, 'magnitude', b))


def test_time_comparisons_work(flight):
    """_trim and regrid_data compare time arrays against datetimes."""
    cutoff = flight.pos[-1][2]
    stamps = np.array(flight.temp[-1])
    assert (stamps > cutoff).sum() == 2


def test_nat_is_refused_rather_than_turned_into_nan():
    from profiles.flight import _as_datetimes
    with pytest.raises(ValueError, match='NaT'):
        _as_datetimes(np.array(['2026-01-01', 'NaT'], dtype='datetime64[ns]'))
