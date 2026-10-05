"""
Import hygiene, one binning rule, single parse per flight, and Meta.

Each of these was a review finding: importing profiles changed the host
process's warning filters and selected a matplotlib backend; lat/lon were
binned [start, end) while everything else was (start, end]; the first flight
of a batch was parsed twice; Meta needed an optional package at import and
stopped reading after the first replaced field.
"""
import subprocess
import sys
from datetime import datetime, timedelta

import numpy as np
import pytest

import profiles.processing as processing
import profiles.utils as utils
from profiles.flight import FlightLog
from profiles.processing import ProcessingConfig, process_flights
from profiles.unit_registry import units
from test.harness import RESOLUTION, RES_UNITS, staged_bin


# --- import hygiene -------------------------------------------------------

def test_import_leaves_process_state_alone():
    code = (
        "import sys, warnings\n"
        # Third-party libraries add their own filters when first imported;
        # only what profiles adds on top is under test.
        "import metpy.units, pint, pandas, scipy.signal, requests, xarray, "
        "netCDF4\n"
        "before = list(warnings.filters)\n"
        "import profiles, profiles.utils, profiles.processing, "
        "profiles.Profile, profiles.Meta\n"
        "assert warnings.filters == before, 'warning filters changed'\n"
        "assert 'matplotlib.pyplot' not in sys.modules, 'pyplot imported'\n")
    done = subprocess.run([sys.executable, '-c', code],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


def test_regrid_data_empty_bin_is_nan_without_leaking_a_filter():
    import warnings
    t0 = datetime(2020, 1, 1)
    data = np.array([1., 2.]) * units.m
    times = [t0 + timedelta(seconds=1), t0 + timedelta(seconds=2)]
    grid = [t0, t0 + timedelta(seconds=2), t0 + timedelta(seconds=10),
            t0 + timedelta(seconds=20)]
    before = list(warnings.filters)
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        out = utils.regrid_data(data, times, grid, units)
    assert np.isnan(out.magnitude[1:]).all()
    assert out.magnitude[0] == 1.5
    assert warnings.filters == before


# --- one binning rule -----------------------------------------------------

def test_group_bins_match_regrid_data_bins():
    t0 = datetime(2020, 1, 1)
    times = [t0 + timedelta(seconds=s) for s in range(0, 12)]
    data = np.arange(12.) * units.m
    # The last gridded time is past every sample, and bin 2 is empty.
    grid = [t0, t0 + timedelta(seconds=3), t0 + timedelta(seconds=3.5),
            t0 + timedelta(seconds=20)]
    groups = list(utils.regrid_data_group(data, times, grid, units))
    gridded = utils.regrid_data(data, times, grid, units)
    assert len(groups) == len(gridded) == len(grid) - 1
    for group, value in zip(groups, gridded.magnitude):
        mags = group['values'].magnitude
        mean = np.nanmean(mags) if len(mags) else np.nan
        assert (np.isnan(mean) and np.isnan(value)) or mean == value
    # (start, end]: the sample at exactly t0+3s belongs to bin 0, not bin 1.
    assert list(groups[0]['values'].magnitude) == [1., 2., 3.]


@pytest.fixture(scope='module')
def bin_path(tmp_path_factory):
    return str(staged_bin(tmp_path_factory.mktemp('hygiene')))


def test_position_has_the_same_length_as_every_other_variable(bin_path):
    config = ProcessingConfig(resolution=RESOLUTION, res_units=RES_UNITS,
                              dev=True, nc_level=None)
    profile = process_flights([bin_path], config, compute=False,
                              on_error='raise')[0].profiles[0]
    n = len(profile.alt)
    assert len(profile.lat) == len(profile.lon) == len(profile.alt_MSL) == n
    assert len(profile.pres) == n


# --- one parse per flight -------------------------------------------------

def test_first_flight_is_parsed_once(bin_path, monkeypatch):
    parsed = []

    class Counting(FlightLog):
        def __init__(self, *args, **kwargs):
            parsed.append(args[0])
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(processing, 'FlightLog', Counting)
    # No profile_start_height, so the first flight is gridded a second time
    # once its own grid is known - that must not re-read the BIN.
    config = ProcessingConfig(resolution=RESOLUTION, res_units=RES_UNITS,
                              dev=True, nc_level=None)
    results = process_flights([bin_path], config, compute=False,
                              on_error='raise')
    assert results[0].profiles
    assert parsed == [bin_path]


# --- Meta -----------------------------------------------------------------

def test_meta_imports_without_dronelogbook():
    code = (
        "import sys\n"
        "sys.modules['dronelogbook'] = None\n"
        "from profiles import Meta\n"
        "Meta.Meta()\n")
    done = subprocess.run([sys.executable, '-c', code],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


def _csv(path, **fields):
    import pandas as pd
    pd.DataFrame({k: [v] for k, v in fields.items()}).to_csv(path,
                                                             index=False)
    return str(path)


def test_meta_reads_every_field_from_the_flight_file(tmp_path, capsys):
    from profiles.Meta import Meta
    header = _csv(tmp_path / 'h.csv', timestamp='20200101_000000',
                  PIC='a', objective='wind', battery_id=7)
    flight = _csv(tmp_path / 'f.csv', PIC='b', objective='thermo',
                  battery_id=8, scoop_id='s')
    meta = Meta(header, flight)
    # Non-string values replaced without a TypeError, and the loop went on
    # past the first replaced field.
    assert meta.all_fields['PIC'] == 'b'
    assert meta.all_fields['objective'] == 'thermo'
    assert meta.all_fields['battery_id'] == 8
    assert meta.all_fields['scoop_id'] == 's'
    assert meta.all_fields['date_utc'] == '20200101'


def test_meta_timestamp_is_not_a_substring_test(tmp_path):
    from profiles.Meta import Meta
    # 'time' and 'stamp' are substrings of 'timestamp' but are not fields;
    # an equality test must still protect only the timestamp itself.
    header = _csv(tmp_path / 'h.csv', timestamp='20200101_000000')
    flight = _csv(tmp_path / 'f.csv', timestamp='20200202_000000')
    meta = Meta(header, flight)
    assert meta.all_fields['timestamp'] == '20200202_000000'
    assert meta.all_fields['date_utc'] == '20200202'


def test_meta_without_timestamp_does_not_raise(tmp_path):
    from profiles.Meta import Meta
    meta = Meta(_csv(tmp_path / 'h.csv', PIC='a'))
    assert meta.all_fields['timestamp'] is None
    assert meta.all_fields['date_utc'] is None
