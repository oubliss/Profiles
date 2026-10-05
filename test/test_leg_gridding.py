"""
Leg detection and gridding: descents, pressure grids, phantom legs, and the
constructor's own leg finder.

All of these were reachable through documented options and all of them
raised (descent, hPa, an out-of-range profile_num) or silently produced a
mostly-empty profile (the phantom leg at the top of flight616).
"""
import builtins
import warnings
from datetime import datetime, timedelta

import numpy as np
import pytest

import profiles.utils as utils
from profiles.flight import FlightLog
from profiles.Profile import Profile
from profiles.processing import (ProcessingConfig, process_flights,
                                 profiles_from_flight)
from profiles.unit_registry import units
from test.harness import staged_bin

GRIDDED = ('alt', 'pres', 'temp', 'rh', 'theta', 'T_d', 'mixing_ratio', 'q',
           'speed', 'dir', 'u', 'v', 'lat', 'lon', 'alt_MSL', 'time')


@pytest.fixture(scope='module')
def bin_path(tmp_path_factory):
    return str(staged_bin(tmp_path_factory.mktemp('legs')))


@pytest.fixture(scope='module')
def flight(bin_path):
    return FlightLog(bin_path, False, nc_level=None)


def _compute(profile):
    profile.compute_thermo()
    profile.compute_wind()
    return profile


def _all_same_length(profile):
    n = len(profile.gridded_centers)
    lengths = {name: len(getattr(profile, name)) for name in GRIDDED}
    assert set(lengths.values()) == {n}, lengths
    assert len(profile.gridded_times) == n + 1
    assert len(profile.gridded_base) == n + 1


# ---- descent ------------------------------------------------------------

def test_descent_grids_a_real_profile(bin_path, flight):
    config = ProcessingConfig(ascent=False)
    found = profiles_from_flight(bin_path, config, flight=flight)

    # flight616's descent is the *second* detected leg; the first "descent"
    # is a 27 m dip below the peak and is below the minimum extent.
    assert len(found) == 1
    profile = _compute(found[0])

    assert len(profile.gridded_centers) > 100
    _all_same_length(profile)

    # Flight order: the first level is the top, altitude falls, time rises.
    alt = profile.alt.magnitude
    assert np.all(np.diff(alt) < 0)
    assert np.all(np.diff(profile.gridded_base.magnitude) < 0)
    assert alt[0] > 1500 and alt[-1] < 400
    times = profile.gridded_times
    assert all(b > a for a, b in zip(times, times[1:]))

    temp = profile.temp.to('degC').magnitude
    assert np.all(np.isfinite(temp))
    assert 15 < temp.min() and temp.max() < 30
    assert profile.n_populated_levels == len(profile.gridded_centers)


def test_descent_shares_the_ascent_lattice(bin_path, flight):
    up = profiles_from_flight(bin_path, ProcessingConfig(), flight=flight)[0]
    down = profiles_from_flight(
        bin_path, ProcessingConfig(ascent=False), flight=flight,
        base_start=up._base_start)[0]
    lattice = (down.gridded_base.magnitude - up._base_start.magnitude) / 10
    np.testing.assert_allclose(lattice, np.round(lattice), atol=1e-9)
    assert down._base_start == up._base_start


def test_regrid_base_descent_on_a_synthetic_leg():
    """Descending 100 m at 1 m/s, sampled at 1 Hz."""
    t0 = datetime(2024, 1, 1)
    times = [t0 + timedelta(seconds=i) for i in range(101)]
    alt = (100.0 - np.arange(101)) * units.m

    centres_t, centres, edge_t, edges = utils.regrid_base(
        base=alt, base_times=times, new_res=10 * units.m, ascent=False,
        units=units, indices=(times[0], times[-1]))

    assert len(centres) == len(centres_t)
    assert len(edge_t) == len(edges) == len(centres) + 1
    assert np.all(np.diff(edges.magnitude) < 0)
    assert all(b > a for a, b in zip(edge_t, edge_t[1:]))
    # An edge at height h is passed at t = 100 - h on the way down.
    for h, when in zip(edges.magnitude, edge_t):
        assert abs((when - t0).total_seconds() - (100 - h)) <= 1

    # regrid_data's (start, end] bins average the right samples.
    gridded = utils.regrid_data(alt, times, edge_t, units)
    np.testing.assert_allclose(gridded.magnitude, centres.magnitude,
                               atol=1.01)


# ---- pressure grids -----------------------------------------------------

@pytest.mark.parametrize('ascent', [True, False])
@pytest.mark.parametrize('res_units,resolution', [('hPa', 10), ('Pa', 1000)])
def test_pressure_grid(bin_path, flight, ascent, res_units, resolution):
    config = ProcessingConfig(resolution=resolution, res_units=res_units,
                              ascent=ascent)
    found = profiles_from_flight(bin_path, config, flight=flight)
    assert len(found) == 1
    profile = _compute(found[0])

    _all_same_length(profile)
    # 10 hPa is ~100 m; the profile is ~1400 m deep.
    assert 12 <= len(profile.gridded_centers) <= 16
    spacing = np.diff(profile.gridded_base.to('hPa').magnitude)
    np.testing.assert_allclose(np.abs(spacing), 10, atol=1e-6)
    # Pressure falls on the way up, rises on the way down.
    sign = -1 if ascent else 1
    assert np.all(sign * spacing > 0)
    assert np.all(np.isfinite(profile.alt.magnitude))
    assert 15 < profile.temp.to('degC').magnitude.min()


def test_pressure_grid_ignores_a_height_start_with_a_warning(bin_path, flight):
    config = ProcessingConfig(resolution=10, res_units='hPa',
                              profile_start_height=350)
    with pytest.warns(UserWarning, match='profile_start_height'):
        found = profiles_from_flight(bin_path, config, flight=flight)
    assert found[0].gridded_base[0].to('hPa').magnitude > 900
    results = process_flights([bin_path], config, on_error='raise')
    assert results[0].ok and len(results[0].profiles) == 1


def test_nearest_index_maps_between_clocks():
    base = datetime(2024, 1, 1)
    times = [base + timedelta(seconds=s) for s in (0, 1, 2, 3)]
    assert utils.nearest_index(times, times[2]) == 2
    assert utils.nearest_index(times, base + timedelta(seconds=1.4)) == 1
    assert utils.nearest_index(times, base + timedelta(seconds=1.6)) == 2
    assert utils.nearest_index(times, base - timedelta(seconds=5)) == 0
    assert utils.nearest_index(times, base + timedelta(seconds=9)) == 3


# ---- phantom legs -------------------------------------------------------

def test_phantom_leg_is_not_detected(flight):
    everything = flight.find_legs(min_extent=0)
    kept = flight.find_legs()
    assert len(everything) == 2
    assert len(kept) == 1
    assert kept[0] == everything[0]

    start, peak, end = everything[1]
    assert (peak - start) < timedelta(seconds=5)


def test_min_extent_is_configurable(bin_path, flight):
    # Detection filter off: the phantom is reported again. On its own grid
    # it is 0 levels deep and is skipped rather than returned...
    config = ProcessingConfig(min_leg_extent=0)
    with pytest.warns(UserWarning, match='skipping profile 2'):
        assert len(profiles_from_flight(bin_path, config, flight=flight)) == 1

    # ...and on a forced grid it is a (mostly empty) profile.
    config = ProcessingConfig(min_leg_extent=0, profile_start_height=350)
    assert len(profiles_from_flight(bin_path, config, flight=flight)) == 2


def test_min_levels_counts_levels_with_data(bin_path):
    """With the detector filter off, min_levels alone must still drop the
    phantom: on a common grid it spans ~140 levels, 119 of them empty."""
    config = ProcessingConfig(min_leg_extent=0, min_levels=20,
                              profile_start_height=350)
    results = process_flights([bin_path], config, on_error='raise')
    assert len(results[0].profiles) == 1

    config = ProcessingConfig(min_leg_extent=0, profile_start_height=350)
    results = process_flights([bin_path], config, compute=False,
                              on_error='raise')
    phantom = results[0].profiles[1]
    assert len(phantom.gridded_centers) > 100
    assert phantom.n_populated_levels < 20


def test_default_config_yields_no_phantom(bin_path):
    config = ProcessingConfig(min_levels=20, profile_start_height=350)
    results = process_flights([bin_path], config, on_error='raise')
    assert len(results[0].profiles) == 1
    str(results[0].profiles[0])


def test_never_returns_a_zero_level_profile(bin_path, flight):
    # A grid whose floor is above the leg's top has no levels.
    config = ProcessingConfig(profile_start_height=5000)
    with pytest.warns(UserWarning, match='skipping profile'):
        found = profiles_from_flight(bin_path, config, flight=flight)
    assert found == []

    with pytest.raises(ValueError, match='no levels'):
        Profile(bin_path, 10, 'm', 1, raw_profile=flight,
                base_start=5000 * units.m)


def test_entry_points_grid_identically(bin_path, flight):
    config = ProcessingConfig(profile_start_height=350, min_levels=20)
    direct = profiles_from_flight(bin_path, config, flight=flight)
    batch = process_flights([bin_path], config, compute=False,
                            on_error='raise')[0].profiles

    assert len(direct) == len(batch) == 1
    assert len(direct[0].gridded_centers) == len(batch[0].gridded_centers) == 142
    np.testing.assert_array_equal(direct[0].gridded_base.magnitude,
                                  batch[0].gridded_base.magnitude)
    assert direct[0].gridded_times == batch[0].gridded_times


# ---- the constructor's own leg finder -----------------------------------

def test_profile_without_index_list_finds_the_same_legs(bin_path, flight):
    direct = Profile(bin_path, 10, 'm', 1, raw_profile=flight,
                     confirm_bounds=False)
    batch = profiles_from_flight(bin_path, ProcessingConfig(), flight=flight)
    assert direct.indices == (batch[0].indices[0], batch[0].indices[1])
    np.testing.assert_array_equal(direct.alt.magnitude, batch[0].alt.magnitude)


@pytest.mark.parametrize('ascent', [True, False])
def test_out_of_range_profile_num_raises_promptly(bin_path, flight, ascent):
    with pytest.raises(IndexError, match=r'profile_num=2 but only 1 '):
        Profile(bin_path, 10, 'm', 2, ascent=ascent, raw_profile=flight)


def test_profile_num_is_one_based(bin_path, flight):
    with pytest.raises(ValueError, match='1-based'):
        Profile(bin_path, 10, 'm', 0, raw_profile=flight)


def test_identify_profile_forwards_its_arguments(monkeypatch):
    """The re-selection path passed to_return positionally into
    confirm_bounds. Drive it with a scripted user: reject once, accept."""
    t0 = datetime(2024, 1, 1)
    n = 400
    times = [t0 + timedelta(seconds=i) for i in range(n)]
    shape = np.concatenate([np.linspace(0, 200, 150), np.full(50, 200.0),
                            np.linspace(200, 0, 150), np.zeros(50)])
    alts = shape * units.m

    asked = []

    def user(prompt=''):
        asked.append(prompt)
        assert len(asked) < 12, 'identify_profile keeps asking'
        if 'start height' in prompt:
            return '20'
        return 'n' if sum('Correct' in a for a in asked) == 1 else 'y'

    monkeypatch.setattr(builtins, 'input', user)

    # Before the fix the reject branch raised: to_return (a list) landed in
    # confirm_bounds and the start height was dropped.
    legs = utils.identify_profile(alts, times, confirm_bounds=True,
                                  profile_start_height=10 * units.m)
    assert len(legs) >= 1
    start, peak, end = legs[0]
    assert start < peak < end
    assert any('start height' in a for a in asked)
