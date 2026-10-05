"""
Grid invariants.

Every gridded variable in a Profile, Thermo_Profile and Wind_Profile must have
exactly as many points as the vertical coordinate it is indexed by. Violating
that is what put a trailing fill row into every published thermo_/wind_ file:
alt came from the N+1 bin edges while temp, rh and pres came from the N bin
averages, and both were written into the same unlimited NetCDF dimension.

The second test pins the coordinate convention itself - bin centres, not bin
edges - which was a systematic -resolution/2 bias.
"""
import warnings

import numpy as np
import pytest

from profiles import Profile_Set
from test.harness import (PROFILE_START_HEIGHT, RESOLUTION, RES_UNITS,
                          staged_bin)

PROFILE_VARS = ('alt', 'pres', 'lat', 'lon', 'alt_MSL', 'time')
THERMO_VARS = ('alt', 'pres', 'temp', 'rh', 'theta', 'T_d', 'mixing_ratio', 'q')
WIND_VARS = ('alt', 'pres', 'speed', 'dir', 'u', 'v')


@pytest.fixture(scope='module')
def processed(tmp_path_factory):
    bin_path = staged_bin(tmp_path_factory.mktemp('grid'))
    # Profile_Set is deprecated; the grid it builds is the Profile's own and
    # is what is checked here, so the warning is expected, not interesting.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        profile_set = Profile_Set.Profile_Set(
            resolution=RESOLUTION, res_units=RES_UNITS, ascent=True,
            dev=True, confirm_bounds=False, nc_level=None,
            profile_start_height=PROFILE_START_HEIGHT)
        profile_set.add_all_profiles(str(bin_path))
    # Thermo and wind now live on the Profile itself; the triple is kept
    # so the per-group assertions below still read clearly.
    return [(p, p.compute_thermo(), p.compute_wind())
            for p in profile_set.profiles]


def _lengths(obj, names):
    return {name: len(getattr(obj, name)) for name in names
            if getattr(obj, name, None) is not None}


def test_all_gridded_variables_share_one_length(processed):
    for i, (profile, thermo, wind) in enumerate(processed):
        n = len(profile.gridded_centers)

        lengths = {}
        lengths.update({f'profile.{k}': v
                        for k, v in _lengths(profile, PROFILE_VARS).items()})
        lengths.update({f'thermo.{k}': v
                        for k, v in _lengths(thermo, THERMO_VARS).items()})
        lengths.update({f'wind.{k}': v
                        for k, v in _lengths(wind, WIND_VARS).items()})

        wrong = {k: v for k, v in lengths.items() if v != n}
        assert not wrong, (
            f'profile {i}: expected every gridded variable to have {n} points '
            f'(len(gridded_centers)); these differ: {wrong}')


def test_edges_are_one_longer_than_centres(processed):
    """gridded_times/gridded_base delimit the bins; there are N+1 of them."""
    for i, (profile, _, _) in enumerate(processed):
        n = len(profile.gridded_centers)
        assert len(profile.gridded_times) == n + 1, f'profile {i}: times'
        assert len(profile.gridded_base) == n + 1, f'profile {i}: base'


def test_child_altitude_is_bin_centres_not_edges(processed):
    """Thermo/Wind alt must match Profile's centres, not its edges.

    Handing them gridded_base put every value half a bin low - 5 m at the
    10 m resolution used here.
    """
    for i, (profile, thermo, wind) in enumerate(processed):
        centres = profile.gridded_centers.magnitude
        edges = profile.gridded_base.magnitude

        np.testing.assert_allclose(
            thermo.alt.magnitude, centres, rtol=0, atol=1e-9,
            err_msg=f'profile {i}: thermo.alt is not the bin centres')
        np.testing.assert_allclose(
            wind.alt.magnitude, centres, rtol=0, atol=1e-9,
            err_msg=f'profile {i}: wind.alt is not the bin centres')

        # And confirm the two conventions really do differ by half a bin,
        # so this test is not silently vacuous.
        np.testing.assert_allclose(centres - edges[:len(centres)],
                                   RESOLUTION / 2.0, rtol=0, atol=1e-9)


def test_saved_netcdf_has_no_fill_row(processed, tmp_path):
    """The written file must not stretch past the data.

    time came from the N+1 bin edges while every variable had N points, so
    the unlimited dimension grew to N+1 and NetCDF padded every variable with
    a masked value. That trailing row is in published thermo_/wind_ files.
    """
    import netCDF4

    profile, _, _ = processed[0]
    writers = {'thermo': profile._save_thermo_netCDF,
               'wind': profile._save_wind_netCDF}
    for label, write in writers.items():
        path = tmp_path / f'{label}.cdf'
        write(str(path))

        with netCDF4.Dataset(path) as handle:
            n_time = len(handle.dimensions['time'])
            assert n_time == len(profile.gridded_centers), (
                f'{label}: time dimension is {n_time}, expected '
                f'{len(profile.gridded_centers)}')

            padded = [name for name, var in handle.variables.items()
                      if var[:].ndim == 1 and var[:].size == n_time
                      and np.ma.is_masked(var[:][-1])]
            assert not padded, f'{label}: fill row in {padded}'
