"""
Writer defects found in the second code review.

README quick start, CF attributes, the q unit, coefficient provenance, and
what a writer is allowed to touch outside the output file (network, OS
APIs, deprecated calls).
"""
import hashlib
import os
import subprocess
import warnings

import netCDF4
import numpy as np
import pytest

from profiles.io import provenance
from profiles.processing import ProcessingConfig, all_profiles, process_flights
from test.harness import (PROFILE_START_HEIGHT, RESOLUTION, RES_UNITS,
                          staged_bin)


@pytest.fixture(scope='module')
def profile(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('writer_review')
    config = ProcessingConfig(
        resolution=RESOLUTION, res_units=RES_UNITS, dev=True, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT, min_levels=50)
    return all_profiles(process_flights([str(staged_bin(tmp))], config))[0]


def _attributes(path):
    with netCDF4.Dataset(path) as handle:
        return {name: handle.getncattr(name) for name in handle.ncattrs()}


# --- README quick start / file names (l9l.12) -----------------------------

def test_cf_writer_works_with_no_metadata_and_no_path(profile):
    """The README quick start: no metadata, no explicit path."""
    assert profile.meta is None
    profile.save_cfnetcdf('N934UA', terrain_elevation=340)
    expected = f'{profile.file_path}.c1.10.cf_ascent.nc'
    try:
        assert os.path.exists(expected)
    finally:
        if os.path.exists(expected):
            os.remove(expected)


def test_the_two_combined_writers_do_not_share_a_file_name(profile):
    profile.save_netcdf()
    profile.save_cfnetcdf('N934UA', terrain_elevation=340)
    plain = f'{profile.file_path}.c1.10.ascent.nc'
    cf = f'{profile.file_path}.c1.10.cf_ascent.nc'
    try:
        assert os.path.exists(plain) and os.path.exists(cf)
        assert 'featureType' not in _attributes(plain)
        assert 'featureType' in _attributes(cf)
    finally:
        for path in (plain, cf):
            if os.path.exists(path):
                os.remove(path)


def test_names_differ_with_metadata_too(profile, monkeypatch):
    monkeypatch.setattr(profile, 'meta',
                        {'location': 'Some Place', 'platform_id': 'P1'})
    plain = profile._c1_output_path('/out/x', profile._ascent_filename_tag)
    cf = profile._c1_output_path('/out/x', 'cf_' + profile._ascent_filename_tag)
    assert plain != cf
    assert os.path.dirname(plain) == '/out'


# --- save_cfnetcdf attributes (l9l.13) ------------------------------------

@pytest.fixture(scope='module')
def cf_file(profile, tmp_path_factory):
    path = tmp_path_factory.mktemp('cf') / 'cf.nc'
    profile.save_cfnetcdf('N934UA', 340, str(path))
    return path


def test_both_references_are_kept(cf_file):
    attributes = _attributes(cf_file)
    assert attributes['Reference1'].startswith('Segales')
    assert attributes['Reference2'].startswith('Bell')


def test_pressure_and_altitude_use_cf_standard_names(cf_file):
    with netCDF4.Dataset(cf_file) as handle:
        pressure = handle.variables['pressure']
        altitude = handle.variables['altitude']
        assert pressure.standard_name == 'air_pressure'
        assert 'positive' not in pressure.ncattrs()
        assert 'axis' not in pressure.ncattrs()
        assert altitude.standard_name == 'altitude'
        assert altitude.positive == 'up'
        assert altitude.axis == 'Z'


def test_exactly_one_vertical_axis(cf_file):
    with netCDF4.Dataset(cf_file) as handle:
        axes = [name for name, variable in handle.variables.items()
                if getattr(variable, 'axis', None) == 'Z']
    assert axes == ['altitude']


def test_trajectory_feature_type_has_a_cf_role_variable(cf_file):
    with netCDF4.Dataset(cf_file) as handle:
        assert handle.featureType == 'trajectory'
        roles = [name for name, variable in handle.variables.items()
                 if getattr(variable, 'cf_role', None) == 'trajectory_id']
        assert roles == ['trajectory_id']
        assert handle.variables['trajectory_id'].dimensions == ()
        assert str(handle.variables['trajectory_id'][...]).startswith('N934UA_')


def test_time_keeps_sub_second_resolution(profile, cf_file):
    with netCDF4.Dataset(cf_file) as handle:
        time = handle.variables['time']
        assert time.dtype == np.dtype('f8')
        expected = netCDF4.date2num(profile.time,
                                    units='seconds since 1970-01-01T00:00:00')
        np.testing.assert_allclose(time[:], expected, atol=1e-6)


def test_unresolved_identity_does_not_raise(profile, tmp_path, monkeypatch):
    """netCDF4 cannot store None; both writers used to raise on it."""
    monkeypatch.setattr(profile, 'tail_number', None)
    monkeypatch.setattr(profile, 'copter_id', None)
    profile.save_cfnetcdf('N934UA', 340, str(tmp_path / 'a.nc'))
    profile.save_netcdf(str(tmp_path / 'b.nc'))
    for name in ('a.nc', 'b.nc'):
        attributes = _attributes(tmp_path / name)
        # the provenance block's own defaults
        assert attributes['copter_id'] == -999
        assert attributes['tail_number'] == 'unknown'


# --- q units (l9l.14) ------------------------------------------------------

def test_q_converts_correctly_between_units(profile):
    kg_per_kg = profile.q.to('kg/kg').magnitude
    g_per_kg = profile.q.to('g/kg').magnitude
    np.testing.assert_allclose(g_per_kg, kg_per_kg * 1e3, rtol=1e-12)
    # Plausible surface humidity: a few g/kg, not thousandths of one.
    assert 1.0 < np.nanmean(g_per_kg) < 40.0
    assert 1e-3 < np.nanmean(kg_per_kg) < 4e-2


def test_q_is_consistent_with_mixing_ratio(profile):
    """q = r / (1 + r) is within r's own magnitude, so the two agree."""
    r = profile.mixing_ratio.to('dimensionless').magnitude
    np.testing.assert_allclose(profile.q.to('kg/kg').magnitude, r / (1 + r),
                               rtol=1e-6)


def test_written_q_is_in_g_per_kg(profile, tmp_path):
    path = tmp_path / 'c.nc'
    profile.save_cfnetcdf('N934UA', 340, str(path))   # no q: CF file omits it
    plain = tmp_path / 'p.nc'
    profile.save_netcdf(str(plain))
    with netCDF4.Dataset(plain) as handle:
        assert handle.variables['q'].units == 'g/kg'
        np.testing.assert_allclose(handle.variables['q'][:],
                                   profile.q.to('g/kg').magnitude)


# --- coefficient provenance (l9l.15) ---------------------------------------

def _run(directory, *args):
    subprocess.run(['git', '-C', str(directory), '-c', 'user.name=t',
                    '-c', 'user.email=t@example.com', *args],
                   check=True, capture_output=True)


@pytest.fixture
def coefficient_repo(tmp_path):
    repo = tmp_path / 'SensorCoefficients'
    repo.mkdir()
    _run(repo, 'init', '-q')
    table = repo / 'MasterCoefList.csv'
    table.write_text('a,b\n1,2\n')
    _run(repo, 'add', 'MasterCoefList.csv')
    _run(repo, 'commit', '-q', '-m', 'table')
    head = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'],
                          capture_output=True, text=True).stdout.strip()
    return repo, table, head


def test_revision_is_the_commit_that_owns_the_table(coefficient_repo):
    _, table, head = coefficient_repo
    assert provenance.coefficient_revision(table) == head


def test_revision_follows_a_symlink_to_the_real_checkout(coefficient_repo,
                                                         tmp_path):
    _, table, head = coefficient_repo
    config_dir = tmp_path / 'wxuas'
    config_dir.mkdir()
    link = config_dir / 'MasterCoefList.csv'
    link.symlink_to(table)
    assert provenance.coefficient_revision(link) == head


def test_uncommitted_edit_is_marked_dirty(coefficient_repo):
    _, table, head = coefficient_repo
    table.write_text('a,b\n1,3\n')
    assert provenance.coefficient_revision(table) == head + '-dirty'


def test_enclosing_repository_is_not_mistaken_for_the_owner(coefficient_repo):
    """A table that merely sits inside some repo has no revision there."""
    repo, _, _ = coefficient_repo
    nested = repo / 'untracked'
    nested.mkdir()
    stray = nested / 'MasterCoefList.csv'
    stray.write_text('a,b\n')
    assert provenance.coefficient_revision(stray) == 'unknown'


def test_revision_of_a_missing_table_is_unknown(tmp_path):
    assert provenance.coefficient_revision(tmp_path / 'nope.csv') == 'unknown'
    assert provenance.coefficient_sha256(tmp_path / 'nope.csv') == 'unknown'


def test_sha256_identifies_the_table(coefficient_repo):
    _, table, _ = coefficient_repo
    assert provenance.coefficient_sha256(table) == \
        hashlib.sha256(table.read_bytes()).hexdigest()


def test_provenance_records_revision_and_hash_of_the_real_table(profile):
    attributes = provenance.provenance_attributes(profile)
    table = os.path.join(attributes['coefficient_directory'],
                         'MasterCoefList.csv')
    assert attributes['coefficient_sha256'] == \
        hashlib.sha256(open(table, 'rb').read()).hexdigest()
    revision = attributes['coefficient_revision']
    # test/data/coefs is tracked here: a commit hash, not 'unknown'.
    assert revision.split('-')[0] != 'unknown'
    assert len(revision.split('-')[0]) == 40


# --- side effects of a write (l9l.21) --------------------------------------

def test_default_save_makes_no_network_call(profile, tmp_path, monkeypatch):
    import profiles.utils as utils

    def refuse(*args, **kwargs):
        raise AssertionError('save_netcdf touched the network')

    monkeypatch.setattr(utils.requests, 'get', refuse)
    profile.save_netcdf(str(tmp_path / 'x.nc'))
    assert 'flight_location' not in _attributes(tmp_path / 'x.nc')


def test_place_lookup_is_opt_in_and_fails_soft(profile, tmp_path, monkeypatch):
    import profiles.utils as utils

    def down(*args, **kwargs):
        raise OSError('offline')

    monkeypatch.setattr(utils.requests, 'get', down)
    with pytest.warns(RuntimeWarning, match='place lookup failed'):
        profile.save_netcdf(str(tmp_path / 'x.nc'), lookup_place=True)
    assert os.path.exists(tmp_path / 'x.nc')


def test_place_lookup_sends_user_agent_and_timeout(monkeypatch):
    import profiles.utils as utils
    seen = {}

    class Reply:
        def raise_for_status(self):
            pass

        def json(self):
            return {'display_name': 'A, B, C'}

    def get(url, params=None, **kwargs):
        seen.update(kwargs)
        return Reply()

    monkeypatch.setattr(utils.requests, 'get', get)
    assert utils.get_place_from_lat_lon(35.0, -97.0) == 'A, B'
    assert seen['timeout'] > 0
    assert seen['headers']['User-Agent'].startswith('profiles-uas/')


def test_writers_do_not_need_os_uname(profile, tmp_path, monkeypatch):
    """os.uname does not exist on Windows."""
    monkeypatch.delattr(os, 'uname')
    profile.save_netcdf(str(tmp_path / 'x.nc'))
    profile.save_cfnetcdf('N934UA', 340, str(tmp_path / 'y.nc'))
    assert _attributes(tmp_path / 'x.nc')['processing_machine']


def test_writers_raise_no_deprecation_warnings(profile, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter('error', DeprecationWarning)
        profile.save_netcdf(str(tmp_path / 'x.nc'))
        profile.save_cfnetcdf('N934UA', 340, str(tmp_path / 'y.nc'))


def test_processing_datetime_is_utc_iso(profile):
    from datetime import datetime
    stamp = provenance.provenance_attributes(profile)['processing_datetime']
    assert stamp.endswith('Z')
    parsed = datetime.fromisoformat(stamp[:-1])
    assert abs((datetime.now() - parsed).total_seconds()) < 86400 * 2


# --- small cleanups (l9l.32) ----------------------------------------------

def test_explicit_file_path_is_not_overwritten_by_the_flights(profile,
                                                              tmp_path):
    from profiles.Profile import Profile
    flight = profile._raw_profile
    elsewhere = str(tmp_path / 'chosen.bin')
    rebuilt = Profile(elsewhere, RESOLUTION, RES_UNITS, 1, dev=True,
                      raw_profile=flight, confirm_bounds=False,
                      profile_start_height=PROFILE_START_HEIGHT)
    assert rebuilt.file_path == str(tmp_path / 'chosen')


def test_profile_set_forwards_coefficient_dir(monkeypatch, tmp_path):
    from profiles import Profile_Set, processing
    seen = {}

    def capture(path, config, **kwargs):
        seen['config'] = config
        raise RuntimeError('stop')

    monkeypatch.setattr(processing, 'profiles_from_flight', capture)
    profile_set = Profile_Set.Profile_Set(coefficient_dir=str(tmp_path))
    with pytest.warns(DeprecationWarning):
        with pytest.raises(RuntimeError, match='stop'):
            profile_set.add_all_profiles('whatever.BIN')
    assert seen['config'].coefficient_dir == str(tmp_path)
