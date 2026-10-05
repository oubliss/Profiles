"""
Batch processing API.

process_flights replaces Profile_Set. The first test is the one that
matters: the new path must produce the same numbers as the old one.
"""
import warnings

import numpy as np
import pytest

from profiles import Profile_Set
from profiles.processing import (ProcessingConfig, all_profiles,
                                 process_flights, profiles_from_flight)
from test.harness import (PROFILE_START_HEIGHT, RESOLUTION, RES_UNITS,
                          staged_bin)


@pytest.fixture(scope='module')
def bin_path(tmp_path_factory):
    return str(staged_bin(tmp_path_factory.mktemp('processing')))


@pytest.fixture(scope='module')
def config():
    return ProcessingConfig(resolution=RESOLUTION, res_units=RES_UNITS,
                            ascent=True, dev=True, confirm_bounds=False,
                            nc_level=None,
                            profile_start_height=PROFILE_START_HEIGHT)


def test_matches_profile_set(bin_path, config):
    """The deprecated container and the new function must agree exactly."""
    results = process_flights([bin_path], config)
    new = all_profiles(results)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        old_set = Profile_Set.Profile_Set(
            resolution=RESOLUTION, res_units=RES_UNITS, ascent=True,
            dev=True, confirm_bounds=False, nc_level=None,
            profile_start_height=PROFILE_START_HEIGHT)
        old_set.add_all_profiles(bin_path)
    old = old_set.profiles
    for profile in old:
        profile.compute_thermo()
        profile.compute_wind()

    assert len(new) == len(old)
    for a, b in zip(new, old):
        for name in ('alt', 'pres', 'temp', 'rh', 'theta', 'speed', 'dir'):
            np.testing.assert_array_equal(
                getattr(a, name).magnitude, getattr(b, name).magnitude,
                err_msg=f'{name} differs between process_flights and Profile_Set')


def test_computes_by_default(bin_path, config):
    profile = all_profiles(process_flights([bin_path], config))[0]
    assert profile.temp is not None
    assert profile.speed is not None


def test_compute_can_be_deferred(bin_path, config):
    profile = all_profiles(process_flights([bin_path], config,
                                           compute=False))[0]
    assert not profile._thermo_computed
    assert not profile._wind_computed


def test_errors_are_collected_not_raised(bin_path, config, tmp_path):
    decoy = tmp_path / 'notalog.BIN'
    decoy.write_text('{"status":3,"description":"Flight Id does not exist."}')

    results = process_flights([str(decoy), bin_path], config)

    assert not results[0].ok
    assert isinstance(results[0].error, ValueError)
    assert results[1].ok
    # One bad file must not lose the good ones.
    assert len(all_profiles(results)) > 0


def test_errors_can_be_raised_instead(config, tmp_path):
    decoy = tmp_path / 'notalog.BIN'
    decoy.write_text('not a log')
    with pytest.raises(ValueError):
        process_flights([str(decoy)], config, on_error='raise')


def test_min_levels_discards_phantom_profiles(bin_path, config):
    """flight616's second 'profile' is a 2 s, 1 m wiggle at the top."""
    everything = all_profiles(process_flights([bin_path], config,
                                              compute=False))
    assert len(everything) == 2

    strict = ProcessingConfig(**{**config.__dict__, 'min_levels': 141})
    kept = all_profiles(process_flights([bin_path], strict, compute=False))
    assert len(kept) == 1


def test_profiles_come_back_sorted(bin_path, config):
    profiles = all_profiles(process_flights([bin_path], config,
                                            compute=False))
    starts = [p.gridded_times[0] for p in profiles]
    assert starts == sorted(starts)


def test_default_config_is_usable():
    config = ProcessingConfig()
    assert config.resolution == 10
    assert config.res_units == 'm'
    assert config.profile_start_height is None


def test_default_config_auto_detects_calibration():
    assert ProcessingConfig().calibration == 'auto'


def test_reference_flight_resolves_to_table_calibration(bin_path):
    """It reports serial numbers, so it predates onboard calibration."""
    from profiles.Coef_Manager import TableCalibration
    from profiles.flight import FlightLog

    flight = FlightLog(bin_path, nc_level=None)
    assert isinstance(flight.calibration_source, TableCalibration)
    assert flight.calibration_source.temperature_from == 'resistance'


def test_calibration_mode_reaches_the_flight(bin_path, config):
    """A forced mode must survive the config -> FlightLog -> Profile hop."""
    import dataclasses

    from profiles.Coef_Manager import OnboardCalibration
    from profiles.processing import profiles_from_flight

    forced = dataclasses.replace(config, calibration='onboard',
                                 min_levels=5)
    profiles = profiles_from_flight(bin_path, forced)
    assert profiles, 'expected at least one profile'

    profile = profiles[0]
    assert isinstance(profile._raw_profile.calibration_source,
                      OnboardCalibration)

    profile.compute_thermo()
    assert 'onboard' in profile.calibration_record['temperature_source']


def test_forcing_onboard_changes_the_temperature(bin_path, config):
    """The two paths must actually differ, or the switch means nothing."""
    import dataclasses

    import numpy as np

    from profiles.processing import profiles_from_flight

    def first_temp(mode):
        cfg = dataclasses.replace(config, calibration=mode, min_levels=5)
        profile = profiles_from_flight(bin_path, cfg)[0]
        return profile.compute_thermo().temp.magnitude

    table, onboard = first_temp('table'), first_temp('onboard')
    assert not np.allclose(table, onboard, equal_nan=True), (
        'table and onboard calibration produced identical temperature')
