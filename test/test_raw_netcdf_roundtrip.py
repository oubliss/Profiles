"""
FlightLog a0 NetCDF round trip.

_save_netCDF / _read_netCDF had no coverage, which is how the rotation block
came to overwrite rot_list[6] three times: PN ended up holding PD's values,
and PE and PD were left as raw netCDF Variable objects rather than pint
quantities (so the next arithmetic on them raised instead of converting).
"""
import numpy as np
import pytest

from profiles.flight import FlightLog
from test.harness import staged_bin

# (index into FlightLog.rotation, variable name in the file)
ROTATION_FIELDS = [(0, 'VE'), (1, 'VN'), (2, 'VD'),
                   (3, 'roll'), (4, 'pitch'), (5, 'yaw'),
                   (6, 'PN'), (7, 'PE'), (8, 'PD')]


@pytest.fixture(scope='module')
def round_tripped(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('roundtrip')
    original = FlightLog(str(staged_bin(tmp)), dev=True, nc_level=None)
    nc_path = tmp / 'flight616_a0.nc'
    original._save_netCDF(str(nc_path))
    return original, FlightLog(str(nc_path), dev=True)


@pytest.mark.parametrize('index,name', ROTATION_FIELDS)
def test_rotation_survives_round_trip(round_tripped, index, name):
    original, reloaded = round_tripped

    assert hasattr(reloaded.rotation[index], 'magnitude'), (
        f'{name} came back as {type(reloaded.rotation[index]).__name__} '
        f'rather than a pint Quantity - it was never unit-converted on read')

    np.testing.assert_array_equal(
        np.asarray(original.rotation[index].magnitude),
        np.asarray(reloaded.rotation[index].magnitude),
        err_msg=f'{name} (rotation[{index}]) changed across the round trip')


def test_position_components_are_distinct(round_tripped):
    """wind_data() read index 6 for all three position components."""
    original, _ = round_tripped
    wind = original.wind_data()

    assert not np.array_equal(wind['pos_n'].magnitude, wind['pos_e'].magnitude)
    assert not np.array_equal(wind['pos_e'].magnitude, wind['pos_d'].magnitude)

    for key, index in (('pos_n', 6), ('pos_e', 7), ('pos_d', 8)):
        np.testing.assert_array_equal(
            wind[key].magnitude, original.rotation[index].magnitude,
            err_msg=f'{key} does not map to rotation[{index}]')


def test_serial_numbers_survive_round_trip(round_tripped):
    original, reloaded = round_tripped
    for key in ('copterID', 'imet1', 'imet2', 'imet3', 'rh1', 'rh2', 'rh3'):
        assert float(reloaded.serial_numbers[key]) == float(
            original.serial_numbers[key]), f'{key} changed across the round trip'


# ---------------------------------------------------------------------------
# Every group, with units
# ---------------------------------------------------------------------------

def _same(first, second, label):
    """Equal values and, for quantities, equal units."""
    if hasattr(first, 'magnitude'):
        assert hasattr(second, 'magnitude'), f'{label} lost its units'
        assert first.units == second.units, (
            f'{label}: {first.units} became {second.units}')
        first, second = first.magnitude, second.magnitude
    np.testing.assert_array_equal(np.asarray(first), np.asarray(second),
                                  err_msg=f'{label} changed')


def _times_equal(first, second, label):
    assert len(first) == len(second), f'{label}: length changed'
    assert all(a == b for a, b in zip(first, second)), f'{label} times changed'


@pytest.mark.parametrize('group', ['temp', 'rh', 'pos', 'pres', 'rotation'])
def test_tuple_groups_round_trip(round_tripped, group):
    """temp (incl. resi), rh (incl. temp_rh), pos, pres (incl. both temps)."""
    original, reloaded = round_tripped
    before, after = getattr(original, group), getattr(reloaded, group)

    assert len(before) == len(after), f'{group} slot count changed'
    for index, (a, b) in enumerate(zip(before[:-1], after[:-1])):
        _same(a, b, f'{group}[{index}]')
    _times_equal(before[-1], after[-1], f'{group} time')


def test_temperatures_are_kelvin_with_resistances(round_tripped):
    """The a0 round trip used to hand back volts/millivolts, ~100 K high."""
    _, reloaded = round_tripped
    data = reloaded.thermo_data()
    for number in (1, 2, 3):
        assert str(data[f'temp{number}'].units) == 'kelvin'
        assert str(data[f'resi{number}'].units) == 'ohm'
        assert 250 < np.nanmean(data[f'temp{number}'].magnitude) < 320
        assert np.nanmean(data[f'resi{number}'].magnitude) > 0
    # RH-sensor and barometer temperatures were labelled "F" (farad).
    assert str(reloaded.rh[1].units) == 'kelvin'
    # BARO.Temp/GndTemp are degC in ArduPilot's log-message documentation.
    assert str(reloaded.pres[1].units) == 'degree_Celsius'
    assert str(reloaded.pres[2].units) == 'degree_Celsius'
    assert 0 < np.nanmean(reloaded.pres[1].magnitude) < 60


def test_wind_events_messages_round_trip(round_tripped):
    original, reloaded = round_tripped

    for index, (a, b) in enumerate(zip(original.wind[:-1], reloaded.wind[:-1])):
        _same(a, b, f'wind[{index}]')
    _times_equal(original.wind[-1], reloaded.wind[-1], 'wind time')

    _same(original.events[0], reloaded.events[0], 'events')
    _times_equal(original.events[1], reloaded.events[1], 'events time')

    assert list(reloaded.messages[0]) == list(original.messages[0])
    _times_equal(original.messages[1], reloaded.messages[1], 'messages time')


def test_every_serial_number_round_trips(round_tripped):
    original, reloaded = round_tripped
    assert 'wind' in original.serial_numbers
    for key, value in original.serial_numbers.items():
        assert reloaded.serial_numbers[key] == value, key


def test_gridded_output_from_a0_matches_the_bin(tmp_path_factory):
    """The point of an a0 file: reprocessing it must give the same profile."""
    from profiles.processing import ProcessingConfig, profiles_from_flight

    tmp = tmp_path_factory.mktemp('reprocess')
    bin_path = str(staged_bin(tmp))
    original = FlightLog(bin_path, dev=True, nc_level=None)
    a0_path = str(tmp / 'a0.nc')
    original._save_netCDF(a0_path)

    config = ProcessingConfig(confirm_bounds=False, profile_start_height=350,
                              nc_level=None, dev=True)
    from_bin = profiles_from_flight(bin_path, config, flight=original)
    from_a0 = profiles_from_flight(a0_path, config,
                                   flight=FlightLog(a0_path, dev=True))

    assert len(from_bin) == len(from_a0) > 0
    for expected, actual in zip(from_bin, from_a0):
        expected.compute_thermo()
        actual.compute_thermo()
        # The derived thermodynamic values, not just the measured ones: they
        # are what a user reprocessing from an a0 file actually wants.
        for name in ('temp', 'rh', 'pres', 'alt', 'theta', 'T_d',
                     'mixing_ratio', 'q'):
            want = getattr(expected, name)
            got = getattr(actual, name)
            assert want is not None and got is not None, name
            assert getattr(want, 'units', None) == getattr(got, 'units', None), \
                f'{name} units differ'
            assert np.isfinite(want.magnitude).any(), f'{name} is all NaN'
            assert np.allclose(want.magnitude, got.magnitude,
                               equal_nan=True), f'gridded {name} differs'
        assert np.array_equal(expected.temp_flags, actual.temp_flags)
        assert np.array_equal(expected.rh_flags, actual.rh_flags)


# ---------------------------------------------------------------------------
# Writer / reader symmetry
# ---------------------------------------------------------------------------

def test_log_without_copter_id_wind_or_events(tmp_path):
    """Pre-2021 logs: no SYSID_THISMAV, no WIND, no EV messages."""
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    del flight.serial_numbers['copterID']
    flight.wind = None
    flight.events = None

    nc_path = str(tmp_path / 'sparse_a0.nc')
    flight._save_netCDF(nc_path)
    reloaded = FlightLog(nc_path, dev=True)

    assert 'copterID' not in reloaded.serial_numbers
    assert reloaded.wind is None
    assert reloaded.events is None
    np.testing.assert_array_equal(flight.temp[0].magnitude,
                                  reloaded.temp[0].magnitude)


def test_message_text_is_not_truncated(tmp_path):
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    long_text = 'EKF3 IMU0 is using GPS and then some more text ' * 8
    flight.messages = ([long_text] + list(flight.messages[0][1:]),
                       flight.messages[1])

    nc_path = str(tmp_path / 'messages_a0.nc')
    flight._save_netCDF(nc_path)

    assert FlightLog(nc_path, dev=True).messages[0][0] == long_text


def test_old_volt_files_read_as_kelvin(tmp_path):
    """a0 files from earlier 1.4.0-dev builds: 'volt<n>' holding kelvin."""
    import netCDF4
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    nc_path = str(tmp_path / 'old_a0.nc')
    # The writer skips slots that are not quantities, which leaves no resi<n>
    # variables - as in the old files.
    full = flight.temp
    flight.temp = tuple(
        value.magnitude if index % 2 else value
        for index, value in enumerate(full[:-2])) + full[-2:]
    flight._save_netCDF(nc_path)
    flight.temp = full

    with netCDF4.Dataset(nc_path, 'a') as dataset:
        group = dataset['temp']
        for number in (1, 2, 3, 4):
            group.renameVariable(f'temp{number}', f'volt{number}')
            group.variables[f'volt{number}'].units = 'mV'
        dataset['rh'].variables['temp1'].units = 'F'

    reloaded = FlightLog(nc_path, dev=True)
    assert str(reloaded.temp[0].units) == 'kelvin'
    np.testing.assert_array_equal(reloaded.temp[0].magnitude,
                                  flight.temp[0].magnitude)
    assert np.isnan(reloaded.temp[1].magnitude).all()
    assert str(reloaded.rh[1].units) == 'kelvin'


# ---------------------------------------------------------------------------
# a0 never overwrites its own input
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', ['FLIGHT.JSON', 'x.Bin', 'x.json', 'Y.BIN'])
def test_fallback_name_is_never_the_input(tmp_path, name):
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    source = tmp_path / name
    source.write_text('original log')
    flight.file_path = str(source)

    flight._save_netCDF(str(source))        # not a .nc path: fallback is used

    assert source.read_text() == 'original log'
    assert (tmp_path / (source.stem + '.nc')).exists()


def test_refuses_to_overwrite_input(tmp_path):
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    source = tmp_path / 'input.nc'
    source.write_text('original log')
    flight.file_path = str(source)

    with pytest.raises(ValueError, match='refusing'):
        flight._save_netCDF(str(source))
    assert source.read_text() == 'original log'


# ---------------------------------------------------------------------------
# Profile file types
# ---------------------------------------------------------------------------

def _profile_from(tmp_path, name):
    from profiles.processing import ProcessingConfig, profiles_from_flight
    flight = FlightLog(str(staged_bin(tmp_path)), dev=True, nc_level=None)
    # Profile takes its path from the flight it is handed.
    flight.file_path = str(tmp_path / name)
    config = ProcessingConfig(confirm_bounds=False, profile_start_height=350,
                              nc_level=None, dev=True)
    return profiles_from_flight(flight.file_path, config, flight=flight)


def test_profile_accepts_cdf(tmp_path):
    profiles = _profile_from(tmp_path, 'flight.cdf')
    assert profiles and profiles[0].file_path == str(tmp_path / 'flight')


def test_profile_rejects_unknown_extension_with_value_error(tmp_path):
    with pytest.raises(ValueError, match='unrecognised extension'):
        _profile_from(tmp_path, 'flight.txt')
