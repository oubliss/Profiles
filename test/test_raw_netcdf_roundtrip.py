"""
Raw_Profile a0 NetCDF round trip.

_save_netCDF / _read_netCDF had no coverage, which is how the rotation block
came to overwrite rot_list[6] three times: PN ended up holding PD's values,
and PE and PD were left as raw netCDF Variable objects rather than pint
quantities (so the next arithmetic on them raised instead of converting).
"""
import numpy as np
import pytest

from profiles.Raw_Profile import Raw_Profile
from test.harness import staged_bin

# (index into Raw_Profile.rotation, variable name in the file)
ROTATION_FIELDS = [(0, 'VE'), (1, 'VN'), (2, 'VD'),
                   (3, 'roll'), (4, 'pitch'), (5, 'yaw'),
                   (6, 'PN'), (7, 'PE'), (8, 'PD')]


@pytest.fixture(scope='module')
def round_tripped(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('roundtrip')
    original = Raw_Profile(str(staged_bin(tmp)), dev=True, nc_level=None)
    nc_path = tmp / 'flight616_a0.nc'
    original._save_netCDF(str(nc_path))
    return original, Raw_Profile(str(nc_path), dev=True)


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
