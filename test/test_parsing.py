"""
Schema-driven parsing.

Each message group becomes an xarray Dataset with named variables on a
shared time coordinate. These tests cover the new interface; that it still
reproduces the old positional tuples exactly is covered by the baseline
snapshots in test_baseline.py.
"""
import numpy as np
import pytest
import xarray as xr

from profiles import schema
from profiles.flight import FlightLog
from profiles.parsing import parse
from profiles.readers import iter_bin
from test.harness import staged_bin


@pytest.fixture(scope='module')
def raw(tmp_path_factory):
    path = staged_bin(tmp_path_factory.mktemp('parse'))
    return FlightLog(str(path), dev=True, nc_level=None)


def test_every_required_group_is_present(raw):
    for name in schema.REQUIRED_GROUPS:
        assert name in raw.data, f'{name} missing'
        assert isinstance(raw.data[name], xr.Dataset)


def test_variables_are_named_not_positional(raw):
    assert set(raw.data['pos'].data_vars) == {
        'lat', 'lon', 'alt_MSL', 'alt_rel_home', 'alt_rel_orig'}
    assert set(raw.data['rotation'].data_vars) == {
        'speed_east', 'speed_north', 'speed_down', 'roll', 'pitch', 'yaw',
        'pos_n', 'pos_e', 'pos_d'}


def test_position_components_are_distinct(raw):
    """rotation[6] used to be read for pos_n, pos_e and pos_d alike."""
    rotation = raw.data['rotation']
    assert not np.array_equal(rotation['pos_n'].values, rotation['pos_e'].values)
    assert not np.array_equal(rotation['pos_e'].values, rotation['pos_d'].values)


def test_every_variable_carries_units_or_is_a_ratio(raw):
    unitless = {'R13', 'R23', 'R33', 'fan_flag',
                'gyr_x', 'gyr_y', 'gyr_z', 'acc_x', 'acc_y', 'acc_z'}
    for name, dataset in raw.data.items():
        for variable in dataset.data_vars:
            has_units = bool(dataset[variable].attrs.get('units'))
            assert has_units or variable in unitless, (
                f'{name}.{variable} has no units and is not a known ratio')


def test_each_group_shares_one_time_coordinate(raw):
    for name, dataset in raw.data.items():
        dimension = f'{name}_time'
        assert dimension in dataset.coords
        length = dataset.sizes[dimension]
        for variable in dataset.data_vars:
            assert dataset[variable].sizes[dimension] == length


def test_bar2_supersedes_baro(raw):
    """flight616 logs both; only the external barometer should be kept."""
    assert raw.data['pres'].attrs['source_message_type'] == 'BAR2'
    assert raw.baro == 'BAR2'


def test_preferred_selection_discards_the_superseded_record():
    """A BARO run before the first BAR2 must not survive into the series."""
    def message(message_type, stamp, press):
        return {'meta': {'type': message_type, 'timestamp': stamp},
                'data': {'Press': press, 'Temp': 0.0, 'GndTemp': 0.0,
                         'Alt': 0.0}}

    stream = [message('BARO', 1.6e9 + i, 100.0 + i) for i in range(5)]
    stream += [message('BAR2', 1.6e9 + 10 + i, 200.0 + i) for i in range(3)]
    stream += [message('BARO', 1.6e9 + 20 + i, 300.0 + i) for i in range(4)]

    pres = parse(stream)['groups']['pres']
    assert pres.attrs['source_message_type'] == 'BAR2'
    np.testing.assert_array_equal(pres['pres'].values, [200.0, 201.0, 202.0])


def test_missing_field_becomes_nan_not_a_short_series():
    stream = [{'meta': {'type': 'POS', 'timestamp': 1.6e9 + i},
               'data': {'Lat': 35.0, 'Lng': -97.0, 'Alt': 300.0 + i}}
              for i in range(4)]

    pos = parse(stream)['groups']['pos']
    assert pos.sizes['pos_time'] == 4
    assert np.all(np.isnan(pos['alt_rel_home'].values))
    np.testing.assert_array_equal(pos['lat'].values, [35.0] * 4)


def test_serial_numbers_are_read_from_parameters():
    stream = [{'meta': {'type': 'PARM', 'timestamp': 1.6e9},
               'data': {'Name': 'SYSID_THISMAV', 'Value': 6.0}},
              {'meta': {'type': 'PARM', 'timestamp': 1.6e9},
               'data': {'Name': 'USER_SENSORS1', 'Value': 62275.0}},
              {'meta': {'type': 'PARM', 'timestamp': 1.6e9},
               'data': {'Name': 'USER_SENSORS5', 'Value': 16.0}}]

    serials = parse(stream)['serial_numbers']
    assert serials['copterID'] == 6.0
    assert serials['imet1'] == 62275
    assert serials['rh1'] == 16


def test_schema_covers_every_type_the_reader_requests():
    from profiles.readers.mavlink import WANTED_TYPES
    assert sorted(WANTED_TYPES) == schema.all_message_types(), (
        'the reader and the schema disagree about which types matter')
