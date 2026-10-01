"""
nc_level must mean the same thing in all three classes.

Raw_Profile tested `nc_level == 'low'` while Thermo_Profile and Wind_Profile
tested `nc_level is not None`. The documented "write nothing" value, the
string 'none', is truthy, so it suppressed the raw file while still writing
thermo_ and wind_ ones. scripts/process_data_from_dlb.py passes exactly that.
"""
import warnings

import pytest

from profiles.utils import writes_netcdf


@pytest.mark.parametrize('value', [None, 'none', 'None', 'NONE', ' none '])
def test_none_like_values_write_nothing(value):
    assert writes_netcdf(value) is False


@pytest.mark.parametrize('value', ['low', 'LOW', ' Low '])
def test_low_writes(value):
    assert writes_netcdf(value) is True


def test_unrecognised_value_warns_and_writes_nothing():
    with pytest.warns(UserWarning, match='unrecognised nc_level'):
        assert writes_netcdf('medium') is False


def test_no_warning_for_documented_values():
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert writes_netcdf('low') is True
        assert writes_netcdf('none') is False
        assert writes_netcdf(None) is False
