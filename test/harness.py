"""
Reference pipeline used to capture and check baseline snapshots.

The point of this module is to pin down what the package currently produces
for one known flight, so that the restructuring work in later stages can prove
it changed nothing by accident. Any intentional change to these numbers has to
be re-captured deliberately with ``python -m test.capture_baseline`` and
explained in the commit message.
"""
import shutil
from collections import OrderedDict

import numpy as np

from test import BASE_TEST_PATH

BIN_NAME = 'flight616_20210609_080023.BIN'
RESOLUTION = 10
RES_UNITS = 'm'
PROFILE_START_HEIGHT = 350

# Variables snapshotted from each object. Missing attributes are recorded as
# absent rather than skipped, so a variable that disappears is also a diff.
PROFILE_VARS = ('time', 'alt', 'pres', 'lat', 'lon', 'alt_MSL',
                'gridded_times', 'gridded_base')
THERMO_VARS = ('temp', 'rh', 'pres', 'alt', 'theta', 'T_d', 'mixing_ratio',
               'q', 'temp_flags', 'rh_flags', 'time', 'gridded_times',
               'lat', 'lon')
WIND_VARS = ('speed', 'dir', 'u', 'v', 'alt', 'pres', 'time',
             'gridded_times', 'lat', 'lon')


def _to_array(value):
    """ Flatten a pint Quantity, datetime sequence or plain sequence to floats.

    Datetimes become epoch seconds so that snapshots are comparable with
    np.allclose and survive a round trip through .npz.
    """
    if value is None:
        return None

    magnitude = getattr(value, 'magnitude', value)
    arr = np.asarray(magnitude)

    if arr.dtype == object or np.issubdtype(arr.dtype, np.datetime64):
        as_datetime = np.asarray(arr, dtype='datetime64[us]')
        return as_datetime.astype('int64') / 1e6

    return arr.astype(float)


def _collect(prefix, obj, names, out):
    for name in names:
        value = getattr(obj, name, None)
        key = f'{prefix}.{name}'
        if value is None:
            out[key + '.__absent__'] = np.array([1.0])
        else:
            out[key] = _to_array(value)


def run_reference_pipeline(bin_path, lowpass=False):
    """ Process the reference flight and return every output array.

    :param pathlib.Path bin_path: the .BIN to process. Pass a copy in a tmp
       directory - reading a .BIN currently writes a large .json next to it.
    :param bool lowpass: apply the filtering the production scripts apply
    :rtype: OrderedDict[str, np.ndarray]
    """
    from profiles import Profile_Set

    profile_set = Profile_Set.Profile_Set(
        resolution=RESOLUTION, res_units=RES_UNITS, ascent=True, dev=True,
        confirm_bounds=False, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT)
    profile_set.add_all_profiles(str(bin_path))

    out = OrderedDict()
    out['n_profiles'] = np.array([float(len(profile_set.profiles))])

    for i, profile in enumerate(profile_set.profiles):
        if lowpass:
            profile.lowpass_filter(wind=True, thermo=False, Fc=.06)

        profile.compute_thermo()
        profile.compute_wind()

        # Keys keep their pre-merge names so snapshots stay comparable;
        # all three groups now read off the one Profile.
        _collect(f'p{i}.profile', profile, PROFILE_VARS, out)
        _collect(f'p{i}.thermo', profile, THERMO_VARS, out)
        _collect(f'p{i}.wind', profile, WIND_VARS, out)

    return out


def staged_bin(tmp_dir):
    """ Copy the reference .BIN into tmp_dir and return the copy's path.

    Reading a .BIN writes an ~83 MB .json beside it; keep that out of the repo.
    """
    source = BASE_TEST_PATH / BIN_NAME
    destination = tmp_dir / BIN_NAME
    shutil.copy(source, destination)
    return destination
