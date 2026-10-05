"""
The v2 entry point, process_flights, against the same p0 snapshot the
deprecated Profile_Set path is held to (test_baseline.py).

harness.run_reference_pipeline drives Profile_Set; this drives the supported
API with the equivalent configuration and compares profile 0 array for array
with the existing snapshot - nothing is re-captured. Like test_baseline, the
snapshot's retired phantom profile (p1.*) is excluded. The lowpass variant
is not repeated: process_flights has no filtering step.
"""
from collections import OrderedDict

import numpy as np
import pytest

from profiles.processing import ProcessingConfig, all_profiles, process_flights
from test import BASELINE_PATH
from test.harness import (PROFILE_START_HEIGHT, PROFILE_VARS, RESOLUTION,
                          RES_UNITS, THERMO_VARS, WIND_VARS, _collect,
                          staged_bin)
from test.test_baseline import _describe

#: Same position tolerances, and for the same reason, as test_baseline.
POSITION_ATOL = {'lat': 1e-7, 'lon': 1e-7, 'alt_MSL': 0.5}


@pytest.fixture(scope='module')
def v2_arrays(tmp_path_factory):
    path = staged_bin(tmp_path_factory.mktemp('v2flight'))
    config = ProcessingConfig(
        resolution=RESOLUTION, res_units=RES_UNITS, ascent=True, dev=True,
        confirm_bounds=False, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT)
    results = process_flights([str(path)], config)
    assert all(result.ok for result in results), \
        [result.error for result in results]
    profiles = all_profiles(results)
    assert len(profiles) == 1

    out = OrderedDict()
    _collect('p0.profile', profiles[0], PROFILE_VARS, out)
    _collect('p0.thermo', profiles[0], THERMO_VARS, out)
    _collect('p0.wind', profiles[0], WIND_VARS, out)
    return out


def test_process_flights_matches_the_p0_snapshot(v2_arrays):
    expected = np.load(BASELINE_PATH / 'flight616_10m.npz')
    names = [n for n in expected.files if n.startswith('p0.')]
    assert names, 'snapshot has no p0 arrays'

    assert sorted(v2_arrays) == sorted(names)

    differences = []
    for name in names:
        want, got = expected[name], np.asarray(v2_arrays[name])
        atol = POSITION_ATOL.get(name.rsplit('.', 1)[-1], 0)
        if want.shape != got.shape or not np.allclose(
                want, got, rtol=1e-9, atol=atol, equal_nan=True):
            differences.append(_describe(name, want, got))

    assert not differences, (
        f'{len(differences)} of {len(names)} arrays differ from the '
        f'snapshot through process_flights:\n  ' + '\n  '.join(differences))
