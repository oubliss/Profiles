"""
Process OK3DM flights with the v2 API.

Kept alongside the original process_ok3dm.py rather than replacing it, so
the two can be compared on the same data before the old one is retired.
"""
import os
from datetime import datetime
from glob import glob

import numpy as np

from profiles.processing import ProcessingConfig, all_profiles, process_flights

DATA_DIR = os.path.expanduser('~/Data/OK3DM')
SITE_ALTITUDE = 397          # KAEFS, m MSL
AIRCRAFT = 'N944UA'

config = ProcessingConfig(
    resolution=10,
    res_units='m',
    ascent=True,
    dev=True,
    profile_start_height=SITE_ALTITUDE,
    nc_level=None,
    # Every OK3DM flight logs SYSID_THISMAV = 1, which maps to four tail
    # numbers with different wind coefficients. Say which one.
    tail_number=AIRCRAFT,
    wind_algorithm='linear',
    # Short wiggles are already dropped at leg detection (min_leg_extent,
    # 50 m by default); this is a backstop against sparse profiles.
    min_levels=20,
)

bin_files = sorted(glob(os.path.join(DATA_DIR, '*.BIN')))
print(f'{len(bin_files)} file(s) in {DATA_DIR}')

results = process_flights(bin_files, config)

for result in results:
    name = os.path.basename(result.path)
    if not result.ok:
        print(f'  SKIP {name}: {type(result.error).__name__}: {result.error}')
        continue
    print(f'  OK   {name}: {len(result.profiles)} profile(s)')

for profile in all_profiles(results):
    stamp = profile.time[0].strftime('%Y%m%d_%H%M%SZ')
    out = os.path.join(DATA_DIR, f'coptersonde_{stamp}.nc')
    profile.save_cfnetcdf(AIRCRAFT, SITE_ALTITUDE, out)
    print(f'wrote {out}')
