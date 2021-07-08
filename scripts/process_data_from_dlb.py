import numpy as np
import os
from datetime import datetime

import dronelogbook as dlb
import profiles
from profiles import Meta
from profiles import Profile, Profile_Set

download_dir = '/Users/tyler.bell/Data/uas'
start = datetime(2020, 8, 27, 0, 0)
end = datetime(2020, 8, 27, 12)

# Connect to DLB
conn = dlb.DLB()

# Get all the flights and the times
flights = np.asarray(sorted(conn.get_flights()))
flight_times = np.array([f.flight_time for f in flights])

# Find the flights that fall within the datetime range we want
ind = np.where((flight_times >= start) & (flight_times <= end))

# Download the data locally
bin_files = []
for f in flights[ind]:
    filename = f"{download_dir}/flight{f.raw_data['flight_number']}_{f.flight_time:%Y%m%d_%H%M%S}.BIN"
    bin_files.append(filename)

    if not os.path.exists(filename):
        print("Downloading " + filename)
        conn.get_binary(f.guid, filename)


a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=False,
                            nc_level="low", profile_start_height=349)

# Process the files
for f, fn in zip(flights[ind], bin_files):
    metadata = Meta.Meta(guid=f.guid)
    a.add_all_profiles(fn, metadata=metadata)

at = []
for p in a.profiles:
    at.append(p.get_thermo_profile())

aw = []
for p in a.profiles:
    aw.append(p.get_wind_profile())


