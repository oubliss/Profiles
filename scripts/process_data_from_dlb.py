import numpy as np
import os
from datetime import datetime

import dronelogbook as dlb
import profiles
from profiles import Meta
from profiles import Profile, Profile_Set, plotting
from my_lib.utils import elevation

download_dir = '/Users/tyler.bell/Data/uas/perils/IOP4'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

start = datetime(2022, 4, 13)
end = datetime(2022, 4, 13, 23)

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
    print(f)
    filename = f"{download_dir}/flight{f.raw_data['flight_number']}_{f.flight_time:%Y%m%d_%H%M%S}.BIN"
    bin_files.append(filename)

    if not os.path.exists(filename):
        print("Downloading " + filename)
        conn.get_binary(f.guid, filename)

# Get all the unique locations as well
flight_locations = list(set([(f.place_name, f.place_guid) for f in flights[ind]]))

for place, place_guid in flight_locations:
    # if place != "Lake Village":
    #     continue

    alt = conn.get_place(place_guid).get_usgs_alt()  # Get the ground elevation from the USGS
    print(f"Processing files from {place} starting at altitude: {round(alt+10)}")
    a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=False,
                                nc_level="low", profile_start_height=round(alt+10))

    # Process the files
    for f, fn in zip(flights[ind], bin_files):
        if f.place_name != place:
            continue
        # print(fn)

        if not os.path.exists(fn):
            continue
        metadata = Meta.Meta(guid=f.guid)

        try:
            a.add_all_profiles(fn, metadata=metadata)
        except Exception as e:
            import traceback
            print("Error on " + fn)
            print(traceback.print_exc())

    at = []
    for p in a.profiles:
        if len(p.gridded_times) > 3:
            at.append(p.get_thermo_profile())

    aw = []
    for p in a.profiles:
        try:
            aw.append(p.get_wind_profile())
        except Exception:
            import traceback
            print("Error on " + fn)
            print(traceback.print_exc())


