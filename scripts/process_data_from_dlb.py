import numpy as np
import os
from datetime import datetime, timedelta

import dronelogbook as dlb
import profiles
from profiles import Meta
from profiles import Profile, Profile_Set, plotting

# download_dir = '/Users/tyler.bell/Data/OUTFLOW/IOP1'
# download_dir = '/Users/tyler.bell/Data/tracer_sgp/'
# download_dir = '/Users/tyler.bell/ATDD/test/'
# download_dir = '/Users/tyler.bell/data/M2HATS/'
# download_dir = '/Users/tyler.bell/Data/uas_wspd/uas/'
download_dir = '/Users/tyler.bell/Data/OUTFLOW/IOP3'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

# start = datetime(2023, 3, 2)
# end = datetime(2023, 3, 4)
# start = datetime(2021, 4, 8)
# end = datetime(2021, 4, 9)
# start = datetime(2022, 3, 30)
# end = datetime(2022, 4, 1)
# start = datetime(2021, 6, 4)
# end = datetime(2021, 7, 8)

start = datetime(2024, 5, 1)
end = datetime(2024, 5, 2)

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
    filename = f"{download_dir}/{f.place_name.replace(' ', '')}_flight{f.raw_data['flight_number']}_{f.flight_time:%Y%m%d_%H%M%S}.BIN"
    bin_files.append(filename)

    if not os.path.exists(filename):
        print("Downloading " + filename)
        conn.get_binary(f.guid, filename)

# Get all the unique locations as well
flight_locations = list(set([(f.place_name, f.place_guid) for f in flights[ind]]))

did_not_process = []

for place, place_guid in flight_locations:
    # if place != "KAEFS":
    #     continue

    alt = conn.get_place(place_guid).get_usgs_alt()  # Get the ground elevation from the USGS
    print(f"Processing files from {place} starting at altitude: {round(alt+10)}")
    a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=False,
                                nc_level=None, profile_start_height=alt+10)

    # Process the files
    for f, fn in zip(flights[ind], bin_files):
        if f.place_name != place:
            continue

        if not os.path.exists(fn):
            continue
        try:

            metadata = Meta.Meta(guid=f.guid)

            # If the json is already there, use that
            if os.path.exists(fn.replace('.BIN', '.json')):
                fn = fn.replace('.BIN', '.json')

            print(fn)

            a.add_all_profiles(fn, metadata=metadata)

            # # Go ahead and process the profile
            p = a.profiles[-1]
            # p.lowpass_filter()
            p.lowpass_filter(wind=True, thermo=False, Fc=.06)
            # p.lowpass_filter(wind=False, thermo=True)
            p.get_thermo_profile()
            p.get_wind_profile()
            p.save_netcdf()
            # p.save_cfnetcdf()

        except Exception as e:
            import traceback
            did_not_process.append((fn, traceback.print_exc()))
            print("Error on " + fn)
            print(traceback.print_exc())




    # Newer way
    # for p in a.profiles:
    #     if len(p.gridded_times) > 3:
    #         p.lowpass_filter(wind=True, thermo=False, Fc=.06)
    #         p.get_thermo_profile()
    #         p.get_wind_profile()
    #
    #         p.save_netcdf()

    # a.save_netCDF(download_dir)

    # older way
    # at = []
    # for p in a.profiles:
    #     if len(p.gridded_times) > 3:
    #         at.append(p.get_thermo_profile())
    #
    # aw = []
    # for p in a.profiles:
    #     try:
    #         aw.append(p.get_wind_profile())
    #     except Exception:
    #         import traceback
    #         print("Error on " + fn)
    #         print(traceback.print_exc())
    #
    # a.save_netCDF(download_dir)
