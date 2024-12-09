import numpy as np
import os
from datetime import datetime, timedelta

import dronelogbook as dlb
import profiles
from profiles import Meta
from profiles import Profile, Profile_Set, plotting

download_dir = '/Users/tyler.bell/Data/SCALES/'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

start = datetime(2024, 9, 8)
end = datetime(2024, 9, 9)

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
    theplace = f.place_name.replace(' ', '').replace(',', '')
    filename = f"{download_dir}/{theplace}/{f.place_name.replace(' ', '')}_flight{f.raw_data['flight_number']}_{f.flight_time:%Y%m%d_%H%M%S}.BIN"
    # filename = f"{download_dir}/ATDD/{f.place_name.replace(' ', '')}_flight{f.raw_data['flight_number']}_{f.flight_time:%Y%m%d_%H%M%S}.BIN"
    bin_files.append(filename)

    if not os.path.exists(filename):
        print("Downloading " + filename)
        conn.get_binary(f.guid, filename)

# Get all the unique locations as well
flight_locations = list(set([(f.place_name, f.place_guid) for f in flights[ind]]))

did_not_process = []

for place, place_guid in flight_locations:

    # if 'Westheimer' not in place and "KAEFS" not in place:
    #     continue

    # if 'Westheimer'  in place or "KAEFS"  in place:
    #     continue

    if place is '':
        continue

    alt = conn.get_place(place_guid).get_usgs_alt()  # Get the ground elevation from the USGS
    print(f"Processing files from {place} starting at altitude: {round(alt+10)}")
    a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=False,
                                nc_level=None, profile_start_height=alt+10, legacy_peak_id=True)

    # Process the files
    for f, fn in zip(flights[ind], bin_files):
        if f.place_name != place:
            continue

        if not os.path.exists(fn):
            continue
        try:
            flight = conn.get_flight(f.guid, recursive=True)
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
            # p.get_wind_profile(algorithm='quadratic')
            # p.save_netcdf()

            serial_number = flight.drone.raw_data['serial_number']
            wmo_id = serial_number[0:2]+serial_number[-3:]
            theplace = flight.place_name.replace(' ', '').replace(',', '')

            filename = os.path.join(download_dir, theplace, f"UASDC_037_{wmo_id}_{p.gridded_times[0]:%Y%m%d%H%M%SZ}.nc")
            # filename = os.path.join(download_dir, 'ATDD', f"UASDC_015_CS2B_{p.gridded_times[0]:%Y%m%d%H%M%SZ}.nc")
            p.save_cfnetcdf(wmo_id, alt, filename)

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
