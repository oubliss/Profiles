import os
from glob import glob

from profiles import Profile, Profile_Set

file_dir = '/Users/tyler.bell/Data/OUTFLOW/IOP3/'
bin_files = glob(os.path.join(file_dir, '*.BIN'))
alt = 340
print(bin_files)
a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=True,
                            nc_level=None, profile_start_height=alt)

# Process the files
for fn in bin_files:

    # If the json is already there, use that
    if os.path.exists(fn.replace('.BIN', '.json')):
        fn = fn.replace('.BIN', '.json')
    try:
        a.add_all_profiles(fn)

        # Go ahead and process the profile

    except Exception as e:
        import traceback

        print("Error on " + fn)
        print(traceback.print_exc())

for p in a.profiles:
    if len(p.gridded_times) > 3:
        p.lowpass_filter(wind=True, thermo=False, Fc=.06)
        p.get_thermo_profile()
        p.get_wind_profile()

        p.save_netcdf(file_dir+f"coptersonde_{p.gridded_times[0].strftime('%Y%m%d_%H%M%S')}.cdf")


