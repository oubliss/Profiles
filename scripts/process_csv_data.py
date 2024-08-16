import os
from glob import glob

from profiles import Profile, Profile_Set

file_dir = '/Users/tyler.bell/Data/OUTFLOW/IOP3/'

operator_id = '037'
aircraft_id = 'BL020'
site_alt = 340

csv_files = glob(os.path.join(file_dir, '*.csv'))
print(csv_files)
a = Profile_Set.Profile_Set(resolution=5, res_units='m', ascent=True, dev=True, confirm_bounds=False,
                            nc_level=None, profile_start_height=10, legacy_peak_id=True, tail_number='N940UA')

# Process the files
for fn in csv_files:

    try:
        a.add_all_profiles(fn)

    except Exception as e:
        import traceback

        print("Error on " + fn)
        print(traceback.print_exc())

for p in a.profiles:
    if len(p.gridded_times) > 3:
        p.lowpass_filter(wind=True, thermo=False, Fc=.06)
        p.get_thermo_profile()
        p.get_wind_profile()

        file_name = file_dir+f"WMODC_{operator_id}_{aircraft_id}_{p.gridded_times[0].strftime('%Y%m%d%H%M%SZ')}.nc"

        p.save_cfnetcdf(aircraft_id, site_alt, file_name)


