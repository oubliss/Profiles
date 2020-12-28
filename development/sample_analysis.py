from profiles.Profile_Set import Profile_Set
from profiles.Profile import Profile
import profiles.plotting as plotting
import matplotlib.pyplot as plt
import numpy as np
from metpy.plots import Hodograph

"""
Read data. The raw data will be read in any case because ".json" was specified,
but pre-processed data can be read by replacing ".json" with ".nc".

You will be asked for a profile start height. This is the altidute above which
the data should be considered valid and part of a distinct profile. If multiple
profiles are flown in one flight, this profile start height MUST be above the
altitude the craft flies down to at the end of each profile.
   No units are needed, simply enter the appropriate y-value from the chart
of altitude vs. time shown, then hit "enter". You will be asked to confirm the
times that the profile selection algorithm has chosen for the start, peak, and
end of the profile. Simply hit "enter" if the times were correctly identified.
If they were not identified correctly, you will be asked to try again using a
different profile start height.
   See the docs for details about parameters.
"""

# Example using Profile_Set
a = Profile_Set(resolution=1, res_units='m', ascent=True, dev=True, profile_start_height=400, confirm_bounds=False, nc_level='low')
a.add_all_profiles("/home/jessicablunt/data_templates/CS_N935UA_KAEFS_flights/20200710_1814.json")
a.add_all_profiles("/home/jessicablunt/data_templates/CS_N935UA_KAEFS_flights/20200710_1832.json")



# Example using Profile
# b = Profile("/home/jessica/GitHub/data_templates/BIN/00000010.BIN", 10, 'm', 1,
#             dev=True, ascent=True, scoop_id='A')

"""
Calculate thermodynamic variables from raw data. The values returned have
already gone through a quality control process based on variance and bias. If
a Thermo_Profile object was already created from the same ".json" file AND the
same vertical resolution AND the same value of ascent (True/False) was used,
the pre-processed netCDF is read to save time and avoid redundant calculations.
"""

# Example using Profile_Set
at = []
for p in a.profiles:
    at.append(p.get_thermo_profile())


aw = []
for p in a.profiles:
    aw.append(p.get_wind_profile())
"""
# bt = b.get_thermo_profile()
# bw = b.get_wind_profile()



Here's an example of one way to view the processed data. See the docs for a
full list of accessible variables.
"""
# Example using Profiles
# plotting.contour_height_time(a.profiles)  # not yet functional
a1 = a.profiles[0]
plotting.plot_skewT(temp=a1.get("temp"), pres=a1.get("pres"), t_d=a1.get("T_d"),
                    u=a1.get("u"), v=a1.get("v"), time=a1.get("gridded_times"),
                    units=a1.get("_units"))
plt.show()
a2 = a.profiles[1]
plotting.plot_skewT(temp=a2.get("temp"), pres=a2.get("pres"), t_d=a2.get("T_d"),
                    u=a2.get("u"), v=a2.get("v"), time=a2.get("gridded_times"),
                    units=a2.get("_units"))
plt.show()


for w in aw:
     # Create a hodograph
     plt.figure()
     h = Hodograph(component_range=10)
     h.plot(w.u, w.v)
     #ws = w.speed
     #wt = w.gridded_times
     #plt.plot(w.speed, w.gridded_times)
     plt.show()


# Example using Profile
# plt.figure()
# plt.plot(bt.temp, bt.alt)
# plt.show()

# Example using Profile
# Create a hodograph
# fig1 = plt.figure()
# ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, theta_direction=-1,
#                     theta_offset=-np.pi/2)
# ax1.set_ylim(0, 20)
# ax1.set_rticks(np.arange(0, 20, 2.5))
# ax1.set_rlabel_position(270)
# ax1.tick_params('y', labelrotation=-45, labelsize='x-small')
# ax1.yaxis.set_label_coords(-0.15, 0.5)
# ax1.plot(bw.dir * np.pi/180, bw.speed, lw=2.5)
# plt.show()

"""
Now look in your data directory. There are .nc files that can be processed
faster than .json for the same result IF you want the same resolution. If you
used binary files as input, there are also new .json files that can be used
to analyze the data at different resolutions faster. Try replacing .json or
.bin with .nc in lines 26 and 31 above.
"""
