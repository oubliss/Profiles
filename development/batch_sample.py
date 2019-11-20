#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 12:58:46 2019

@author: jessica
"""
from profiles.Profile_Set import Profile_Set
import matplotlib.pyplot as plt
from metpy.plots import Hodograph
import os
from profiles import plotting

datadir = '/home/jessica/GitHub/data_templates/20191014'
# Example using Profiles
a = Profile_Set(resolution=15, res_units='Pa', ascent=True, dev=True,
                confirm_bounds=False, profile_start_height=365,
                nc_level="none")
for file_name in os.listdir(datadir):
    if ".json" in file_name:
        a.add_all_profiles(os.path.join(datadir, file_name),
                           meta_path_flight="/home/jessica/GitHub/Utilities/Logs/20191115_1212N963UA (CS3D)_flight_log.csv",
                           meta_path_header="/home/jessica/GitHub/Utilities/Logs/20191115N963UA (CS3D)_log_header.csv")

a.save_netCDF("pancake.nc")
'''
# Example using Profile (singular)
b = Profile("/home/jessica/GitHub/data_templates/00000003.JSON", 10, 'm', 1,
            dev=True, ascent=True)
'''


"""
Calculate thermodynamic variables from raw data. The values returned have
already gone through a quality control process based on variance and bias. If
a Thermo_Profile object was already created from the same ".json" file AND the
same vertical resolution AND the same value of ascent (True/False) was used,
the pre-processed netCDF is read to save time and avoid redundant calculations.
"""

# Example using Profiles
at = []
for p in a.profiles:
    at.append(p.get_thermo_profile())

aw = []
for p in a.profiles:
    aw.append(p.get_wind_profile())

fig = plotting.contour_height_time(a.profiles,
                                     var=['theta', 'q'],
                                     use_pres=True)
fig.show()
