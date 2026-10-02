"""
Manages data from a collection of flights or profiles at a specific location
"""
import os
import netCDF4
import datetime as dt
import numpy as np
from profiles.unit_registry import units
from profiles.Profile import Profile
from profiles.flight import FlightLog
import profiles.utils as utils
import profiles.processing as processing
import warnings
from copy import deepcopy


class Profile_Set():
    """ This class manages data (in the form of Profile objects) from one or
    many flights.

    :var list<Profile> profiles: list of Profile objects at this location
    :var bool ascent: True if data from the ascending leg of the profile is \
       to be used. If False, the descending leg will be processed instead
    :var bool dev: True if data from developmental flights is to be uploaded
    :var int resolution: the vertical resolution desired
    :var str res_units: the units in which the vertical resolution is given
    :var bool confirm_bounds: if True, the user will be asked to verify the \
       automatically-determined start, peak, and end times of each profile
    :var int profile_start_height: either passed to the constructor or \
       provided by the user during processing
    :var Meta meta: reads and processes metadata from oucass-checklist
    """

    def __init__(self, resolution=10, res_units='m', ascent=True,
                 dev=False, confirm_bounds=True, profile_start_height=None,
                 nc_level=None, legacy_peak_id=False, tail_number=None):
        """ Creates a Profiles object.

        :param int resolution: resolution to which data should be
           calculated in units of altitude or pressure
        :param str res_units: units of resolution in a format which can \
           be parsed by pint
        :param bool ascent: True to use ascending leg of flight, False to use \
           descending leg
        :param bool dev: True if data is from a developmental flight
        :param confirm_bounds: False to bypass user confirmation of \
           automatically identified start, peak, and end times
        :param int profile_start_height: if provided, the user will not be \
           prompted to enter the starting height for each profile separately.\
           This can be usefull when processing many profiles from the same \
           mission, but at least one profile should be processed without this \
           parameter to determine its correct value.
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'. To generate a single, Profile_Set-level file, call \
           Profile_Set.save_netCDF where you are done adding data.
        """
        self.resolution = resolution
        self.res_units = res_units
        self.ascent = ascent
        self.dev = dev
        self.confirm_bounds = confirm_bounds
        self.profiles = []
        self.legacy_peaks = legacy_peak_id
        self.tail_number = tail_number

        if profile_start_height is not None:
            self.profile_start_height = profile_start_height * units.m
            self._base_start = profile_start_height * units.m
        else:
            self.profile_start_height = None
            self._base_start = None

        #self.meta = None
        self._nc_level = nc_level
        self._root_dir = ""

    def add_all_profiles(self, file_path, metadata=None):
        """ Read a file, split it into profiles, and add them all.

        Deprecated. This now delegates to profiles.processing, which is
        where the leg detection and grid standardisation actually live.

        :param str file_path: the data file
        :param metadata: a Meta object, or None
        :rtype: int
        :return: the number of profiles held after adding
        """
        warnings.warn(
            'Profile_Set is deprecated; use '
            'profiles.processing.process_flights() with a ProcessingConfig.',
            DeprecationWarning, stacklevel=2)

        file_path = os.path.abspath(file_path)

        config = processing.ProcessingConfig(
            resolution=self.resolution, res_units=self.res_units,
            ascent=self.ascent, dev=self.dev,
            confirm_bounds=self.confirm_bounds,
            profile_start_height=(self.profile_start_height.magnitude
                                  if self.profile_start_height is not None
                                  else None),
            nc_level=self._nc_level, legacy_peak_id=self.legacy_peaks,
            tail_number=self.tail_number)

        profiles = processing.profiles_from_flight(
            file_path, config, metadata=metadata, base_start=self._base_start)

        if self._base_start is None and profiles:
            self._base_start = profiles[0]._base_start
            profiles = processing.profiles_from_flight(
                file_path, config, metadata=metadata,
                base_start=self._base_start)

        self.profiles.extend(profiles)
        self.profiles.sort()
        print(len(self.profiles), "profile(s) including those added from file",
              file_path)
        return len(self.profiles)

    def merge(self, to_add):
        """ Loads all Profile objects from a pre-existing Profiles into this
        Profiles. All flights must be from the same location.

        :param Profiles to_add: the Profiles object to be merged in
        """

        if to_add.resolution != self.resolution or \
           to_add.res_units != self.res_units:
            print("NOTICE: All future Profiles added will have resolution "
                  + str(self.resolution*self.res_units))

        if to_add.ascent != self.ascent:
            if self.ascent:
                print("NOTICE: All future Profiles added will be treated as \
                      ascending")
            else:
                print("NOTICE: All future Profiles added will be treated as \
                      descending")

        if to_add.profile_start_height != self.profile_start_height:
            print("NOTICE: All future Profiles added will start at height "
                  + str(self.profile_start_height))

        if to_add.dev != self.dev:
            if self.dev:
                print("NOTICE: All future Profiles added will be considered \
                      developmental")
            else:
                print("NOTICE: All future Profiles added will be considered \
                      operational")

        if len(self.profiles) > 0:
            units = self.profiles[0]._units

        for new_profile in to_add.profiles:

            self.profiles.append(deepcopy(new_profile))  # this doesn't include units

    def save_netCDF(self, file_path):
        """
        Stores all attributes of this Profile_Set object as a NetCDF

        :param string file_path: the file name to which attributes should be
           saved

        """

        file_name = str(self.profiles[0]._thermo_profile._meta.get("location"))[:5] + \
                    str(self.profiles[0]._thermo_profile._meta.get("platform_id")) + ".c1." + \
                    self.profiles[0]._thermo_profile._meta.get('timestamp').split('_')[0] + ".cdf"
        f_name = os.path.join(os.path.dirname(file_path), file_name)
        main_file = netCDF4.Dataset(f_name,
                                    "w", format="NETCDF4", mmap=False)
        # if self._meta is not None:
        #     self._meta.write_public_meta(
        #         os.path.join(self._root_dir,"processed", file_path)[:-3]
        #         + "_meta.txt")
        main_file.dev = str(self.dev)
        main_file.resolution = self.resolution
        main_file.res_units = self.res_units
        main_file.ascent = str(self.ascent)

        for i in range(len(self.profiles)):
            profile_group = main_file.createGroup("Profile" + str(i))
            profile_group.createDimension("time", None)

            #
            # Thermo
            #

            thermo = self.profiles[i]._thermo_profile

            if thermo is not None:

                # RH
                try:
                    rh_var = profile_group.createVariable("rh", "f8",
                                                          ("time",))
                    rh_var[:] = thermo.rh.magnitude
                    rh_var.units = str(thermo.rh.units)
                except Exception:
                    continue
                # TEMP
                try:
                    temp_var = profile_group.createVariable("temp", "f8",
                                                            ("time",))
                    temp_var[:] = thermo.temp.magnitude
                    temp_var.units = str(thermo.temp.units)
                except Exception:
                    continue
                # MIXING RATIO
                try:
                    mr_var = profile_group.createVariable("mr", "f8",
                                                          ("time",))
                    mr_var[:] = thermo.mixing_ratio.magnitude
                    mr_var.units = str(thermo.mixing_ratio.units)
                except Exception:
                    continue
                # POTENTIAL TEMPERATURE
                try:
                    theta_var = profile_group.createVariable("theta", "f8",
                                                          ("time",))
                    theta_var[:] = thermo.theta.magnitude
                    theta_var.units = str(thermo.theta.units)
                except Exception:
                    continue
                # DEWPOINT TEMPERATURE
                try:
                    dewp_var = profile_group.createVariable("T_d", "f8",
                                                          ("time",))
                    dewp_var[:] = thermo.T_d.magnitude
                    dewp_var.units = str(thermo.T_d.units)
                except Exception:
                    continue
                # TIME
                try:
                    time_var = profile_group.createVariable("time", "f8",
                                                            ("time",))
                    time_var[:] = netCDF4.date2num(thermo.gridded_times,
                                                   units='microseconds since \
                                                   2010-01-01 00:00:00:00')
                    time_var.units = \
                        'microseconds since 2010-01-01 00:00:00:00'
                except Exception:
                    continue
                # ALT
                try:
                    alt_var = profile_group.createVariable("alt", "f8",
                                                           ("time",))
                    alt_var[:] = thermo.alt.magnitude
                    alt_var.units = str(thermo.alt.units)
                except Exception:
                    continue
                # PRES
                try:
                    pres_var = profile_group.createVariable("pres", "f8",
                                                            ("time",))
                    pres_var[:] = thermo.pres.magnitude
                    pres_var.units = str(thermo.pres.units)
                except Exception:
                    continue
                #LAT
                try:
                    lat_var = profile_group.createVariable("lat", "f8", ("time",))
                    lat_var[:] = thermo.lat.magnitude
                    lat_var.units = str(thermo.lat.units)
                except Exception:
                    continue
                # LON
                try:
                    lon_var = profile_group.createVariable("lon", "f8", ("time",))
                    lon_var[:] = thermo.lon.magnitude
                    lon_var.units = str(thermo.lon.units)
                except Exception:
                    continue

            #
            # Wind
            #

            wind = self.profiles[i]._wind_profile

            if wind is not None:

                # DIRECTION
                try:
                    dir_var = profile_group.createVariable("dir", "f8",
                                                           ("time",))
                    dir_var[:] = wind.dir.magnitude
                    dir_var.units = str(wind.dir.units)
                except Exception:
                    continue
                # SPEED
                try:
                    spd_var = profile_group.createVariable("speed", "f8",
                                                           ("time",))
                    spd_var[:] = wind.speed.magnitude
                    spd_var.units = str(wind.speed.units)
                except Exception:
                    continue
                # U
                try:
                    u_var = profile_group.createVariable("u", "f8",
                                                         ("time",))
                    u_var[:] = wind.u.magnitude
                    u_var.units = str(wind.u.units)
                except Exception:
                    continue
                # V
                try:
                    v_var = profile_group.createVariable("v", "f8", ("time",))
                    v_var[:] = wind.v.magnitude
                    v_var.units = str(wind.v.units)
                except Exception:
                    continue
                # PRES
                try:
                    pres_var = profile_group.createVariable("pres", "f8",
                                                            ("time",))
                    pres_var[:] = wind.pres.magnitude
                    pres_var.units = str(wind.pres.units)
                except Exception:
                    continue
                # TIME
                try:
                    time_var = profile_group.createVariable("time", "f8",
                                                            ("time",))
                    time_var[:] = netCDF4.date2num(wind.gridded_times,
                                                   units='microseconds since \
                                                   2010-01-01 00:00:00:00')
                    time_var.units = \
                        'microseconds since 2010-01-01 00:00:00:00'
                except Exception:
                    continue
                # ALT
                try:
                    alt_var = profile_group.createVariable("alt", "f8",
                                                           ("time",))
                    alt_var[:] = wind.alt.magnitude
                    alt_var.units = str(wind.alt.units)
                except Exception:
                    continue
                #LAT
                try:
                    lat_var = profile_group.createVariable("lat", "f8", ("time",))
                    lat_var[:] = wind.lat.magnitude
                    lat_var.units = str(wind.lat.units)
                except Exception:
                    continue
                # LON
                try:
                    lon_var = profile_group.createVariable("lon", "f8", ("time",))
                    lon_var[:] = wind.lon.magnitude
                    print(lon_var)
                    lon_var.units = str(wind.lon.units)
                except Exception:
                    continue


        #
        # META
        #
        # if self.meta is not None:
        #     # print("\n\n" + str(self.meta.public_fields) + "\n\n")
        #     for key in np.unique(self.meta.public_fields):
        #         if self.meta.get(key) is not None:
        #             # print(key)
        #             main_file.key = self.meta.get(key)
        #             main_file.renameAttribute("key", key)

        main_file.close()
    """

    def __str__(self):
        to_return = "=====================================================\n" \
                    + "Profile Set with " + str(len(self.profiles)) + \
                    " Profiles\n"
        for profile in self.profiles:
            to_return = to_return + "\t" + str(profile) + "\n"
        return to_return
    """
