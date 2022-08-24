"""
Manages data from a single flight or profile
"""
from metpy.units import units
import profiles.utils as utils
import sys
import os
from profiles.Raw_Profile import Raw_Profile
from profiles.Thermo_Profile import Thermo_Profile
from profiles.Wind_Profile import Wind_Profile
from copy import deepcopy, copy
import numpy as np


class Profile():
    """ A Profile object contains data from a profile (if altitude or pressure
    is specified under resolution) or flight (if the resolution is in units
    of time)

    :var bool dev: True if data is from developmental flights
    :var Quantity resolution: resolution of the data in units of altitude or \
       pressure
    :var tuple indices: the bounds of the profile to be processed as\
       (start_time, end_time)
    :var bool ascent: True if data from the ascending leg should be processed,\
       otherwise the descending leg will be processed instead
    :var String file_path: the path to you .bin, .json, or .nc data file
    :var np.Array<Datetime> gridded_times: the times at which data points are \
       generated
    :var np.Array<Quantity> gridded_base: the value of the vertical coordinate\
       at each data point
    """

    # def __init__(self, *args, **kwargs):
    #     if len([*args]) > 0:
    #         self._init2(*args, **kwargs)

    def __init__(self, file_path, resolution, res_units, profile_num,
               ascent=True, dev=False, confirm_bounds=True,
               index_list=None, scoop_id=None, raw_profile=None,
               profile_start_height=None, nc_level='low', base_start=None,
               metadata=None):
        """ Creates a Profile object.

        :param string file_path: data file
        :param int resolution: resolution to which data should be
           calculated in units of altitude or pressure
        :param str res_units: units of resolution in a format which can \
           be parsed by pint
        :param int profile_num: 1 or greater. Identifies profile when file \
           contains multiple
        :param bool ascent: True to use ascending leg of flight, False to use \
           descending leg
        :param bool dev: True if data is from a developmental flight
        :param confirm_bounds: False to bypass user confirmation of \
           automatically identified start, peak, and end times
        :param list<tuple> index_list: Profile start, peak, and end indices if\
           known - leave as None in most cases
        :param str scoop_id: the sensor package used, if known
        :param Raw_Profile raw_profile: the partially-processed file - use \
           this if you have it, there's no need to make the computer do extra \
           work.
        :param int profile_start_height: if given, replaces prompt to user \
           asking for starting height of a profile. Recommended value is None\
           if you're only processing one profile.
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        :param Quantity base_start: lowest altitude value of the base after gridding
        :param profiles.Meta metadata: Meta object with metadata included
        """


        self._nc_level = nc_level

        # self.jms_path = coefs_path
        if raw_profile is not None:
            self._raw_profile = raw_profile
        else:
            self._raw_profile = Raw_Profile(file_path, dev, scoop_id,
                                            metadata=metadata)
        self._units = self._raw_profile.get_units()
        self._pos = self._raw_profile.pos_data()
        self._pres = (self._raw_profile.pres[0], self._raw_profile.pres[-1])
        self.meta = self._raw_profile.meta
        file_path = self._raw_profile.file_path

        if profile_start_height is not None:
            profile_start_height = profile_start_height
        try:
            if index_list is None:
                index_list = \
                 utils.identify_profile(self._pos["alt_MSL"],
                                        self._pos["time"], confirm_bounds,
                                        profile_start_height=\
                                        profile_start_height)
            indices = index_list[profile_num - 1]
        except IndexError:
            print("Analysis shows that the given file has fewer than " +
                  str(profile_num) + " profiles. If you are certain the file "
                  + "does contain more profiles than we have found, try again "
                  + "with a different starting height. \n\n")
            return self.__init__(file_path, resolution, res_units, profile_num,
                                 ascent=True, dev=False, confirm_bounds=True)

        if ascent:
            self.indices = (indices[0], indices[1])
        else:
            self.indices = (indices[1], indices[2])
        self._wind_profile = None
        self._thermo_profile = None
        self.dev = dev  # TODO this is not used
        self.resolution = resolution * self._units.parse_expression(res_units)
        self.ascent = ascent
        if ".nc" in file_path or ".NC" in file_path:
            self.file_path = file_path[:-3]
        elif ".json" in file_path or ".JSON" in file_path:
            self.file_path = file_path[:-5]
        elif ".bin" in file_path or ".BIN" in file_path:
            self.file_path = file_path[:-4]
        else:
            print("File type not recognized")
            sys.exit(0)

        if(self.resolution.dimensionality ==
           self._units.get_dimensionality('m')):
            base = self._pos['alt_MSL']
            base_time = self._pos['time']
        elif(self.resolution.dimensionality ==
             self._units.get_dimensionality('Pa')):
            base = self._pres[0]
            base_time = self._pres[1]
        self.gridded_times, self.gridded_base \
                = utils.regrid_base(base=base, base_times=base_time,
                                    new_res=self.resolution, ascent=ascent,
                                    units=self._units, indices=self.indices,
                                    base_start=base_start)
        self._base_start = self.gridded_base[0]

        self.__load_pos__()


    def __load_pos__(self):
        self.lat = []
        self.lon = []
        self.alt_MSL = []
        for item in utils.regrid_data_group(data=list(zip(self._pos['lat'], self._pos['lon'], self._pos['alt_MSL'])), data_times=self._pos['time'], gridded_times=self.gridded_times):
            self.lat.append(
                np.nanmean(list(map(lambda latlon: latlon[0].magnitude, item['values'])))
            )
            self.lon.append(
                np.nanmean(list(map(lambda latlon: latlon[1].magnitude, item['values'])))
            )
            self.alt_MSL.append(
                np.nanmean(list(map(lambda latlon: latlon[2].magnitude, item['values'])))
            )
        self.lat =  np.array(self.lat) * units.deg
        self.lon =  np.array(self.lon) * units.deg
        self.alt_MSL =  np.array(self.alt_MSL) * units.m 

    def get(self, varname):
        """
        Returns the requested variable, which may be in Profile or one of its
        attributes (ex. temp is in thermo_profile)

        :param str varname: the name of the requested variable
        :return: the requested variable
        """

        if varname in ['lat', 'lon', 'alt_MSL']:
            return self.__getattribute__(varname)

        try:
            return self.__getattribute__(varname)
        except AttributeError:
            pass
        if self._thermo_profile is not None:
            try:
                return self._thermo_profile.__getattribute__(varname)
            except AttributeError:
                pass
        if self._wind_profile is not None:
            try:
                return self._wind_profile.__getattribute__(varname)
            except AttributeError:
                pass
        try:
            return self._raw_profile.__getattribute__(varname)
        except AttributeError:
            print("The requested variable " + varname + " does not exist. Call "
                  "get_thermo_profile and get_wind_profile before trying "
                  "again.")

    def get_wind_profile(self, file_path=None):
        """ If a Wind_Profile object does not already exist, it is created when
        this method is called.

        :return: the Wind_Profile object
        :rtype: Wind_Profile
        """

        if file_path is None:
            file_path = self.file_path

        if self._wind_profile is None:
            wind_data = self._raw_profile.wind_data()
            self._wind_profile = \
                Wind_Profile(wind_data, self.resolution,
                             gridded_times=self.gridded_times,
                             gridded_base=self.gridded_base,
                             indices=self.indices, ascent=self.ascent,
                             units=self._units, file_path=file_path,
                             pos=self._pos,
                             meta=self.meta,
                             nc_level=self._nc_level)
            if len(self._wind_profile.gridded_times) > len(self.gridded_times):
                new_len = len(self.gridded_times)
                self._wind_profile.trucate_to(new_len)
            elif len(self._wind_profile.gridded_times) < \
                    len(self.gridded_times):
                new_len = len(self._wind_profile.gridded_times)
                self.gridded_times = self.gridded_times[:new_len]
                self.gridded_base = self.gridded_base[:new_len]
            if self._thermo_profile is not None:
                new_len = len(self._wind_profile.gridded_times)
                self._thermo_profile.truncate_to(new_len)
        return self._wind_profile

    def get_thermo_profile(self, file_path=None):
        """ If a Thermo_Profile object does not already exist, it is created
        when this method is called.

        :return: the Thermo_Profile object
        :rtype: Thermo_Profile
        """
        if file_path is None:
            file_path = self.file_path

        if self._thermo_profile is None:
            thermo_data = self._raw_profile.thermo_data()
            self._thermo_profile = \
                Thermo_Profile(thermo_data, self.resolution,
                               gridded_times=self.gridded_times,
                               gridded_base=self.gridded_base,
                               indices=self.indices, ascent=self.ascent,
                               units=self._units, file_path=file_path,
                               pos=self._pos,
                               meta=self.meta,
                               nc_level=self._nc_level)
            if len(self._thermo_profile.gridded_times) > \
                    len(self.gridded_times):
                new_len = len(self.gridded_times)
                self._thermo_profile.trucate_to(new_len)
            elif len(self._thermo_profile.gridded_times) < \
                    len(self.gridded_times):
                new_len = len(self._thermo_profile.gridded_times)
                self.gridded_times = self.gridded_times[:new_len]
                self.gridded_base = self.gridded_base[:new_len]
            if self._wind_profile is not None:
                new_len = len(self._thermo_profile.gridded_times)
                self._wind_profile.truncate_to(new_len)
        return self._thermo_profile

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            print(key)
            if key in "_units":
                continue
            if key in "_pos":
                setattr(result, key, copy(value))
                continue
            try:
                value = value.magnitude
            except AttributeError:
                value = value
            setattr(result, key, deepcopy(value, memo))
        return result

    def __str__(self):
        if self._wind_profile is not None:
            wind_str = str(self._wind_profile)
        else:
            wind_str = ""
        if self._thermo_profile is not None:
            thermo_str = str(self._thermo_profile)
        else:
            thermo_str = ""
        return "Profile object:\n\t\tLocation data: " + str(type(self._pos)) \
               + "\n" + wind_str + thermo_str

    def __gt__(self, other):
        if self.gridded_times[0] > other.gridded_times[0]:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.gridded_times[0] < other.gridded_times[0]:
            return True
        else:
            return False

    def __eq__(self, other):
        def __lt__(self, other):
            if self.gridded_times[0] == other.gridded_times[0]:
                return True
            else:
                return False
