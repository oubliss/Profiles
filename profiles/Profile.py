
"""
Manages data from a single flight or profile
"""
from datetime import datetime, timedelta
from metpy.units import units
import profiles.utils as utils
import profiles
import sys
import os
from profiles.Raw_Profile import Raw_Profile
from profiles.Thermo_Profile import Thermo_Profile
from profiles.Wind_Profile import Wind_Profile
from profiles.Coef_Manager import Coef_Manager
from copy import deepcopy, copy
import numpy as np
import netCDF4


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
        self._wind_data = self._raw_profile.wind_data().copy()
        self._thermo_data = self._raw_profile.thermo_data().copy()
        self.meta = self._raw_profile.meta
        file_path = self._raw_profile.file_path

        if profile_start_height is not None:
            profile_start_height = profile_start_height * self._units.m
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
        self._ascent_filename_tag = 'ascent' if ascent else 'descent'

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

            # self.time/self.alt is the mean time/alt between two values in gridded_times
            # whereas gridded_times/gridded_base are the edges to use in the averaging in
            # get_wind_profile and get_thermo_profile
            self.time, self.alt, self.gridded_times, self.gridded_base \
                = utils.regrid_base(base=base, base_times=base_time,
                                    new_res=self.resolution, ascent=ascent,
                                    units=self._units, indices=self.indices,
                                    base_start=base_start)

            self.pres = utils.regrid_data(data=self._raw_profile.pres[0],
                                          data_times=self._raw_profile.pres[-1],
                                          gridded_times=self.gridded_times,
                                          units=self._units)

        elif(self.resolution.dimensionality ==
             self._units.get_dimensionality('Pa')):
            base = self._pres[0]
            base_time = self._pres[1]

            self.time, self.pres, self.gridded_times, self.gridded_base \
                    = utils.regrid_base(base=base, base_times=base_time,
                                        new_res=self.resolution, ascent=ascent,
                                        units=self._units, indices=self.indices,
                                        base_start=base_start)

            self.alt = utils.regrid_data(data=self._pos['alt_MSL'],
                                          data_times=self._pos['time'],
                                          gridded_times=self.gridded_times,
                                          units=self._units)


        self._base_start = self.gridded_base[0]
        self.copter_id = self._raw_profile.serial_numbers['copterID']
        self.tail_number = Coef_Manager().get_tail_n(self.copter_id)

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

    def lowpass_filter(self, wind=True, thermo=True, Fs=10., Fc=0.1, n=501):
        """
        Based on contributions from Dr. Brian Green (OU SoM)

        Lowpass zero-phase FIR filter the raw data read from JSON files
        Also filter the thrust vectors here for better speed and direction
        Finally, apply response time correction convolution to the keys in tau
        Fs: sampling frequency in Hz
        Fc: cutoff frequency in Hz
        n: number of filter coefficients to include
        """
        from scipy import signal

        # Make a copy of the data from raw profile. ALWAYS want to filter from the raw data
        wind_data = self._raw_profile.wind_data().copy()
        thermo_data = self._raw_profile.thermo_data().copy()

        # Set up my filter
        # F = Fc/Fs  # Don't need this since I'm keeping things in physical units
        fir_coef = signal.remez(n,                       # Number of coeffs
                                [0., Fc/10, Fc, 0.5*Fs], # Bands to specify
                                [1, 0],                  # Desired gain of the bands
                                fs=Fs)                   # Sampling Frequency
        # apply filter to elements
        N = len(fir_coef)
        N0 = int((N-1)/2)

        # Determine the wind data we want to filter (manually set here)
        vars_to_filter = ['roll', 'pitch', 'yaw']
        if wind:
            for var in wind_data.keys():
                if var not in vars_to_filter:
                    continue
                # print(f"    Filtering {var}")

                # # zero pad the signal to the left and right
                # aux = np.pad(wind_data[var].magnitude, N0, mode="constant")
                #
                # # Pre-allocate memory
                # filt1 = np.zeros(len(wind_data[var]), dtype=float)
                #
                # # Run filter over input signal
                # i0 = int((N-1)/2 + 1)
                # i1 = len(filt1)+1
                # for i in range(i0, i1):
                #     filt1[i-N0] = np.sum(fir_coef * aux[(i-N0):(i+N0+1)])

                filt2 = signal.filtfilt(fir_coef, 1, wind_data[var].magnitude, padlen=N0, padtype="constant")

                self._wind_data[var] = filt2 * wind_data[var].units

        # Determine the wind data we want to filter (manually set here)
        vars_to_filter = ['temp1', 'temp2', 'temp3', 'temp4',
                  'rh1', 'rh2', 'rh3', 'rh4',
                  'resi1', 'resi2', 'resi3', 'resi4']
        if thermo:
            for var in thermo_data.keys():
                if var not in vars_to_filter:
                    continue

                # zero pad the signal to the left and right
                # aux = np.pad(thermo_data[var].magnitude, N0, mode="edge")
                #
                # # Pre-allocate memory
                # filt1 = np.zeros(len(thermo_data[var]), dtype=float)
                #
                # # Run filter over input signal
                # i0 = int((N - 1) / 2 + 1)
                # i1 = len(filt1) + 1
                # for i in range(i0, i1):
                #     filt1[i - N0] = np.sum(fir_coef * aux[(i - N0):(i + N0 + 1)])

                filt2 = signal.filtfilt(fir_coef, 1, thermo_data[var].magnitude, padlen=N0, padtype="constant")
                self._thermo_data[var] = filt2 * thermo_data[var].units

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

    def get_wind_profile(self, file_path=None, algorithm='linear'):
        """ If a Wind_Profile object does not already exist, it is created when
        this method is called.

        :return: the Wind_Profile object
        :rtype: Wind_Profile
        """

        if file_path is None:
            file_path = self.file_path

        if self._wind_profile is None:
            wind_data = self._wind_data
            self._wind_profile = \
                Wind_Profile(wind_data, self.resolution,
                             algorithm=algorithm,
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
            thermo_data = self._thermo_data
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

    def save_netcdf(self, file_path=None):
        if file_path is None:
            file_path = self.file_path

        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self.meta is not None:
            file_name = str(self.meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                        str(self.meta.get("platform_id")) + "CMT"  + self._ascent_filename_tag + ".c1." + \
                        self.time[0].strftime("%Y%m%d.%H%M%S") + ".cdf"
            file_name = os.path.join(os.path.dirname(file_path), file_name)

        else:
            raise IOError("Please specify a file name or include metadata when saving Profile netcdfs")

        if self._wind_profile is None and self._thermo_profile is None:
            print("No wind or thermo profile to save")
            return

        main_file = netCDF4.Dataset(file_name, "w", format="NETCDF4", mmap=False)

        # Vital ncattrs
        main_file.setncattr("conventions", "NC-1.8")
        main_file.setncattr("processing_version", profiles.__version__)
        main_file.setncattr("copter_id", self.copter_id)
        main_file.setncattr("tail_number", self.tail_number)
        
        try:
            main_file.setncattr("flight_location", utils.get_place_from_lat_lon(self.lat[0].magnitude, self.lon[0].magnitude))
        except Exception:
            pass 
        
        main_file.setncattr('datafile_created_on_date', datetime.utcnow().isoformat())
        main_file.setncattr('datafile_created_on_machine',  os.uname().nodename)

        main_file.setncattr("reference1", "Segales, A. R., B. R. Greene, T. M. Bell, W. Doyle, J. J. Martin, "
                                          "E. A. Pillar-Little, and P. B. Chilson, 2020: The CopterSonde: an insight"
                                          " into the development of a smart unmanned aircraft system for atmospheric "
                                          "boundary layer research. Atmospheric Measurement Techniques, 13, 2833–2848, "
                                          "https://doi.org/10.5194/amt-13-2833-2020.")
        main_file.setncattr("reference2", "Bell, T. M., B. R. Greene, P. M. Klein, M. Carney, and "
                                          "P. B. Chilson, 2020: Confronting the boundary layer data gap: evaluating new "
                                          "and existing methodologies of probing the lower atmosphere. Atmospheric "
                                          "Measurement Techniques, 13, 3855–3872, https://doi.org/10.5194/amt-13-3855-2020.")

        # Create the dimensions
        main_file.createDimension("time", None)

        # TIME
        # Be sure to use self.time instead of self.gridded_times since we want to store the mean time between two levels
        # of gridded times
        time_var = main_file.createVariable("time", "f8", ("time",))
        time_var[:] = netCDF4.date2num(self.time,
                                       units='microseconds since \
                                                       2010-01-01 00:00:00:00')
        time_var.units = 'microseconds since 2010-01-01 00:00:00:00'

        # Do base_time and time_offset like ARM
        bt = abs((self.time[0] - datetime(1970, 1, 1)).total_seconds())
        bt_var = main_file.createVariable('base_time', 'i8')
        bt_var.setncattr('long_name', 'Base time in Epoch')
        bt_var.setncattr('ancillary_variables', 'time_offset')
        bt_var.setncattr('units', 'seconds since 1970-01-01 00:00:00 UTC')
        bt_var[:] = bt

        to = netCDF4.date2num(self.time,
                              units=f'seconds since {self.time[0]:%Y-%m-%d %H:%M:%S UTC}')
        to_var = main_file.createVariable('time_offset', 'f4', dimensions=('time',))
        to_var.setncattr('long_name', 'Time offset from base_time')
        to_var.setncattr('units', f'seconds since {self.time[0]:%Y-%m-%d %H:%M:%S UTC}')
        to_var.setncattr('ancillary_variables', 'base_time')
        to_var[:] = to

        # ALT
        alt_var = main_file.createVariable("alt", "f8", ("time",))
        alt_var[:] = self.alt_MSL.magnitude
        alt_var.units = str(self.alt_MSL.units)
        # PRES
        pres_var = main_file.createVariable("pres", "f8", ("time",))
        pres_var[:] = self.pres.magnitude
        pres_var.units = str(self.pres.units)
        # LAT
        lat_var = main_file.createVariable("lat", "f8", ("time",))
        lat_var[:] = self.lat.magnitude
        lat_var.units = str(self.lat.units)
        # LON
        lon_var = main_file.createVariable("lon", "f8", ("time",))
        lon_var[:] = self.lon.magnitude
        lon_var.units = str(self.lon.units)

        if self._thermo_profile is not None:
            # TEMP
            temp_var = main_file.createVariable("tdry", "f8", ("time",))
            temp_var[:] = self._thermo_profile.temp.magnitude
            temp_var.units = str(self._thermo_profile.temp.units)
            temp_var.long_name = "Dry bulb temperature"

            # MIXING RATIO
            mr_var = main_file.createVariable("mr", "f8", ("time",))
            mr_var[:] = self._thermo_profile.mixing_ratio.magnitude
            mr_var.units = str(self._thermo_profile.mixing_ratio.units)
            mr_var.long_name = "Water vapor mixing ratio"
            # THETA
            theta_var = main_file.createVariable("theta", "f8", ("time",))
            theta_var[:] = self._thermo_profile.theta.magnitude
            theta_var.units = str(self._thermo_profile.theta.units)
            theta_var.long_name = "Potential temperature"
            # T_D
            Td_var = main_file.createVariable("Td", "f8", ("time",))
            Td_var[:] = self._thermo_profile.T_d.magnitude
            Td_var.units = str(self._thermo_profile.T_d.units)
            Td_var.long_name = "Dew point temperature"
            # Q
            q_var = main_file.createVariable("q", "f8", ("time",))
            q_var[:] = self._thermo_profile.q.magnitude * 1e3
            q_var.units = str(self._thermo_profile.q.units)
            q_var.long_name = "Specific humidity"

        if self._wind_profile is not None:
            # DIRECTION
            dir_var = main_file.createVariable("dir", "f8", ("time",))
            dir_var[:] = self._wind_profile.dir.magnitude
            dir_var.units = str(self._wind_profile.dir.units)
            dir_var.long_name = "Wind direction"
            # SPEED
            spd_var = main_file.createVariable("wspd", "f8", ("time",))
            spd_var[:] = self._wind_profile.speed.magnitude
            spd_var.units = str(self._wind_profile.speed.units)
            spd_var.long_name = "Wind speed"
            # U
            u_var = main_file.createVariable("wind_u", "f8", ("time",))
            u_var[:] = self._wind_profile.u.magnitude
            u_var.units = str(self._wind_profile.u.units)
            u_var.long_name = "westward wind component"
            # V
            v_var = main_file.createVariable("wind_v", "f8", ("time",))
            v_var[:] = self._wind_profile.v.magnitude
            v_var.units = str(self._wind_profile.v.units)
            v_var.long_name = "northward wind component"

        # Close the netCDF file
        main_file.close()


        return None

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
