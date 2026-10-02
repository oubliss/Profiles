
"""
Manages data from a single flight or profile
"""
from datetime import datetime, timedelta
from profiles.unit_registry import units
import profiles.utils as utils
import profiles.qc as qc
import profiles.calibration as calibration
import metpy.calc
import warnings
from profiles.retrievals import thermo as thermo_retrieval
from profiles.retrievals import wind as wind_retrieval
import profiles
import sys
import os
from profiles.flight import FlightLog
from profiles.Coef_Manager import Coef_Manager
from copy import deepcopy, copy
import numpy as np
import netCDF4


#: Which time coordinate each thermodynamic variable is sampled on.
#: Fragments are tested in order, so the more specific ones come first:
#: "temp_rh1" must land on the humidity clock, not the temperature one.
#: fan_flag is deliberately absent - the original left it untrimmed, and
#: nothing downstream of the trim reads it.
_THERMO_TIME_BASES = {
    'resi': 'time_temp',
    'temp_rh': 'time_rh',
    'temp_pres': 'time_pres',
    'alt_pres': 'time_pres',
    'pres': 'time_pres',
    'rh': 'time_rh',
    'temp': 'time_temp',
}

#: Same for wind. Only the attitude and velocity series were ever trimmed;
#: pressure and altitude come off the Profile, which has already gridded
#: them, so they pass through.
_WIND_TIME_BASES = {
    'roll': 'time', 'pitch': 'time', 'yaw': 'time',
    'speed_east': 'time', 'speed_north': 'time', 'speed_down': 'time',
}


def _time_base_for(key, selectors):
    """ Which time array a variable is sampled against, or None to pass through.

    A time array is itself trimmed, by its own mask - that is what keeps
    each series the same length as the clock it is indexed by.

    :param str key: variable name
    :param dict selectors: fragment -> time key, most specific first
    :rtype: str or None
    """
    if key in set(selectors.values()):
        return key

    if 'time' in key or key in ('serial_numbers', 'units'):
        return None

    for fragment, time_key in selectors.items():
        if fragment in key:
            return time_key
    return None


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
               metadata=None, **kwargs):
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
        :param FlightLog raw_profile: the partially-processed file - use \
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
            # scoop_id is not a FlightLog parameter; passing it
            # positionally put it in the nc_level slot, which silently
            # suppressed the a0 file and the coefficient application.
            self._raw_profile = FlightLog(file_path, dev,
                                            nc_level=nc_level,
                                            metadata=metadata, **kwargs)
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

                try:
                    # Find the window where the scoop fan was active
                    foo = np.where(np.array(self._thermo_data['fan_flag']) > 0)[0]
                    fan_start = self._thermo_data['time_temp'][foo.min()]
                    fan_stop = self._thermo_data['time_temp'][foo.max()]

                    index_list = utils.identify_profile_peaks(self._pos["alt_MSL"].magnitude, self._pos['time'],
                                                              window=(fan_start, fan_stop))
                except ValueError:
                    index_list = utils.identify_profile(self._pos["alt_MSL"],
                                                        self._pos["time"], confirm_bounds,
                                                        profile_start_height=profile_start_height)

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
        self._wind_computed = False
        self._thermo_computed = False
        self.dev = dev  # TODO this is not used
        self.resolution = resolution * self._units.parse_expression(res_units)
        self.ascent = ascent
        self._ascent_filename_tag = 'ascent' if ascent else 'descent'

        if ".nc" in file_path or ".NC" in file_path:
            self.file_path = file_path[:-3]
        elif ".json" in file_path or ".JSON" in file_path:
            self.file_path = file_path[:-5]
        elif ".bin" in file_path or ".BIN" in file_path or ".csv" in file_path:
            self.file_path = file_path[:-4]
        else:
            print("File type not recognized")
            sys.exit(0)

        if(self.resolution.dimensionality ==
           self._units.get_dimensionality('m')):
            base = self._pos['alt_MSL']
            base_time = self._pos['time']

            # Two grids come back from regrid_base and they are NOT the same
            # length. gridded_times/gridded_base are the N+1 bin *edges* used
            # to delimit the averaging; self.time/self.alt are the N bin
            # *centres*, which is what every regridded variable lines up with.
            # Handing the edges to Thermo_Profile/Wind_Profile as their
            # vertical coordinate was a half-resolution low bias plus a
            # trailing fill row - see gridded_centers below.
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


        # The vertical coordinate at bin centres, whichever coordinate was
        # chosen. This is what the gridded variables are aligned to.
        if (self.resolution.dimensionality ==
                self._units.get_dimensionality('m')):
            self.gridded_centers = self.alt
        else:
            self.gridded_centers = self.pres

        self._base_start = self.gridded_base[0]
        try:
            self.copter_id = self._raw_profile.serial_numbers['copterID']
            self.tail_number = Coef_Manager().get_tail_n(self.copter_id)
        except Exception:
            self.copter_id = -999
            self.tail_number = self._raw_profile.tail_number


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
        if self._thermo_computed:
            try:
                return self.__getattribute__(varname)
            except AttributeError:
                pass
        if self._wind_computed:
            try:
                return self.__getattribute__(varname)
            except AttributeError:
                pass
        try:
            return self._raw_profile.__getattribute__(varname)
        except AttributeError:
            print("The requested variable " + varname + " does not exist. Call "
                  "compute_thermo and compute_wind before trying "
                  "again.")

    #: QC thresholds: (max spread of sensor means, max spread of sensor
    #: standard deviations), in each variable's own units.
    #: TODO move into a config object - Stage 5.
    QC_THRESHOLDS = {'temp': (0.25, 0.1), 'rh': (0.4, 0.2)}

    def _trim(self, data, selectors):
        """ Restrict a data dict to this profile's time bounds.

        Returns a new dict. The old Thermo_Profile and Wind_Profile trimmed
        the caller's dictionary in place, so constructing either one
        silently shortened Profile._thermo_data / _wind_data from 7286
        samples to 4050.

        :param dict data: keys to arrays, as from FlightLog.thermo_data()
        :param dict selectors: variable-name predicate -> time key, deciding
           which time base each variable is sampled on
        :rtype: dict
        """
        if self.indices[0] is None:
            return dict(data)

        masks = {}
        for time_key in set(selectors.values()):
            stamps = np.array(data[time_key])
            masks[time_key] = np.where(stamps > self.indices[0],
                                       stamps < self.indices[1], False)

        trimmed = {}
        for key, value in data.items():
            time_key = _time_base_for(key, selectors)
            if time_key is None:
                trimmed[key] = value
                continue

            keep = np.where(masks[time_key])
            if hasattr(value, 'magnitude'):
                trimmed[key] = value.magnitude[keep] * value.units
            else:
                trimmed[key] = np.array(value)[keep]

        return trimmed

    def _average_ensemble(self, series, flags):
        """ Mean across sensors, excluding any the flags reject.

        :param list series: one array per sensor
        :param list flags: one flag per sensor, 0 meaning good
        :rtype: np.ndarray
        """
        kept = [values for values, flag in zip(series, flags) if flag == qc.GOOD]
        if not kept:
            # Every sensor rejected: report NaN rather than an empty mean,
            # which is what the previous NaN-fill-then-nanmean produced.
            return np.full(len(series[0]), np.nan)
        return np.nanmean(np.vstack(kept), axis=0)

    def compute_thermo(self):
        """ Calibrate, QC and grid the thermodynamic variables.

        Populates temp, rh, mixing_ratio, theta, T_d, q and the per-sensor
        flag arrays temp_flags and rh_flags. alt, pres, time, lat and lon
        are already on the Profile - they were previously recomputed here
        from the same inputs, giving identical numbers three times over.

        :rtype: Profile
        """
        if self._thermo_computed:
            return self

        data = self._trim(self._thermo_data, _THERMO_TIME_BASES)
        serial_numbers = data['serial_numbers']

        temp_raw = calibration.calibrate_temperature(data, serial_numbers)
        rh_raw = calibration.calibrate_humidity(data, serial_numbers)

        self.temp_flags = qc.qc(temp_raw, *self.QC_THRESHOLDS['temp'])
        self.rh_flags = qc.qc(rh_raw, *self.QC_THRESHOLDS['rh'])

        temp = self._average_ensemble(temp_raw, self.temp_flags) \
            * self._units.kelvin
        rh = self._average_ensemble(rh_raw, self.rh_flags) \
            * self._units.percent

        self.temp = utils.regrid_data(data=temp, data_times=data['time_temp'],
                                      gridded_times=self.gridded_times,
                                      units=self._units)
        self.rh = utils.regrid_data(data=rh, data_times=data['time_rh'],
                                    gridded_times=self.gridded_times,
                                    units=self._units)

        derived = thermo_retrieval.derive(self.pres, self.temp, self.rh)
        self.mixing_ratio = derived['mixing_ratio']
        self.theta = derived['theta']
        self.T_d = derived['T_d']
        self.q = derived['q']

        self._thermo_computed = True
        if utils.writes_netcdf(self._nc_level):
            self._save_thermo_netCDF(self.file_path)
        return self

    def compute_wind(self, algorithm='linear'):
        """ Retrieve wind from airframe tilt and grid it.

        :param str algorithm: 'linear' or 'quadratic'
        :rtype: Profile
        """
        if self._wind_computed:
            return self

        if algorithm not in wind_retrieval.ALGORITHM_EQUATIONS:
            raise ValueError(
                f'unknown wind algorithm {algorithm!r}; available: '
                f'{sorted(wind_retrieval.ALGORITHM_EQUATIONS)}')

        data = self._trim(self._wind_data, _WIND_TIME_BASES)
        equation_name = wind_retrieval.ALGORITHM_EQUATIONS[algorithm]
        coefficients = utils.coef_manager.get_coefs(
            'Wind', self.tail_number, equation_name)

        direction, speed = wind_retrieval.retrieve(
            data['roll'], data['pitch'], data['yaw'],
            coefficients, equation_name)
        direction = direction % (2 * np.pi)

        self.dir = utils.regrid_data(data=direction, data_times=data['time'],
                                     gridded_times=self.gridded_times,
                                     units=self._units)
        self.speed = utils.regrid_data(data=speed, data_times=data['time'],
                                       gridded_times=self.gridded_times,
                                       units=self._units)
        self.u, self.v = metpy.calc.wind_components(self.speed, self.dir)

        self._wind_computed = True
        if utils.writes_netcdf(self._nc_level):
            self._save_wind_netCDF(self.file_path)
        return self

    def get_wind_profile(self, file_path=None, algorithm='linear'):
        """ Deprecated. Use compute_wind(); the wind lives on the Profile.

        :rtype: Profile
        """
        warnings.warn(
            'get_wind_profile() is deprecated; wind is computed onto the '
            'Profile itself by compute_wind().',
            DeprecationWarning, stacklevel=2)
        return self.compute_wind(algorithm=algorithm)

    def get_thermo_profile(self, file_path=None):
        """ Deprecated. Use compute_thermo(); the values live on the Profile.

        :rtype: Profile
        """
        warnings.warn(
            'get_thermo_profile() is deprecated; thermodynamic variables '
            'are computed onto the Profile itself by compute_thermo().',
            DeprecationWarning, stacklevel=2)
        return self.compute_thermo()

    # The per-variable thermo_/wind_ writers, moved here verbatim from the
    # classes they used to live on. Stage 5 replaces all four writers in
    # this file with one that takes a processing level.
    def _save_thermo_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """
        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self.meta is not None:
            file_name = str(self.meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                        str(self.meta.get("platform_id")) + "CMT" + \
                        "thermo_" + self._ascent_filename_tag + ".c1." + \
                        self.meta.get("timestamp").replace("_", ".") + ".cdf"
            file_name = os.path.join(os.path.dirname(file_path), file_name)

        else:
            raise IOError("Please specify a file name or include metadata when saving Profile netcdfs")


        main_file = netCDF4.Dataset(file_name, "w",
                                    format="NETCDF4", mmap=False)
        # File NC compliant to version 1.8
        main_file.setncattr("Conventions", "NC-1.8")

        #
        # Get the flags in
        #
        flag_dict = {0: "good",
                     2: "bias",
                     3: "lag",
                     4: "empty"}
        rh_flags = main_file.createGroup("rh_flags")
        for i in range(len(self.rh_flags)):
            rh_flags.setncattr("sensor" + str(i+1), flag_dict[self.rh_flags[i]])
        temp_flags = main_file.createGroup("temp_flags")
        for i in range(len(self.temp_flags)):
            temp_flags.setncattr("sensor" + str(i+1), flag_dict[self.temp_flags[i]])

        main_file.createDimension("time", None)
        # TIME
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
        to_var.setncattr('units', f'seconds since {self.gridded_times[0]:%Y-%m-%d %H:%M:%S UTC}')
        to_var.setncattr('ancillary_variables', 'base_time')
        to_var[:] = to

        # PRES
        pres_var = main_file.createVariable("pres", "f8", ("time",))
        pres_var[:] = self.pres.magnitude
        pres_var.units = str(self.pres.units)
        # RH
        rh_var = main_file.createVariable("rh", "f8", ("time",))
        rh_var[:] = self.rh.magnitude
        rh_var.units = str(self.rh.units)
        # ALT
        alt_var = main_file.createVariable("alt", "f8", ("time",))
        alt_var[:] = self.alt.magnitude
        alt_var.units = str(self.alt.units)
        # TEMP
        temp_var = main_file.createVariable("temp", "f8", ("time",))
        temp_var[:] = self.temp.magnitude
        temp_var.units = str(self.temp.units)
        # MIXING RATIO
        mr_var = main_file.createVariable("mr", "f8", ("time",))
        mr_var[:] = self.mixing_ratio.magnitude
        mr_var.units = str(self.mixing_ratio.units)
        # THETA
        theta_var = main_file.createVariable("theta", "f8", ("time",))
        theta_var[:] = self.theta.magnitude
        theta_var.units = str(self.theta.units)
        # T_D
        Td_var = main_file.createVariable("Td", "f8", ("time",))
        Td_var[:] = self.T_d.magnitude
        Td_var.units = str(self.T_d.units)
        # Q
        q_var = main_file.createVariable("q", "f8", ("time",))
        q_var[:] = self.q.magnitude
        q_var.units = str(self.q.units)
        # LAT
        lat_var = main_file.createVariable("lat", "f8", ("time",))
        lat_var[:] = self.lat.magnitude
        lat_var.units = str(self.lat.units)
        # LON
        lon_var = main_file.createVariable("lon", "f8", ("time",))
        lon_var[:] = self.lon.magnitude
        lon_var.units = str(self.lon.units)

        main_file.close()

    def _save_wind_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """

        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self.meta is not None:
            file_name = str(self.meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                    str(self.meta.get("platform_id")) + "CMT" + \
                    "wind_" + self._ascent_filename_tag + ".c1." + \
                    self.meta.get("timestamp").replace("_", ".") + ".cdf"
            file_name = os.path.join(os.path.dirname(file_path), file_name)

        else:
            raise IOError("Please specify a file name or include metadata when saving Profile netcdfs")


        main_file = netCDF4.Dataset(file_name, "w",
                                    format="NETCDF4", mmap=False)
        # File NC compliant to version 1.8
        main_file.setncattr("Conventions", "NC-1.8")
        
        main_file.createDimension("time", None)
        # DIRECTION
        dir_var = main_file.createVariable("dir", "f8", ("time",))
        dir_var[:] = self.dir.magnitude
        dir_var.units = str(self.dir.units)
        # SPEED
        spd_var = main_file.createVariable("speed", "f8", ("time",))
        spd_var[:] = self.speed.magnitude
        spd_var.units = str(self.speed.units)
        # U
        u_var = main_file.createVariable("u", "f8", ("time",))
        u_var[:] = self.u.magnitude
        u_var.units = str(self.u.units)
        # V
        v_var = main_file.createVariable("v", "f8", ("time",))
        v_var[:] = self.v.magnitude
        v_var.units = str(self.v.units)
        # ALT
        alt_var = main_file.createVariable("alt", "f8", ("time",))
        alt_var[:] = self.alt.magnitude
        alt_var.units = str(self.alt.units)
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

        # TIME
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
        to_var.setncattr('units', f'seconds since {self.gridded_times[0]:%Y-%m-%d %H:%M:%S UTC}')
        to_var.setncattr('ancillary_variables', 'base_time')
        to_var[:] = to

        main_file.close()

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

        if not (self._wind_computed or self._thermo_computed):
            print("No wind or thermo data to save; call compute_thermo() "
                  "and/or compute_wind() first")
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

        if self._thermo_computed:
            # TEMP
            temp_var = main_file.createVariable("tdry", "f8", ("time",))
            temp_var[:] = self.temp.magnitude
            temp_var.units = str(self.temp.units)
            temp_var.long_name = "Dry bulb temperature"

            # MIXING RATIO
            mr_var = main_file.createVariable("mr", "f8", ("time",))
            mr_var[:] = self.mixing_ratio.magnitude
            mr_var.units = str(self.mixing_ratio.units)
            mr_var.long_name = "Water vapor mixing ratio"
            # THETA
            theta_var = main_file.createVariable("theta", "f8", ("time",))
            theta_var[:] = self.theta.magnitude
            theta_var.units = str(self.theta.units)
            theta_var.long_name = "Potential temperature"
            # T_D
            Td_var = main_file.createVariable("Td", "f8", ("time",))
            Td_var[:] = self.T_d.magnitude
            Td_var.units = str(self.T_d.units)
            Td_var.long_name = "Dew point temperature"
            # Q
            q_var = main_file.createVariable("q", "f8", ("time",))
            q_var[:] = self.q.magnitude * 1e3
            q_var.units = str(self.q.units)
            q_var.long_name = "Specific humidity"

        if self._wind_computed:
            # DIRECTION
            dir_var = main_file.createVariable("dir", "f8", ("time",))
            dir_var[:] = self.dir.magnitude
            dir_var.units = str(self.dir.units)
            dir_var.long_name = "Wind direction"
            # SPEED
            spd_var = main_file.createVariable("wspd", "f8", ("time",))
            spd_var[:] = self.speed.magnitude
            spd_var.units = str(self.speed.units)
            spd_var.long_name = "Wind speed"
            # U
            u_var = main_file.createVariable("wind_u", "f8", ("time",))
            u_var[:] = self.u.magnitude
            u_var.units = str(self.u.units)
            u_var.long_name = "westward wind component"
            # V
            v_var = main_file.createVariable("wind_v", "f8", ("time",))
            v_var[:] = self.v.magnitude
            v_var.units = str(self.v.units)
            v_var.long_name = "northward wind component"

        # Close the netCDF file
        main_file.close()


        return None

    def save_cfnetcdf(self, platform_name, terrain_elevation, file_path=None):
        if file_path is None:
            file_path = self.file_path

        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self.meta is not None:
            file_name = str(self.meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                        str(self.meta.get("platform_id")) + "CMT" + self._ascent_filename_tag + ".c1." + \
                        self.time[0].strftime("%Y%m%d.%H%M%S") + ".cdf"
            file_name = os.path.join(os.path.dirname(file_path), file_name)

        else:
            raise IOError("Please specify a file name or include metadata when saving Profile netcdfs")

        if not (self._wind_computed or self._thermo_computed):
            print("No wind or thermo data to save; call compute_thermo() "
                  "and/or compute_wind() first")
            return

        main_file = netCDF4.Dataset(file_name, "w", format="NETCDF4")

        # Vital ncattrs
        main_file.setncattr("Conventions", "CF-1.8, WMO-CF-1.0")
        main_file.setncattr("wmo__cf_profile", "FM 303-2024")
        main_file.setncattr("featureType", "trajectory")
        main_file.setncattr("platform_name", f"{platform_name}")
        main_file.setncattr("flight_id", f"{platform_name}_{self.gridded_times[0]}")
        main_file.setncattr("site_terrain_elevation_height", terrain_elevation)
        main_file.setncattr("processing_level", 'c1')
        main_file.setncattr("processing_version", profiles.__version__)
        main_file.setncattr("copter_id", self.copter_id)
        main_file.setncattr("tail_number", self.tail_number)

        main_file.setncattr('datafile_created_on_date', datetime.utcnow().isoformat())
        main_file.setncattr('datafile_created_on_machine', os.uname().nodename)

        main_file.setncattr("Reference1", "Segales, A. R., B. R. Greene, T. M. Bell, W. Doyle, J. J. Martin, "
                                          "E. A. Pillar-Little, and P. B. Chilson, 2020: The CopterSonde: an insight"
                                          " into the development of a smart unmanned aircraft system for atmospheric "
                                          "boundary layer research. Atmospheric Measurement Techniques, 13, 2833–2848, "
                                          "https://doi.org/10.5194/amt-13-2833-2020.")
        main_file.setncattr("Reference1", "Bell, T. M., B. R. Greene, P. M. Klein, M. Carney, and "
                                          "P. B. Chilson, 2020: Confronting the boundary layer data gap: evaluating new "
                                          "and existing methodologies of probing the lower atmosphere. Atmospheric "
                                          "Measurement Techniques, 13, 3855–3872, https://doi.org/10.5194/amt-13-3855-2020.")

        # Create the dimensions
        main_file.createDimension("obs", None)

        # TIME
        # Be sure to use self.time instead of self.gridded_times since we want to store the mean time between two levels
        # of gridded times
        time_var = main_file.createVariable("time", "i8", ("obs",))
        time_var[:] = netCDF4.date2num(self.time,
                                       units='seconds since 1970-01-01T00:00:00')
        time_var.units = 'seconds since 1970-01-01T00:00:00'
        time_var.setncattr('standard_name', 'time')
        time_var.setncattr('long_name', 'time')
        time_var.setncattr('axis', 'T')

        # # Do base_time and time_offset like ARM
        # bt = abs((self.time[0] - datetime(1970, 1, 1)).total_seconds())
        # bt_var = main_file.createVariable('base_time', 'i8')
        # bt_var.setncattr('long_name', 'Base time in Epoch')
        # bt_var.setncattr('ancillary_variables', 'time_offset')
        # bt_var.setncattr('units', 'seconds since 1970-01-01 00:00:00 UTC')
        # bt_var[:] = bt
        #
        # to = netCDF4.date2num(self.time,
        #                       units=f'seconds since {self.time[0]:%Y-%m-%d %H:%M:%S UTC}')
        # to_var = main_file.createVariable('time_offset', 'f4', dimensions=('time',))
        # to_var.setncattr('long_name', 'Time offset from base_time')
        # to_var.setncattr('units', f'seconds since {self.time[0]:%Y-%m-%d %H:%M:%S UTC}')
        # to_var.setncattr('ancillary_variables', 'base_time')
        # to_var[:] = to

        # ALT
        alt_var = main_file.createVariable("altitude", "f8", ("obs",))
        alt_var.units = str(self.alt.units)
        alt_var.setncattr('long_name', 'altitude_above_sea_level')
        alt_var.setncattr('standard_name', 'altitude_above_sea_level')
        alt_var.setncattr('positive', 'up')
        alt_var.setncattr('axis', 'Z')
        alt_var[:] = self.alt.magnitude

        # PRES
        pres_var = main_file.createVariable("pressure", "f8", ("obs",))
        pres_var.units = str(self.pres.units)
        pres_var.setncattr('long_name', 'air_pressure')
        pres_var.setncattr('standard_name', 'pressure')
        pres_var.setncattr('positive', 'up')
        pres_var.setncattr('axis', 'Z')
        pres_var[:] = self.pres.magnitude

        # LAT
        lat_var = main_file.createVariable("lat", "f8", ("obs",))
        lat_var.units = 'degrees_north'
        lat_var.setncattr('standard_name', 'latitude')
        lat_var.setncattr('axis', 'Y')
        lat_var[:] = self.lat.magnitude

        # LON
        lon_var = main_file.createVariable("lon", "f8", ("obs",))
        lon_var.units = 'degrees_east'
        lon_var.setncattr('standard_name', 'longitude')
        lon_var.setncattr('axis', 'X')
        lon_var[:] = self.lon.magnitude

        if self._thermo_computed:
            # TEMP
            temp_var = main_file.createVariable("air_temperature", "f8", ("obs",))
            temp_var[:] = self.temp.magnitude
            temp_var.units = str(self.temp.units)
            temp_var.standard_name = "air_temperature"
            temp_var.long_name = "bulk temperature of the air"

            temp_var = main_file.createVariable("relative_humidity", "f8", ("obs",))
            temp_var[:] = self.rh.magnitude
            temp_var.units = str(self.rh.units)
            temp_var.standard_name = "relative_humidity"
            temp_var.long_name = "Relative humidity of the air"

            # MIXING RATIO
            mr_var = main_file.createVariable("humidity_mixing_ratio", "f8", ("obs",))
            mr_var[:] = self.mixing_ratio.magnitude
            mr_var.units = 'kg/kg'
            mr_var.standard_name = "humidity_mixing_ratio"
            mr_var.long_name = "Humidity Mixing Ratio"
            # # THETA
            # theta_var = main_file.createVariable("theta", "f8", ("time",))
            # theta_var[:] = self.theta.magnitude
            # theta_var.units = str(self.theta.units)
            # theta_var.long_name = "Potential temperature"
            # # T_D
            Td_var = main_file.createVariable("dewpoint", "f8", ("obs",))
            Td_var[:] = self.T_d.magnitude
            Td_var.units = str(self.T_d.units)
            Td_var.long_name = "Dew point temperature"
            Td_var.standard_name = "dew_point_temperature"
            # # Q
            # q_var = main_file.createVariable("q", "f8", ("time",))
            # q_var[:] = self.q.magnitude * 1e3
            # q_var.units = str(self.q.units)
            # q_var.long_name = "Specific humidity"

        if self._wind_computed:
            # DIRECTION
            dir_var = main_file.createVariable("wind_direction", "f8", ("obs",))
            dir_var[:] = self.dir.magnitude
            dir_var.units = str(self.dir.units)
            dir_var.standard_name = 'wind_from_direction'
            dir_var.long_name = "Wind direction"
            # SPEED
            spd_var = main_file.createVariable("wind_speed", "f8", ("obs",))
            spd_var[:] = self.speed.magnitude
            spd_var.units = str(self.speed.units)
            spd_var.standard_name = "wind_speed"
            spd_var.long_name = "Wind speed"
            # # U
            # u_var = main_file.createVariable("wind_u", "f8", ("time",))
            # u_var[:] = self.u.magnitude
            # u_var.units = str(self.u.units)
            # u_var.long_name = "westward wind component"
            # # V
            # v_var = main_file.createVariable("wind_v", "f8", ("time",))
            # v_var[:] = self.v.magnitude
            # v_var.units = str(self.v.units)
            # v_var.long_name = "northward wind component"

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
        computed = [name for name, done in
                    (('thermo', self._thermo_computed),
                     ('wind', self._wind_computed)) if done]
        return (f'Profile({self._ascent_filename_tag}, '
                f'{len(self.gridded_centers)} levels at {self.resolution:~P}, '
                f'starting {self.time[0]:%Y-%m-%d %H:%M:%S}Z, '
                f'computed: {", ".join(computed) or "none"})')

    def _sort_key(self):
        """Profiles order by when they started."""
        return self.gridded_times[0]

    def __lt__(self, other):
        return self._sort_key() < other._sort_key()

    def __gt__(self, other):
        return self._sort_key() > other._sort_key()

    def __le__(self, other):
        return self._sort_key() <= other._sort_key()

    def __ge__(self, other):
        return self._sort_key() >= other._sort_key()

    def __eq__(self, other):
        # The previous __eq__ defined a nested __lt__ and then fell off the
        # end, so every Profile compared unequal to every other including
        # itself.
        if not isinstance(other, Profile):
            return NotImplemented
        return self._sort_key() == other._sort_key()

    def __hash__(self):
        return hash(self._sort_key())
