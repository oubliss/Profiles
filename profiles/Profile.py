
"""
Manages data from a single flight or profile
"""
from datetime import datetime, timedelta, timezone
import platform
from profiles.unit_registry import units
import profiles.utils as utils
import profiles.qc as qc
from profiles import io as profile_io
import profiles.calibration as calibration
import metpy.calc
import warnings
from profiles.retrievals import thermo as thermo_retrieval
from profiles.retrievals import wind as wind_retrieval
import profiles
import os
from profiles.flight import FlightLog
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


class EmptyProfileError(ValueError):
    """ A leg gridded to no levels, so there is no Profile to return. """


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
        if file_path is None:
            file_path = self._raw_profile.file_path

        if profile_start_height is not None:
            profile_start_height = profile_start_height * self._units.m
        if profile_num < 1:
            raise ValueError(
                f'profile_num is 1-based; got {profile_num}')

        if index_list is None:
            # One leg finder for every entry point. This used to be a
            # private copy with the raw fan window and no settle time, so a
            # directly built Profile and process_flights disagreed.
            # confirm_bounds is not forwarded: its default here is True,
            # and peak detection has never plotted from this constructor
            # (a blocking window per Profile). Plot via find_legs.
            index_list = self._raw_profile.find_legs(
                confirm_bounds=False, ascent=ascent,
                profile_start_height=(None if profile_start_height is None
                                      else profile_start_height.magnitude))

        if profile_num > len(index_list):
            # Raise rather than re-run the constructor: that re-parsed the
            # file, found the same legs and recursed until RecursionError.
            raise IndexError(
                f'profile_num={profile_num} but only {len(index_list)} '
                f'{"ascending" if ascent else "descending"} profile(s) were '
                f'found in {file_path}. If the file does contain more, try '
                f'a different profile_start_height or legacy_peak_id.')
        indices = index_list[profile_num - 1]

        if ascent:
            self.indices = (indices[0], indices[1])
        else:
            self.indices = (indices[1], indices[2])
        self._wind_computed = False
        self._thermo_computed = False
        self.qc_thresholds = dict(self.DEFAULT_QC_THRESHOLDS)
        #: Coefficient rows actually applied, recorded for the output.
        self.calibration_record = {}
        #: A bias correction applied on top, if any.
        self.bias_correction = None
        self.dev = dev  # TODO this is not used
        self.resolution = resolution * self._units.parse_expression(res_units)
        self.ascent = ascent
        self._ascent_filename_tag = 'ascent' if ascent else 'descent'

        # Match the extension itself: a substring test accepted any path
        # containing ".nc" (e.g. "run.nc_old/x.bin"), had no ".cdf" branch
        # although FlightLog reads it, and exited the interpreter on a miss.
        root, extension = os.path.splitext(file_path)
        if extension.lower() in ('.nc', '.cdf', '.json', '.bin', '.csv'):
            self.file_path = root
        else:
            raise ValueError(
                f'{file_path!r}: unrecognised extension {extension!r} '
                f'(expected .bin, .json, .nc, .cdf or .csv)')

        on_altitude = (self.resolution.dimensionality ==
                       self._units.get_dimensionality('m'))

        # profile_start_height is a height in metres. It used to take effect
        # only through process_flights (which converted it to base_start),
        # so profiles_from_flight and process_flights gridded the same
        # flight differently.
        self._base_start_given = base_start is not None
        if base_start is None and profile_start_height is not None:
            if on_altitude:
                base_start = profile_start_height
                self._base_start_given = True
            else:
                warnings.warn(
                    f'profile_start_height is a height in metres and '
                    f'cannot start a grid in {res_units}; ignoring it. '
                    f'Pass base_start as a pressure to fix the first '
                    f'edge.', UserWarning, stacklevel=2)

        if on_altitude:
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
            # The leg times are on the GPS clock; regrid_base maps them
            # onto the barometer's own.
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


        else:
            raise ValueError(
                f'res_units {res_units!r} is neither a length nor a '
                f'pressure')

        if len(self.time) == 0 or len(self.gridded_times) < 2:
            raise EmptyProfileError(
                f'the leg {self.indices[0]} to {self.indices[1]} grids to '
                f'no levels at {self.resolution}')

        # What the leg physically covers on the grid's vertical coordinate.
        # n_populated_levels needs it: a short leg on a long common grid
        # still has a full-length coordinate and, because the walk in
        # regrid_base hands successive levels successive samples, bins that
        # each hold a sample or two from a spot the leg never left.
        first = utils.nearest_index(base_time, self.indices[0])
        last = utils.nearest_index(base_time, self.indices[1])
        covered = base[first:last + 1].magnitude
        self._leg_range = (np.nanmin(covered) * base.units,
                           np.nanmax(covered) * base.units)

        # The vertical coordinate at bin centres, whichever coordinate was
        # chosen. This is what the gridded variables are aligned to.
        if (self.resolution.dimensionality ==
                self._units.get_dimensionality('m')):
            self.gridded_centers = self.alt
        else:
            self.gridded_centers = self.pres

        # The lowest edge of the grid, which is what another flight needs to
        # share these levels. An ascent is gridded bottom-up, a descent top
        # down; on a pressure grid the lowest edge is the highest pressure.
        edges = self.gridded_base.magnitude
        self._base_start = self.gridded_base[int(np.nanargmin(edges)
                                                 if on_altitude
                                                 else np.nanargmax(edges))]
        try:
            self.copter_id = self._raw_profile.serial_numbers['copterID']
        except KeyError:
            self.copter_id = -999
        try:
            # An explicit tail number wins; the registry is only consulted
            # when none was given, through the flight's own source and date.
            self.tail_number = self._raw_profile.resolve_tail_number()
        except Exception:
            # The copterID read above is still valid; only the registry
            # lookup failed.
            self.tail_number = self._raw_profile.tail_number


        self.__load_pos__()


    @property
    def n_populated_levels(self):
        """ Levels the leg actually reached and that hold data, as opposed
        to levels that only exist because the grid was extended.

        A short leg gridded onto a long common grid has a full-length
        vertical coordinate, so len(gridded_centers) says nothing about
        whether there is anything in the profile. A level counts when its
        centre is within half a bin of the vertical range the leg covered
        and at least one temperature sample was taken in its time bin.

        :rtype: int
        """
        centers = self.gridded_centers
        half = abs(self.resolution.to(centers.units).magnitude) / 2
        low, high = (q.to(centers.units).magnitude for q in self._leg_range)
        values = centers.magnitude
        reached = (values >= low - half) & (values <= high + half)

        stamps = np.asarray(self._thermo_data['time_temp'],
                            dtype='datetime64[us]')
        edges = np.asarray(self.gridded_times, dtype='datetime64[us]')
        in_bin = (np.searchsorted(stamps, edges[1:], side='right')
                  - np.searchsorted(stamps, edges[:-1], side='right'))
        return int(np.sum(reached & (in_bin > 0)))

    def __load_pos__(self):
        # Same (start, end] bins as every other gridded variable.
        def _grid(key, unit):
            return utils.regrid_data(
                data=self._pos[key], data_times=self._pos['time'],
                gridded_times=self.gridded_times,
                units=self._units).to(unit)

        self.lat = _grid('lat', units.deg)
        self.lon = _grid('lon', units.deg)
        self.alt_MSL = _grid('alt_MSL', units.m)

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

    #: Default QC thresholds: (max spread of sensor means, max spread of
    #: sensor standard deviations), in each variable's own units. Override
    #: per-Profile via qc_thresholds, or for a batch via ProcessingConfig.
    DEFAULT_QC_THRESHOLDS = {'temp': (0.25, 0.1), 'rh': (0.4, 0.2)}

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

        temp_raw = calibration.calibrate_temperature(
            data, serial_numbers, record=self.calibration_record,
            source=self._raw_profile.calibration_source,
            when=self._raw_profile.start_time)
        rh_raw = calibration.calibrate_humidity(data, serial_numbers)

        self.temp_flags = qc.qc(temp_raw, *self.qc_thresholds['temp'])
        self.rh_flags = qc.qc(rh_raw, *self.qc_thresholds['rh'])

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
        coefficients = self._raw_profile.calibration_source.get_coefs(
            'Wind', self.tail_number, equation_name,
            when=self._raw_profile.start_time)
        self.calibration_record['wind'] = coefficients

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
        file_name = profile_io.resolve(
            file_path, self.meta, 'c1',
            self.meta.get("timestamp").replace("_", ".") if self.meta else '',
            resolution=self.resolution.magnitude,
            tag='thermo_' + self._ascent_filename_tag)


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
        # q is held in kg/kg; the files report g/kg.
        q_var[:] = self.q.to('g/kg').magnitude
        q_var.units = 'g/kg'
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

        file_name = profile_io.resolve(
            file_path, self.meta, 'c1',
            self.meta.get("timestamp").replace("_", ".") if self.meta else '',
            resolution=self.resolution.magnitude,
            tag='wind_' + self._ascent_filename_tag)


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

    def _set_identity(self, handle):
        """ copter_id and tail_number, skipping whichever is unresolved
        (netCDF4 cannot store None)."""
        for name in ('copter_id', 'tail_number'):
            value = getattr(self, name, None)
            if value is not None:
                handle.setncattr(name, value)

    @staticmethod
    def _set_created(handle):
        handle.setncattr('datafile_created_on_date',
                         datetime.now(timezone.utc).replace(tzinfo=None)
                         .isoformat())
        handle.setncattr('datafile_created_on_machine', platform.node())

    def _c1_output_path(self, file_path, tag):
        """ Where a combined c1 file goes.

        An explicit .nc/.cdf path wins; otherwise the name comes from
        metadata, or - with none - from the input file's path, as the a0
        writer does. The two combined writers pass different tags so that
        calling both does not overwrite one with the other.

        :param str file_path: the caller's path, or None for the input's
        :param str tag: product tag, e.g. 'ascent' or 'cf_ascent'
        :rtype: str
        """
        if file_path is None:
            file_path = self.file_path
        resolution = self.resolution.magnitude
        return profile_io.resolve(
            file_path, self.meta, 'c1',
            self.time[0].strftime("%Y%m%d.%H%M%S"),
            resolution=resolution, tag=tag,
            fallback=f'{file_path}.c1.{resolution}.{tag}.nc')

    def save_netcdf(self, file_path=None, lookup_place=False):
        """ Write the combined c1 file.

        :param str file_path: an explicit .nc/.cdf path, or None to name
           the file from metadata (or the input file when there is none)
        :param bool lookup_place: also record a place name for the first
           fix as ``flight_location``. This asks a public web service
           (OpenStreetMap Nominatim), so it is off by default: a writer
           should not need a network, and a failed lookup only warns.
        """
        file_name = self._c1_output_path(file_path, self._ascent_filename_tag)

        if not (self._wind_computed or self._thermo_computed):
            print("No wind or thermo data to save; call compute_thermo() "
                  "and/or compute_wind() first")
            return

        main_file = netCDF4.Dataset(file_name, "w", format="NETCDF4", mmap=False)

        # Vital ncattrs
        main_file.setncattr("conventions", "NC-1.8")
        main_file.setncattr("processing_version", profiles.__version__)
        self._set_identity(main_file)

        if lookup_place:
            try:
                main_file.setncattr("flight_location",
                                    utils.get_place_from_lat_lon(
                                        self.lat[0].magnitude,
                                        self.lon[0].magnitude))
            except Exception as error:
                warnings.warn(f'place lookup failed, flight_location not '
                              f'written: {error}', RuntimeWarning,
                              stacklevel=2)

        self._set_created(main_file)
        main_file.setncattr('processing_level', 'c1')

        # Which coefficients, which thresholds, which table revision.
        for name, value in profile_io.provenance_attributes(self).items():
            main_file.setncattr(name, value)

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
            q_var[:] = self.q.to('g/kg').magnitude
            q_var.units = 'g/kg'
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

        profile_io.write_qc_variables(main_file, self)

        # Close the netCDF file
        main_file.close()


        return None

    def save_cfnetcdf(self, platform_name, terrain_elevation, file_path=None):
        """ Write the combined c1 file with CF-1.8 / WMO-CF attributes.

        :param str platform_name: platform identifier, e.g. the tail number
        :param terrain_elevation: site terrain elevation, in metres
        :param str file_path: an explicit .nc/.cdf path, or None to name
           the file from metadata (or the input file when there is none)
        """
        file_name = self._c1_output_path(
            file_path, 'cf_' + self._ascent_filename_tag)

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
        self._set_identity(main_file)
        self._set_created(main_file)

        for name, value in profile_io.provenance_attributes(self).items():
            main_file.setncattr(name, value)

        main_file.setncattr("Reference1", "Segales, A. R., B. R. Greene, T. M. Bell, W. Doyle, J. J. Martin, "
                                          "E. A. Pillar-Little, and P. B. Chilson, 2020: The CopterSonde: an insight"
                                          " into the development of a smart unmanned aircraft system for atmospheric "
                                          "boundary layer research. Atmospheric Measurement Techniques, 13, 2833–2848, "
                                          "https://doi.org/10.5194/amt-13-2833-2020.")
        main_file.setncattr("Reference2", "Bell, T. M., B. R. Greene, P. M. Klein, M. Carney, and "
                                          "P. B. Chilson, 2020: Confronting the boundary layer data gap: evaluating new "
                                          "and existing methodologies of probing the lower atmosphere. Atmospheric "
                                          "Measurement Techniques, 13, 3855–3872, https://doi.org/10.5194/amt-13-3855-2020.")

        # Create the dimensions
        main_file.createDimension("obs", None)

        # featureType trajectory needs a variable carrying cf_role; with one
        # flight per file it is a scalar.
        trajectory = main_file.createVariable('trajectory_id', str, ())
        trajectory.setncattr('cf_role', 'trajectory_id')
        trajectory.setncattr('long_name', 'flight identifier')
        trajectory[...] = f'{platform_name}_{self.gridded_times[0]}'

        # TIME
        # Be sure to use self.time instead of self.gridded_times since we want to store the mean time between two levels
        # of gridded times
        # f8: bin centres fall between whole seconds, which i8 truncated.
        time_var = main_file.createVariable("time", "f8", ("obs",))
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
        alt_var.setncattr('long_name', 'altitude above mean sea level')
        alt_var.setncattr('standard_name', 'altitude')
        alt_var.setncattr('positive', 'up')
        alt_var.setncattr('axis', 'Z')
        alt_var[:] = self.alt.magnitude

        # PRES
        pres_var = main_file.createVariable("pressure", "f8", ("obs",))
        pres_var.units = str(self.pres.units)
        pres_var.setncattr('long_name', 'air pressure')
        pres_var.setncattr('standard_name', 'air_pressure')
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
            # q_var[:] = self.q.to('g/kg').magnitude
            # q_var.units = 'g/kg'
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

        profile_io.write_qc_variables(main_file, self)
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
