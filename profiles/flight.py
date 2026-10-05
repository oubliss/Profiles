"""
One flight log, at the rate each sensor reported.

Holds the data as parsed - native sample rates, several independent time
coordinates, no gridding and no profile selection. Those belong to Profile.
"""
import netCDF4
import numpy as np
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
from profiles.unit_registry import units  # shared pint registry
import profiles.utils as utils
import profiles.readers as readers
import profiles.parsing as parsing
from profiles.io import naming
import profiles.calibration as calibration
from profiles import Coef_Manager
from profiles.retrievals import wind as wind_retrieval
from profiles import schema
import os

from .utils import event_IDs

def _as_datetimes(values):
    """ numpy datetime64 back to Python datetimes.

    Everything downstream - netCDF4.date2num, list().index() on a time
    array, strftime, comparison against datetimes - expects datetime
    objects. Samples stamped before the GPS clock was set are dropped by
    profiles.parsing, so there is no NaT to represent; one reaching here
    would otherwise become NaN and fail much later, in the a0 writer or a
    time comparison.

    :param values: datetime64 array or sequence
    :rtype: list
    :raises ValueError: a value is NaT
    """
    stamps = np.asarray(values, dtype='datetime64[us]')
    if np.isnat(stamps).any():
        raise ValueError('time array contains NaT; samples without a valid '
                         'timestamp must be dropped, not carried')
    return [stamp.astype(object) for stamp in stamps]


def _times(dataset):
    """Time coordinate of a group Dataset as Python datetimes."""
    dimension = f'{dataset.attrs["group"]}_time'
    return _as_datetimes(dataset[dimension].values)


def _arrays_equal(first, second):
    """ np.array_equal that also works on pint Quantities (which do not
    implement it): units must match and values must be equal, NaN included
    where the same slots are NaN in both."""
    first_unit = getattr(first, 'units', None)
    second_unit = getattr(second, 'units', None)
    if first_unit != second_unit:
        return False
    if first_unit is not None:
        first, second = first.magnitude, second.magnitude
    return np.array_equal(first, second, equal_nan=_is_float(first))


def _is_float(array):
    return np.asarray(array).dtype.kind == 'f'


def _quantity(dataset, name):
    """A Dataset variable as a pint Quantity, using its recorded units."""
    variable = dataset[name]
    unit = variable.attrs.get('units')
    values = np.asarray(variable.values)
    return values * units.parse_expression(unit) if unit else values


#: Time reference of every a0 time coordinate.
_A0_EPOCH = "microseconds since 2010-01-01 00:00:00:00"


def _nc_values(group, name):
    """A float variable of an a0 group as a plain array, masked values NaN.

    A netCDF Variable is incompatible with pint, and a masked array would
    carry its fill value into the arithmetic.
    """
    return np.asarray(np.ma.filled(group.variables[name][:], np.nan))


def _nc_times(group):
    """The ``time`` variable of an a0 group as Python datetimes."""
    return list(netCDF4.num2date(
        group.variables["time"][:], units=_A0_EPOCH,
        only_use_cftime_datetimes=False, only_use_python_datetimes=True))


def _slot(series, number):
    """ Sensor ``number``'s calibrated series, or None if the slot is empty.

    ``series`` is slot-aligned (see calibration._sensor_series): index
    number - 1 is sensor number, and an absent sensor is all NaN.
    """
    if series is None or number > len(series):
        return None
    values = np.asarray(series[number - 1], dtype=float)
    return values if np.isfinite(values).any() else None


def _read_calibrated(group, prefix):
    """ Slot-aligned calibrated series from an a0 group, or None if none.

    Sensors the writer skipped come back as NaN so that list index and
    sensor number agree, whichever sensors were present.
    """
    found = [number for number in range(1, schema.N_SENSORS + 1)
             if f'{prefix}{number}' in group.variables]
    if not found:
        return None
    length = len(group.variables[f'{prefix}{found[0]}'])
    return [_nc_values(group, f'{prefix}{number}')
            if number in found else np.full(length, np.nan)
            for number in range(1, schema.N_SENSORS + 1)]


def _nc_units(variable, default):
    """ Units recorded on an a0 variable, or ``default``.

    Earlier builds wrote "F" for temperatures that were really kelvin (RH
    sensors) or degrees Celsius (barometer, once mislabelled Fahrenheit);
    pint reads "F" as farad, so that label, or none, falls back to what the
    values actually are.
    """
    text = getattr(variable, 'units', None)
    if not text or text == 'F':
        text = default
    return units.parse_expression(text)


def _same_file(first, second):
    """True if two paths name the same file, whether or not it exists yet."""
    try:
        return os.path.samefile(first, second)
    except OSError:
        return (os.path.normcase(os.path.realpath(first))
                == os.path.normcase(os.path.realpath(second)))


class FlightLog():
    """ Contains data from one flight log.

    :var tuple temp: temperature as (Temp1, Resi1, Temp2, Resi2, ..., time)
    :var tuple rh: relative humidity as (rh1, T1, rh2, T2, ..., time)
    :var tuple pos: GPS data as (lat, lon, alt_MSL, alt_rel_home,
                                 alt_rel_orig, time)
    :var tuple pres: barometer data as (pres, temp, ground_temp, alt_AGL,
                                        time)
    :var tuple rotation: UAS position data as (VE, VN, VD, roll, pitch, yaw,
                                               time)
    :var bool dev: True if the data is from a developmental flight
    :var str baro: contains 4-letter code for the type of barom sensor used
    :var dict serial_numbers: Contains serial number or 0 for each sensor
    :var Meta meta: processes metadata
    """

    def __init__(self, file_path, dev=False, nc_level='low', metadata=None,
                 tail_number=None, calibration='auto', coefficient_dir=None,
                 baro_instance=1, ekf_core=0):
        """ Creates a FlightLog and reads in data in the appropriate
        format. *If meta_path_flight or meta_path_header includes scoop_id,
        the scoop_id constructor parameter will be overwritten*

        :param string file_path: file name
        :param bool dev: True if the flight was developmental, false otherwise
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        :param profiles.Meta metadata: Meta Object
        :param str calibration: 'auto' reads the log to decide whether
           the thermistors were calibrated onboard or need the
           coefficient tables; 'table' and 'onboard' force one.
        :param coefficient_dir: directory of coefficient tables for this
           flight's calibration source; resolved through profiles.config
           when omitted
        :param int baro_instance: which barometer to use, the ``I`` field
           of current firmware's BARO messages. 1 is the scoop barometer
           on the OK3DM/CopterSonde airframes; confirm it for others.
           Raises ValueError if the log has barometers but not this one.
           Old logs with BARO/BAR2 messages and no instance field always
           use BAR2 and ignore this setting.
        :param int ekf_core: which EKF core to use, the ``C`` field of
           XKF1 messages
        """
        self.meta = None
        if metadata is not None:
            self.meta = metadata
        self.temp = None
        self.rh = None
        self.pos = None
        self.pres = None
        self.rotation = None
        self.wind = None
        self.rpm = None
        self.imu = None
        # Parsed Datasets by group; only populated when read from a log.
        self.data = {}
        self.dev = dev
        self.baro = "BARO"
        self.baro_instance = None     # instance actually used, if known
        self._instances = {'pres': baro_instance, 'rotation': ekf_core}
        self.serial_numbers = {}
        self.file_path = file_path
        self.calib_temp = None
        self.calib_rh = None
        self.calib_speed = None
        self.calib_dir = None
        self.file_type = None
        self.tail_number = tail_number
        self._calibration_mode = calibration
        self._coefficient_dir = coefficient_dir
        self._calibration_source = None

        # Set dummy serial numbers - these will allow the file 
        # to be processed even if the JSON and checklist files 
        # do not provide serial numbers
        # IMET
        for sensor_number in np.add(range(4), 1):
            self.serial_numbers["imet" + str(sensor_number)] = 0
        # RH
        for sensor_number in np.add(range(4), 1):
            self.serial_numbers["rh" + str(sensor_number)] = 0
        # WIND
        self.serial_numbers["wind"] = 0

        extension = os.path.splitext(file_path)[1].lower()

        if extension == '.json':
            self.file_type = 'json'
            self._read_messages(readers.iter_json(file_path),
                                nc_level=nc_level)
        elif extension == '.bin':
            # Read straight from the log. This used to dump the .BIN to a
            # newline-delimited .json beside the original and re-read that;
            # for an 18 MB OK3DM flight the intermediate was 104 MB.
            self.file_type = 'bin'
            self._read_messages(readers.iter_bin(file_path),
                                nc_level=nc_level)
        elif extension == '.csv':
            self.file_type = 'csv'
            self._read_csv(file_path)
        elif extension in ('.nc', '.cdf'):
            self.file_type = 'nc'
            self._read_netCDF(file_path)
        else:
            raise ValueError(
                f'{file_path!r}: unrecognised extension {extension!r} '
                f'(expected .bin, .json, .nc, .cdf or .csv)')



    @property
    def calibration_source(self):
        """ Where this flight's calibrated values come from.

        Resolved once from the log's serial numbers, or forced by the
        ``calibration`` constructor argument. Assigning to it overrides
        both.

        :rtype: profiles.Coef_Manager.CalibrationSource
        """
        if self._calibration_source is None:
            self._calibration_source = Coef_Manager.source_for_flight(
                self.serial_numbers, directory=self._coefficient_dir,
                mode=self._calibration_mode)
        return self._calibration_source

    @calibration_source.setter
    def calibration_source(self, source):
        self._calibration_source = source

    @property
    def start_time(self):
        """ The flight's first valid timestamp, or None if it has none.

        Coefficient rows can be dated (a recalibrated sensor has one row per
        calibration), and the flight's start is the date they are selected
        by. Taken as the earliest first timestamp across the logged streams
        so it does not depend on which message type arrived first.

        :rtype: pandas.Timestamp or None
        """
        starts = []
        for stream in (self.temp, self.rh, self.pos, self.pres):
            if stream is None or len(stream) == 0:
                continue
            for value in stream[-1]:
                try:
                    stamp = pd.Timestamp(value)
                except (TypeError, ValueError):
                    continue
                if not pd.isna(stamp):
                    starts.append(stamp)
                    break
        return min(starts) if starts else None

    def resolve_tail_number(self):
        """ The airframe this flight was flown with.

        An explicitly given tail number wins. The registry is consulted only
        when none was given - the log's vehicle ID is a short number that can
        map to several aircraft, so it is the weaker evidence.

        :rtype: str
        :raises KeyError: no tail number given and the log has no copterID
        """
        if self.tail_number is not None:
            return self.tail_number
        return self.calibration_source.get_tail_n(
            self.serial_numbers['copterID'], when=self.start_time)

    def apply_thermo_coeffs(self):
        """ Calibrate every temperature and humidity sensor individually.

        Stores per-sensor arrays on self.calib_temp and self.calib_rh for
        the a0 file. The gridded profile calibrates from the same functions.
        """
        thermo_data = self.thermo_data()
        serial_numbers = thermo_data['serial_numbers']

        self.calib_temp = calibration.calibrate_temperature(
            thermo_data, serial_numbers, source=self.calibration_source,
            when=self.start_time)
        self.calib_rh = calibration.calibrate_humidity(
            thermo_data, serial_numbers)

    def apply_wind_coeffs(self, equation_name='E1'):
        """ Retrieve wind from airframe tilt over the whole flight.

        :param str equation_name: calibration equation for this airframe
        """
        wind_data = self.wind_data()

        try:
            tail_number = self.resolve_tail_number()
        except KeyError:
            print("No CopterID found. Please specify a tail number upon "
                  "FlightLog creation to calc winds (needed for CSV reads)")
            return

        coefficients = self.calibration_source.get_coefs(
            'Wind', tail_number, equation_name, when=self.start_time)
        self.calib_dir, self.calib_speed = wind_retrieval.retrieve(
            wind_data['roll'], wind_data['pitch'], wind_data['yaw'],
            coefficients, equation_name)

    def sensor_window(self, settle_seconds=5):
        """ The interval over which the scoop fan was aspirating the sensors.

        Outside it the sensors are not ventilated and their readings are not
        trustworthy, so leg detection is restricted to this window.

        :param float settle_seconds: delay added to the start to let the
           sensors reach equilibrium once the fan spins up
        :rtype: tuple
        :return: (start, stop) datetimes. Falls back to the full record when
           the fan flag is unusable.
        """
        thermo = self.thermo_data()
        times = thermo['time_temp']

        running = np.where(np.array(thermo['fan_flag']) > 0)[0]
        if running.size == 0:
            print("Error with the fan_flag.... Just using the start and end "
                  "times of the file...")
            return times[0], times[-1]

        return (times[running.min()] + timedelta(seconds=settle_seconds),
                times[running.max()])

    #: Smallest climb (or descent), in metres, that find_legs reports as a
    #: profile. See find_legs for where 20 comes from.
    MIN_LEG_EXTENT = 50.0

    def find_legs(self, legacy=False, confirm_bounds=False,
                  profile_start_height=None, min_extent=MIN_LEG_EXTENT,
                  ascent=True):
        """ Locate each vertical profile flown during this flight.

        This lived in two places - Profile_Set.add_all_profiles and
        Profile.__init__ - which had drifted apart: one passed the fan
        window and the other did not.

        :param bool legacy: use the pre-2021 altitude-threshold finder
           instead of peak detection
        :param bool confirm_bounds: plot what was found for a sanity check
        :param profile_start_height: starting height for the legacy finder
        :param float min_extent: peak detection only. Legs that do not climb
           (descend, for ``ascent=False``) at least this many metres are not
           reported. 0 or None disables the check. The legacy finder is not
           filtered: its caller already chose a start height.
        :param bool ascent: which direction ``min_extent`` is measured in -
           True for the climb, False for the descent, None for either
        :rtype: list[tuple]
        :return: (start, peak, end) times for each profile found

        Peak detection uses a one-metre prominence, so a small altitude
        wiggle at the top of a real profile used to be reported as a profile
        of its own (flight616: 08:07:51-08:07:53, 1.1 m). The same flight's
        descent also starts with a 27 m dip below that peak, which is a
        "profile" when processing descents. The default extent of 50 m
        clears both with margin while staying below every real profile in
        the logs checked: flight616 and the 2026 OK3DM logs that parse have
        profiles of 76 m to 1420 m. The only other legs found there were
        13-17 m hover/ground-test wiggles (flights 2858 and 2861), which
        grid to one or two 10 m levels, so nothing of value is lost.
        """
        pos = self.pos_data()

        if legacy:
            return utils.identify_profile(
                pos['alt_MSL'], pos['time'], confirm_bounds, to_return=[],
                profile_start_height=profile_start_height)

        legs = utils.identify_profile_peaks(
            pos['alt_MSL'].magnitude, pos['time'],
            window=self.sensor_window(), confirm_bounds=confirm_bounds)

        return utils.filter_legs_by_extent(
            legs, pos['alt_MSL'].to('m').magnitude, pos['time'],
            min_extent, ascent=ascent)

    def pos_data(self):
        """ Gets data needed by the Profile constructor.

        rtype: dict
        return: {"lat":, "lon":, "alt_MSL":, "time":, "units"}
        """

        to_return = {}

        to_return["lat"] = self.pos[0]
        to_return["lon"] = self.pos[1]
        to_return["alt_MSL"] = self.pos[2]
        to_return["time"] = self.pos[-1]
        to_return["units"] = units

        return to_return

    def thermo_data(self):
        """ Gets data needed by the Thermo_Profile constructor.

        rtype: dict
        return: {"temp1":, "temp2":, ..., "tempj":, \
                 "resi1":, "resi2":, ..., "resij": , "time_temp": \
                 "rh1":, "rh2":, ..., "rhk":, "time_rh":, \
                 "temp_rh1":, "temp_rh2":, ..., "temp_rhk":, \
                 "pres":, "temp_pres":, "ground_temp_pres":, \
                 "alt_pres":, "time_pres"}
        """
        to_return = {}
        for sensor_number in [a + 1 for a in
                              range(int(len(self.temp) / 2) - 1)]:
            to_return["temp" + str(sensor_number)] \
                = self.temp[sensor_number*2 - 2]
            to_return["resi" + str(sensor_number)] \
                = self.temp[sensor_number*2 - 1]

        to_return["fan_flag"] = self.temp[-2]
        to_return["time_temp"] = np.array(self.temp[-1])

        for sensor_number in [a + 1 for a in range((len(self.rh)-1) // 2)]:
            to_return["rh" + str(sensor_number)] = self.rh[sensor_number * 2 - 2]
            to_return["temp_rh" + str(sensor_number)] = self.rh[sensor_number*2 - 1]

        to_return["time_rh"] = self.rh[-1]

        to_return["pres"] = self.pres[0]
        to_return["temp_pres"] = self.pres[1]
        to_return["ground_temp_pres"] = self.pres[2]
        to_return["alt_pres"] = self.pres[3]
        to_return["time_pres"] = self.pres[-1]

        to_return["serial_numbers"] = self.serial_numbers
        return to_return

    def wind_data(self):
        """ Gets data needed by the Wind_Profile constructor.

        rtype: list
        return: {"speed_east":, "speed_north":, "speed_down":, \
                 "roll":, "pitch":, "yaw":, "time":}
        """
        to_return = {}

        # rotation is formatted: (VE, VN, VD, roll, pitch, yaw, time)
        to_return["speed_east"] = self.rotation[0]
        to_return["speed_north"] = self.rotation[1]
        to_return["speed_down"] = self.rotation[2]
        to_return["roll"] = self.rotation[3]  # These are in radians
        to_return["pitch"] = self.rotation[4]
        to_return["yaw"] = self.rotation[5]
        to_return["pos_n"] = self.rotation[6]
        to_return["pos_e"] = self.rotation[7]
        to_return["pos_d"] = self.rotation[8]
        to_return["time"] = self.rotation[-1]

        to_return["alt"] = self.pres[3]
        to_return["pres"] = self.pres[0]
        to_return['time_pres'] = self.pres[-1]


        to_return["serial_numbers"] = self.serial_numbers

        return to_return

    def _read_csv(self, file_path):

        csv_header = ["date", 'lat', 'lon', 'alt', 'pressure',
                      'roll', 'pitch', 'yaw',
                      'gyry', 'gyrx', 'gyrz',
                      'vx', 'vy', 'vz',
                      'accx', 'accy', 'accz',
                      'temp1', 'temp2', 'temp3', 'temp4', 'temp5',
                      'rh1', 'rh2', 'rh3', 'rh4', 'rh5',
                      'gpsboottime', 'pressureboottime', 'attitudeboottime', 'imetboottime',
                      'temp_r1', 'temp_r2', 'temp_r3', 'temp_r4', 'temp_r5']
        data_types = {}
        for name in csv_header:
            data_types[name] = float  # Everything should be floats
        data_types['date'] = str  # Except the date string

        sensor_names = {}

        # Read in the CSV
        data = pd.read_csv(file_path, names=csv_header, dtype=data_types)

        # Convert to a dict for ease of not working with a pandas dataframe
        data = data.to_dict('list')
        data_len = len(data['date'])

        ######
        # Format into temp list following the format used for the full logs
        ######
        sensor_names["IMET"] = {}
        temp_list = [[] for x in range(10)]  # Ignoring the 5th sensor spot in the csvs so we don't break things...
        sensor_numbers = np.add(range(int((len(temp_list)-2) / 2)), 1)

        for num in sensor_numbers:
            sensor_names["IMET"]["temp" + str(num)] = 2 * num - 2
            sensor_names["IMET"]["temp_r" + str(num)] = 2 * num - 1

        sensor_names["IMET"]["Fan"] = -2
        sensor_names["IMET"]["Time"] = -1

        # Read fields into temp_list
        for key, value in sensor_names["IMET"].items():
            try:
                if 'Time' in key:
                    temp_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data['date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                else:
                    temp_list[value] = data[key]

            except KeyError:
                # Any expected variable that was not logged will show
                # as a list of NaN.
                temp_list[value] += [np.nan for foo in range(data_len)]

        ######
        # Format into rhum list following the format used for the full logs
        ######
        sensor_names["RHUM"] = {}
        rh_list = [[] for x in range(10)]  # Ignoring the 5th sensor spot in the csvs so we don't break things...
        sensor_numbers = np.add(range(int((len(rh_list) - 2) / 2)), 1)

        for num in sensor_numbers:
            sensor_names["RHUM"]["rh" + str(num)] = 2 * num - 2
            sensor_names["RHUM"]["rh_t" + str(num)] = 2 * num - 1
        sensor_names["RHUM"]["Time"] = -1

        # Read fields into rh_list
        for key, value in sensor_names["RHUM"].items():
            try:
                if 'Time' in key:
                    rh_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data[
                        'date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                else:
                    rh_list[value] = data[key]

            except KeyError:
                # Any expected variable that was not logged will show
                # as a list of NaN.
                rh_list[value] += [np.nan for foo in range(data_len)]
            except IndexError:
                print("Error in Raw_Profile - 227")


        ######
        # Read in the GPS data
        ######
        pos_list = [[] for x in range(6)]

        sensor_names["POS"] = {}

        # Determine field names
        sensor_names["POS"]["Lat"] = 0
        sensor_names["POS"]["Lng"] = 1
        sensor_names["POS"]["Alt"] = 2
        sensor_names["POS"]["RelHomeAlt"] = 3
        sensor_names["POS"]["RelOriginAlt"] = 4
        sensor_names["POS"]["TimeUS"] = -1

        # Read fields into gps_list, including TimeUS
        for key, value in sensor_names["POS"].items():
            try:
                if 'Time' in key:
                    pos_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data['date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                else:
                    if "Rel" in key:
                        pos_list[value] = data['alt']
                    elif "Lat" in key:
                        pos_list[value] = data['lat']
                    elif "Lng" in key:
                        pos_list[value] = data['lon']
                    elif 'Alt' in key:
                        pos_list[value] = data['alt']
                    else:
                        pos_list[value] += [np.nan for foo in range(data_len)]
            except KeyError:
                pos_list[value] += [np.nan for foo in range(data_len)]

        ######
        # Read in Pressure data
        ######
        pres_list = [[] for x in range(5)]

        sensor_names[self.baro] = {}

        # Determine field names
        sensor_names[self.baro]["Press"] = 0
        sensor_names[self.baro]["Temp"] = 1
        sensor_names[self.baro]["GndTemp"] = 2
        sensor_names[self.baro]["Alt"] = 3
        sensor_names[self.baro]["TimeUS"] = 4

        # Read fields into gps_list, including TimeUS
        for key, value in sensor_names[self.baro].items():
            try:
                if 'Time' in key:
                    pres_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data['date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                else:
                    if "Press" in key:
                        pres_list[value] = data['pressure']
                    elif "Alt" in key:
                        pres_list[value] = data['alt']
                    else:
                        pres_list[value] += [np.nan for foo in range(data_len)]
            except KeyError:
                pres_list[value] += [np.nan for foo in range(data_len)]

        ######
        # Read in rotation data
        ######
        rotation_list = [[] for x in range(10)]

        sensor_names["NKF1"] = {}

        # Determine field names
        sensor_names["NKF1"]["VE"] = 0
        sensor_names["NKF1"]["VN"] = 1
        sensor_names["NKF1"]["VD"] = 2
        sensor_names["NKF1"]["Roll"] = 3
        sensor_names["NKF1"]["Pitch"] = 4
        sensor_names["NKF1"]["Yaw"] = 5
        sensor_names["NKF1"]["PN"] = 6
        sensor_names["NKF1"]["PE"] = 7
        sensor_names["NKF1"]["PD"] = 8
        sensor_names["NKF1"]["TimeUS"] = -1

        for key, value in sensor_names["NKF1"].items():
            try:
                if 'Time' in key:
                    rotation_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data['date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                elif 'VE' in key:
                    rotation_list[value] = data['vx']
                elif 'VN' in key:
                    rotation_list[value] = data['vy']
                elif 'VD' in key:
                    rotation_list[value] = data['vz']
                else:
                    rotation_list[value] = data[key.lower()]
            except KeyError:
                rotation_list[value] += [np.nan for foo in range(data_len)]

        ######
        # Read in IMU data
        ######

        imu_list = [[] for x in range(7)]

        sensor_names['IMU'] = {}

        sensor_names['IMU']['GyrX'] = 0
        sensor_names['IMU']['GyrY'] = 1
        sensor_names['IMU']['GyrZ'] = 2
        sensor_names['IMU']['AccX'] = 3
        sensor_names['IMU']['AccY'] = 4
        sensor_names['IMU']['AccZ'] = 5
        sensor_names['IMU']['TimeUS'] = -1

        for key, value in sensor_names['IMU'].items():
            try:
                if 'Time' in key:
                    imu_list[value] = [dt.strptime(d[:-2], '%Y-%m-%dT%H:%M:%S.%f') for d in data['date']]  # Need the [:-2] because the microsecond string is 7 chars long but python can only deode 6

                else:
                    imu_list[value] = data[key.lower()]

            except KeyError:
                imu_list[value] += [np.nan for forr in range(data_len)]


        #####
        # Add in the units
        #####

        # Temperature
        for i in range(int((len(temp_list) - 1) / 2)):
            try:
                temp_list[2*i] = np.array(temp_list[2*i]) * units.K
                temp_list[2*i + 1] = np.array(temp_list[2*i + 1]) * units.ohm
            except IndexError:
                # print("No data for sensor ", i + 1)
                continue

        # RH
        for i in range(len(rh_list) - 1):
            # rh
            if i % 2 == 0:
                rh_list[i] = np.array(rh_list[i]) * units.percent
            # temp
            else:
                rh_list[i] = np.array(rh_list[i]) * units.kelvin

        # POS
        ground_alt = 0  # Hard coded since we don't have MSL alt in these files for some reason...
        # Profiles have not yet been separated.
        pos_list[0] = np.array(pos_list[0]) * units.deg  # lat
        pos_list[1] = np.array(pos_list[1]) * units.deg  # lng
        pos_list[2] = np.array(pos_list[2]) * units.m  # alt
        pos_list[3] = np.array(pos_list[3]) * units.m  # relHomeAlt
        pos_list[4] = np.array(pos_list[4]) * units.m  # relOrigAlt

        # PRES
        pres_list[0] = np.array(pres_list[0]) * units.Pa
        pres_list[1] = np.array(pres_list[1]) * units.degC
        pres_list[2] = np.array(pres_list[2]) * units.degC
        pres_list[3] = np.array(np.add(pres_list[3], ground_alt)) * units.m

        # ROTATION
        for i in range(len(rotation_list) - 1):
            if i < 3:
                rotation_list[i] = np.array(rotation_list[i]) \
                                            * units.m / units.s
            elif i >= 6:
                rotation_list[i] = np.array(rotation_list[i]) \
                                   * units.m
            else:
                rotation_list[i] = np.rad2deg(np.array(rotation_list[i])) * units.deg

        # IMU List
        if imu_list is not None:
            self.imu = tuple(imu_list)

        #
        # Convert to tuple
        #
        self.temp = tuple(temp_list)
        self.rh = tuple(rh_list)
        self.pos = tuple(pos_list)
        self.pres = tuple(pres_list)
        self.rotation = tuple(rotation_list)


    def _read_messages(self, messages, nc_level='low'):
        """ Build the raw arrays from a stream of log messages.

        Called by the constructor for both .BIN and .json input; the reader
        for each normalises to the same {"meta": ..., "data": ...} shape.

        Parsing itself is driven by profiles.schema and produces one xarray
        Dataset per group, kept on self.data. The positional tuples
        (self.temp, self.rh, ...) are then derived from those Datasets so
        that existing consumers keep working; they are the legacy interface
        and will go away once those consumers read by name.

        :param iterable messages: normalised log messages, in file order
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """
        parsed = parsing.parse(messages, instances=self._instances)
        self.data = parsed['groups']

        self.serial_numbers.update(parsed['serial_numbers'])

        missing = [name for name in schema.REQUIRED_GROUPS
                   if name not in self.data]
        if missing:
            raise ValueError(
                f'{self.file_path!r} contains no usable data: no '
                + ', '.join(f'{name} messages' for name in missing)
                + ' were found.')

        self.baro = self.data['pres'].attrs['source_message_type']
        self.baro_instance = self.data['pres'].attrs.get('source_instance')

        self._build_legacy_tuples(parsed)

        if utils.writes_netcdf(nc_level):
            self.apply_thermo_coeffs()
            self.apply_wind_coeffs()
            self._save_netCDF(self.file_path)

    def _build_legacy_tuples(self, parsed):
        """ Populate the positional attributes from the parsed Datasets.

        Interleaving (temp1, resi1, temp2, resi2, ...) and the exact slot
        order are what thermo_data(), wind_data(), _save_netCDF() and
        is_equal() still expect.
        """
        temp = self.data['temp']
        self.temp = tuple(
            [_quantity(temp, f'{kind}{n}')
             for n in range(1, schema.N_SENSORS + 1)
             for kind in ('temp', 'resi')]
            + [np.asarray(temp['fan_flag'].values), _times(temp)])

        rh = self.data['rh']
        self.rh = tuple(
            [_quantity(rh, f'{kind}{n}')
             for n in range(1, schema.N_SENSORS + 1)
             for kind in ('rh', 'temp_rh')]
            + [_times(rh)])

        pos = self.data['pos']
        self.pos = tuple(
            [_quantity(pos, name) for name in
             ('lat', 'lon', 'alt_MSL', 'alt_rel_home', 'alt_rel_orig')]
            + [_times(pos)])

        # Barometric altitude is logged relative to the first MSL fix.
        ground_alt = float(pos['alt_MSL'].values[0])
        pres = self.data['pres']
        self.pres = (_quantity(pres, 'pres'),
                     _quantity(pres, 'temp'),
                     _quantity(pres, 'ground_temp'),
                     np.add(pres['alt'].values, ground_alt) * units.m,
                     _times(pres))

        rotation = self.data['rotation']
        self.rotation = tuple(
            [_quantity(rotation, name) for name in
             ('speed_east', 'speed_north', 'speed_down', 'roll', 'pitch',
              'yaw', 'pos_n', 'pos_e', 'pos_d')]
            + [_times(rotation)])

        # The groups below carry no units in the legacy interface.
        if 'wind' in self.data:
            wind = self.data['wind']
            self.wind = tuple(
                [np.asarray(wind[name].values) for name in
                 ('wdir', 'wspeed', 'R13', 'R23', 'R33')] + [_times(wind)])
        else:
            self.wind = None

        if 'imu' in self.data:
            imu = self.data['imu']
            self.imu = tuple(
                [np.asarray(imu[name].values) for name in
                 ('gyr_x', 'gyr_y', 'gyr_z', 'acc_x', 'acc_y', 'acc_z')]
                + [_times(imu)])

        events = parsed['events']
        self.events = (events[0], _as_datetimes(events[1])) if events else None

        texts = parsed['messages']
        self.messages = ((texts[0], _as_datetimes(texts[1])) if texts
                         else ([], []))

        rpm = parsed['rpm']
        if rpm is not None:
            self.rpm = tuple(rpm[:-1] + [_as_datetimes(rpm[-1])])

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        Rebuilds the same positional tuples _build_legacy_tuples produces
        from a log, so a FlightLog read from an a0 file is interchangeable
        with one read from the .BIN it was written from. Units come from
        each variable's own ``units`` attribute where the writer recorded
        them (temperatures and resistances); a0 files written by earlier
        1.4.0-dev builds labelled temperatures in K as ``volt<n>`` / "mV"
        and never stored resistances, and are read as temperatures in K
        with the resistance slots left NaN.

        Groups the writer makes optional (wind, events, rpm, calibrated
        values) are optional here too.

        :param string file_path: file name
        """

        main_file = netCDF4.Dataset(file_path, "r", format="NETCDF4",
                                    mmap=False)

        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        # SERIAL NUMBERS
        # Start from the dummy serials set in the constructor and overlay
        # whatever the file recorded: copterID is absent from logs without
        # SYSID_THISMAV, and 'wind' is carried alongside the sensor serials.
        if "serial_numbers" in main_file.groups:
            sn_grp = main_file.groups["serial_numbers"]
            for name in sn_grp.ncattrs():
                value = sn_grp.getncattr(name)
                self.serial_numbers[name] = (value.item()
                                             if isinstance(value, np.generic)
                                             else value)

        #
        # POSITION - this should be first
        #
        pos = main_file.groups["pos"]
        self.pos = (_nc_values(pos, "lat") * units.deg,
                    _nc_values(pos, "lng") * units.deg,
                    _nc_values(pos, "alt") * units.m,
                    _nc_values(pos, "alt_rel_home") * units.m,
                    _nc_values(pos, "alt_rel_orig") * units.m,
                    _nc_times(pos))

        #
        # TEMPERATURE
        #
        temp = main_file.groups["temp"]
        temp_times = _nc_times(temp)
        found = [int(name[4:]) for name in temp.variables
                 if name[:4] in ('temp', 'volt', 'resi') and name[4:].isdigit()]
        # Throughout the file it is assumed that there are 4 sensors of
        # each type; a sensor that did not report is a NaN slot so the
        # (temp, resi) pairing thermo_data() relies on is never shifted.
        temp_list = []
        for number in range(1, max(schema.N_SENSORS, max(found, default=0)) + 1):
            for kind, default in (('temp', 'kelvin'), ('resi', 'ohm')):
                name = f'{kind}{number}'
                legacy = kind == 'temp' and name not in temp.variables
                if legacy:
                    # Earlier 1.4.0-dev files: temperature in K under 'volt'.
                    name = f'volt{number}'
                if name in temp.variables:
                    unit = (units.kelvin if legacy
                            else _nc_units(temp.variables[name], default))
                    temp_list.append(_nc_values(temp, name) * unit)
                else:
                    temp_list.append(np.full(len(temp_times), np.nan)
                                     * units.parse_expression(default))
        temp_list.append(_nc_values(temp, "fan_flag"))
        temp_list.append(temp_times)
        self.temp = tuple(temp_list)

        #
        # RELATIVE HUMIDITY
        #
        rh = main_file.groups["rh"]
        rh_times = _nc_times(rh)
        found = [int(name[2:]) for name in rh.variables
                 if name[:2] == 'rh' and name[2:].isdigit()]
        rh_list = []
        for number in range(1, max(schema.N_SENSORS, max(found, default=0)) + 1):
            if f'rh{number}' in rh.variables:
                rh_list.append(_nc_values(rh, f'rh{number}') * units.percent)
                # Earlier builds labelled this "F" (farad to pint) though
                # the values are kelvin.
                rh_list.append(_nc_values(rh, f'temp{number}') * _nc_units(
                    rh.variables[f'temp{number}'], 'kelvin'))
            else:
                rh_list.append(np.full(len(rh_times), np.nan) * units.percent)
                rh_list.append(np.full(len(rh_times), np.nan) * units.kelvin)
        rh_list.append(rh_times)
        self.rh = tuple(rh_list)

        #
        # PRESSURE
        #
        pres = main_file.groups["pres"]
        self.pres = (_nc_values(pres, "pres") * units.Pa,
                     _nc_values(pres, "temp")
                     * _nc_units(pres.variables["temp"], 'degC'),
                     _nc_values(pres, "temp_ground")
                     * _nc_units(pres.variables["temp_ground"], 'degC'),
                     _nc_values(pres, "alt") * units.m,
                     _nc_times(pres))

        #
        # ROTATION
        #
        rotation = main_file.groups["rotation"]
        self.rotation = tuple(
            [_nc_values(rotation, name) * unit for name, unit in (
                ("VE", units.m / units.s), ("VN", units.m / units.s),
                ("VD", units.m / units.s), ("roll", units.deg),
                ("pitch", units.deg), ("yaw", units.deg),
                # Estimated distance from origin (N, E, Down components)
                ("PN", units.m), ("PE", units.m), ("PD", units.m))]
            + [_nc_times(rotation)])

        #
        # WIND - absent from logs that carried no WIND messages. Like the
        # tuple built from a log, it carries no units.
        #
        if "wind" in main_file.groups:
            wind = main_file.groups["wind"]
            self.wind = tuple(
                [_nc_values(wind, name) for name in
                 ("wdir", "wspd", "R13", "R23", "R33")]
                + [_nc_times(wind)])
        else:
            self.wind = None

        #
        # RPM - absent from logs with no motor messages
        #
        if "rpm" in main_file.groups:
            rpm = main_file.groups["rpm"]
            motors = []
            while f'rpm{len(motors) + 1}' in rpm.variables:
                motors.append(_nc_values(rpm, f'rpm{len(motors) + 1}'))
            self.rpm = tuple(motors + [_nc_times(rpm)])

        #
        # Calibrated values, where the writer had them
        #
        self.calib_temp = _read_calibrated(temp, 'calib_temp')
        self.calib_rh = _read_calibrated(rh, 'calib_rh')

        if "calib_wind" in main_file.groups:
            calib_wind = main_file.groups["calib_wind"]
            self.calib_speed = (_nc_values(calib_wind, "calib_wspd")
                                * units.m / units.s)
            self.calib_dir = _nc_values(calib_wind, "calib_wdir") * units.deg

        #
        # Other Attributes
        #
        self.baro = getattr(main_file, "baro", self.baro)
        instance = getattr(main_file, "baro_instance", None)
        self.baro_instance = None if instance is None else int(instance)
        # if main_file.dev contains the string "True", then this is a
        # developmental flight.
        self.dev = "True" in getattr(main_file, "dev", "")

        #
        # Events - absent from logs before August 2021
        #
        if "events" in main_file.groups:
            events = main_file.groups["events"]
            self.events = (_nc_values(events, "events"), _nc_times(events))
        else:
            self.events = None

        #
        # Messages
        #
        if "messages" in main_file.groups:
            messages = main_file.groups["messages"]
            texts = messages.variables["messages"][:]
            self.messages = ([str(text) for text in texts],
                             _nc_times(messages))
        else:
            self.messages = ([], [])

        main_file.close()

    def _save_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        Temperatures are written as ``temp<n>`` and resistances as
        ``resi<n>`` with the units the quantity actually carries, so
        _read_netCDF can rebuild the same (temp, resi, ...) tuple. They used
        to be written as ``volt<n>`` labelled "mV", which read back about
        100 K high once treated as millivolts and left no resistances for
        the pairs thermo_data() expects.

        :param string file_path: file name
        :raises ValueError: if the output path would be the input log
        """

        # Not .replace(): 'FLIGHT.JSON' or 'x.Bin' matched none of the
        # patterns and the "output" was the input, which was then truncated.
        file_name = naming.resolve(
            file_path, self.meta, 'a0',
            self.meta.get("timestamp").replace("_", ".") if self.meta else '',
            resolution=None,
            fallback=os.path.splitext(self.file_path)[0] + '.nc')

        if _same_file(file_name, self.file_path):
            raise ValueError(
                f'refusing to write the a0 file over its own input '
                f'{self.file_path!r}; pass a different output path')

        main_file = netCDF4.Dataset(file_name, "w",
                                    format="NETCDF4", mmap=False)

        # File NC compliant to version 1.8
        main_file.setncattr("Conventions", "NC-1.8")

        # SERIAL NUMBERS
        # Whatever the log provided: copterID needs SYSID_THISMAV, which
        # not every log carries, and the wind serial sits beside the sensors.
        sn_grp = main_file.createGroup("/serial_numbers")
        for name, value in self.serial_numbers.items():
            if value is not None:
                sn_grp.setncattr(name, value)

        # EVENTS
        if self.events is not None:  # This maintains compatability for files pre August 2021
            events_grp = main_file.createGroup("/events")
            events_grp.createDimension("event_time", None)
            new_var = events_grp.createVariable("time", "f8", ("event_time",))
            new_var[:] = netCDF4.date2num(self.events[-1], units=_A0_EPOCH)
            new_var.units = _A0_EPOCH

            new_var = events_grp.createVariable("events", "f8", ("event_time",))
            new_var[:] = self.events[0]
            new_var.comment1 = "Event IDs last updated Aug 2021"
            new_var.units = event_IDs

        # MESSAGES
        message_grp = main_file.createGroup("/messages")
        message_grp.createDimension("message_time", None)
        new_var = message_grp.createVariable("time", "f8", ("message_time",))
        new_var.units = _A0_EPOCH
        # Variable-length strings: a fixed '<U13' silently cut longer text.
        text_var = message_grp.createVariable("messages", str, ("message_time",))
        text_var.units = 'string'
        if len(self.messages[0]):
            new_var[:] = netCDF4.date2num(self.messages[-1], units=_A0_EPOCH)
            text_var[:] = np.array([str(text) for text in self.messages[0]],
                                   dtype=object)

        # TEMP
        temp_grp = main_file.createGroup("/temp")
        temp_grp.createDimension("temp_time", None)
        temp_sensor_numbers = np.add(range(int((len(self.temp)-1)/2)), 1)
        for num in temp_sensor_numbers:
            for kind, slot in (("temp", 2*num-2), ("resi", 2*num-1)):
                values = self.temp[slot]
                if not hasattr(values, 'magnitude'):
                    # This sensor didn't report
                    continue
                new_var = temp_grp.createVariable(f"{kind}{num}", "f8",
                                                  ("temp_time",))
                new_var[:] = values.magnitude
                new_var.units = str(values.units)
        new_var = temp_grp.createVariable("time", "f8", ("temp_time",))
        new_var[:] = netCDF4.date2num(self.temp[-1], units=_A0_EPOCH)
        new_var.units = _A0_EPOCH

        # Scoop fan flag
        new_var = temp_grp.createVariable("fan_flag", "f8", ("temp_time",))
        new_var[:] = self.temp[-2]
        new_var.units = "unitless"
        new_var.comment1 = "0 -> Scoop fan off (no active sensor aspiration)"
        new_var.comment2 = "1 -> Scoop fan on (active sensor aspiration)"

        if self.calib_temp is not None:
            for num in temp_sensor_numbers:
                values = _slot(self.calib_temp, num)
                if values is None:
                    # This sensor didn't report
                    continue
                new_var = temp_grp.createVariable("calib_temp" + str(num), "f8",
                                                  ("temp_time",))
                new_var[:] = values
                new_var.units = 'K'

        # RH
        rh_grp = main_file.createGroup("/rh")
        rh_grp.createDimension("rh_time", None)
        rh_sensor_numbers = np.add(range(int((len(self.rh)-1)/2)), 1)
        for num in rh_sensor_numbers:
            if not (hasattr(self.rh[2*num-2], 'magnitude')
                    and hasattr(self.rh[2*num-1], 'magnitude')):
                # This sensor didn't report
                continue
            new_rh = rh_grp.createVariable("rh" + str(num),
                                           "f8", ("rh_time", ))
            new_temp = rh_grp.createVariable("temp" + str(num),
                                             "f8", ("rh_time", ))
            new_rh[:] = self.rh[2*num-2].magnitude
            new_temp[:] = self.rh[2*num-1].magnitude
            new_rh.units = "%"
            # Was "F", which pint reads as farad; these are kelvin.
            new_temp.units = str(self.rh[2*num-1].units)
        new_var = rh_grp.createVariable("time", "i8", ("rh_time",))
        new_var[:] = netCDF4.date2num(self.rh[-1], units=_A0_EPOCH)
        new_var.units = _A0_EPOCH

        if self.calib_rh is not None:
            for num in rh_sensor_numbers:
                values = _slot(self.calib_rh, num)
                if values is None:
                    # This sensor didn't report
                    continue
                new_var = rh_grp.createVariable("calib_rh" + str(num), "f8",
                                                ("rh_time",))
                new_var[:] = values
                new_var.units = '%'

        # POS
        pos_grp = main_file.createGroup("/pos")
        pos_grp.createDimension("pos_time", None)
        lat = pos_grp.createVariable("lat", "f8", ("pos_time", ))
        lng = pos_grp.createVariable("lng", "f8", ("pos_time", ))
        alt = pos_grp.createVariable("alt", "f8", ("pos_time", ))
        alt_rel_home = pos_grp.createVariable("alt_rel_home", "f8",
                                              ("pos_time", ))
        alt_rel_orig = pos_grp.createVariable("alt_rel_orig", "f8",
                                              ("pos_time", ))
        time = pos_grp.createVariable("time", "i8", ("pos_time",))

        lat[:] = self.pos[0].magnitude
        lng[:] = self.pos[1].magnitude
        alt[:] = self.pos[2].magnitude
        alt_rel_home[:] = self.pos[3].magnitude
        alt_rel_orig[:] = self.pos[4].magnitude
        time[:] = netCDF4.date2num(self.pos[-1], units=_A0_EPOCH)

        lat.units = "deg"
        lng.units = "deg"
        alt.units = "m MSL"
        alt_rel_home.units = "m"
        alt_rel_orig.units = "m"
        time.units = _A0_EPOCH

        # PRES
        pres_grp = main_file.createGroup("/pres")
        pres_grp.createDimension("pres_time", None)
        pres = pres_grp.createVariable("pres", "f8", ("pres_time", ))
        temp = pres_grp.createVariable("temp", "f8", ("pres_time", ))
        temp_gnd = pres_grp.createVariable("temp_ground", "f8",
                                           ("pres_time", ))
        alt = pres_grp.createVariable("alt", "f8", ("pres_time", ))
        time = pres_grp.createVariable("time", "i8", ("pres_time", ))

        pres[:] = self.pres[0].magnitude
        temp[:] = self.pres[1].magnitude
        temp_gnd[:] = self.pres[2].magnitude
        alt[:] = self.pres[3].magnitude
        time[:] = netCDF4.date2num(self.pres[-1], units=_A0_EPOCH)

        pres.units = "Pa"
        # Was "F", which pint reads as farad. BARO.Temp is degrees Celsius.
        temp.units = str(self.pres[1].units)
        temp_gnd.units = str(self.pres[2].units)
        alt.units = "m (MSL)"
        time.units = _A0_EPOCH

        # ROTATION
        rot_grp = main_file.createGroup("/rotation")
        rot_grp.createDimension("rot_time", None)
        ve = rot_grp.createVariable("VE", "f8", ("rot_time", ))
        vn = rot_grp.createVariable("VN", "f8", ("rot_time", ))
        vd = rot_grp.createVariable("VD", "f8", ("rot_time", ))
        roll = rot_grp.createVariable("roll", "f8", ("rot_time", ))
        pitch = rot_grp.createVariable("pitch", "f8", ("rot_time", ))
        yaw = rot_grp.createVariable("yaw", "f8", ("rot_time", ))
        pn = rot_grp.createVariable("PN", "f8", ("rot_time", ))
        pe = rot_grp.createVariable("PE", "f8", ("rot_time", ))
        pd = rot_grp.createVariable("PD", "f8", ("rot_time", ))
        time = rot_grp.createVariable("time", "i8", ("rot_time", ))

        ve[:] = self.rotation[0].magnitude
        vn[:] = self.rotation[1].magnitude
        vd[:] = self.rotation[2].magnitude
        roll[:] = self.rotation[3].magnitude
        pitch[:] = self.rotation[4].magnitude
        yaw[:] = self.rotation[5].magnitude
        pn[:] = self.rotation[6].magnitude
        pe[:] = self.rotation[7].magnitude
        pd[:] = self.rotation[8].magnitude
        time[:] = netCDF4.date2num(self.rotation[-1], units=_A0_EPOCH)

        ve.units = "m/s"
        vn.units = "m/s"
        vd.units = "m/s"
        roll.units = "deg"
        pitch.units = "deg"
        yaw.units = "deg"
        pn.units = 'meters'
        pe.units = 'meters'
        pd.units = 'meters'
        time.units = _A0_EPOCH

        # WIND
        if self.wind is not None:
            wind_grp = main_file.createGroup("/wind")
            wind_grp.createDimension('wind_time', None)

            time_var = wind_grp.createVariable("time", 'f8', ('wind_time',))
            wdir_var = wind_grp.createVariable("wdir", 'f8', ('wind_time',))
            wspd_var = wind_grp.createVariable("wspd", 'f8', ('wind_time',))
            r13_var  = wind_grp.createVariable("R13", 'f8', ('wind_time',))
            r23_var  = wind_grp.createVariable("R23", 'f8', ('wind_time',))
            r33_var  = wind_grp.createVariable("R33", 'f8', ('wind_time',))

            time_var[:] = netCDF4.date2num(self.wind[-1], units=_A0_EPOCH)
            wdir_var[:] = self.wind[0]
            wspd_var[:] = self.wind[1]
            r13_var[:] = self.wind[2]
            r23_var[:] = self.wind[3]
            r33_var[:] = self.wind[4]

            time_var.units = _A0_EPOCH
            wdir_var.units = "degrees"
            wspd_var.units = "m/s"
            r13_var.units = "None"
            r23_var.units = "None"
            r33_var.units = "None"

        if self.calib_speed is not None:
            wind_grp = main_file.createGroup("/calib_wind")
            wind_grp.createDimension('wind_time', None)

            time = wind_grp.createVariable("calib_time", 'i8', ('wind_time',))
            calib_wspd_var = wind_grp.createVariable("calib_wspd", 'f8', ('wind_time',))
            calib_wdir_var = wind_grp.createVariable("calib_wdir", 'f8', ('wind_time',))

            time[:] = netCDF4.date2num(self.rotation[-1], units=_A0_EPOCH)
            calib_wspd_var[:] = self.calib_speed.magnitude
            calib_wdir_var[:] = self.calib_dir.magnitude

            time.units = _A0_EPOCH
            calib_wspd_var.units = 'm/s'
            calib_wspd_var.comment = "NOTE: These values are only valid for ascending portions of the profile"
            calib_wdir_var.units = 'degrees'

        # RPM
        if self.rpm is not None:
            rpm_grp = main_file.createGroup("/rpm")
            rpm_grp.createDimension('rpm_time', None)

            time_var = rpm_grp.createVariable("time", 'f8', ('rpm_time',))
            time_var[:] = netCDF4.date2num(self.rpm[-1], units=_A0_EPOCH)
            time_var.units = _A0_EPOCH

            for motor_num in range(len(self.rpm) - 1):
                rpm_var = rpm_grp.createVariable(f"rpm{motor_num+1}", 'f8', ('rpm_time',))
                rpm_var[:] = self.rpm[motor_num]
                rpm_var.units = 'rpm'


        # Assign global attributes and close the file
        main_file.baro = self.baro
        if self.baro_instance is not None:
            main_file.baro_instance = int(self.baro_instance)
        main_file.dev = str(self.dev)

        main_file.close()

    def is_equal(self, other):
        """ Checks if two FlightLogs are the same.

        :param FlightLog other: profile with which to compare this one
        """
        # temps
        for i in range(len(self.temp)):
            if not _arrays_equal(self.temp[i], other.temp[i]):
                print("temp not equal at " + str(i))
                return False

        # rhs
        for i in range(len(self.rh)):
            if not _arrays_equal(self.rh[i], other.rh[i]):
                print("rh not equal at " + str(i))
                return False

        # pos
        for i in range(len(self.pos)):
            if not _arrays_equal(self.pos[i], other.pos[i]):
                print("pos not equal at " + str(i))
                return False

        # rotations
        for i in range(len(self.rotation)):
            if not _arrays_equal(self.rotation[i], other.rotation[i]):
                print("rotation not equal at " + str(i))
                return False

        # pres
        for i in range(len(self.pres)):
            if not _arrays_equal(self.pres[i], other.pres[i]):
                print("pres not equal at " + str(i))
                return False

        # misc
        if self.baro != other.baro:
            print("baro not equal")
            return False
        if self.dev != other.dev:
            print("dev not equal")
            return False

        return True

    def get_units(self):
        """
        :return: units
        """
        return units


#: Previous name for this class. Kept so that code written against 1.3.x
#: keeps importing, including the `from profiles.Raw_Profile import
#: Raw_Profile` form via the shim module of the same name.
Raw_Profile = FlightLog
