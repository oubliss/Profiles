"""
Reads data file (JSON or netCDF) and stores the raw data
"""
import json
import netCDF4
import numpy as np
from datetime import datetime as dt
from metpy.units import units  # this is a pint UnitRegistry
import profiles.mavlogdump_Profiles as mavlogdump_Profiles
import profiles.utils as utils
from profiles.Meta import Meta
import pandas as pd
import os

units.define('percent = 0.01*count = %')
units.define('gPerKg = 0.001*count = g/Kg')


class Raw_Profile():
    """ Contains data from one file. Data is stored as a pandas DataFrame.

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

    def __init__(self, file_path, dev=False, scoop_id=None, nc_level='low', metadata=None):
        """ Creates a Raw_Profile object and reads in data in the appropriate
        format. *If meta_path_flight or meta_path_header includes scoop_id,
        the scoop_id constructor parameter will be overwritten*

        :param string file_path: file name
        :param bool dev: True if the flight was developmental, false otherwise
        :param char scoop_id: The set of sensors flown
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        :param profiles.Meta metadata: Meta Object
        """
        self.meta = None
        if metadata is not None:
            self.meta = metadata
        self.temp = None
        self.rh = None
        self.pos = None
        self.pres = None
        self.rotation = None
        self.dev = dev
        self.baro = "BARO"
        self.serial_numbers = {}
        self.file_path = file_path
        self.calib_temp = None
        self.calib_rh = None

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

        if "json" in file_path or "JSON" in file_path:
            if os.path.basename(file_path)[:-5] + ".nc" in \
               os.listdir(os.path.dirname(file_path)):
                self._read_netCDF(file_path[:-5] + ".nc")
            else:
                self._read_JSON(file_path, nc_level=nc_level)
        elif ".nc" in file_path or ".NC" in file_path:
            self._read_netCDF(file_path)
        elif ".bin" in file_path or ".BIN" in file_path:
            self.file_path = mavlogdump_Profiles.with_args(fmt="json",
                                                      file_name=file_path)
            self._read_JSON(self.file_path, nc_level=nc_level)

        if self.meta is not None:
            scoop_id = self.meta.get("scoop_id")


    def apply_thermo_coeffs(self):

        temp_dict = self.thermo_data()

        temp = []
        rh = []

        temp_raw = []  # List of lists, each containing data from a sensor

        # Fill temp_raw
        use_resistance = False
        use_temp = False
        for key in temp_dict.keys():
            if "resi" in key:
                use_resistance = True
                if use_temp:
                    use_temp = False
                    temp_raw = []
                temp_raw.append(temp_dict[key].magnitude)
            if "temp" in key and "_" not in key and not use_resistance:
                use_temp = True
                temp_raw.append(temp_dict[key].magnitude)

        # Process resistance if needed
        serial_numbers = temp_dict["serial_numbers"]
        if use_resistance:
            for i in range(len(temp_raw)):
                temp_raw[i] = utils.temp_calib(temp_raw[i],
                                               serial_numbers["imet" + str(i + 1)])

        rh_raw = []
        # Fill rh_raw
        for key in temp_dict.keys():
            # Ensure only humidity is processed here
            if "rh" in key and "temp" not in key and "time" not in key:
                rh_raw.append(temp_dict[key].magnitude)
        for i in range(len(rh_raw)):
            rh_raw[i] = utils.rh_calib(rh_raw[i], serial_numbers["rh" + str(i + 1)])


        self.calib_temp = temp_raw
        self.calib_rh = rh_raw

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
        to_return["alt"] = self.pres[3]
        to_return["pres"] = self.pres[0]
        to_return["time"] = self.rotation[6]

        to_return["serial_numbers"] = self.serial_numbers

        return to_return

    def _read_JSON(self, file_path, nc_level='low'):
        """ Reads data from a .JSON file. Called by the constructor.

        :param string file_path: file name
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """

        # Read the file into a list which pandas can normalize and read
        full_data = []
        for line in open(file_path, 'r'):
            full_data.append(json.loads(line))

        """
        Now full_data is a list of JSON element with 2 dictionaries each. If
        we refer to one JSON element as "tweet", the structure can be described
        as follows:

        tweet["meta"] contains "timestamp" and "type".

        tweet["data"] depends on tweet["meta"]["type"]. IMET, for example,
        could contain
            Temp1, float, 1, 286.2941589355469
            Temp2, float, 1, 286.27020263671875
            Temp3, float, 1, 286.0711364746094
            Temp4, float, 1, 0.0
            Time, int, 1, 60898080
            Volt1, float, 1, 4441.3125
            Volt2, float, 1, 4431.5625
            Volt3, float, 1, 4429.875
            Volt4, float, 1, 0.0

        Next, we iterate through full_data and identify line containing the
        types we want to keep. We then extract the data from each element
        using a different code for each type.
        """
        temp_list = None
        rh_list = None
        pos_list = None
        pres_list = None
        rotation_list = None
        # sensor_names will be dictionary of dictionaries formatted
        # {
        #     "valid_from": ,
        #     "IMET": {name: index, name: index, ...},
        #     "RHUM": {name: index, name: index, ...},
        #     ...
        # }
        sensor_names = {}

        for elem in full_data:

            if elem["meta"]["type"] == "PARM" and "SYSID_THISMAV" in elem["data"]["Name"]:

                self.serial_numbers['copterID'] = elem['data']['Value']

            if elem["meta"]["type"] == "PARM" and "USER_SENSORS" in elem["data"]["Name"]:
                index = int(elem['data']['Name'][-1])
                if index <= 4:
                    self.serial_numbers['imet' + str(index)] = int(elem['data']['Value'])
                elif index > 4 and index <= 8:
                    self.serial_numbers['rh' + str(index-4)] = int(elem['data']['Value'])

            if self.baro == "BARO" and elem["meta"]["type"] == "BAR2":
                # remove BARO structure and switch to using BAR2
                self.baro = "BAR2"
                pres_list = None
                sensor_names["BARO"] = None

            # IMET -> Temperature
            if elem["meta"]["type"] == "IMET":

                # First time only - setup temp_list
                if temp_list is None:

                    # Create array of lists with two lists per temperature
                    # sensor reported in the data file - one for temperature
                    # and one for resistance - plus one for times
                    temp_list = [[] for x in
                                 range(9)]
                    sensor_names["IMET"] = {}
                    # Determine field names
                    sensor_numbers = np.add(range(int((len(temp_list)-2) / 2)),
                                            1)
                    for num in sensor_numbers:
                        sensor_names["IMET"]["T"+str(num)] = 2*num - 2
                        sensor_names["IMET"]["R"+str(num)] = 2*num - 1
                    sensor_names["IMET"]["Time"] = -1

                # Read fields into temp_list, including Time
                for key, value in sensor_names["IMET"].items():
                    try:
                        if 'Time' in key:
                            time = dt.utcfromtimestamp(elem["meta"]
                                                       ["timestamp"])
                            if time.year < 2000:
                                raise KeyError("Time formatted incorrectly")
                            else:
                                temp_list[value].append(time)
                        else:
                            temp_list[value].append(elem["data"][key])
                    except KeyError:
                        # Any expected variable that was not logged will show
                        # as a list of NaN.
                        temp_list[value].append(np.nan)
                    except IndexError:
                        print("Error in Raw_Profile - 227")

            # Humidity
            elif elem["meta"]["type"] == "RHUM":

                # First time only - setup rh_list and temp_rh_list
                if rh_list is None:
                    # Create array of lists with one list per RH
                    # sensor reported in the data file, plus one for times
                    rh_list = [[] for x in range(sum(('H' in s and
                                                      'th' not in s)
                               for s in elem["data"].keys()) * 2 + 1)]

                    sensor_names["RHUM"] = {}
                    # Determine field names
                    sensor_numbers = np.add(range(int((len(rh_list)-1)/2)), 1)
                    for num in sensor_numbers:
                        sensor_names["RHUM"]["H"+str(num)] = 2*num - 2
                        sensor_names["RHUM"]["T"+str(num)] = 2*num - 1
                    sensor_names["RHUM"]["Time"] = -1

                # Read fields into rh_list, including Time
                for key, value in sensor_names["RHUM"].items():
                    try:
                        if 'Time' in key:
                            time = dt.utcfromtimestamp(elem["meta"]
                                                           ["timestamp"])
                            if time.year < 2000:
                                raise KeyError("Time formatted incorrectly")
                            else:
                                rh_list[value].append(time)
                        elif 'H' in key:
                            rh_list[value].append(elem["data"][key])
                        elif 'T' in key:
                            rh_list[value].append(elem["data"][key])
                    except KeyError:
                        rh_list[value].append(np.nan)

            # POS
            elif elem["meta"]["type"] == "POS":

                # First time only - setup gps_list
                if pos_list is None:
                    # Create array of lists with one list per [lat, lon, alt,
                    # time]
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
                            time = dt.utcfromtimestamp(elem["meta"]
                                                       ["timestamp"])
                            if time.year < 2000:
                                raise KeyError("Time formatted incorrectly")
                            else:
                                pos_list[value].append(time)
                        else:
                            if "Rel" in key:
                                pos_list[value].append(elem["data"][key])
                            elif 'Alt' in key:
                                pos_list[value].append(elem["data"][key])
                            else:
                                pos_list[value].append(elem["data"][key])
                    except KeyError:
                        pos_list[value].append(np.nan)

            # BARO or BAR2-> Pressure
            elif elem["meta"]["type"] == self.baro:
                # First time only - setup gps_list
                if pres_list is None:
                    # Create array of lists with one list per [pres, temp,
                    # ground_temp, alt, time]
                    pres_list = [[] for x in range(5)]

                    sensor_names[self.baro] = {}

                    # Determine field names
                    sensor_names[self.baro]["Press"] = 0
                    sensor_names[self.baro]["Temp"] = 1
                    sensor_names[self.baro]["GndTemp"] = 2
                    sensor_names[self.baro]["Alt"] = 3
                    sensor_names[self.baro]["TimeUS"] = 4

                # Read fields into pres_list, including TimeUS
                for key, value in sensor_names[self.baro].items():
                    try:
                        if 'Time' in key:
                            time = dt.utcfromtimestamp(elem["meta"]
                                                       ["timestamp"])
                            if time.year < 2000:
                                raise KeyError("Time formatted incorrectly")
                            else:
                                pres_list[value].append(time)
                        elif 'Alt' in key:
                            pres_list[value].append(elem["data"][key])
                        elif 'Temp' in key or 'GndTemp' in key:
                            pres_list[value].append(elem["data"][key])
                        elif 'Press' in key:
                            pres_list[value].append(elem["data"][key])
                        else:
                            print('undefined BARO key: ' + key)
                    except KeyError:
                        pres_list[value].append(np.nan)

            # NKF1 -> Rotation
            elif (elem["meta"]["type"] == "NKF1") | (elem["meta"]["type"] == "XKF1"):

                # First time only - setup gps_list
                if rotation_list is None:
                    # Create array of lists with one list per [ve, vn, vd,
                    # roll, pitch, yaw, time]
                    rotation_list = [[] for x in range(7)]

                    sensor_names["NKF1"] = {}

                    # Determine field names
                    sensor_names["NKF1"]["VE"] = 0
                    sensor_names["NKF1"]["VN"] = 1
                    sensor_names["NKF1"]["VD"] = 2
                    sensor_names["NKF1"]["Roll"] = 3
                    sensor_names["NKF1"]["Pitch"] = 4
                    sensor_names["NKF1"]["Yaw"] = 5
                    sensor_names["NKF1"]["TimeUS"] = -1

                # Read fields into rotation_list, including TimeUS
                for key, value in sensor_names["NKF1"].items():
                    try:
                        if 'Time' in key:
                            time = dt.utcfromtimestamp(elem["meta"]
                                                       ["timestamp"])
                            if time.year < 2000:
                                raise KeyError("Time formatted incorrectly")
                            else:
                                rotation_list[value].append(time)
                        elif 'VE' in key or 'VN' in key or 'VD' in key:
                            rotation_list[value].append(elem["data"][key])
                        else:  # Roll, pitch, yaw
                            rotation_list[value].append(elem["data"][key])
                    except KeyError:
                        rotation_list[value].append(np.nan)
        #
        # Add the units
        #

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
        ground_alt = pos_list[2][0]  # This is the first in the file.
        # Profiles have not yet been separated.
        pos_list[0] = np.array(pos_list[0]) * units.deg  # lat
        pos_list[1] = np.array(pos_list[1]) * units.deg  # lng
        pos_list[2] = np.array(pos_list[2]) * units.m  # alt
        pos_list[3] = np.array(pos_list[3]) * units.m  # relHomeAlt
        pos_list[4] = np.array(pos_list[4]) * units.m  # relOrigAlt

        # PRES
        pres_list[0] = np.array(pres_list[0]) * units.Pa
        pres_list[1] = np.array(pres_list[1]) * units.fahrenheit
        pres_list[2] = np.array(pres_list[2]) * units.fahrenheit
        pres_list[3] = np.array(np.add(pres_list[3], ground_alt)) * units.m

        # ROTATION
        for i in range(len(rotation_list) - 1):
            if i < 3:
                rotation_list[i] = np.array(rotation_list[i]) \
                                            * units.m / units.s
            else:
                rotation_list[i] = np.array(rotation_list[i]) * units.deg

        #
        # Convert to tuple
        #
        self.temp = tuple(temp_list)
        self.rh = tuple(rh_list)
        self.pos = tuple(pos_list)
        self.pres = tuple(pres_list)
        self.rotation = tuple(rotation_list)

        if nc_level in 'low':
            self.apply_thermo_coeffs()
            self._save_netCDF(file_path)

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        :param string file_path: file name
        """
        
        main_file = netCDF4.Dataset(file_path, "r", format="NETCDF4",
                                    mmap=False)

        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        # SERIAL NUMBERS
        self.serial_numbers = {}
        self.serial_numbers["copterID"] = main_file.groups["serial_numbers"].getncattr("copterID")
        for i in range(4):  # Throughout the file, it is assumed that there are 4 sensors of each type
            self.serial_numbers["rh" + str(i+1)] = main_file.groups["serial_numbers"].getncattr("rh" + str(i+1))
            self.serial_numbers["imet" + str(i+1)] = main_file.groups["serial_numbers"].getncattr("imet" + str(i+1))

        #
        # POSITION - this should be first
        #
        pos_list = []
        # lat
        pos_list.append(main_file.groups["pos"].variables["lat"])
        pos_list[0] = np.array(pos_list[0]) * units.deg
        # lng
        pos_list.append(main_file["pos"].variables["lng"])
        pos_list[1] = np.array(pos_list[1]) * units.deg
        # alt
        pos_list.append(main_file["pos"].variables["alt"])
        pos_list[2] = np.array(pos_list[2]) * units.m
        # altitude relative to home
        pos_list.append(main_file["pos"].variables["alt_rel_home"])
        pos_list[3] = np.array(pos_list[3]) * units.m
        # altitude relative to origin
        pos_list.append(main_file["pos"].variables["alt_rel_orig"])
        pos_list[4] = np.array(pos_list[4]) * units.m
        # time
        pos_list.append(netCDF4.num2date(main_file["pos"].
                                         variables["time"][:],
                                         units="microseconds since \
                                         2010-01-01 00:00:00:00"))
        # convert to tuple and save
        self.pos = tuple(pos_list)

        #
        # TEMPERATURE
        #
        temp_list = []
        i = 1
        while(True):
            try:
                temp_list.append(main_file["temp"].variables["volt" + str(i)])
                temp_list[i-1] = np.array(temp_list[i-1]) * units.mV
                i += 1
            except KeyError:
                break
        temp_list.append(netCDF4.num2date(main_file["temp"].
                                          variables["time"][:],
                                          units="microseconds since \
                                          2010-01-01 00:00:00:00"))
        self.temp = tuple(temp_list)

        #
        # RELATIVE HUMIDITY
        #
        rh_list = []
        i = 1
        while(True):
            try:
                rh_list.append(main_file["rh"].variables["rh" + str(i)])
                rh_list[-1] = np.array(rh_list[-1]) * units.percent

                rh_list.append(main_file["rh"].variables["temp" + str(i)])
                rh_list[-1] = np.array(rh_list[-1]) * units.F

                i += 1
            except KeyError:
                break
        rh_list.append(netCDF4.num2date(main_file["rh"].
                                        variables["time"][:],
                                        units="microseconds since \
                                        2010-01-01 00:00:00:00"))
        self.rh = tuple(rh_list)

        #
        # PRESSURE
        #
        pres_list = []
        # pres
        pres_list.append(main_file["pres"].variables["pres"])
        pres_list[0] = np.array(pres_list[0]) * units.Pa
        # temp
        pres_list.append(main_file["pres"].variables["temp"])
        pres_list[1] = np.array(pres_list[1]) * units.F
        # temp_gnd
        pres_list.append(main_file["pres"].variables["temp_ground"])
        pres_list[2] = np.array(pres_list[2]) * units.F
        # alt
        pres_list.append(main_file["pres"].variables["alt"])
        pres_list[3] = np.array(pres_list[3]) * units.m
        # time
        pres_list.append(netCDF4.num2date(main_file["pres"].
                                          variables["time"][:],
                                          units="microseconds since \
                                          2010-01-01 00:00:00:00"))
        self.pres = tuple(pres_list)

        #
        # ROTATION
        #
        rot_list = []
        # VE
        rot_list.append(main_file["rotation"].variables["VE"])
        rot_list[0] = np.array(rot_list[0]) * units.m / units.s
        # VN
        rot_list.append(main_file["rotation"].variables["VN"])
        rot_list[1] = np.array(rot_list[1]) * units.m / units.s
        # VD
        rot_list.append(main_file["rotation"].variables["VD"])
        rot_list[2] = np.array(rot_list[2]) * units.m / units.s
        # roll
        rot_list.append(main_file["rotation"].variables["roll"])
        rot_list[3] = np.array(rot_list[3]) * units.deg
        # pitch
        rot_list.append(main_file["rotation"].variables["pitch"])
        rot_list[4] = np.array(rot_list[4]) * units.deg
        # yaw
        rot_list.append(main_file["rotation"].variables["yaw"])
        rot_list[5] = np.array(rot_list[5]) * units.deg
        # time
        rot_list.append(netCDF4.num2date(main_file["rotation"].
                                         variables["time"][:],
                                         units="microseconds since \
                                         2010-01-01 00:00:00:00"))
        self.rotation = tuple(rot_list)

        #
        # Other Attributes
        #
        self.baro = main_file.baro
        self.dev = "True" in main_file.dev  # if main_file.dev contains the
        # string "True", then this is a developmental flight.

        main_file.close()

    def _save_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """
        file_name = str(self.meta.get("location")) + \
                    str(self.meta.get("platform_id")) + "CMT" + \
                    ".a0." + self.meta.get("timestamp").replace("_", ".") + ".cdf"
        file_name = os.path.join(os.path.dirname(file_path), file_name)
        main_file = netCDF4.Dataset(file_name, "w",
                                    format="NETCDF4", mmap=False)

        # File NC compliant to version 1.8
        main_file.setncattr("Conventions", "NC-1.8")

        # SERIAL NUMBERS
        sn_grp = main_file.createGroup("/serial_numbers")
        sn_grp.setncattr("copterID", self.serial_numbers["copterID"])
        for i in range(4):  # Throughout the file, it is assumed that there are 4 sensors of each type
            sn_grp.setncattr("rh" + str(i+1), self.serial_numbers['rh' + str(i+1)])
            sn_grp.setncattr("imet" + str(i+1), self.serial_numbers['imet' + str(i+1)])

        # TEMP
        temp_grp = main_file.createGroup("/temp")
        temp_grp.createDimension("temp_time", None)
        # temp_grp.base_time = date2num(self.temp[-1][0])
        temp_sensor_numbers = np.add(range(int((len(self.temp)-1)/2)), 1)
        for num in temp_sensor_numbers:
            new_var = temp_grp.createVariable("volt" + str(num), "f8",
                                              ("temp_time",))
            try:
                new_var[:] = self.temp[2*num-2].magnitude
            except AttributeError:
                # This sensor didn't report
                continue
            new_var.units = "mV"
        new_var = temp_grp.createVariable("time", "f8", ("temp_time",))
        new_var[:] = netCDF4.date2num(self.temp[-1],
                                      units="microseconds since \
                                      2010-01-01 00:00:00:00")
        new_var.units = "microseconds since 2010-01-01 00:00:00:00"

        if self.calib_temp is not None:
            for num in temp_sensor_numbers:
                new_var = temp_grp.createVariable("calib_temp" + str(num), "f8",
                                                  ("temp_time",))
                try:
                    new_var[:] = self.calib_temp[num-1]
                except Exception:
                    # This sensor didn't report
                    continue

                new_var.units = 'K'

        # RH
        rh_grp = main_file.createGroup("/rh")
        rh_grp.createDimension("rh_time", None)
        rh_sensor_numbers = np.add(range(int((len(self.rh)-1)/2)), 1)
        for num in rh_sensor_numbers:
            new_rh = rh_grp.createVariable("rh" + str(num),
                                           "f8", ("rh_time", ))
            new_temp = rh_grp.createVariable("temp" + str(num),
                                             "f8", ("rh_time", ))
            new_rh[:] = self.rh[2*num-2].magnitude
            new_temp[:] = self.rh[2*num-1].magnitude
            new_rh.units = "%"
            new_temp.units = "F"
        new_var = rh_grp.createVariable("time", "i8", ("rh_time",))
        new_var[:] = netCDF4.date2num(self.rh[-1],
                                      units="microseconds since \
                                      2010-01-01 00:00:00:00")
        new_var.units = "microseconds since 2010-01-01 00:00:00:00"

        if self.calib_rh is not None:
            for num in rh_sensor_numbers:
                new_var = rh_grp.createVariable("calib_rh" + str(num), "f8",
                                                  ("rh_time",))
                try:
                    new_var[:] = self.calib_rh[num-1]
                except Exception:
                    # This sensor didn't report
                    continue

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
        time[:] = netCDF4.date2num(self.pos[-1], units="microseconds since \
                                   2010-01-01 00:00:00:00")

        lat.units = "deg"
        lng.units = "deg"
        alt.units = "m MSL"
        alt_rel_home.units = "m"
        alt_rel_orig.units = "m"
        time.units = "microseconds since 2010-01-01 00:00:00:00"

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
        time[:] = netCDF4.date2num(self.pres[-1], units="microseconds since \
                                   2010-01-01 00:00:00:00")

        pres.units = "Pa"
        temp.units = "F"
        temp_gnd.units = "F"
        alt.units = "m (MSL)"
        time.units = "microseconds since 2010-01-01 00:00:00:00"

        # ROTATION
        rot_grp = main_file.createGroup("/rotation")
        rot_grp.createDimension("rot_time", None)
        ve = rot_grp.createVariable("VE", "f8", ("rot_time", ))
        vn = rot_grp.createVariable("VN", "f8", ("rot_time", ))
        vd = rot_grp.createVariable("VD", "f8", ("rot_time", ))
        roll = rot_grp.createVariable("roll", "f8", ("rot_time", ))
        pitch = rot_grp.createVariable("pitch", "f8", ("rot_time", ))
        yaw = rot_grp.createVariable("yaw", "f8", ("rot_time", ))
        time = rot_grp.createVariable("time", "i8", ("rot_time", ))

        ve[:] = self.rotation[0].magnitude
        vn[:] = self.rotation[1].magnitude
        vd[:] = self.rotation[2].magnitude
        roll[:] = self.rotation[3].magnitude
        pitch[:] = self.rotation[4].magnitude
        yaw[:] = self.rotation[5].magnitude
        time[:] = netCDF4.date2num(self.rotation[-1], units="microseconds \
                                   since 2010-01-01 00:00:00:00")

        ve.units = "m/s"
        vn.units = "m/s"
        vd.units = "m/s"
        roll.units = "deg"
        pitch.units = "deg"
        yaw.units = "deg"
        time.units = "microseconds since 2010-01-01 00:00:00:00"

        # Assign global attributes and close the file
        main_file.baro = self.baro
        main_file.dev = str(self.dev)

        main_file.close()

    def is_equal(self, other):
        """ Checks if two Raw_Profiles are the same.

        :param Raw_Profile other: profile with which to compare this one
        """
        # temps
        for i in range(len(self.temp)):
            if not np.array_equal(self.temp[i], other.temp[i]):
                print("temp not equal at " + str(i))
                return False

        # rhs
        for i in range(len(self.rh)):
            if not np.array_equal(self.rh[i], other.rh[i]):
                print("rh not equal at " + str(i))
                return False

        # pos
        for i in range(len(self.pos)):
            if not np.array_equal(self.pos[i], other.pos[i]):
                print("pos not equal at " + str(i))
                return False

        # rotations
        for i in range(len(self.rotation)):
            if not np.array_equal(self.rotation[i], other.rotation[i]):
                print("rotation not equal at " + str(i))
                return False

        # pres
        for i in range(len(self.pres)):
            if not np.array_equal(self.pres[i], other.pres[i]):
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
