"""
Calculates and stores basic thermodynamic parameters
"""
import datetime
from metpy import calc
import profiles.utils as utils
import numpy as np
import netCDF4
import os
import datetime as dt
from copy import deepcopy, copy


class Thermo_Profile():
    """ Contains data from one file.

    :var np.array<Quantity> temp: QC'd and averaged temperature
    :var np.array<Quantity> mixing_ratio: calculated mixing ratio
    :var np.array<Quantity> theta: calculated potential temperature
    :var np.array<Quantity> T_d: calculated dewpoint temperature
    :var np.array<Quantity> q: calculated mixing ratio
    :var np.array<Quantity> rh: QC'd and averaged relative humidity
    :var np.array<Quantity> pres: QC'd pressure
    :var np.array<Quantity> alt: altitude
    :var np.array<Datetime> gridded_times: times at which processed data exists
    :var Quantity resolution: vertical resolution in units of time,
           altitude, or pressure to which the data is calculated
    """

    def __init__(self, *args, **kwargs):
        if len([*args]) > 0:
            self._init2(*args, **kwargs)

    def _init2(self, temp_dict, resolution, file_path=None,
               gridded_times=None, gridded_base=None, indices=(None, None),
               ascent=True, units=None, pos=None, meta=None, nc_level='low'):
        """ Creates Thermo_Profile object from raw data at the specified
        resolution.

        :param dict temp_dict: A dictionary of the format \
           {"temp1":, "temp2":, ..., "tempj":, \
            "resi1":, "resi2":, ..., "resij", "time_temp":, \
            "rh1":, "rh2":, ..., "rhk":, "time_rh":, \
            "temp_rh1":, "temp_rh2":, ..., "temp_rhk":, \
            "pres":, "temp_pres":, "ground_temp_pres":, \
            "alt_pres":, "time_pres":, "serial_numbers":}, \
            which is returned by \
            Raw_Profile.thermo_data
        :param Quantity resoltion: vertical resolution in units of altitude \
           or pressure to which the data should be calculated
        :param str file_path: the path to the original data file WITHOUT the \
           suffix .nc or .json
        :param np.Array<Datetime> gridded_times: times at which data points \
           should be calculated
        :param np.Array<Quantity> gridded_base: base values corresponding to \
           gridded_times
        :param bool ascent: True if data should be processed for the ascending\
           leg of the flight, False if descending
        :param metpy.Units units: the unit registry created by Profile
        :param Meta meta: the parent Profile's Meta object
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """
        self._meta = meta
        self._units = units
        if ascent:
            self._ascent_filename_tag = "Ascending"
        else:
            self._ascent_filename_tag = "Descending"

        try:
            self._read_netCDF(file_path + "thermo_" +
                            str(resolution.magnitude) +
                            str(resolution.units) +
                            self._ascent_filename_tag + ".nc")
            return
        except Exception:
            self.resolution = resolution
            self.gridded_times = gridded_times
            self.rh = None
            self.pres = None
            self.temp = None
            self.alt = None
            self.rh_flags = None
            self.temp_flags = None
            self._units = units
            self._datadir = os.path.dirname(file_path + ".json")
        
        if not indices[0] is None:
            # trim profile
            selection_temp = np.where(np.array(temp_dict["time_temp"]) > indices[0],
                                      np.array(temp_dict["time_temp"]) < indices[1],
                                      False)
            selection_rh = np.where(np.array(temp_dict["time_rh"]) > indices[0],
                                    np.array(temp_dict["time_rh"]) < indices[1],
                                    False)
            selection_pres = np.where(np.array(temp_dict["time_pres"]) > indices[0],
                                      np.array(temp_dict["time_pres"]) < indices[1],
                                      False)

            for key in temp_dict.keys():
                if "time" not in key:
                    if ("temp" in key and "rh" not in key and "pres" not in key)\
                            or "resi" in key:
                        temp_dict[key] = temp_dict[key].magnitude[np.where(selection_temp)] * temp_dict[key].units
                    elif "pres" in key:
                        temp_dict[key] = temp_dict[key].magnitude[np.where(selection_pres)] * temp_dict[key].units
                    elif "rh" in key:
                        temp_dict[key] = temp_dict[key].magnitude[np.where(selection_rh)] * temp_dict[key].units
                else:
                     if ("temp" in key and "rh" not in key) \
                             or "resi" in key:
                         temp_dict[key] = np.array(temp_dict[key])[np.where(selection_temp)]
                     elif "pres" in key:
                         temp_dict[key] = np.array(temp_dict[key])[np.where(selection_pres)]
                     elif "rh" in key:
                        temp_dict[key] = np.array(temp_dict[key])[np.where(selection_rh)]

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
                                               serial_numbers["imet"+str(i+1)])
        # End if-else blocks

        rh_raw = []
        # Fill rh_raw
        for key in temp_dict.keys():
            # Ensure only humidity is processed here
            if "rh" in key and "temp" not in key and "time" not in key:
                rh_raw.append(temp_dict[key].magnitude)
        for i in range(len(rh_raw)):
            rh_raw[i] = utils.rh_calib(rh_raw[i], serial_numbers["rh"+str(i+1)])
        alts = np.array(temp_dict["alt_pres"].magnitude)\
            * temp_dict["alt_pres"].units
        pres = np.array(temp_dict["pres"].magnitude)\
            * temp_dict["pres"].units

        time_rh = temp_dict["time_rh"]
        time_pres = temp_dict["time_pres"]
        time_temp = temp_dict["time_temp"]
        # Determine bad sensors
        self.rh_flags = utils.qc(rh_raw, 0.4, 0.2)  # TODO read these from file

        # Remove bad sensors
        for flags_ind in range(len(self.rh_flags)):
            if self.rh_flags[flags_ind] != 0:
                rh_raw[flags_ind] = np.full(len(rh_raw[flags_ind]), np.NaN)

        # Average the sensors
        for i in range(len(rh_raw[0])):
            rh.append(np.nanmean([rh_raw[a][i] for a in
                                  range(len(rh_raw))]))

        rh = np.array(rh) * units.percent
        # Determine which sensors are "bad"
        self.temp_flags = utils.qc(temp_raw, 0.25, 0.1)

        # Remove bad sensors
        temp_ind = 0  # track index in temp_raw  after items are removed.
        for flags_ind in range(len(self.temp_flags)):
            if self.temp_flags[flags_ind] != 0:

                if self.temp_flags[flags_ind] != 4:  # if this isn't the missing data flag...
                    print("Temperature sensor", flags_ind + 1, "removed")

                temp_raw[flags_ind] = np.full(len(temp_raw[flags_ind]), np.nan)

        # Average the sensors
        for i in range(len(temp_raw[0])):

            temp.append(np.nanmean([temp_raw[a][i] for a in
                                   range(len(temp_raw))]))

        temp = np.array(temp) * units.kelvin

        #
        # Regrid to match times specified by Profile
        #

        # grid alt and pres
        if (self.resolution.dimensionality ==
                self._units.get_dimensionality('m')):
            self.alt = gridded_base
            self.pres = utils.regrid_data(data=pres, data_times=time_pres,
                                          gridded_times=self.gridded_times,
                                          units=self._units)
        elif (self.resolution.dimensionality ==
              self._units.get_dimensionality('Pa')):
            self.pres = gridded_base
            self.alt = utils.regrid_data(data=alts, data_times=time_pres,
                                         gridded_times=self.gridded_times,
                                         units=self._units)

        # grid RH
        self.rh = utils.regrid_data(data=rh, data_times=time_rh,
                                    gridded_times=self.gridded_times,
                                    units=self._units)

        # grid temp
        self.temp = utils.regrid_data(data=temp, data_times=time_temp,
                                      gridded_times=self.gridded_times,
                                      units=self._units)

        if pos is not None:
            # grid lat
            self.lat = utils.regrid_data(data=pos['lat'], data_times=pos['time'],
                                         gridded_times=self.gridded_times,
                                         units=self._units)

            self.lon = utils.regrid_data(data=pos['lon'], data_times=pos['time'],
                                         gridded_times=self.gridded_times,
                                         units=self._units)

        else:
            self.lat = np.full_like(self.gridded_times, -999.) * self._units('degrees')
            self.lon = np.full_like(self.gridded_times, -999.) * self._units('degrees')

        """ 
        TB --  I don't think this is needed... 
        and it causes unexpected down stream issues...
        may bite me later... 
        If causing issues later, do 'minlen+1'
        """
        # minlen = min(len(self.alt), len(self.gridded_times), len(self.rh),
        #              len(self.pres), len(self.temp))
        # self.pres = self.pres[0:minlen]
        # self.temp = self.temp[0:minlen]
        # self.rh = self.rh[0:minlen]
        # self.alt = self.alt[0:minlen]
        # self.gridded_times = self.gridded_times[0:minlen]

        # Calculate mixing ratio
        self.mixing_ratio = calc.mixing_ratio_from_relative_humidity(
                            self.pres, self.temp,
                            np.divide(self.rh.magnitude, 100))

        self.theta = calc.potential_temperature(self.pres, self.temp)
        self.T_d = calc.dewpoint_from_relative_humidity(self.temp, self.rh)
        self.q = calc.specific_humidity_from_mixing_ratio(self.mixing_ratio) * \
                 units.gPerKg

        if nc_level in 'low':
            self._save_netCDF(file_path)

    def truncate_to(self, new_len):
        """ Shortens arrays to have no more than new_len data points

        :param new_len: The new, shorter length
        :return: None
        """

        self.pres = self.pres[:new_len]
        self.temp = self.temp[:new_len]
        self.rh = self.rh[:new_len]
        self.alt = self.alt[:new_len]
        self.T_d = self.T_d[:new_len]
        self.theta = self.theta[:new_len]
        self.mixing_ratio = self.mixing_ratio[:new_len]
        self.q = self.q[:new_len]
        self.gridded_times = self.gridded_times[:new_len]

    def _save_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """
        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self._meta is not None:
            file_name = str(self._meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                        str(self._meta.get("platform_id")) + "CMT" + \
                        "thermo_" + self._ascent_filename_tag + ".c1." + \
                        self._meta.get("timestamp").replace("_", ".") + ".cdf"
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
        time_var[:] = netCDF4.date2num(self.gridded_times,
                                       units='microseconds since \
                                       2010-01-01 00:00:00:00')
        time_var.units = 'microseconds since 2010-01-01 00:00:00:00'

        # Do base_time and time_offset like ARM
        bt = abs((self.gridded_times[0] - dt.datetime(1970, 1, 1)).total_seconds())
        bt_var = main_file.createVariable('base_time', 'i8')
        bt_var.setncattr('long_name', 'Base time in Epoch')
        bt_var.setncattr('ancillary_variables', 'time_offset')
        bt_var.setncattr('units', 'seconds since 1970-01-01 00:00:00 UTC')
        bt_var[:] = bt

        to = netCDF4.date2num(self.gridded_times,
                              units=f'seconds since {self.gridded_times[0]:%Y-%m-%d %H:%M:%S UTC}')
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

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        :param string file_path: file name
        """
        file_name = str(self._meta.get("location")) + str(self.resolution.magnitude) + \
                    str(self._meta.get("platform_id")) + "CMT" + \
                    "thermo_" + self._ascent_filename_tag + ".c1." + \
                    self._meta.get("date_utc").replace("_", ".") + ".cdf"
        file_name = os.path.join(os.path.dirname(file_path), file_name)
        main_file = netCDF4.Dataset(file_path, "r",
                                    format="NETCDF4", mmap=False)

        self.temp_flags = [main_file["temp_flags"].getncattr("sensor"+str(i+1)) for i in 
                           range(len(main_file["temp_flags"].ncattrs()))] 
        self.rh_flags = [main_file["rh_flags"].getncattr("sensor"+str(i+1)) for i in 
                           range(len(main_file["rh_flags"].ncattrs()))]
        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        self.alt = np.array(main_file.variables["alt"])* \
            self._units.parse_expression(main_file.variables["alt"].units)
        self.pres = np.array(main_file.variables["pres"]) * \
            self._units.parse_expression(main_file.variables["pres"].units)
        self.rh = np.array(main_file.variables["rh"]) * \
            self._units.parse_expression(main_file.variables["rh"].units)
        self.temp = np.array(main_file.variables["temp"]) * \
            self._units.parse_expression(main_file.variables["temp"].units)
        self.mixing_ratio = np.array(main_file.variables["mr"]) * \
            self._units.parse_expression(main_file.variables["mr"].units)
        self.theta = np.array(main_file.variables["theta"]) * \
            self._units.parse_expression(main_file.variables["theta"].units)
        self.T_d = np.array(main_file.variables["Td"]) * \
            self._units.parse_expression(main_file.variables["Td"].units)
        self.q = np.array(main_file.variables["q"]) * \
            self._units.parse_expression(main_file.variables["q"].units)
        base_time = dt.datetime(2010, 1, 1, 0, 0, 0, 0)
        self.gridded_times = []
        for i in range(len(main_file.variables["time"][:])):
            self.gridded_times.append(base_time + dt.timedelta(microseconds=
                                                               int(main_file.variables
                                                                   ["time"][i])))
            # Hardcoded to microseconds since 2010-1-1
        main_file.close()

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            if key in "_units":
                setattr(result, key, copy(value))
            else:
                setattr(result, key, deepcopy(value, memo))
        return result

    def __str__(self):
        to_return = "\t\tThermo_Profile" \
                    + "\n\t\t\talt:          " + str(type(self.alt)) \
                    + "\n\t\t\tpres:         " + str(type(self.pres)) \
                    + "\n\t\t\trh:           " + str(type(self.rh)) \
                    + "\n\t\t\ttemp:         " + str(type(self.temp)) \
                    + "\n\t\t\tmixing_ratio: " + str(type(self.mixing_ratio))
        return to_return
