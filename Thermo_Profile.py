"""
Calculates and stores basic thermodynamic parameters

Authors Brian Greene, Jessica Blunt, Tyler Bell, Gus Azevedo \n
Copyright University of Oklahoma Center for Autonomous Sensing and Sampling
2019

Component of Profiles v1.0.0
"""
from metpy import calc
import utils
import numpy as np
import netCDF4
import os
from copy import deepcopy, copy

# T1 [0x48]	T2 [0x49] T3 [0x4a]


class Thermo_Profile():
    """ Contains data from one file.

    :var np.array<Quantity> temp: QC'd and averaged temperature
    :var np.array<Quantity> mixing_ratio: calculated mixing ratio
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
               gridded_times=None, ascent=True, units=None, nc_level='low'):
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
        :param str filepath: the path to the original data file WITHOUT the \
           suffix .nc or .json
        :param np.Array<Datetime> gridded_times: times at which data points \
           should be calculated
        :param bool ascent: True if data should be processed for the ascending\
           leg of the flight, False if descending
        :param metpy.Units units: the unit registry created by Profile
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """

        self._units = units

        try:
            self._read_netCDF(file_path)
            return
        except Exception:
            self.resolution = resolution
            self.gridded_times = gridded_times
            self.rh = None
            self.pres = None
            self.temp = None
            self.alt = None
            self._units = units
            self._datadir = os.path.dirname(file_path + ".json")

        if ascent:
            self._ascent_filename_tag = "Ascending"
        else:
            self._ascent_filename_tag = "Descending"

        # If a .json or .bin was given, checks for corrosponding .nc
        if os.path.basename(file_path + "thermo_" +
                            str(self.resolution.magnitude) +
                            str(self.resolution.units) +
                            self._ascent_filename_tag + ".nc") \
           in os.listdir(self._datadir):
            print("Reading thermo_profile from pre-processed netCDF")
            self._read_netCDF(file_path + "thermo_" +
                            str(self.resolution.magnitude) +
                            str(self.resolution.units) +
                            self._ascent_filename_tag + ".nc")
            return

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
                                               serial_numbers
                                               ["imet" + str(i+1)])
        # End if-else blocks

        rh_raw = []
        # Fill rh_raw
        for key in temp_dict.keys():
            # Ensure only humidity is processed here
            if "rh" in key and "temp" not in key and "time" not in key:
                rh_raw.append(temp_dict[key].magnitude)

        alts = np.array(temp_dict["alt_pres"].magnitude)\
            * temp_dict["alt_pres"].units
        pres = np.array(temp_dict["pres"].magnitude)\
            * temp_dict["pres"].units

        time_rh = temp_dict["time_rh"]
        time_pres = temp_dict["time_pres"]
        time_temp = temp_dict["time_temp"]
        # Determine bad sensors
        rh_flags = utils.qc(rh_raw, 0.4, 0.2)  # TODO read these from file

        # Remove bad sensors
        for flags_ind in range(len(rh_flags)):
            if rh_flags[flags_ind] != 0:
                rh_raw[flags_ind] = np.full(len(rh_raw[flags_ind]), np.NaN)

        # Average the sensors
        for i in range(len(rh_raw[0])):
            rh.append(np.nanmean([rh_raw[a][i] for a in
                                  range(len(rh_raw))]))

        rh = np.array(rh) * units.percent

        # Determine which sensors are "bad"
        temp_flags = utils.qc(temp_raw, 0.25, 0.1)

        # Remove bad sensors
        temp_ind = 0  # track index in temp_raw  after items are removed.
        for flags_ind in range(len(temp_flags)):
            if temp_flags[flags_ind] != 0:
                print("Temperature sensor", temp_ind + 1, "removed")
                temp_raw[temp_ind] = \
                    [np.nan]*len(temp_raw[temp_ind])
            else:
                temp_ind += 1

        # Average the sensors
        for i in range(len(temp_raw[0])):

            temp.append(np.nanmean([temp_raw[a][i] for a in
                                   range(len(temp_raw))]))

        temp = np.array(temp) * units.kelvin

        #
        # Regrid to match times specified by Profile
        #

        # grid alt
        self.alt = utils.regrid_data(data=alts, data_times=time_pres,
                                     gridded_times=self.gridded_times,
                                     units=self._units)

        # grid pres
        self.pres = utils.regrid_data(data=pres, data_times=time_pres,
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

        minlen = min(len(self.alt), len(self.gridded_times), len(self.rh),
                     len(self.pres), len(self.temp))
        self.pres = self.pres[0:minlen-1]
        self.temp = self.temp[0:minlen-1]
        self.rh = self.rh[0:minlen-1]
        self.alt = self.alt[0:minlen-1]

        # Calculate mixing ratio
        self.mixing_ratio = calc.mixing_ratio_from_relative_humidity(
                            np.divide(self.rh.magnitude, 100), self.temp,
                            self.pres)

        if nc_level in 'low':
            self._save_netCDF(file_path + "thermo_" +
                              str(self.resolution.magnitude) +
                              str(self.resolution.units) +
                              self._ascent_filename_tag + ".nc")

    def _save_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """
        main_file = netCDF4.Dataset(file_path, "w",
                                    format="NETCDF4", mmap=False)

        main_file.createDimension("time", None)
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
        # TIME
        time_var = main_file.createVariable("time", "f8", ("time",))
        time_var[:] = netCDF4.date2num(self.gridded_times,
                                       units='microseconds since \
                                       2010-01-01 00:00:00:00')
        time_var.units = 'microseconds since 2010-01-01 00:00:00:00'

        main_file.close()

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        :param string file_path: file name
        """
        main_file = netCDF4.Dataset(file_path, "r",
                                    format="NETCDF4", mmap=False)
        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        self.alt = np.array(main_file.variables["alt"])[:-2] * \
            self._units.parse_expression(main_file.variables["alt"].units)
        self.pres = np.array(main_file.variables["pres"])[:-2] * \
            self._units.parse_expression(main_file.variables["pres"].units)
        self.rh = np.array(main_file.variables["rh"])[:-2] * \
            self._units.parse_expression(main_file.variables["rh"].units)
        self.temp = np.array(main_file.variables["temp"])[:-2] * \
            self._units.parse_expression(main_file.variables["temp"].units)
        self.mixing_ratio = np.array(main_file.variables["mr"])[:-2] * \
            self._units.parse_expression(main_file.variables["mr"].units)
        self.gridded_times = \
            np.array(netCDF4.num2date(main_file.variables["time"][:-2],
                                      units=main_file.variables["time"].units))

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
