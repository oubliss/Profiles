"""
Calculates and stores wind parameters

Authors Brian Greene, Jessica Blunt, Tyler Bell, Gus Azevedo \n
Copyright University of Oklahoma Center for Autonomous Sensing and Sampling
2019

Component of Profiles v1.0.0
"""
import numpy as np
import pandas as pd
import os
import utils
import metpy.calc
import netCDF4
from copy import deepcopy, copy


class Wind_Profile():
    """ Processes and holds wind data from one vertical profile

    :var list<Quantity> u: U component of wind
    :var list<Quantity> v: V-component of wind
    :var list<Quantity> dir: wind direction
    :var list<Quantity> speed: wind speed
    :var list<Quantity> pres: air pressure
    :var list<Quantity> alt: altitude
    :var list<Datetime> gridded_times: time of each point
    :var Quantity resolution: the vertical resolution of the processed data
    :var bool ascent: is data from the ascending leg of the flight processed?\
       If not, False.
    """

    def __init__(self, *args, **kwargs):
        if len([*args]) > 0:
            self._init2(*args, **kwargs)

    def _init2(self, wind_dict, resolution, gridded_times=None,
                indices=(None, None), ascent=True, units=None, file_path='',
                nc_level='low'):
        """ Creates Wind_Profile object based on rotation data at the specified
        resolution

        :param dict wind_dict: the dictionary produced by \
           Raw_Profile.get_wind_data()
        :param Quantity resolution: vertical resolution of the processed data
        :param List<Datetime> gridded_times: times for which Profile has \
           requested wind data
        :param tuple<int> indices: if applicable, the user-defined bounds of \
           the profile
        :param bool ascent: is data from the ascending leg of the flight \
           processed? If not, False.
        :param metpy.units units: the units defined by Raw_Profile
        :param str file_path: the original file passed to the package
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """

        self.resolution = resolution
        self.gridded_times = gridded_times
        self.ascent = ascent
        self.pres = wind_dict["pres"]
        self.alt = wind_dict["alt"]
        self._indices = indices
        self._units = units
        if ascent:
            self._ascent_filename_tag = "Ascending"
        else:
            self._ascent_filename_tag = "Descending"

        self._datadir = os.path.dirname(file_path + ".json")
        filepath_nc = file_path + "wind_" + str(resolution.magnitude) + \
            str(resolution.units) + self._ascent_filename_tag + ".nc"

        if os.path.basename(filepath_nc) in os.listdir(self._datadir):
            print("Reading wind_profile from pre-processed netCDF")
            self._read_netCDF(file_path)
            return
        else:
            # If no indices given, use entire file
            if not indices[0] is None:
                # trim profile
                selection = np.where(np.array(wind_dict["time"]) > indices[0],
                                     np.array(wind_dict["time"]) < indices[1],
                                     False)

                wind_dict["roll"] = \
                    np.array(wind_dict["roll"].magnitude)[selection] * \
                    wind_dict["roll"].units
                wind_dict["pitch"] = \
                    np.array(wind_dict["pitch"].magnitude)[selection] * \
                    wind_dict["pitch"].units
                wind_dict["yaw"] = \
                    np.array(wind_dict["yaw"].magnitude)[selection] * \
                    wind_dict["yaw"].units
                wind_dict["speed_east"] = \
                    np.array(wind_dict["speed_east"].magnitude)[selection] * \
                    wind_dict["speed_east"].units
                wind_dict["speed_north"] = \
                    np.array(wind_dict["speed_north"].magnitude)[selection] \
                    * wind_dict["speed_north"].units
                wind_dict["speed_down"] = \
                    np.array(wind_dict["speed_down"].magnitude)[selection] * \
                    wind_dict["speed_down"].units
                wind_dict["time"] = np.array(wind_dict["time"])[selection]

        direction, speed, time = self._calc_winds(wind_dict)

        direction = direction % (2*np.pi)

        #
        # Regrid to res
        #
        self.dir = utils.regrid_data(data=direction, data_times=time,
                                     gridded_times=self.gridded_times,
                                     units=self._units)
        self.speed = utils.regrid_data(data=speed, data_times=time,
                                       gridded_times=self.gridded_times,
                                       units=self._units)
        self.u, self.v = metpy.calc.wind_components(self.speed, self.dir)

        self.alt = utils.regrid_data(data=self.alt, data_times=time,
                                     gridded_times=self.gridded_times,
                                     units=self._units)

        #
        # save NC
        #
        if nc_level in 'low':
            self._save_netCDF(file_path)

    def _calc_winds(self, wind_data):
        """ Calculate wind direction, speed, u, and v. Currently, this only
        works when the craft is HORIZONTALLY STATIONARY.
        :param dict wind_data: dictionary from Raw_Profile.get_wind_data()
        :param bool isCopter: True if rotor-wing, false if fixed-wing
        :rtype: tuple<list>
        :return: (direction, speed)
        """

        # TODO account for moving platform
        tail_num = wind_data["serial_numbers"]["wind"]

        # psi and az represent the copter's direction in spherical coordinates
        psi = np.zeros(len(wind_data["roll"])) * self._units.rad
        az = np.zeros(len(wind_data["roll"])) * self._units.rad

        for i in range(len(wind_data["roll"])):
            # croll is cos(roll), sroll is sin(roll)...
            croll = np.cos(wind_data["roll"][i])
            sroll = np.sin(wind_data["roll"][i])
            cpitch = np.cos(wind_data["pitch"][i])
            spitch = np.sin(wind_data["pitch"][i])
            cyaw = np.cos(wind_data["yaw"][i])
            syaw = np.sin(wind_data["yaw"][i])

            Rx = np.matrix([[1, 0, 0],
                            [0, croll, sroll],
                            [0, -sroll, croll]])
            Ry = np.matrix([[cpitch, 0, -spitch],
                            [0, 1, 0],
                            [spitch, 0, cpitch]])
            Rz = np.matrix([[cyaw, -syaw, 0],
                            [syaw, cyaw, 0],
                            [0, 0, 1]])
            R = Rz * Ry * Rx

            psi[i] = np.arccos(R[2, 2])
            az[i] = np.arctan2(R[1, 2], R[0, 2])

        coefs = pd.read_csv('coefs/MasterCoefList.csv')
        a_spd = float(coefs.A[coefs.SerialNumber == tail_num]
                      [coefs.SensorType == "Wind"])
        b_spd = float(coefs.B[coefs.SerialNumber == tail_num]
                      [coefs.SensorType == "Wind"])

        speed = a_spd * np.sqrt(np.tan(psi)) + b_spd

        speed = speed * self._units.m / self._units.s
        # Throw out negative speeds
        speed[speed.magnitude < 0.] = np.nan

        # Fix negative angles
        az = az.to(self._units.deg)
        iNeg = np.squeeze(np.where(az.magnitude < 0.))
        az[iNeg] = az[iNeg] + 360. * self._units.deg

        # az is the wind direction, speed is the wind speed
        return (az, speed, wind_data["time"])

    def _save_netCDF(self, file_path):
        """ Save a NetCDF file to facilitate future processing if a .JSON was
        read.

        :param string file_path: file name
        """
        minlen = min([len(self.u), len(self.v), len(self.dir),
                      len(self.speed), len(self.alt), len(self.pres),
                      len(self.gridded_times)])

        main_file = netCDF4.Dataset(file_path + "wind_" +
                                    str(self.resolution.magnitude) +
                                    str(self.resolution.units) +
                                    self._ascent_filename_tag + ".nc", "w",
                                    format="NETCDF4", mmap=False)

        main_file.createDimension("time", None)
        # DIRECTION
        dir_var = main_file.createVariable("dir", "f8", ("time",))
        dir_var[:] = self.dir[:minlen-1].magnitude
        dir_var.units = str(self.dir.units)
        # SPEED
        spd_var = main_file.createVariable("speed", "f8", ("time",))
        spd_var[:] = self.speed[:minlen-1].magnitude
        spd_var.units = str(self.speed.units)
        # U
        u_var = main_file.createVariable("u", "f8", ("time",))
        u_var[:] = self.u[:minlen-1].magnitude
        u_var.units = str(self.u.units)
        # V
        v_var = main_file.createVariable("v", "f8", ("time",))
        v_var[:] = self.v[:minlen-1].magnitude
        v_var.units = str(self.v.units)
        # ALT
        alt_var = main_file.createVariable("alt", "f8", ("time",))
        alt_var[:] = self.alt[:minlen-1].magnitude
        alt_var.units = str(self.alt.units)
        # PRES
        pres_var = main_file.createVariable("pres", "f8", ("time",))
        pres_var[:] = self.pres[:minlen-1].magnitude
        pres_var.units = str(self.pres.units)

        # TIME
        time_var = main_file.createVariable("time", "f8", ("time",))
        time_var[:] = netCDF4.date2num(self.gridded_times[:minlen-1],
                                       units='microseconds since \
                                       2010-01-01 00:00:00:00')
        time_var.units = 'microseconds since 2010-01-01 00:00:00:00'

        main_file.close()

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        :param string file_path: file name
        """
        main_file = netCDF4.Dataset(file_path + "wind_" +
                                    str(self.resolution.magnitude) +
                                    str(self.resolution.units) +
                                    self._ascent_filename_tag + ".nc", "r",
                                    format="NETCDF4", mmap=False)
        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        self.dir = np.array(main_file.variables["dir"])[:-2] * \
            self._units.parse_expression(main_file.variables["dir"].units)
        self.speed = np.array(main_file.variables["speed"])[:-2] * \
            self._units.parse_expression(main_file.variables["speed"].units)
        self.u = np.array(main_file.variables["u"])[:-2] * \
            self._units.parse_expression(main_file.variables["u"].units)
        self.v = np.array(main_file.variables["v"])[:-2] * \
            self._units.parse_expression(main_file.variables["v"].units)
        self.alt = np.array(main_file.variables["alt"])[:-2] * \
            self._units.parse_expression(main_file.variables["alt"].units)
        self.pres = np.array(main_file.variables["pres"])[:-2] * \
            self._units.parse_expression(main_file.variables["pres"].units)
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
        to_return = "\t\tWind_Profile" \
                    + "\n\t\t\tdir:   " + str(type(self.dir)) \
                    + "\n\t\t\tspeed: " + str(type(self.speed)) \
                    + "\n\t\t\tu:     " + str(type(self.u)) \
                    + "\n\t\t\tv:     " + str(type(self.v)) \
                    + "\n\t\t\talt:   " + str(type(self.alt)) \
                    + "\n\t\t\tpres:  " + str(type(self.pres)) + "\n"
        return to_return
