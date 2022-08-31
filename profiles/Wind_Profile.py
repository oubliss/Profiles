"""
Calculates and stores wind parameters
"""
import numpy as np
import pandas as pd
import datetime as dt
import os
import profiles.utils as utils
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

    def _init2(self, wind_dict, resolution, file_path=None,
               gridded_times=None, gridded_base=None, indices=(None, None),
               ascent=True, units=None, pos=None, nc_level='low', meta=None):
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
        :param Meta meta: the parent Profile's Meta object
        :param str nc_level: either 'low', or 'none'. This parameter \
           is used when processing non-NetCDF files to determine which types \
           of NetCDF files will be generated. For individual files for each \
           Raw, Thermo, \
           and Wind Profile, specify 'low'. For no NetCDF files, specify \
           'none'.
        """

        self._meta = meta
        if ascent:
            self._ascent_filename_tag = "Ascending"
        else:
            self._ascent_filename_tag = "Descending"

        try:
            self._read_netCDF(file_path + "wind_" +
                              str(resolution.magnitude) +
                              str(resolution.units) +
                              self._ascent_filename_tag + ".nc")
            return

        except Exception:
            self.resolution = resolution
            self.gridded_times = gridded_times
            self.ascent = ascent
            self.pres = wind_dict["pres"]
            self.alt = wind_dict["alt"]
            self._indices = indices
            self._units = units
            self._datadir = os.path.dirname(file_path + ".json")

        
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

        # grid alt and pres
        if (self.resolution.dimensionality ==
                self._units.get_dimensionality('m')):
            self.alt = gridded_base
            self.pres = utils.regrid_data(data=self.pres, data_times=time,
                                          gridded_times=self.gridded_times,
                                          units=self._units)
        elif (self.resolution.dimensionality ==
              self._units.get_dimensionality('Pa')):
            self.pres = gridded_base
            self.alt = utils.regrid_data(data=self.alt, data_times=time,
                                         gridded_times=self.gridded_times,
                                         units=self._units)

        self.dir = utils.regrid_data(data=direction, data_times=time,
                                     gridded_times=self.gridded_times,
                                     units=self._units)
        self.speed = utils.regrid_data(data=speed, data_times=time,
                                       gridded_times=self.gridded_times,
                                       units=self._units)
        self.u, self.v = metpy.calc.wind_components(self.speed, self.dir)

        if pos is not None:
            # grid lat
            self.lat = utils.regrid_data(data=pos['lat'], data_times=pos['time'],
                                         gridded_times=self.gridded_times,
                                         units=self._units)

            self.lon = utils.regrid_data(data=pos['lon'], data_times=pos['time'],
                                         gridded_times=self.gridded_times,
                                         units=self._units)

        else:
            self.lat = np.full_like(self.gridded_times, -999.)
            self.lon = np.full_like(self.gridded_times, -999.)

        """ 
        TB --  I don't think this is needed... 
        and it causes unexpected down stream issues...
        may bite me later... 
        If causing issues later, do 'minlen+1'
        """
        # minlen = min([len(self.u), len(self.v), len(self.dir),
        #               len(self.speed), len(self.alt), len(self.pres),
        #               len(self.gridded_times)])
        # self.u = self.u[0:minlen]
        # self.v = self.v[0:minlen]
        # self.dir = self.dir[0:minlen]
        # self.speed = self.speed[0:minlen]
        # self.alt = self.alt[0:minlen]
        # self.pres = self.pres[0:minlen]
        # self.gridded_times = self.gridded_times[0:minlen]

        #
        # save NC
        #
        if nc_level in 'low':
            self._save_netCDF(file_path)

    def truncate_to(self, new_len):
        """ Shortens arrays to have no more than new_len data points

        :param new_len: The new, shorter length
        :return: None
        """

        self.u = self.u[:new_len]
        self.v = self.v[:new_len]
        self.dir = self.dir[:new_len]
        self.speed = self.speed[:new_len]
        self.alt = self.alt[:new_len]
        self.pres = self.pres[:new_len]
        self.gridded_times = self.gridded_times[:new_len]

    def _calc_winds(self, wind_data):
        """ Calculate wind direction, speed, u, and v. Currently, this only
        works when the craft is HORIZONTALLY STATIONARY.
        :param dict wind_data: dictionary from Raw_Profile.get_wind_data()
        :param bool isCopter: True if rotor-wing, false if fixed-wing
        :rtype: tuple<list>
        :return: (direction, speed)
        """

        # TODO account for moving platform
        tail_num = utils.coef_manager.get_tail_n(wind_data['serial_numbers']['copterID'])

        # psi and az represent the copter's direction in spherical coordinates
        psi = np.zeros(len(wind_data["roll"])) * self._units.rad
        az = np.zeros(len(wind_data["roll"])) * self._units.rad

        for i in range(len(wind_data["roll"])):
            # croll is cos(roll), sroll is sin(roll)...
            croll = np.cos(wind_data["roll"][i]).magnitude
            sroll = np.sin(wind_data["roll"][i]).magnitude
            cpitch = np.cos(wind_data["pitch"][i]).magnitude
            spitch = np.sin(wind_data["pitch"][i]).magnitude
            cyaw = np.cos(wind_data["yaw"][i]).magnitude
            syaw = np.sin(wind_data["yaw"][i]).magnitude

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


        coefs = utils.coef_manager.get_coefs('Wind', tail_num)
        speed = float(coefs['A']) * np.sqrt(np.tan(psi)).magnitude + float(coefs['B'])

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

        if '.nc' in file_path or '.cdf' in file_path:
            file_name = file_path
        elif self._meta is not None:
            file_name = str(self._meta.get("location")).replace(' ', '') + str(self.resolution.magnitude) + \
                    str(self._meta.get("platform_id")) + "CMT" + \
                    "wind_" + self._ascent_filename_tag + ".c1." + \
                    self._meta.get("timestamp").replace("_", ".") + ".cdf"
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

        main_file.close()

    def _read_netCDF(self, file_path):
        """ Reads data from a NetCDF file. Called by the constructor.

        :param string file_path: file name
        """
        main_file = netCDF4.Dataset(file_path, "r",
                                    format="NETCDF4", mmap=False)
        # Note: each data chunk is converted to an np array. This is not a
        # superfluous conversion; a Variable object is incompatible with pint.

        self.dir = np.array(main_file.variables["dir"]) * \
            self._units.parse_expression(main_file.variables["dir"].units)
        self.speed = np.array(main_file.variables["speed"]) * \
            self._units.parse_expression(main_file.variables["speed"].units)
        self.u = np.array(main_file.variables["u"]) * \
            self._units.parse_expression(main_file.variables["u"].units)
        self.v = np.array(main_file.variables["v"]) * \
            self._units.parse_expression(main_file.variables["v"].units)
        self.alt = np.array(main_file.variables["alt"]) * \
            self._units.parse_expression(main_file.variables["alt"].units)
        self.pres = np.array(main_file.variables["pres"]) * \
            self._units.parse_expression(main_file.variables["pres"].units)
        base_time = dt.datetime(2010, 1, 1, 0, 0, 0, 0)
        self.gridded_times = []
        for i in range(len(main_file.variables["time"][:])):
            self.gridded_times.append(base_time + dt.timedelta(microseconds=
                                                               int(main_file.variables
                                                                   ["time"][i])))

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
