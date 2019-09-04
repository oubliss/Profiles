"""
Calculates and stores wind parameters

Authors Brian Greene, Jessica Wiedemeier, Tyler Bell, Gus Azevedo \n
Copyright University of Oklahoma Center for Autonomous Sensing and Sampling
2019
"""
import numpy as np
import pandas as pd
import os
import utils
import metpy.calc
import netCDF4


class Wind_Profile():
    """ Contains data from one file.

    :var list<Quantity> u: U component of wind
    :var list<Quantity> v: V-component of wind
    :var list<Datetime> gridded_times: time of each point
    """

    def __init__(self, wind_dict, resolution, gridded_times=None,
                 indices=(None, None), ascent=True, units=None, file_path=''):
        """ Creates Wind_Profile object based on rotation data at the specified
        resolution.
        """
        # TODO add tail_number to wind_dict

        self.gridded_times = gridded_times
        self.resolution = resolution
        self.gridded_times = gridded_times
        self.indices = indices
        self.ascent = ascent
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
        print("80 speed is: ", speed)
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

        #
        # save NC
        #
        self._save_netCDF(file_path)

    def _calc_winds(self, wind_data):
        """ Calculate wind direction, speed, u, and v. Currently, this only
        works when the craft is HORIZONTALLY STATIONARY.
        :param dict wind_data: dictionary from Raw_Profile.get_wind_data()
        :param bool isCopter: True if rotor-wing, false if fixed-wing
        :rtype: tuple<list>
        :return: (direction, speed)
        """

        # TODO find a way to do this with a fixed-wing
        # TODO account for moving platform
        # TODO get tail_num from wind_dict
        tail_num = 0

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
        print("147 speed is: ", speed)
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
        main_file = netCDF4.Dataset(file_path + "wind_" +
                                    str(self.resolution.magnitude) +
                                    str(self.resolution.units) +
                                    self._ascent_filename_tag + ".nc", "w",
                                    format="NETCDF4", mmap=False)

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
        self.gridded_times = \
            np.array(netCDF4.num2date(main_file.variables["time"][:-2],
                                      units=main_file.variables["time"].units))

        main_file.close()
