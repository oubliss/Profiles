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


    def __init__(self, wind_dict, resolution, algorithm='linear', file_path=None,
               gridded_times=None, gridded_base=None, indices=(None, None),
               ascent=True, units=None, pos=None, nc_level=None, tail_number=None, meta=None):
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
            self.tail_number = tail_number
            self._indices = indices
            self._units = units
            self._datadir = os.path.dirname(file_path + ".json")

        # Extract the proper time arrays
        time_pres = wind_dict['time_pres']

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

        if algorithm == 'linear':
            direction, speed, time = self._calc_winds_linear(wind_dict)
        elif algorithm == 'quadratic':
            print("this wind algorithm isn't implmented")
            return
        elif algorithm == 'leso':
            print("This wind algorithm isn't implmented")
            return

        direction = direction % (2*np.pi)

        #
        # Regrid to res

        # grid alt and pres
        if (self.resolution.dimensionality ==
                self._units.get_dimensionality('m')):
            self.alt = gridded_base

            self.pres = utils.regrid_data(data=self.pres, data_times=time_pres,
                                          gridded_times=self.gridded_times,
                                          units=self._units)
        elif (self.resolution.dimensionality ==
              self._units.get_dimensionality('Pa')):
            self.pres = gridded_base
            self.alt = utils.regrid_data(data=self.alt, data_times=time_pres,
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


        # save NC

        if nc_level is not None:
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


    def ext_leso(self, rpma, rho_air, pn, pe, pd, vn, ve, vd, roll, pitch, yaw, h, order):
        from scipy.signal import StateSpace
        from control.matlab import ss, c2d


        # Dynamic Model with extended state:
        # x_3_dot = h(t) = 0, assuming generalized distrubance to be constant or slowly changing
        # For this study, the generalized disturbance is mainly drag caused by wind

        m = 2.28  # Mass (kg)
        g = 9.81  # Gravity (m/s2)
        dt = 0.1  # Sampling period (must be equal to CS logging rate)
        L = len(pn)
        wind_i = np.zeros((3, L))

        # Propeller model parameters
        # From Gill & D'Andrea 2019

        # APCE-10x7
        R = 0.14
        cl0 = 0.75;   cla = 6.6
        cd0 = 0.12;   cda = 2.8;        Nb = 2
        delta = 0.12; theta_tip = 0.21; c_tip = 5.7e-3

        sigma = Nb * c_tip / (np.pi * R)

        # State space model of mosition and velocity of the drone
        A = np.matrix([[0, 1, 0  ],
                       [0, 0, 1/m],
                       [0, 0, 0  ]])
        B = np.matrix([0, 1, 0]).T
        E = np.matrix([0, 0, 1]).T
        BF = np.concatenate((B, E), axis=1)
        Cv = np.matrix([[1, 0, 0],
                        [0, 1, 0]])

        # LESO gains
        lxp = 4; lyp = 3; lzp = 4  # Position error gains
        lxv = 4; lyv = 3; lzv = 4  # Velocity error gains
        alpha = 1.8   # Gain boost factor for LESOs of higher order

        # Observer System creation
        xcell = {}
        ycell = {}
        zcell = {}

        for i in range(order):
            Lxp = np.matrix([3*lxp, 3*lxp*lxp, lxp*lxp*lxp]).T
            Lxv = np.matrix([3*lxv, 3*lxv*lxv, lxv*lxv*lxv])
            aa = A - np.concatenate((Lxp, Lxv), axis=1)*cv
            ab = np.concatenate((BF, Lxp, Lxv), axis=1)
            ac = np.identity(3)
            ad = 0 * np.concatenate((BF, Lxp, Lxv), axis=1)
            aux = ss(aa, ab, ac, ad)
            xcell[i] = c2d(aux, dt)

            Lyp = np.matrix([3 * lyp, 3 * lyp * lyp, lyp * lyp * lyp]).T
            Lyv = np.matrix([3 * lyv, 3 * lyv * lyv, lyv * lyv * lyv])

            ba = A - np.concatenate((Lyp, Lyv), axis=1) * cv
            bb = np.concatenate((BF, Lyp, Lyv), axis=1)
            bc = np.identity(3)
            bd = 0 * np.concatenate((BF, Lyp, Lyv), axis=1)
            aux = ss(ba, bb, bc, bd)
            ycell[i] = c2d(aux, dt)

            Lzp = np.matrix([3 * lzp, 3 * lzp * lzp, lzp * lzp * lzp]).T
            Lzv = np.matrix([3 * lzv, 3 * lzv * lzv, lzv * lzv * lzv])

            ca = A - np.concatenate((Lzp, Lzv), axis=1) * cv
            cb = np.concatenate((BF, Lzp, Lzv), axis=1)
            cc = np.identity(3)
            cd = 0 * np.concatenate((BF, Lzp, Lzv), axis=1)
            aux = ss(ca, cb, cc, cd)
            zcell[i] = c2d(aux, dt)

        # Preallocate arrays
        Nest = np.zeros((3, L, order))
        Eest = np.zeros((3, L, order))
        Dest = np.zeros((3, L, order))
        Dn = np.zeros((L))
        De = np.zeros((L))
        Dd = np.zeros((L))
        F_h = np.zeros((L))
        F_t = np.zeros((L))
        Tq = np.zeros((L))

        for k in range(L-1):

            # Transform wind from NED to Body
            wind_b = utils.ned2body(wind_i[0, k], wind_i[1, k], wind_i[2, k], roll[k], pitch[k], yaw[k])

            # compute propeller H-forces
            omega = 2*np.pi*rpms[k] / 60  # Average angular velocity of the propellers [rad/s]
            lambda_c = wind_b[2] / (omega * R)  # Perpendicular air velocity (normalized)
            mu = np.abs(wind_b[0]) / (omega * R)  # Tangential air velocity (normalized)

            # Get thrust from extended prop model
            lambda_i = 1/8 * (
                    -4*lambda_c +
                    cla*sigma*(delta-1) +
                    np.sqrt(
                        16*lambda_c**2 +
                        8*cla*(delta-1)*lambda_c*sigma +
                        (1-1/delta)*sigma*(-8*cl0*delta*(1+delta) +
                                           cla*(cla*(delta-1)*delta*sigma - 8*(2*delta+mu**2)*theta_tip)) -
                        8*cl0*mu**2*sigma*np.log(delta)))

            lambda_f = lambda_c + lambda_i

            # Thrust coefficient
            c_ft = (0.5 * sigma / delta) * (
                    (1 - delta) * (cl0 * delta * (1 + delta) - 2 * cla * delta * (lambda_f-theta_tip) +
                                   cla * mu**2 * theta_tip) - cl0 * delta * mu * mu * log(delta))

            # Horizontal drag coeff
            c_fh = (0.5*mu*sigma/delta)*(
                    (1-delta)*(2*cd0*delta + theta_tip*((cla-2*cda)*lambda_f + 2*cda*theta_tip)) -
                    cl0*delta*lambda_f*np.log(delta))

            # Compute forces
            F_t[k] = np.cos(np.deg2rad(4)) * c_ft * 0.5 * rho_air[k] * np.pi * omega**2 * R**4
            F_h[k] = 4*np.sign(wind_b[0]) * c_fh * 0.5 * rho_air[k] * np.pi * omega**2 * R**4

            # Compute thrust using prop model, rpms, and intake air velocity
            Tq[k] = F_t[k]
            if np.isnan(Tq[k]):
                Tq[k] = 0

            # Transform propeller forces from body to NED
            F_H = utils.body2ned(np.matrix([F_h[k], 0, 0]).T, roll[k], pitch[k], yaw[k])
            T_corr = utils.body2ned(np.matrix([0, 0, -Tq[k]]).T, roll[k], pitch[k], yaw[k])
            U = np.matrix([0, 0, g]).T + (T_corr + F_H)/m

            # LESO Core function
            for j in range(order):
                if j == 0:  # 1st order loop estimation
                    Nest[:, k+1, j] = xcell[j].a * Nest[:, k, j] + \
                                      xcell[j].b * np.matrix([U[0], pn[k], vn[k]])
                    Eest[:, k+1, j] = ycell[j].a * Eest[:, k, j] + \
                                      ycell[j].b * np.matrix([U[1], pe[k], ve[k]])
                    Dest[:, k+1, j] = zcell[j].a * Dest[:, k, j] + \
                                      zcell[j].b * np.matrix([U[2], pd[k], vd[k]])
                else:  # 2nd order loop. Computes additional drag forces at higher gains
                    Nest[:, k+1, j] = xcell[j].a * Nest[:, k, j] + \
                                      xcell[j].b * np.matrix([U[0]+Dn[k+1]/m, Nest[0, k+1, j-1], Nest[1, k+1, j-1]]).T
                    Eest[:, k+1, j] = ycell[j].a * Eest[:, k, j] + \
                                      ycell[j].b * np.matrix([U[1]+De[k+1]/m, Eest[0, k+1, j-1], Eest[1, k+1, j-1]]).T
                    Dest[:, k+1, j] = zcell[j].a * Dest[:, k, j] + \
                                      zcell[j].b * np.matrix([U[2]+Dd[k+1]/m, Dest[0, k+1, j-1], Dest[1, k+1, j-1]]).T

                # In case the 2nd order or higher, add the additional drag forces to the 1st order.
                Dn[k+1] = Dn[k+1] + Nest[2, k+1, j]
                De[k+1] = De[k+1] + Nest[2, k+1, j]
                Dd[k+1] = Dd[k+1] + Nest[2, k+1, j]

            # wind velocity estimation
            # Extract drag forces and velocity estimates. Transform to body frame
            D_b = utils.ned2body(np.matrix([Dn[k+1], De[k+1], Dd[k+1]]).T, roll[k], pitch[k], yaw[k])
            V_b = utils.ned2body(np.matrix([Nest[1, k+1, -1], Eest[1, k+1, -1], Dest[1, k+1, -1]]).T,
                                 roll[k], pitch[k], yaw[k])

            # Second order friction equation for each up/down leg
            # Work in progress, not final
            cf_h = [0.037, -.12]; cf_v = -0.14388*Dd[k+1]/9.81 + 0.352; # These shouldn't be hard coded here, they must be adjustable params outside this function.
            wind_b[0] = V_b[0] + np.sign(D_b[0])*np.sqrt(np.abs(D_b[0]))/(cf_h[0]*rho_air[k])
            wind_b[1] = V_b[1] + np.sign(D_b[1])*np.sqrt(np.abs(D_b[1]))/(cf_h[0]*rho_air[k])
            wind_b[2] = V_b[2] + np.sign(D_b[2])*np.sqrt(np.abs(D_b[2]))/(cf_v*rho_air[k])

            # Transform wind vector from body to inertial frame
            wind_i[:, k+1] = utils.body2ned(wind_b, roll[k], pitch[k], yaw[k])

            return Dn, De, Dd, Nest, Eest, Dest, wind_i, Tq, F_h




    def _calc_winds_linear(self, wind_data):
        """ Calculate wind direction, speed, u, and v. Currently, this only
        works when the craft is HORIZONTALLY STATIONARY.
        :param dict wind_data: dictionary from Raw_Profile.get_wind_data()
        :param bool isCopter: True if rotor-wing, false if fixed-wing
        :rtype: tuple<list>
        :return: (direction, speed)
        """

        # TODO account for moving platform
        tail_num = self.tail_number

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


        coefs = utils.coef_manager.get_coefs('Wind', tail_num, 'E1')
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
