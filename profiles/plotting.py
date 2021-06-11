import os
import cmocean
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as datenum
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from scipy.interpolate import interp1d
from profiles.UnitFormatter import UnitFormatter
import metpy.plots as mpplots


vars = {'theta': ["Potential Temperature", 'theta', 'K', cmocean.cm.thermal,
                  1.0],
        'temp': ["Temperature", 'temp', '$^\circ$C', cmocean.cm.thermal, 1.0],
        'T_d': ["Dewpoint Temperature", 'T_d', '$^\circ$C', cmocean.cm.haline,
                1.0],
        'dewp': ["Dewpoint Temperature", 'T_d', '$^\circ$C', cmocean.cm.haline,
                 1.0],
        'r': ["Mixing Ratio", 'mixing_ratio', 'g Kg$^{-1}$', cmocean.cm.haline,
              0.5],
        'mr': ["Mixing Ratio", 'mixing_ratio', 'g Kg$^{-1}$', cmocean.cm.haline,
               0.5],
        'q': ["Specific Humidity", 'q', 'g Kg$^{-1}$', cmocean.cm.haline, 0.5],
        'rh': ["Relative Humidity", 'rh', '%', cmocean.cm.haline, 5.0],
        'speed': ["Wind Speed", 'speed', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'ws': ["Wind Speed", 'speed', 'm s$^{-1}$', 'gist_stern_r', 5.0],
        'u': ["U", 'u', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'v': ["V", 'v', 'm s$^{-1}$', cmocean.cm.speed, 5.0],
        'dir': ["Wind Direction", 'dir', '$^\circ$', cmocean.cm.phase, 360.,
                'wind'],
        'pres': ["Pressure", 'pres', 'Pa', cmocean.cm.haline, 15.0],
        'p': ["Pressure", 'pres', 'Pa', cmocean.cm.haline, 15.0],
        'alt': ["Altitude", 'alt', 'm', cmocean.cm.haline, 10.0]}


# File path to logos added to plots
fpath_logos = os.path.join(os.getcwd(), 'resources', 'CircleLogos.png')


def contour_height_time(profiles, var=['temp'], use_pres=False):
    """ contourHeightTime creates a filled contour plot of the first element of
       var in a time-height coordinate system. If len(var) > 1, it also
       overlays unfilled contours of the remaining elements. No more than 4
       variables can be plotted at once.
       Accepted variable names are:

       * 'theta'
       * 'temp'
       * 'T_d'
       * 'dewp'
       * 'r'
       * 'mr'
       * 'q'
       * 'rh'
       * 'speed':
       * 'u'
       * 'v'
       * 'dir'
       * 'pres'
       * 'p'
       * 'alt'

    :param list profiles: a list of all profiles to be included in the plot
    :param list<str> var: names of the variable to be plotted
    :rtype: matplotlib.figure.Figure
    :return: the contoured plot
    """

    plt.figure()  # Don't append this to an existing figure

    legend_handles = []
    for var_i in var:
        if var_i not in vars.keys():
            print(var_i + " was not recognized. Try one of these:\n" +
                  str(vars.keys()))

    times = []  # datenum.date2num(list)
    z = []  # unitless
    data = {}  # also unitless
    data_units = {}
    for var_i in var:
        data[var_i] = []

    linestyles = ['solid', 'dashed', 'dashdot']
    style_ind = 0

    for i in range(len(profiles)):
        # Get data from Profile objects
        times.append(list(profiles[i].gridded_times))
        z.append(profiles[i].get("gridded_base").magnitude)
        for var_i in var:
            data[var_i].append(list(profiles[i].get(vars[var_i][1]).magnitude))
            if var_i not in data_units.keys():
                data_units[var_i] = profiles[i].get(vars[var_i][1]).units

    # Now there are 3 parallel lists for each profile.
    # Force them to share z
    max_len = 0
    which_i = -1
    for i in range(len(z)):
        if len(z[i]) > max_len:
            max_len = len(z[i])
            which_i = i

    z = z[which_i]
    # There is now only one list for z - all profiles have to share

    time_flat = np.array(times[0])
    for p in range(len(times)):
        if p > 0:
            time_flat = np.concatenate((time_flat, times[p]))
    timerange = datenum.drange(np.nanmin(time_flat),
                               np.nanmax(time_flat),
                               (np.nanmax(time_flat)
                                -np.nanmin(time_flat))/100)
    q = (np.nanmax(time_flat)-np.nanmin(time_flat))/100
    for i in range(len(times)):
        diff = max_len - len(times[i])
        for j in range(diff):
            times[i].append(None)
            for var_i in var:
                data[var_i][i].append(np.nan)

    # Convert datetime to datenum
    for p in range(len(times)):
        for i in range(len(times[p])):
            try:
                times[p][i] = datenum.date2num(times[p][i])
            except AttributeError:
                times[p][i] = np.nan

    # Switch to arrays to make indexing easier
    times = np.array(times, dtype=float)
    z = np.array(z)
    if use_pres:
        z *= 0.01
    XX, YY = np.meshgrid(timerange, z)

    # Prepare for interpolation
    data_grid = {}
    fig = None
    for var_i in var:
        data_grid[var_i] = np.full_like(XX, np.nan)
        data[var_i] = np.array(data[var_i])

    # For z in z interp1d with time as x and data as y
    # Force back into one grid
    for var_i in var:
        for i in range(len(z)):
            a = list(np.array(times[:, i]).ravel())
            j = 0
            while j < len(a):
                if np.isnan((a[j])):
                    a.remove(a[j])
                else:
                    j += 1
            if len(a) < 2:
                continue
            interp_fun = interp1d(np.array(times[:, i]).ravel(), data[var_i][:, i],
                                  fill_value='extrapolate', kind='cubic')
            data_grid[var_i][i, :] = interp_fun(XX[i, :])

        # Set up figure
        if fig is None:
            start = 0
            end = -1
            for r in range(len(data_grid[var_i])):
                if not np.isnan(data_grid[var_i][r][0]):
                    start = r
                    break
            for r in range(len(data_grid[var_i])):
                if r > start and np.isnan(data_grid[var_i][r][0]):
                    end = r
                    break
            fig, ax = plt.subplots(1, figsize=(16, 9))
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
            ax.xaxis.set_major_formatter(datenum.DateFormatter('%H:%M:%S'))
            plt.ylim((z[start], z[end]))
            if use_pres:
                plt.ylabel("Pressure (hPa)", fontsize=15)
            else:
                plt.ylabel("Altitude (m MSL)", fontsize=15)
            plt.xlabel("Time (UTC)", fontsize=15)
            ax.tick_params(labelsize=14)

            # Make filled contour plot
            cfax = ax.pcolormesh(XX[start:end], YY[start:end],
                                 data_grid[var_i][start:end],
                                 cmap=vars[var_i][3])
            cbar = plt.colorbar(cfax, ax=ax, pad=0.01, )
            cbar.set_label(vars[var_i][0] + " (" + str(data_units[var_i]) + ")",
                           rotation=270, fontsize=20, labelpad=30)
        else:
            # Make unfilled contour plot
            cfax = ax.contour(XX[start:end], YY[start:end],
                              data_grid[var_i][start:end],
                              np.linspace(np.nanmin(data[var_i]),
                                          np.nanmax(data[var_i]), 10),
                              colors='black', linestyles=linestyles[style_ind])
            legend_handles.append(
                mlines.Line2D([], [], color='black', label=vars[var_i][0],
                              linestyle=linestyles[style_ind], marker='.',
                              markersize=1))
            style_ind += 1
            plt.clabel(cfax, fontsize=12,
                       fmt=UnitFormatter(unit=vars[var_i][2],
                                               places=1))


    for p_times in times:
        ax.scatter(p_times.astype(float), z, c='black', s=0.5)
    legend_handles.append(mlines.Line2D([], [], color='black',
                          label="Data collection points",
                          linestyle='dotted'))
    ax.legend(handles=legend_handles, fontsize=14, framealpha=1.0, loc=4)

    return fig

# TODO plot all variables
# TODO determine x_lim
# TODO Test with big flight data
def plot_skewT(profiles, wind_barbs=False, barb_density=10):
    r""" Plots a SkewT diagram.
    :param list<number> profiles: profiles which contain T_d, press, and temp data
    :param bool wind_barbs: if True, plot wind barbs. Requires that profiles contain u, v data.
    :param int barb_density: n for which every nth barb is plotted
    :rtype: matplotlib.figure.Figure
    :return: fig containing a SkewT diagram of the data
    """
    # Create plot
    fig = mpplots.SkewT(rotation=30, aspect=350)
    fig._fig.set_figheight(9)
    fig._fig.set_figwidth(9)
    fig.ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    fig.plot_dry_adiabats(linewidth=0.5, label="Dry Adiabats")
    fig.plot_moist_adiabats(linewidth=0.5, label="Moist Adiabats")
    fig.plot_mixing_lines(linewidth=0.5, label="Mixing Ratio")
    units = profiles[0]._units

    for profile in profiles:

        # Pressure is used several times, so copy here to minimize passing
        pres = profile.get("pres").to(units.hectopascal)
        # Add data to plot
        fig.plot(pres, profile.get("temp"), 'r', label="Temperature")
        fig.plot(pres, profile.get("T_d"), 'g', label="Dewpoint")
        
        
        u = profile.get("u").to(units.knots)
        v = profile.get("v").to(units.knots)

        if wind_barbs:
            fig.plot_barbs(pres.astype('float64')[::barb_density],
                           u.astype('float64')[::barb_density],
                           v.astype('float64')[::barb_density])

    plt.legend(loc='upper left')
    fig.ax.set_ylim(np.nanmax(pres.to(units.hPa).magnitude) + 10,
                    np.nanmin(pres.to(units.hPa).magnitude) - 20)
    fig.ax.set_xlim(np.nanmin(profiles[0].get("T_d").to(units.degC).magnitude) - 5,
                    np.nanmax(profiles[0].get("temp").to(units.degC).magnitude) + 10)
    return fig
'''
def meteogram(fpath):
    """ Graphically displays Mesonet data.

    Four subplots are created, each with time on the horizontal axis. The top
    plot is of T and Td, the second of P, the third of wind speed and
    direction, and the bottom of solar radiation.

    :param string fpath: the file path for the Mesonet timeseries data
    :rtype: matplotlib.figure.Figure
    :return: figure containing four horizontal subplots
    """

    return


def plot_var_time(var=None, t=None, times=None):
    """ Plots var vs time.

    :param list<list<Quantity>> var: the variables to be plotted
    :param list<datetime> t: times corresponding to the data
    :param tuple<datetime> times: start and end times to highlight
    :rtype: matplotlib.figure.Figure
    :return: plot of var vs. time
    """

    return





def summary(temp=None, pres=None, t_d=None, u=None, v=None, dt=None,
            loc=(None, None), flight_name=None):
    """ Creates a figure with a SkewT, hodograph, map, and logos. The plot
    title will be <flight_name> at <loc>, <time>. Logos can be
    changed by altering the contents of directory resources.

    :param list<number> temp: Temperatures to plot in C
    :param list<number> pres: Pressures to plot in ?
    :param list<number> t_d: Dewpoints to plot in C
    :param list<number> u: U-component of wind in kts
    :param list<number> v: V-component of wind in kts
    :param datetime dt: start time of flight (for title)
    :param tuple<number> loc: lat, lon pair of location data
    :param string flight_name: name of the flight (for title)
    :rtype: matplotlib.figure.Figure
    :return: the summary figure (described above)
    """

    # Create hodograph
    h = Hodograph()

    # Create SkewT
    s = plot_skewT()

    # Create map

    # Create logos
    logos = mpimg.imread()

    # Draw all

    return h, s, logos  # change to subplots format


def rh_comp_co2(rh):
    """ Plot averaged, QC'd relative humidity against time for sensors inside
       and sensors outside of the CO2 box.

    :param tuple rh: relative humidity as (rh1, rh2, ..., time: ms)
    """

# class plotSkewT():
#     def __init__(self, T=None, pres=None, Td=None, u=None, v=None,
#                  dt_start=None, **kwargs):
#         self.T = T
#         self.pres = pres
#         self.Td = Td
#         self.u = u
#         self.v = v
#         self.dt_start = dt_start
#         argDict = {'parcel': None,
#                    'lclpres': None,
#                    'lcltemp': None,
#                    'SBCAPE': None,
#                    'pmeso': None,
#                    'T2meso': None,
#                    'T9meso': None,
#                    'RH2meso': None,
#                    'Td2meso': None,
#                    'umeso': None,
#                    'vmeso': None,
#                    'radmeso': None,
#                    'tmeso': None,
#                    'loc': None,
#                    'copnum': None,
#                    'lat': None,
#                    'lon': None}
#
#         argDict.update(kwargs)
#         self.argDict = argDict
#         self.windkts = np.sqrt(u ** 2 + v ** 2)
#
#         # convert RH to Td if not defined
#         if (self.argDict['Td2meso'] is None) & \
#                 (self.argDict['RH2meso'] is not None):
#             self.argDict['Td2meso'] = np.array(
#                 mcalc.dewpoint_rh(self.argDict['T2meso'] * units.degC,
#                                   self.argDict['RH2meso'] / 100.))
#
#     def plot(self):
#         fig = plt.figure(figsize=(9, 8))
#         gs0 = gridspec.GridSpec(1, 2)
#         gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0],
#                                                wspace=0)
#         gs2 = gridspec.GridSpecFromSubplotSpec(6, 1, subplot_spec=gs0[1],
#                                                wspace=0, hspace=0.4)
#         # fig.subplots_adjust(wspace=0.1, hspace=0.4)
#         skew = SkewT(fig, rotation=20, subplot=gs1[:, :])
#
#         skew.plot(self.pres, self.T, 'r', linewidth=2)
#         skew.plot(self.pres, self.Td, 'g', linewidth=2)
#         skew.plot_barbs(self.pres[0::4], self.u[0::4], self.v[0::4],
#                         x_clip_radius=0.12, y_clip_radius=0.12)
#
#         # if mesonet data exist, plot
#         if self.argDict['pmeso'] is not None:
#             # plot mesonet data
#             p9meso = self.argDict['pmeso'] * (1. - (7. * 9.8) / (
#                 287. * (self.argDict['T2meso'] + 273.15)))
#             skew.plot(self.argDict['pmeso'], self.argDict['T2meso'], 'k*',
#                       linewidth=2, label='Mesonet 2 m T')
#             skew.plot(p9meso, self.argDict['T9meso'], 'r*',
#                       linewidth=2, label='Mesonet 9 m T')
#             skew.plot(self.argDict['pmeso'], self.argDict['Td2meso'], 'g*',
#                       linewidth=2, label='Mesonet 2 m Td')
#             skew.plot_barbs(self.pres[0], self.argDict['umeso'],
#                             self.argDict['vmeso'], barbcolor='r')  # , label='Mesonet 10 m Wind')
#
#         hand, lab = skew.ax.get_legend_handles_labels()
#
#         # if parcel data exist, plot with lcl
#         if self.argDict['lclpres'] is not None:
#             skew.plot(self.argDict['lclpres'],
#                       self.argDict['lcltemp'], 'ko',
#                       markerfacecolor='black')
#             skew.plot(self.pres, self.argDict['parcel'], 'k',
#                       linewidth=2)
#
#         # axis limits
#         # check temperature range - if below 0C, use range -20 to +10
#         # if all above 0C, use range 0 to +30
#         # if straddling 0C, use range -10 to +20
#         if (np.nanmax(self.T) <= 0.):
#             self.xmin = -30.
#             self.xmax = 10.
#         elif (np.nanmin(self.Td) <= 0.) & (np.nanmax(self.T) >= 0.):
#             self.xmin = -10.
#             self.xmax = 30.
#         elif (np.nanmin(self.Td) >= 0.):
#             self.xmin = 0.
#             self.xmax = 40.
#         else:
#             print '>>x axis limits error '
#             self.xmin = -10.
#             self.xmax = 30.
#         # y limits - use LCL as max if higher than profile
#         if (self.argDict['lclpres']) < np.nanmin(self.pres):
#             self.ymin = round((self.argDict['lclpres']), -1) - 10
#         else:
#             self.ymin = round(np.nanmin(self.pres), -1) - 10
#         self.ymax = round(np.nanmax(self.pres), -1) + 10
#
#         skew.ax.set_ylim(self.ymax, self.ymin)
#         skew.ax.set_xlim(self.xmin, self.xmax)
#         skew.ax.set_yticks(np.arange(self.ymin, self.ymax + 10, 10))
#         skew.ax.set_xlabel('Temperature ($^\circ$C)')
#         skew.ax.set_ylabel('Pressure (hPa)')
#         self.titleName = '{0} {1} UTC - {2}'.format(
#             self.argDict['copnum'],
#             self.dt_start.strftime('%d-%b-%Y %H:%M:%S'),
#             self.argDict['loc'])
#         skew.ax.set_title(self.titleName)
#
#         skew.plot_dry_adiabats(linewidth=0.75)
#         skew.plot_moist_adiabats(linewidth=0.75)
#         skew.plot_mixing_lines(linewidth=0.75)
#
#         ## Hodograph
#         ax_hod = fig.add_subplot(gs2[:2, 0])
#         if np.nanmax(self.windkts) > 18.:
#             comprange = 35
#         else:
#             comprange = 20
#
#         h = Hodograph(ax_hod, component_range=comprange)
#         h.add_grid(increment=5)
#         h.plot_colormapped(self.u, self.v, self.pres, cmap=cmocean.cm.deep_r)
#         ax_hod.set_title('Hodograph (kts)')
#         ax_hod.yaxis.set_ticklabels([])
#
#         ## Map
#         if self.argDict['loc'] == 'HAIL':
#             lllat = 55.34
#             urlat = 70.54
#             lat_0 = 65.
#             lon_0 = 24.
#             ax_map = fig.add_subplot(gs2[3:5, 0])
#             m = Basemap(width=1600000, height=900000, projection='lcc',
#                         resolution='l', lat_1=lllat, lat_2=urlat, lat_0=lat_0,
#                         lon_0=lon_0)
#             m.drawcountries()
#             m.shadedrelief()
#             x, y = m(24.555, 65.038)
#             plt.plot(x, y, 'b.')
#
#         elif self.argDict['loc'] in ['K04V', 'CRES', 'MOFF', 'SAGF', 'K1V8']:
#             lllat = 30.
#             urlat = 50.
#             lat_0 = 37.9
#             lon_0 = -105.7
#             ax_map = fig.add_subplot(gs2[3:5, 0])
#             m = Basemap(width=1600000, height=900000, projection='lcc',
#                         resolution='l', lat_1=lllat, lat_2=urlat, lat_0=lat_0,
#                         lon_0=lon_0)
#             m.drawstates()
#             m.shadedrelief()
#             x, y = m(self.argDict['lon'], self.argDict['lat'])
#             plt.plot(x, y, 'b.')
#
#         else:
#             lllat = 33.6
#             urlat = 37.3
#             lat_0 = 35.45
#             lon_0 = -97.5
#             ax_map = fig.add_subplot(gs2[3:5, 0])
#             m = Basemap(width=1600000, height=900000, projection='lcc',
#                         resolution='l', lat_1=lllat, lat_2=urlat, lat_0=lat_0,
#                         lon_0=lon_0)
#             m.drawstates()
#             #	m.drawcounties()
#             m.shadedrelief()
#             x, y = m(self.argDict['lon'], self.argDict['lat'])
#             plt.plot(x, y, 'b.')
#
#         # Data readings
#         ax_data = fig.add_subplot(gs2[2, 0])
#         plt.axis('off')
#         datastr = ('LCL: %.0f hPa, %.0f$^\circ$C\n' + \
#                    'Parcel Buoyancy: %.0f J kg$^{-1}$\n' + \
#                    # '0-%.0f m bulk shear: %.0f kts\n' + \
#                    '10 m T: %.0f$^\circ$C, Td: %.0f$^\circ$C') % \
#                   (self.argDict['lclpres'], self.argDict['lcltemp'],
#                    self.argDict['SBCAPE'].magnitude,
#                    # sampleHeights_m[-3], bulkshear,
#                    self.T[0], self.Td[0])
#         boxprops = dict(boxstyle='round', facecolor='none')
#         ax_data.text(0.5, 0.95, datastr, transform=ax_data.transAxes,
#                      fontsize=14,
#                      va='top', ha='center', bbox=boxprops)
#         # Legend for mesonet data
#         if self.argDict['loc'] not in ['K04V', 'CRES', 'MOFF', 'SAGF', 'K1V8']:
#             ax_data.legend(hand, lab, loc='upper center', \
#                            bbox_to_anchor=(0.5, 0.15), ncol=2, frameon=False)
#
#         ## Logos
#         ax_png = fig.add_subplot(gs2[5, 0])
#         img = mpimg.imread(fpath_logos, format='png')
#         plt.axis('off')
#         plt.imshow(img, aspect='equal')
#
#         fig.tight_layout()
#         return fig
#
# class meteogram():
#     def __init__(self, fmeso, tstart, tend, tsunrise):
#         mesodata = np.genfromtxt(fmeso, delimiter=',')
#         mesotimes_str = np.genfromtxt(fmeso, delimiter=',',
#                                       dtype=str, usecols=2)
#         self.fname = fmeso.split(os.sep)[-1]
#
#         self.p = (mesodata[:, 9] + 700.)
#         self.RH2m = mesodata[:, 4]
#         self.T2m = mesodata[:, 5]
#         self.u2m = mesodata[:, 15]
#         self.T9m = mesodata[:, 14]
#         self.u10m = mesodata[:, 12]
#         self.dir10m = mesodata[:, 13]
#         self.srad = mesodata[:, 6]
#
#         self.Rd = 287.
#         self.rho = self.p / ((self.T2m + 273.15) * self.Rd)
#         self.Td2m = np.array(mcalc.dewpoint_rh(self.T2m * units.degC,
#                                                self.RH2m / 100.))
#
#         self.mesotimes_dt = [datetime.strptime(t,
#                                                '"%Y-%m-%d %H:%M:%S"') for t in mesotimes_str]
#         self.mesotimes_t = mpdates.date2num(self.mesotimes_dt)
#
#         # self.tstart = tstart
#         # self.tend = tend
#         self.tsunrise = tsunrise
#
#         count = 0
#         for i in self.mesotimes_dt:
#             if (i.hour == int(tstart[0][:2])) & (i.minute == int(tstart[0][2:])):
#                 istart = count
#             elif (i.hour == int(tend[0][:2])) & (i.minute == int(tstart[0][2:])):
#                 iend = count
#             count += 1
#         self.irange = range(istart, iend)
#
#     def plot(self):
#         fig, axarr = plt.subplots(nrows=2, ncols=2, sharex='col',
#                                   figsize=(16, 9))
#         figtitle = '{0:s} Meteogram {1:s}'.format(self.fname.split('.')[1],
#                                                   self.fname.split('.')[0])
#         plt.suptitle(figtitle, fontsize=20)
#
#         # T & Td
#         axarr[0, 0].plot(self.mesotimes_t[self.irange], self.T2m[self.irange],
#                          color=(213./255, 94./255, 0), label='Temperature 2m',
#                          linewidth=2)
#         axarr[0, 0].plot(self.mesotimes_t[self.irange], self.T9m[self.irange],
#                          color=(204./255, 121./255, 167./255), label='Temperature 9m',
#                          linewidth=2)
#         axarr[0, 0].plot(self.mesotimes_t[self.irange], self.Td2m[self.irange],
#                          color=(0, 114./255, 178./255), label='Dewpoint Temperature',
#                          linewidth=2)
#         axarr[0, 0].set_title('Temperature and Dewpoint', fontsize=20)
#         axarr[0, 0].tick_params(labeltop=False, right=True, labelright=True, labelsize=16)
#         axarr[0, 0].set_ylabel('Temperature [$^\circ$C]', fontsize=18)
#         axarr[0, 0].grid(axis='y')
#         axarr[0, 0].legend(loc=0, fontsize=14)
#
#         # p
#         axarr[0, 1].plot(self.mesotimes_t[self.irange], self.p[self.irange], 'k', linewidth=2)
#         axarr[0, 1].set_title('Pressure', fontsize=20)
#         axarr[0, 1].tick_params(labeltop=False, right=True, labelright=True, labelsize=16)
#         axarr[0, 1].set_ylabel('Pressure [hPa]', fontsize=18)
#         axarr[0, 1].grid(axis='y')
#
#         # wind speed and direction
#         axarr_2 = axarr[1, 0].twinx()
#         axarr[1, 0].plot(self.mesotimes_t[self.irange], self.u10m[self.irange],
#                          color=(0, 158./255, 115./255), linewidth=2)
#         axarr_2.plot(self.mesotimes_t[self.irange], self.dir10m[self.irange],
#                      color=(213./255, 94./255, 0), marker='o', markersize=3,
#                      linestyle='')
#         axarr[1, 0].set_title('Wind Speed and Direction', fontsize=20)
#         axarr[1, 0].set_ylabel('Wind Speed [m s$^{-1}$]', color=(0, 158./255, 115./255), fontsize=18)
#         axarr_2.set_ylabel('Wind Direction [$^\circ$]', color=(213./255, 94./255, 0), fontsize=18)
#         axarr_2.set_yticks(range(0, 405, 45))
#         axarr_2.set_yticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N'])
#         axarr_2.tick_params(labelsize=16)
#         axarr[1, 0].grid(axis='y')
#         axarr[1, 0].xaxis.set_major_locator(mpdates.HourLocator(interval=1))
#         axarr[1, 0].xaxis.set_major_formatter(mpdates.DateFormatter('%H'))
#         axarr[1, 0].set_xlabel('Time [UTC]', fontsize=18)
#         axarr[1, 0].tick_params(labelsize=16)
#
#         # solar radiation
#         axarr[1, 1].plot(self.mesotimes_t[self.irange], self.srad[self.irange],
#                          color=(230./255, 159./255, 0), linewidth=2)
#         axarr[1, 1].set_title('Solar Radiation', fontsize=20)
#         axarr[1, 1].set_ylabel('Solar Radiation [W m$^{-2}$]', fontsize=18)
#         axarr[1, 1].xaxis.set_major_locator(mpdates.HourLocator(interval=1))
#         axarr[1, 1].xaxis.set_major_formatter(mpdates.DateFormatter('%H'))
#         axarr[1, 1].set_xlabel('Time [UTC]', fontsize=18)
#         axarr[1, 1].grid(axis='y')
#         axarr[1, 1].tick_params(labeltop=False, right=True, labelright=True, labelsize=16)
#
#         # vertical lines
#         # for i in range(len(axarr)):
#         #     for j in range(len(axarr)):
#         #         axarr[j, i].axvline(self.tstart, color='k', linestyle='--')
#         #         axarr[j, i].axvline(self.tend, color='k', linestyle='--')
#         #         axarr[j, i].axvline(self.tsunrise, color='r', linestyle='-')
#
#        return fig, axarr
'''
