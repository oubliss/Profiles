"""
Utils contains misc. functions to aid in data analysis.
"""
import sys
import os
import warnings
import requests
import numpy as np
from datetime import timedelta
from metpy.units import units as u
from scipy.signal import find_peaks

from .Coef_Manager import Coef_Manager
# QC moved to profiles.qc; re-exported so utils.qc keeps working.
from .qc import qc, _bias, _s_dev, _reject_outliers  # noqa: F401


package_path = os.path.dirname(os.path.abspath(__file__))

_coef_manager = None


def get_coef_manager():
    """ Return the shared Coef_Manager, constructing it on first use.

    Deprecated. Processing no longer goes through this: every lookup is made
    on the flight's own calibration source (FlightLog.calibration_source),
    so the coefficients applied follow the flight rather than process-global
    state. It remains for scripts that reach for it directly.

    Construction reads the coefficient tables off disk, so it is deferred
    until something actually needs coefficients. Importing this package must
    not require ~/.wxuas to exist, and callers (notably the test suite) must
    be able to point profiles.conf.coef_info at a different directory after
    import but before the first lookup.

    :rtype: profiles.Coef_Manager.Coef_Manager
    :return: the process-wide Coef_Manager
    """
    warnings.warn(
        'utils.get_coef_manager()/utils.coef_manager are deprecated; use '
        'the flight\'s calibration_source, or TableCalibration(directory).',
        DeprecationWarning, stacklevel=2)
    global _coef_manager
    if _coef_manager is None:
        _coef_manager = Coef_Manager()
    return _coef_manager


def reset_coef_manager():
    """ Discard the cached Coef_Manager so the next lookup rebuilds it.

    Needed when coef_info is repointed at a different coefficient directory
    part-way through a process, which the tests do.
    """
    global _coef_manager
    _coef_manager = None


def __getattr__(name):
    """ Keep ``utils.coef_manager`` working as a lazily-built module attribute.

    Several modules reach for ``utils.coef_manager`` directly. PEP 562 module
    __getattr__ lets that keep working without constructing the manager at
    import time.
    """
    if name == "coef_manager":
        return get_coef_manager()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

event_IDs = """
// DATA - event logging
#define DATA_AP_STATE                       7
// 8 was DATA_SYSTEM_TIME_SET
#define DATA_INIT_SIMPLE_BEARING            9
#define DATA_ARMED                          10
#define DATA_DISARMED                       11
#define DATA_AUTO_ARMED                     15
#define DATA_LAND_COMPLETE_MAYBE            17
#define DATA_LAND_COMPLETE                  18
#define DATA_NOT_LANDED                     28
#define DATA_LOST_GPS                       19
#define DATA_FLIP_START                     21
#define DATA_FLIP_END                       22
#define DATA_SET_HOME                       25
#define DATA_SET_SIMPLE_ON                  26
#define DATA_SET_SIMPLE_OFF                 27
#define DATA_SET_SUPERSIMPLE_ON             29
#define DATA_AUTOTUNE_INITIALISED           30
#define DATA_AUTOTUNE_OFF                   31
#define DATA_AUTOTUNE_RESTART               32
#define DATA_AUTOTUNE_SUCCESS               33
#define DATA_AUTOTUNE_FAILED                34
#define DATA_AUTOTUNE_REACHED_LIMIT         35
#define DATA_AUTOTUNE_PILOT_TESTING         36
#define DATA_AUTOTUNE_SAVEDGAINS            37
#define DATA_SAVE_TRIM                      38
#define DATA_SAVEWP_ADD_WP                  39
#define DATA_FENCE_ENABLE                   41
#define DATA_FENCE_DISABLE                  42
#define DATA_ACRO_TRAINER_DISABLED          43
#define DATA_ACRO_TRAINER_LEVELING          44
#define DATA_ACRO_TRAINER_LIMITED           45
#define DATA_GRIPPER_GRAB                   46
#define DATA_GRIPPER_RELEASE                47
#define DATA_PARACHUTE_DISABLED             49
#define DATA_PARACHUTE_ENABLED              50
#define DATA_PARACHUTE_RELEASED             51
#define DATA_LANDING_GEAR_DEPLOYED          52
#define DATA_LANDING_GEAR_RETRACTED         53
#define DATA_MOTORS_EMERGENCY_STOPPED       54
#define DATA_MOTORS_EMERGENCY_STOP_CLEARED  55
#define DATA_MOTORS_INTERLOCK_DISABLED      56
#define DATA_MOTORS_INTERLOCK_ENABLED       57
#define DATA_ROTOR_RUNUP_COMPLETE           58  // Heli only
#define DATA_ROTOR_SPEED_BELOW_CRITICAL     59  // Heli only
#define DATA_EKF_ALT_RESET                  60
#define DATA_LAND_CANCELLED_BY_PILOT        61
#define DATA_EKF_YAW_RESET                  62
#define DATA_AVOIDANCE_ADSB_ENABLE          63
#define DATA_AVOIDANCE_ADSB_DISABLE         64
#define DATA_AVOIDANCE_PROXIMITY_ENABLE     65
#define DATA_AVOIDANCE_PROXIMITY_DISABLE    66
#define DATA_GPS_PRIMARY_CHANGED            67
#define DATA_WINCH_RELAXED                  68
#define DATA_WINCH_LENGTH_CONTROL           69
#define DATA_WINCH_RATE_CONTROL             70
"""


NC_LEVELS = ('low', 'none')


def writes_netcdf(nc_level):
    """ Should this nc_level write per-object NetCDF files?

    The three classes used to disagree: Raw_Profile tested ``== 'low'`` while
    Thermo_Profile and Wind_Profile tested ``is not None``. The documented
    "write nothing" value, the string 'none', is truthy, so passing it wrote
    thermo_ and wind_ files while suppressing the raw one.

    :param nc_level: 'low' to write, None or 'none' to write nothing
    :rtype: bool
    """
    if nc_level is None:
        return False

    normalised = str(nc_level).strip().lower()
    if normalised not in NC_LEVELS:
        warnings.warn(
            f"unrecognised nc_level {nc_level!r}; expected one of "
            f"{NC_LEVELS} or None. Treating it as 'none' and writing no "
            f"NetCDF files.", stacklevel=2)
        return False

    return normalised == 'low'


def nearest_index(times, target):
    """ Index of the sample in ``times`` closest in time to ``target``.

    Leg times come off the GPS clock but the barometer, thermistors and
    humidity sensors each keep their own, so an exact ``list.index`` lookup
    only ever worked for the clock the leg was found on.

    :param sequence<datetime> times: sample times, non-decreasing
    :param datetime target: the time to find
    :rtype: int
    """
    stamps = np.asarray(times, dtype='datetime64[us]')
    want = np.datetime64(target, 'us')
    after = int(np.searchsorted(stamps, want, side='left'))
    if after == 0:
        return 0
    if after >= len(stamps):
        return len(stamps) - 1
    # side='left' lands on the first sample at or after the target, so an
    # exact match is returned as such (the old .index() behaviour).
    if want - stamps[after - 1] < stamps[after] - want:
        return after - 1
    return after


def leg_extents(legs, alts, alt_times):
    """ Vertical extent of each leg in each direction.

    :param list<tuple> legs: (start, peak, end) times
    :param np.Array<float> alts: altitudes in metres
    :param np.Array<Datetime> alt_times: times corresponding to alts
    :rtype: list<tuple>
    :return: (rise, fall) per leg, metres: peak minus start, peak minus end
    """
    alts = np.asarray(alts, dtype=float)
    out = []
    for start, peak, end in legs:
        top = alts[nearest_index(alt_times, peak)]
        out.append((top - alts[nearest_index(alt_times, start)],
                    top - alts[nearest_index(alt_times, end)]))
    return out


def filter_legs_by_extent(legs, alts, alt_times, min_extent, ascent=True):
    """ Drop legs that do not climb (or descend) at least ``min_extent``.

    Peak detection with a one-metre prominence reports every wiggle at the
    top of a real profile as a profile of its own. Because the legs are
    valley-peak-valley triples, the wiggle is a *rise* of about a metre but
    its "end" is the bottom of the real descent, so it is a perfectly good
    descent. The extent therefore has to be measured in the direction being
    processed.

    :param float min_extent: metres; 0 or None keeps everything
    :param bool ascent: measure the rise (start to peak) if True, the fall
       (peak to end) if False, or either if None
    :rtype: list<tuple>
    """
    if not min_extent or not legs:
        return list(legs)

    kept = []
    for leg, (rise, fall) in zip(legs, leg_extents(legs, alts, alt_times)):
        if ascent is None:
            extent = max(rise, fall)
        else:
            extent = rise if ascent else fall
        if extent >= min_extent:
            kept.append(leg)
    return kept


def _times_at_levels(values, times, first, last, levels):
    """ Time at which a series first reaches each level, walking forward.

    The walk always advances one sample past each hit, so successive levels
    get distinct samples even where the series is flat or noisy.
    """
    found = []
    i = first
    for elem in levels:
        # Bound check first: i can be advanced past last by the increment
        # below, and `and` does not short-circuit a subscript written on
        # its left, so the original order raised IndexError when a profile
        # ran to the last sample in the file.
        while i < last and values[i] < elem:
            i += 1
        # i can only run past the end of a leg that is shorter than the
        # grid asked of it (a common base_start on a short leg); those
        # levels get the last sample, which leaves their bins empty.
        found.append(times[min(i, len(times) - 1)])
        i += 1
    return found


def regrid_base(base=None, base_times=None, new_res=None, ascent=True,
                units=None, indices=(None, None), base_start=None):
    """ Calculates times at which data means should be calculated.

    The vertical coordinate is altitude, or negative pressure so that it
    still increases with height. For an ascent everything comes back in
    ascending order. For a descent everything comes back in *flight order*:
    the first level is the top and altitude decreases, so that the edge times
    still increase along the leg and the (start, end] time bins used by
    :func:`regrid_data` keep working. Levels sit on the same
    ``base_start + n*res`` lattice either way, so an ascent and a descent can
    share a grid.

    :param np.Array<Quantity> base: Measurements of the variable serving as \
       the vertical coordinate
    :param np.Array<Datetime> base_times: Times coresponding to base
    :param Quantity new_res: The resolution to which base should be gridded. \
       This must have the same dimension (i.e. both length or both pressure) \
       as base.
    :param bool ascent: True if data from ascending leg of profile is to be \
       analyzed, false if descending
    :param pint.UnitRegistry units: The unit registry defined in Profile
    :param tuple indices: (start, end) times of the leg being gridded - \
       (start, peak) for an ascent, (peak, end) for a descent. They are \
       matched to the nearest sample on base's own clock, which need not be \
       the clock they were detected on. A (start, peak, end) triple is also \
       accepted and the leg picked by ``ascent``.
    :param Quantity base_start: the lowest edge of the grid, in the units of \
       base. For a pressure grid this is the highest pressure.
    :rtype: tuple(np.Array<Datetime>, np.Array<Quantity>)
    :return: times at which the craft is at vertical points n*res above \
       the profile starting height and the corrosponding base values
    """
    is_pressure = new_res.dimensionality == units.Pa.dimensionality

    if base_start is not None and (base_start.dimensionality
                                   != new_res.dimensionality):
        raise ValueError(
            f'base_start {base_start} is not in the dimension of the grid '
            f'resolution {new_res}')

    # np.arange below steps in bare magnitudes, so everything has to be in
    # base's units first - a resolution of 5 hPa against a barometer in Pa
    # stepped 5 Pa.
    new_res = new_res.to(base.units)
    if base_start is not None:
        base_start = base_start.to(base.units)

    # Change indices to a 2-tuple with sample indices instead of times
    if indices[0] is None:
        indices = (0, len(base) - 1)
    else:
        if len(indices) == 3:
            indices = (indices[0], indices[1]) if ascent \
                else (indices[1], indices[2])
        indices = (nearest_index(base_times, indices[0]),
                   nearest_index(base_times, indices[1]))

    # Use negative pressure so that the vertical coordinate increases upward
    if is_pressure:
        base = -1*base
        if base_start is not None:
            base_start = -1*base_start

    # bottom and top of the leg, in sample indices
    lowest, highest = (indices[0], indices[1]) if ascent \
        else (indices[1], indices[0])

    # Regrid base
    if base_start is None:
        floor = base[lowest]
    else:
        floor = base_start
    ceiling = base[highest]

    if base_start is None:
        new_base = np.arange((floor + 0.5*new_res).magnitude,
                             (ceiling - 0.5*new_res).magnitude,
                             new_res.magnitude)
        base_edges = np.arange(floor.magnitude, ceiling.magnitude,
                               new_res.magnitude)
    else:
        new_base = np.arange((base_start + 0.5*new_res).magnitude,
                             (ceiling - 0.5 * new_res).magnitude,
                             new_res.magnitude)
        base_edges = np.arange(base_start.magnitude, ceiling.magnitude,
                               new_res.magnitude)

    new_base = np.array(new_base) * base.units
    base_edges = np.array(base_edges) * base.units

    # Find the times where the levels and edges occur in the profile
    if ascent:
        new_times = _times_at_levels(base, base_times, indices[0],
                                     indices[1], new_base)
        time_edges = _times_at_levels(base, base_times, indices[0],
                                      indices[1], base_edges)
    else:
        # Walk the leg backwards in time, which is an ascent in the
        # vertical coordinate, then put everything back in flight order.
        span = slice(indices[0], indices[1] + 1)
        flipped_base = base[span][::-1]
        flipped_times = list(base_times[span])[::-1]
        last = len(flipped_times) - 1
        new_times = _times_at_levels(flipped_base, flipped_times, 0, last,
                                     new_base)[::-1]
        time_edges = _times_at_levels(flipped_base, flipped_times, 0, last,
                                      base_edges)[::-1]
        new_base = new_base[::-1]
        base_edges = base_edges[::-1]

    if is_pressure:
        new_base = -1*new_base
        base_edges = -1*base_edges

    # Remove duplicates:
    # new_times, indices = np.unique(new_times, return_index=True)
    # new_base = new_base[indices]
    return (new_times, new_base, time_edges, base_edges)


def _bin_indices(data_times, gridded_times):
    """ Yield, for each bin between consecutive gridded_times, the indices of
    the samples in it.

    This is the one binning rule used for every gridded variable, so they
    all have ``len(gridded_times) - 1`` bins with the same boundaries: bin i
    is the half-open interval (gridded_times[i], gridded_times[i+1]]. A bin
    with no samples yields an empty index array rather than being skipped.

    :param np.Array<Datetime> data_times: Times coresponding to the data
    :param np.Array<Datetime> gridded_times: The times returned by regrid_base
    :rtype: generator of np.Array<int>
    """
    # Built once, not once per bin.
    times = np.array(data_times)
    for i in range(len(gridded_times) - 1):
        in_bin = (times > gridded_times[i]) & (times <= gridded_times[i + 1])
        yield np.where(in_bin)[0]


def regrid_data(data=None, data_times=None, gridded_times=None, units=None):
    """ Returns data interpolated to an evenly spaced array based on
    gridded_times.

    Bins are (gridded_times[i], gridded_times[i+1]]; a bin containing no
    samples is NaN.

    :param np.Array<Quantity> data: a non-base variable (i.e. not yor chosen \
       vertical coordinate)
    :param np.Array<Datetime> data_times: Times coresponding to data
    :param pint.UnitRegistry units: The unit registry defined in Profile
    :param np.Array<Datetime> gridded_times: The times returned by regrid_base
    :rtype: np.Array<Quantity>
    :return: gridded_data
    """
    magnitudes = data.magnitude
    gridded_data = []
    # An all-NaN or empty bin is NaN by design; the mean of nothing warns.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for foo in _bin_indices(data_times, gridded_times):
            try:
                gridded_data.append(np.nanmean(magnitudes[foo]))
            except IndexError:
                raise ValueError(
                    "The data time array should be the same length as the "
                    "data itself") from None

    return np.array(gridded_data) * data.units


def regrid_data_group(data=None, data_times=None, gridded_times=None, units=None):
    """ Yield the samples in each bin between consecutive gridded_times.

    Uses the same bins as :func:`regrid_data` - (start, end], one per pair
    of gridded times, empty bins included - so a variable gridded from these
    groups has the same length and boundaries as every other variable.

    :param np.Array<Quantity> data: a non-base variable (i.e. not yor chosen \
       vertical coordinate)
    :param np.Array<Datetime> data_times: Times coresponding to data
    :param pint.UnitRegistry units: The unit registry defined in Profile
    :param np.Array<Datetime> gridded_times: The times returned by regrid_base
    :rtype: generator of dict
    :return: dicts with start_time, end_time and the bin's values
    """
    for i, foo in enumerate(_bin_indices(data_times, gridded_times)):
        yield {
            'start_time': gridded_times[i],
            'end_time': gridded_times[i + 1],
            'values': data[foo]
        }


def regrid_data_group_interp(data=None, data_times=None, gridded_times=None, units=None):
    """ Returns data interpolated to an evenly spaced array based on
    gridded_times.

    :param np.Array<Quantity> data: a non-base variable (i.e. not yor chosen \
       vertical coordinate)
    :param np.Array<Datetime> data_times: Times coresponding to data
    :param pint.UnitRegistry units: The unit registry defined in Profile
    :param np.Array<Datetime> gridded_times: The times returned by regrid_base
    :rtype: np.Array<Quantity>
    :return: gridded_data
    """

    #
    # Average around selected points
    #
    data_index = 0  # This tracks the most recent data element processed

    for i in range(len(gridded_times)):
        #
        # Find the data indices in the specified time range
        #
        val = gridded_times[i]
        data_seg_start_ind = None
        data_seg_end_ind = None

        while data_index < len(data):
            if data_times[data_index] <= val:
                data_seg_start_ind = data_index
            else:
                data_seg_end_ind = data_index
                break
            data_index += 1

        # Calculate and store the segment mean
        if data_seg_start_ind is not None and data_seg_end_ind is not None:
            yield {
                'key': val,
                'values': data[data_seg_start_ind:data_seg_end_ind+1]
            }


def steinhart_hart(resistance, coefs):
    """ Resistance to temperature with one coefficient row.

    :param list<Quantity> resistance: resistances recorded by a thermistor
    :param dict coefs: a MasterCoefList row with A, B and C
    :rtype: list<Quantity>
    :return: list of temperatures in K
    """
    a = float(coefs["A"])
    b = float(coefs["B"])
    c = float(coefs["C"])

    # A zero or negative resistance (a dead channel) is inf/NaN by design.
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.power(np.add(np.add(b * np.log(resistance), a),
                        c * np.power(np.log(resistance), 3)), -1)


def temp_calib(resistance, sn, source=None, when=None):
    """ Converts resistance to temperature using the coefficients for the \
       sensor specified.

    Prefer profiles.calibration.calibrate_temperature, which takes the
    flight's calibration source. When no ``source`` is given this falls back
    to the deprecated process-wide manager, which is how it always behaved.

    :param list<Quantity> resistance: resistances recorded by temperature \
       sensors
    :param int sn: the serial number of the sensor reporting
    :param source: a CalibrationSource to look the coefficients up in
    :param when: flight time, to select among dated coefficient rows
    :rtype: list<Quantity>
    :return: list of temperatures in K
    """
    if source is None:
        source = get_coef_manager()
    return steinhart_hart(resistance, source.get_coefs("Imet", sn, when=when))


def rh_calib(raw, sn):
    """ Return RH as reported. No per-sensor correction is applied.

    This is deliberately a pass-through, and has been in effect since well
    before 1.4.0 - the previous implementation looked up the sensor's 'A'
    coefficient, divided it by 1000, and then unconditionally overwrote the
    result with 0 on the next line, so no offset ever reached the data. That
    has been made explicit here rather than left looking like a live
    calculation.

    Reinstating a correction is a scientific decision, not a cleanup, and
    needs three things settled first:

    * RH rows in MasterCoefList carry equation E3 with *two* coefficients
      (A and B); a single additive offset cannot be the whole model.
    * The /1000 scaling does not match the stored magnitudes. For the
      sensors on the reference flight A is 9.10E-02, 1.80E-01 and 9.30E-02,
      which after dividing by 1000 would shift RH by ~1e-4 %, i.e. nothing.
    * Whatever is chosen has to be recorded in the output, or files
      processed before and after become indistinguishable.

    :param list<Quantity> raw: raw RH
    :param int sn: serial number of the humidity sensor (currently unused)
    :rtype: list<Quantity>
    :return: raw, unchanged
    """
    return raw


def get_place_from_lat_lon(lat, lon, zoom=16, timeout=5):
    """
    Pings the nominatim API to get a place name from a lat/lon pair.

    Example zoom levels:
     - 16: Street level
     - 12: Town Level
     - 8: County Level

    :param lat:
    :param lon:
    :param zoom:
    :param timeout: seconds to wait before giving up
    :return:
    """

    url = "https://nominatim.openstreetmap.org/reverse.php"

    params = {
        'lat': lat,
        'lon': lon,
        'zoom': zoom,
        'format': 'jsonv2'
    }

    # format query string and return query value
    # Nominatim's usage policy requires an identifying User-Agent.
    from profiles import __version__
    headers = {'User-Agent': f'profiles-uas/{__version__}'}
    result = requests.get(url, params, headers=headers, timeout=timeout)
    result.raise_for_status()

    return ','.join(result.json()['display_name'].split(',')[:-1])


def ned2body(xned, yned, zned, roll, pitch, yaw):
    croll = np.cos(roll)
    sroll = np.sin(roll)
    cpitch = np.cos(pitch)
    spitch = np.sin(pitch)
    cyaw = np.cos(yaw)
    syaw = np.sin(yaw)

    Rx = np.array([[1, 0, 0],
                    [0, croll, sroll],
                    [0, -sroll, croll]])
    Ry = np.array([[cpitch, 0, -spitch],
                    [0, 1, 0],
                    [spitch, 0, cpitch]])
    Rz = np.array([[cyaw, -syaw, 0],
                    [syaw, cyaw, 0],
                    [0, 0, 1]])

    R = Rz @ Ry @ Rx

    a = R @ np.array([xned, yned, zned])

    return float(a[0]), float(a[1]), float(a[2])


def body2ned(xb, yb, zb, roll, pitch, yaw):
    croll = np.cos(roll)
    sroll = np.sin(roll)
    cpitch = np.cos(pitch)
    spitch = np.sin(pitch)
    cyaw = np.cos(yaw)
    syaw = np.sin(yaw)

    Rx = np.array([[croll, -sroll, 0],
                    [sroll, croll, 0],
                    [0, 0, 1]])
    Ry = np.array([[cpitch, 0, spitch],
                    [0, 1, 0],
                    [-spitch, 0, cpitch]])
    Rz = np.array([[1, 0, 0],
                    [0, cyaw, -syaw],
                    [0, syaw, cyaw]])

    R = Rx @ Ry @ Rz

    a = R @ np.array([xb, yb, zb])

    return float(a[0]), float(a[1]), float(a[2])


def identify_profile_peaks(alts, alt_times, window=None,
                           confirm_bounds=False, **kwargs):
    """
    New method for identifying the individual flight legs of a CopterSonde flight without prior knowledge
    of the surface altitude. This uses scipy.signal.find_peaks and defaults to a prominence of 1 as an input.
    Other inputs to find_peaks are available via the kwargs.

    :param np.Array<float> alts: Altitude of the UAS
    :param np.Array<datetime> alt_times: UTC time at each altitude
    :param tuple<datetime> window: Tuple of two times that should be considered the min max of the data to be considered
    :param bool confirm_bounds: Make a quick figure that shows the identified peaks and valleys as a sanity check
    :param kwargs: Kwargs for scipy.signal.find_peaks
    :return: List of tuples with the start, peak, and end time of each identified profile
    """

    # Run find peaks on the altitudes to find the peaks and valleys
    peaks, foo = find_peaks(alts, prominence=1, **kwargs)
    valleys, foo = find_peaks(-1*alts, prominence=1, **kwargs)

    if window is not None:
        # Convert the window from times to indices to make things easy later on
        window_min = np.argmin(np.abs(window[0] - np.array(alt_times)))
        window_max = np.argmin(np.abs(window[1] - np.array(alt_times)))
        window = (window_min, window_max)

        # Concat these arrays and sort.
        p_and_v = np.sort(np.concatenate((window, peaks, valleys)))

        # Remove any indices outside of the specified window, if provide
        foo = np.where((p_and_v >= window[0]) & (p_and_v <= window[1]))
        p_and_v = p_and_v[foo]

    else:
        p_and_v = np.sort(np.concatenate((peaks, valleys)))

    # Check the bounds if desired
    if confirm_bounds:
        import matplotlib.pyplot as plt
        # User verifies selection
        fig2 = plt.figure()
        plt.plot(alts, figure=fig2)
        plt.grid(axis="y", which="both", figure=fig2)
        plt.vlines(p_and_v, min(alts) - 10, max(alts) + 10)

        plt.show(block=True)

    # Function to format the profiles into the expected data structure for Profiles
    def _format_profiles(start, inds):
        try:
            foo = (alt_times[inds[start]], alt_times[inds[start+1]], alt_times[inds[start+2]])
            result = _format_profiles(start+2, inds)
            return [foo] + result
        except IndexError:
            return []

    return _format_profiles(0, p_and_v)


def identify_profile(alts, alt_times, confirm_bounds=True,
                     profile_start_height=None, to_return=None, ind=0):
    """ Identifies the temporal bounds of all profiles in the data file. These
    assumptions must be valid:
    * The craft starts and ends each profile below profile_start_height
    * The craft does not go above profile_start_height until the first
    profile is started
    * The craft does not go above profile_start_height after the last
    profile is ended.

    :param np.Array<Quantity> alts: recorded altitudes; units don't matter
    :param np.Array<Datetime> alt_times: times coresponding to alts
    :param bool confirm_bounds: if True, will ask user for verification that \
       the start, peak, and end times of the profile have been properly \
       identified
    :param Quantity profile_start_height: if this is given, the user will not be \
       prompted to enter a start height for each profile. This is recommended \
       when processing many profiles from the same mission. At least one \
       profile should be processed without this option to determine the correct\
       value.
    :param int ind: used privately for recurrsion - leave this alone
    :param list to_return: used privately for recurrsion - leave this alone
    :rtype: list<tuple>
    :return: a list of times defining the profiles in the format \
       (time_start, time_max_height, time_end)
    """

    if to_return is None:
        to_return = []

    isDone = False
    # Get the starting height from the user
    if profile_start_height is None:
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from pandas.plotting import register_matplotlib_converters
        register_matplotlib_converters()
        fig1 = plt.figure()
        plt.plot(alt_times, alts, figure=fig1)
        plt.grid(axis="y", which="both", figure=fig1)

        myFmt = mdates.DateFormatter('%M')
        fig1.gca().xaxis.set_major_formatter(myFmt)

        plt.show(block=False)

        try:
            profile_start_height = int(input('Wrong file? Enter "q" to quit. '
                                             + '\nProfile start height: ')) * \
                                             alts.units
        except ValueError:
            sys.exit(0)
        plt.close()

    # If no profiles exist after ind, return an empty index list.
    if max(alts[ind:]) < profile_start_height:
        return []

    # Declare variables used in loop
    start_ind_asc = None
    end_ind_des = len(alts) - 12
    peak_ind = None

    # Check through valid alts for start_ind_asc, peak_ind, end_ind_des in
    # that order
    while ind < len(alts) - 10:
            if(start_ind_asc is None):  # should be ascending
                # finds at what time the altitude range is reached going up

                # Error if starts on a descent
                if(alts[ind] > profile_start_height and
                   alts[ind + 10] < alts[ind]):
                    print("Error separating profiles: start height is first reached on a descent")
                    break

                # Set start_ind_asc when the craft is first above
                # profile_start_height
                if alts[ind] > profile_start_height:
                    start_ind_asc = ind

                ind += 1

            elif(end_ind_des == len(alts) - 12):  # should be descending
                if end_ind_des == ind:
                    end_ind_des += 1

                # Set end_ind_des when the craft is again below
                # profile_start_height for the first time since start_ind_asc
                if alts[ind] < profile_start_height:
                    end_ind_des = ind

                    # Now that the bounds of the profile have been found, we
                    # find the index of the maximum altitude.
                    peak_ind = list(alts).index(np.nanmax(alts[start_ind_asc:end_ind_des]),
                                                start_ind_asc, end_ind_des)

                ind += 1

            # The current profile has been processed; we just need to check
            # if there are more profiles in the file
            else:
                if peak_ind is None:
                    peak_ind = list(alts).index(np.nanmax(alts[start_ind_asc:end_ind_des]),
                                                start_ind_asc, end_ind_des)
                if confirm_bounds:
                    import matplotlib.pyplot as plt
                    # User verifies selection
                    fig2 = plt.figure()
                    plt.plot(range(len(alt_times)), alts, figure=fig2)
                    plt.grid(axis="y", which="both", figure=fig2)
                    plt.vlines([start_ind_asc, peak_ind, end_ind_des],
                               min(alts).magnitude - 50, max(alts.magnitude) + 50)

                    plt.show(block=False)

                    # Get user opinion
                    valid = input('Correct? (Y/n): ')
                    # If good, wrap up the profile
                    if valid in "yYyesYes" or valid == "":
                        plt.close()
                        if end_ind_des is None:
                            print("Could not find end time des (LineTag B)")
                            return
                        isDone = True
                        break
                    elif valid in "nNnoNo":
                        plt.close()
                        # Re-ask for the start height: that is the one
                        # thing the user can change to get a different
                        # answer. to_return used to land positionally in
                        # confirm_bounds.
                        to_return = identify_profile(
                            alts, alt_times, confirm_bounds=confirm_bounds,
                            profile_start_height=None, to_return=to_return)
                    else:
                        print("Invalid choice. Re-selecting profile...")
                        plt.close()
                        # Re-ask for the start height: that is the one
                        # thing the user can change to get a different
                        # answer. to_return used to land positionally in
                        # confirm_bounds.
                        to_return = identify_profile(
                            alts, alt_times, confirm_bounds=confirm_bounds,
                            profile_start_height=None, to_return=to_return)
                else:
                    isDone = True
                    break

                ind += 1

    if(isDone):
        # Add the profile if it is not already in to_return
        pending_profile = (alt_times[start_ind_asc],
                           alt_times[peak_ind],
                           alt_times[end_ind_des])
        if not _profile_in(pending_profile, to_return):
            to_return.append((alt_times[start_ind_asc],
                              alt_times[peak_ind],
                              alt_times[end_ind_des]))

            print("Profile from ", alt_times[start_ind_asc],
                  "to", alt_times[end_ind_des], "added")
        # Check if more profiles in file
        if ind + 100 < len(alts) \
           and max(alts[ind + 100::]) > profile_start_height:
            # There is another profile before the end of the
            # file - find it.
            a = profile_start_height
            to_return =\
                identify_profile(alts, alt_times,
                                 confirm_bounds=confirm_bounds,
                                 profile_start_height=a,
                                 to_return=to_return,
                                 ind=ind+100)
    return to_return


def _profile_in(indices, all_indices):
    """ Helper function for identify_profile to ensure similar or overlapping
    profiles not included

    :param tuple indices: the identifying tuple for the profile to look for
    :param list<tuple> all_indices: identifying tuples for all included profiles
    """
    for profile_n in all_indices:
        # Check for similar
        if ((profile_n[0]-indices[0] <= timedelta(seconds=5)
           and indices[0]-profile_n[0] <= timedelta(seconds=5)) or
           (profile_n[1]-indices[1] <= timedelta(seconds=5)
           and indices[1]-profile_n[1] <= timedelta(seconds=5)) or
           (profile_n[2]-indices[2] <= timedelta(seconds=5)
           and indices[2]-profile_n[2] <= timedelta(seconds=5))):
            return True
        # Check for overlapping
        if indices[1] > profile_n[0] and indices[1] < profile_n[2]:
            return True
    return False
