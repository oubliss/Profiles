"""
Utils contains misc. functions to aid in data analysis.
"""
import sys
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
from pandas.plotting import register_matplotlib_converters
from pint import UnitStrippedWarning
from metpy.units import units as u

from .Coef_Manager import Coef_Manager


package_path = os.path.dirname(os.path.abspath(__file__))
coef_manager = Coef_Manager()  # All required input is given in __init__.py

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("error", category=UnitStrippedWarning)
register_matplotlib_converters()

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

def regrid_base(base=None, base_times=None, new_res=None, ascent=True,
                units=None, indices=(None, None), base_start=None):
    """ Calculates times at which data means should be calculated.

    :param np.Array<Quantity> base: Measurements of the variable serving as \
       the vertical coordinate
    :param np.Array<Datetime> base_times: Times coresponding to base
    :param Quantity new_res: The resolution to which base should be gridded. \
       This must have the same dimension (i.e. both length or both pressure) \
       as base.
    :param bool ascent: True if data from ascending leg of profile is to be \
       analyzed, false if descending
    :param pint.UnitRegistry units: The unit registry defined in Profile
    :param tuple indices: start and end times
    :param Quantity base_start: lowest altitude value of gridded_base
    :rtype: tuple(np.Array<Datetime>, np.Array<Quantity>)
    :return: times at which the craft is at vertical points n*res above \
       the profile starting height and the corrosponding base values
    """
    # Use negative pressure so that the max of the data list is the peak

    # Change indices to a 2-tuple with indices instead of times, start and end
    if indices[0] is None:
        indices = (0, len(base))
    else:
        a = list(base_times).index(indices[0])
        if ascent:
            b = list(base_times).index(indices[1])
        else:
            b = list(base_times).index(indices[2])
        indices = (a, b)

    if new_res.dimensionality == units.Pa.dimensionality:
        base = -1*base
        if base_start is not None:
            base_start = -1*base_start

    # Regrid base
    if base_start is None:
        new_base = np.arange((base[indices[0]] + 0.5*new_res).magnitude,
                             (base[indices[1]] - 0.5*new_res).magnitude,
                             new_res.magnitude)
    else:
        new_base = np.arange(base_start.magnitude,
                             (base[indices[1]] - 0.5 * new_res).magnitude,
                             new_res.magnitude)

    new_base = np.array(new_base) * base.units

    ind_in_grid = []
    i = indices[0]
    for elem in new_base:
        while base[i] < elem and i < indices[1]:
            i += 1
        ind_in_grid.append(i)
        i += 1

    new_times = [base_times[i] for i in ind_in_grid]


    if new_res.dimensionality == units.Pa.dimensionality:
        new_base = -1*new_base

    # Remove duplicates:
    # new_times, indices = np.unique(new_times, return_index=True)
    # new_base = new_base[indices]
    return (new_times, new_base)


def regrid_data(data=None, data_times=None, gridded_times=None, units=None):
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

    gridded_data = []
    for i in range(len(gridded_times)-1):
        #
        # Find the data indices in the specified time range
        #
        start_time = gridded_times[i]
        end_time = gridded_times[i+1]
        data_seg_start_ind = None
        data_seg_end_ind = None

        while data_index < len(data):
            if data_times[data_index] >= start_time:
                data_seg_start_ind = data_index
                break
            data_index += 1

        while data_index < len(data):
            if data_times[data_index] >= end_time:
                data_seg_end_ind = data_index
                break
            data_index += 1

        # Calculate and store the segment mean
        if data_seg_start_ind is not None and data_seg_end_ind is not None:
            gridded_data.append(np.nanmean(data.magnitude[data_seg_start_ind:
                                           data_seg_end_ind]))

    gridded_data = np.array(gridded_data) * data.units

    return (gridded_data)



def temp_calib(resistance, sn):
    """ Converts resistance to temperature using the coefficients for the \
       sensor specified OR generalized coefficients if the serial number (sn)\
       is not recognized.

    :param list<Quantity> resistance: resistances recorded by temperature \
       sensors
    :param int sn: the serial number of the sensor reporting
    :rtype: list<Quantity>
    :return: list of temperatures in K
    """
    coefs = coef_manager.get_coefs("Imet", sn)
    a = float(coefs["A"])
    b = float(coefs["B"])
    c = float(coefs["C"])

    return np.power(np.add(np.add(b * np.log(resistance), a),
                    c * np.power(np.log(resistance), 3)), -1)


def rh_calib(raw, sn):
    """ Adds the sensor offsets

    :param list<Quanitity> raw: raw RH
    :param int sn: serial number of the humidity sensor
    :rtype: list<Quantity>
    :return: list of calibrated rh
    """
    offset = coef_manager.get_coefs('RH', sn)['A']
    try:
        offset = float(coef_manager.get_coefs('RH', sn)['A']) / 1000
    except Exception:
        offset = 0

    return np.add(raw, offset)


def qc(data, max_bias, max_variance):
    """ Determines which sensors are not reliable from a given set. Be sure
       to only include like sensors (not both temperature inside and outside
                                     the CO2 sensor) in Data.

    :param list<Quantity> data: a list containing one list for each sensor
       in the ensemble, i.e. all external RH sensors
    :param Quantity max_bias: the maximum absolute difference between the \
       mean of one sensor and the mean of all sensors of that type. This \
       should be determined experimentally for each type of sensor.
    :param Quantity max_variance: the maximum absolute difference between the \
       standard deviation of one sensor and the standard deviation of all \
       sensors of that type. This should be determined experimentally for \
       each type of sensor.
    :rtype: list<int> of length len(data)
    :return: list containing 0 in the position of each "good" sensor, 2 in the
       position of each sensor flagged for bias, 3 in the position of each
       sensor flagged for response time, and 4 in the position of each flagged as empty
    """

    if isinstance(data, u.Quantity):
        data = data.magnitude

    good_nonempty = [1] * len(data)
    for i in range(len(data)):
        if np.nanmean(data[i]) == 0:
            good_nonempty[i] = 4
        else:
            good_nonempty[i] = 0

    # _bias: returns list of length number of sensors; 0 means data is good
    good_means = _bias(data, max_bias)
    # _s_dev: returns list of length number of sensors; 0 means data is good
    good_sdevs = _s_dev(data, max_variance)

    combined_sensor_flags = [1] * len(data)

    # Combine good_means and good_sdevs, leaving 0 only where the sensor
    # passed both tests.
    for i in range(len(data)):
        combined_sensor_flags[i] = max([good_means[i], good_sdevs[i], good_nonempty[i]])

    return combined_sensor_flags


def _bias(data, max_abs_error):
    """ This method identifies sensors with excessive biases and returns a
    list flagging sensors determined to be questionable.

    :param np.Array<Quantity> data: a list containing one list for each sensor
       in the ensemble, i.e. all external RH sensors
    :param Quantity max_abs_error: sensors with means more than
       max_abs_error from the mean of sensor means will be flagged
    :rtype: list of length len(data)
    :return: list containing 0s by default and 2 in the position of each sensor
       flagged for bias.
    """

    to_return = np.zeros(len(data))
    # Calculate the mean of each sensor
    means = np.zeros(len(data))
    for i in range(len(data)):
        means[i] = np.nanmean(data[i])

    while(True):
        # Identify the sensor with the mean farthest from the mean of means
        max_diff = 0
        furthest_from_mean = 0  # index of sensor furthest from mean

        for j in range(len(data)):

            if(np.abs(np.nanmean(means)-means[j]) >
               np.abs(np.nanmean(means)-means[furthest_from_mean])):
                furthest_from_mean = j

            for k in range(len(data)):
                if(np.abs(means[j]-means[k]) > max_diff):
                    max_diff = np.abs(means[j]-means[k])

            # If the furthest sensor is farther than max_abs_error from the
            # mean of means, eliminate it and perform analysis again.
            if(max_diff > max_abs_error):
                to_return[furthest_from_mean] = 2
                means[furthest_from_mean] = np.NaN
            else:
                return to_return


def _s_dev(data, max_abs_error):
    """ This method identifies sensors with excessively low or high
    variabilities and returns a list flagging sensors determined to be
    questionable.

    :param np.Array<Quantity> data: a list containing one list for each sensor
       in the ensemble, i.e. all external RH sensors
    :param Quantity max_abs_error: sensors with standard deviations farther \
       from the average standard deviation will be flagged.
    :rtype: list of length len(data)
    :return: list containing 0s by default and 3 in the position of each sensor
       flagged for variability.
    """
    to_return = np.zeros(len(data))
    # Calculate the mean of each sensor
    sdevs = np.zeros(len(data))
    for i in range(len(data)):
        sdevs[i] = np.nanstd(data[i])

    while(True):
        # Identify the sensor with the mean farthest from the mean of means
        max_diff = 0
        furthest_from_mean = 0  # index of sensor furthest from mean

        for j in range(len(data)):

            if(np.abs(np.nanmean(sdevs)-sdevs[j]) >
               np.abs(np.nanmean(sdevs)-sdevs[furthest_from_mean])):
                furthest_from_mean = j

            for k in range(len(data)):
                if(np.abs(sdevs[j]-sdevs[k]) > max_diff):
                    max_diff = np.abs(sdevs[j]-sdevs[k])

            # If the furthest sensor is farther than max_abs_error from the
            # mean of means, eliminate it and perform analysis again.
            if(max_diff > max_abs_error):
                to_return[furthest_from_mean] = 3
                sdevs[furthest_from_mean] = np.NaN
            else:
                return to_return


def identify_profile(alts, alt_times, confirm_bounds=True,
                     profile_start_height=None, to_return=[], ind=0):
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

    isDone = False
    # Get the starting height from the user
    if profile_start_height is None:
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
    end_ind_des = None
    peak_ind = None

    # Check through valid alts for start_ind_asc, peak_ind, end_ind_des in
    # that order
    while ind < len(alts) - 10:
            if(start_ind_asc is None):  # should be ascending
                # finds at what time the altitude range is reached going up

                # Error if starts on a descent
                if(alts[ind] > profile_start_height and
                   alts[ind + 10] < alts[ind]):
                    print("Error separating profiles: start height is first \
                          reached on a descent")
                    break

                # Set start_ind_asc when the craft is first above
                # profile_start_height
                if alts[ind] > profile_start_height:
                    start_ind_asc = ind

                ind += 1

            elif(end_ind_des is None):  # should be descending

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
                if confirm_bounds:
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
                    if valid in "yYyesYes" or valid is "":
                        plt.close()
                        if end_ind_des is None:
                            print("Could not find end time des (LineTag B)")
                            return
                        isDone = True
                        break
                    elif valid in "nNnoNo":
                        plt.close()
                        to_return = identify_profile(alts, alt_times,
                                                     to_return)
                    else:
                        print("Invalid choice. Re-selecting profile...")
                        plt.close()
                        to_return = identify_profile(alts, alt_times,
                                                     to_return)
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
