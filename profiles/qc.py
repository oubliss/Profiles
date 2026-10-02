"""
Ensemble QC for redundant sensors.

The CopterSonde carries several nominally identical temperature and
humidity sensors. Where one disagrees with the rest it is flagged, and the
remaining ensemble is averaged.

Flags are integers so they can be written straight into a NetCDF variable:

    0  good
    2  rejected for bias (mean too far from the ensemble)
    3  rejected for variability (standard deviation too far from the ensemble)
    4  empty - the sensor position is not populated

Nothing here touches a file, a unit registry or a Profile.
"""
import warnings

import numpy as np
from metpy.units import units as u

GOOD = 0
BIAS = 2
VARIABILITY = 3
EMPTY = 4

#: Human-readable flag names, for writing into output attributes.
FLAG_MEANINGS = {GOOD: 'good', BIAS: 'bias', VARIABILITY: 'lag',
                 EMPTY: 'empty'}


def _reject_outliers(statistics, max_abs_error, flag):
    """ Iteratively flag the sensor furthest from the ensemble consensus.

    While the spread across the still-accepted sensors exceeds
    max_abs_error, drop the one furthest from their mean and re-test.

    :param np.ndarray statistics: one summary statistic per sensor (mean for
       bias, standard deviation for variability)
    :param float max_abs_error: largest spread tolerated across the ensemble
    :param int flag: flag value to record for each rejected sensor
    :rtype: np.ndarray
    :return: array of length len(statistics), 0 where accepted
    """
    statistics = np.array(statistics, dtype=float)
    flags = np.zeros(len(statistics))

    while True:
        accepted = ~np.isnan(statistics)

        # With fewer than two sensors left there is no consensus to compare
        # against, so no further rejection is defensible.
        if accepted.sum() < 2:
            return flags

        spread = statistics[accepted].max() - statistics[accepted].min()
        if spread <= max_abs_error:
            return flags

        deviation = np.abs(statistics - np.mean(statistics[accepted]))
        deviation[~accepted] = -np.inf

        worst = int(np.argmax(deviation))
        flags[worst] = flag
        statistics[worst] = np.nan


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
    with warnings.catch_warnings():
        # An all-NaN sensor is legitimate here; it simply never participates.
        warnings.simplefilter("ignore", category=RuntimeWarning)
        means = np.array([np.nanmean(sensor) for sensor in data], dtype=float)

    return _reject_outliers(means, max_abs_error, flag=2)


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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sdevs = np.array([np.nanstd(sensor) for sensor in data], dtype=float)

    return _reject_outliers(sdevs, max_abs_error, flag=3)


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
