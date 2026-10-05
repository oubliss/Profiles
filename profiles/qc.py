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


def is_populated(series):
    """ Does this sensor position carry real measurements?

    An unfitted position is logged as zeros, and a sensor that never
    reported is all NaN. Either way it has nothing to contribute to an
    ensemble comparison.

    :param series: one sensor's readings
    :rtype: bool
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(series)
    return bool(np.isfinite(mean)) and mean != 0


def qc(data, max_bias, max_variance):
    """ Flag sensors in an ensemble that disagree with the rest.

    Only populated sensors take part in the comparison. Including the
    unfitted positions was catastrophic: the CopterSonde logs three
    thermistors in four slots, so every flight had an all-zero sensor in
    the ensemble. _bias then saw a ~283 K spread and _s_dev a ~0.9 K one,
    and because rejection stops at two survivors it threw out the two real
    sensors and kept the two empty ones. Temperature came out entirely NaN
    on every flight where both remaining positions were empty.

    :param list<Quantity> data: one series per sensor, all the same type
    :param max_bias: largest tolerated spread across the sensors' means
    :param max_variance: largest tolerated spread across their standard
       deviations
    :rtype: list<int> of length len(data)
    :return: GOOD, BIAS, VARIABILITY or EMPTY for each sensor
    """
    if isinstance(data, u.Quantity):
        data = data.magnitude

    flags = [EMPTY if not is_populated(series) else GOOD for series in data]

    populated = [i for i, flag in enumerate(flags) if flag == GOOD]
    if len(populated) < 2:
        # Nothing to compare against; a lone sensor is accepted as-is.
        return flags

    subset = [data[i] for i in populated]
    bias_flags = _bias(subset, max_bias)
    sdev_flags = _s_dev(subset, max_variance)

    for position, index in enumerate(populated):
        flags[index] = int(max(bias_flags[position], sdev_flags[position]))

    return flags
