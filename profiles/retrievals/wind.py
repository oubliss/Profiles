"""
Wind from airframe tilt.

A multirotor holding position leans into the wind, and the lean angle maps
to wind speed through a per-airframe calibration. The lean direction gives
the wind direction.

This lived in three places: FlightLog.apply_wind_coeffs, and both
Wind_Profile._calc_winds_linear and _calc_winds_quadratic. All three built
the same rotation matrices in the same order; they differed only in which
calibration equation they applied at the very end.

References
----------
Segales et al. 2020, https://doi.org/10.5194/amt-13-2833-2020
"""
import numpy as np
from profiles.unit_registry import units

#: Calibration equations, keyed by the Equation column of MasterCoefList.
#: Each maps the square root of tan(tilt) to wind speed in m/s. Registering
#: them here means the coefficient row selects its own maths, instead of the
#: call site hardcoding 'E1' or 'E5' and assuming the shape.
EQUATIONS = {}


def equation(name):
    """Register a calibration equation under its MasterCoefList name."""
    def register(function):
        EQUATIONS[name] = function
        return function
    return register


@equation('E1')
def _linear(root_tan_psi, coefficients):
    """speed = A * sqrt(tan(psi)) + B"""
    return (float(coefficients['A']) * root_tan_psi
            + float(coefficients['B']))


@equation('E5')
def _quadratic(root_tan_psi, coefficients):
    """speed = A * tan(psi) + B * sqrt(tan(psi))"""
    return (float(coefficients['A']) * root_tan_psi ** 2.
            + float(coefficients['B']) * root_tan_psi)


def tilt_and_azimuth(roll, pitch, yaw):
    """ Tilt from vertical and the direction of lean, from attitude.

    :param roll: roll angles, a pint Quantity in angular units
    :param pitch: pitch angles
    :param yaw: yaw angles
    :rtype: tuple
    :return: (psi, azimuth), both pint Quantities in radians. psi is the
       angle between the airframe's thrust axis and vertical; azimuth is the
       compass direction it leans toward.

    Implementation note: this builds a 3x3 rotation matrix per sample, which
    is how it has always been done and is preserved for bit-exactness. The
    closed form - psi = arccos(cos(pitch) cos(roll)), azimuth from the third
    column of the product - is 92x faster and gives a bit-identical psi, but
    azimuth differs by up to 4.4e-16 rad from the matmul's summation order.
    That is far below anything that matters physically; it is just not worth
    breaking a bit-exact comparison for a loop that costs ~40 ms per flight.
    """
    n_samples = len(roll)
    psi = np.zeros(n_samples) * units.rad
    azimuth = np.zeros(n_samples) * units.rad

    for i in range(n_samples):
        croll = np.cos(roll[i]).magnitude
        sroll = np.sin(roll[i]).magnitude
        cpitch = np.cos(pitch[i]).magnitude
        spitch = np.sin(pitch[i]).magnitude
        cyaw = np.cos(yaw[i]).magnitude
        syaw = np.sin(yaw[i]).magnitude

        rotate_x = np.array([[1, 0, 0],
                             [0, croll, sroll],
                             [0, -sroll, croll]])
        rotate_y = np.array([[cpitch, 0, -spitch],
                             [0, 1, 0],
                             [spitch, 0, cpitch]])
        rotate_z = np.array([[cyaw, -syaw, 0],
                             [syaw, cyaw, 0],
                             [0, 0, 1]])
        rotation = rotate_z @ rotate_y @ rotate_x

        psi[i] = np.arccos(rotation[2, 2])
        azimuth[i] = np.arctan2(rotation[1, 2], rotation[0, 2])

    return psi, azimuth


def speed_from_tilt(psi, coefficients, equation_name):
    """ Apply an airframe's calibration to tilt angle.

    :param psi: tilt from vertical, a pint Quantity
    :param dict coefficients: a MasterCoefList row
    :param str equation_name: which registered equation to use
    :rtype: pint.Quantity
    :return: wind speed in m/s, with physically impossible negatives as NaN
    :raises KeyError: if the equation is not registered
    """
    if equation_name not in EQUATIONS:
        raise KeyError(
            f'unknown wind calibration equation {equation_name!r}; '
            f'registered: {sorted(EQUATIONS)}')

    root_tan_psi = np.sqrt(np.tan(psi)).magnitude
    speed = EQUATIONS[equation_name](root_tan_psi, coefficients)
    speed = speed * units.m / units.s

    # A calibration fitted over a limited tilt range can extrapolate below
    # zero near vertical; there is no such thing as a negative wind speed.
    speed[speed.magnitude < 0.] = np.nan
    return speed


def to_compass(azimuth):
    """ Wrap an azimuth onto [0, 360) degrees.

    :param azimuth: a pint Quantity in angular units
    :rtype: pint.Quantity
    """
    azimuth = azimuth.to(units.deg)
    negative = np.squeeze(np.where(azimuth.magnitude < 0.))
    azimuth[negative] = azimuth[negative] + 360. * units.deg
    return azimuth


def retrieve(roll, pitch, yaw, coefficients, equation_name='E1'):
    """ Wind speed and direction from a time series of attitude.

    Assumes the airframe is holding station: a translating aircraft leans
    for reasons other than wind and this will read that lean as wind.

    :param roll: roll angles, a pint Quantity
    :param pitch: pitch angles
    :param yaw: yaw angles
    :param dict coefficients: the airframe's MasterCoefList row
    :param str equation_name: which calibration equation the row uses
    :rtype: tuple
    :return: (direction, speed) as pint Quantities
    """
    psi, azimuth = tilt_and_azimuth(roll, pitch, yaw)
    return to_compass(azimuth), speed_from_tilt(psi, coefficients,
                                                equation_name)
