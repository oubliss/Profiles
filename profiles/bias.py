"""
Bias corrections, applied on top of whatever produced the calibrated value.

A third layer, versioned separately from the sensor coefficients: the
thermistor calibration says what the sensor read, the bias correction says
how that reading departs from truth across the measurement range. It
applies whether the underlying value came from a table lookup or from the
CopterSonde's onboard calibration, which is why it is not part of either.

Both corrections the lab has are MATLAB surface fits over relative
humidity and temperature, so one evaluator covers them. A coefficient
named ``pij`` multiplies ``rh**i * temp**j``, which is exactly MATLAB's
convention, so coefficients can be pasted straight out of the fit report.

    poly22:  p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
    poly41:  p00 + p10*x + p01*y + p20*x^2 + p11*x*y
             + p30*x^3 + p21*x^2*y + p40*x^4 + p31*x^3*y

Provenance matters here as much as the numbers. Each correction records
where it came from and over what range it was fitted, and
``applied_outside_range`` reports extrapolation rather than hiding it.
"""
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

#: Coefficient names look like p<i><j>: i is the rh power, j the temp power.
_TERM = re.compile(r'^p(\d)(\d)$')


@dataclass(frozen=True)
class SurfaceCorrection:
    """ A polynomial surface correction in relative humidity and temperature.

    :var name: identifier recorded in the output
    :var coefficients: {'p00': ..., 'p10': ...}, MATLAB naming
    :var source: where the fit came from, for provenance
    :var rh_range: (low, high) relative humidity the fit was made over
    :var temp_range: (low, high) temperature in C the fit was made over
    :var notes: anything a reader needs to apply it correctly
    """
    name: str
    coefficients: dict
    source: str = ''
    rh_range: Optional[tuple] = None
    temp_range: Optional[tuple] = None
    notes: str = ''
    _terms: tuple = field(default=(), repr=False, compare=False)

    def __post_init__(self):
        terms = []
        for key, value in self.coefficients.items():
            match = _TERM.match(key)
            if match is None:
                raise ValueError(
                    f'{self.name}: coefficient {key!r} is not of the form '
                    f'p<i><j>, e.g. p21 for rh^2 * temp')
            terms.append((int(match.group(1)), int(match.group(2)),
                          float(value)))
        object.__setattr__(self, '_terms', tuple(sorted(terms)))

    def __call__(self, rh, temp):
        """ Corrected relative humidity.

        :param rh: relative humidity in percent
        :param temp: temperature in degrees Celsius
        :rtype: np.ndarray
        """
        rh = np.asarray(rh, dtype=float)
        temp = np.asarray(temp, dtype=float)

        total = np.zeros(np.broadcast(rh, temp).shape, dtype=float)
        for i, j, coefficient in self._terms:
            total = total + coefficient * rh ** i * temp ** j
        return total

    def applied_outside_range(self, rh, temp):
        """ Fraction of samples falling outside the fitted range.

        :rtype: float
        """
        outside = np.zeros(np.broadcast(rh, temp).shape, dtype=bool)
        for values, bounds in ((rh, self.rh_range), (temp, self.temp_range)):
            if bounds is None:
                continue
            values = np.asarray(values, dtype=float)
            with np.errstate(invalid='ignore'):
                outside |= (values < bounds[0]) | (values > bounds[1])

        total = outside.size
        return float(outside.sum()) / total if total else 0.0

    def provenance(self):
        """ Attributes describing this correction, for the output file.

        :rtype: dict
        """
        attributes = {
            'rh_bias_correction': self.name,
            'rh_bias_correction_source': self.source,
            'rh_bias_correction_form':
                ' + '.join(f'p{i}{j}*rh^{i}*temp^{j}'
                           for i, j, _ in self._terms),
        }
        for i, j, coefficient in self._terms:
            attributes[f'rh_bias_p{i}{j}'] = coefficient
        if self.rh_range:
            attributes['rh_bias_fitted_rh_range'] = list(self.rh_range)
        if self.temp_range:
            attributes['rh_bias_fitted_temp_range'] = list(self.temp_range)
        if self.notes:
            attributes['rh_bias_correction_notes'] = self.notes
        return attributes


#: Chamber fit over RH 20-95% at 10, 22 and 35 C, for aircraft whose
#: built-in temperature sensor was NOT separately corrected.
RH_POLY22_UNCORRECTED_T = SurfaceCorrection(
    name='rh_poly22_uncorrected_t',
    coefficients={'p00': 19.9512, 'p10': 1.1893, 'p01': -0.1269,
                  'p20': -6.5933e-04, 'p11': -6.2418e-05, 'p02': 1.6144e-04},
    source='RH_bias_correction - 20 to 95.txt',
    rh_range=(20., 95.),
    temp_range=(10., 35.),
    notes='Apply when no correction was applied to the built-in '
          'temperature sensor.')

#: The same fit for aircraft whose built-in temperature sensor WAS corrected.
RH_POLY22_CORRECTED_T = SurfaceCorrection(
    name='rh_poly22_corrected_t',
    coefficients={'p00': 20.1526, 'p10': 1.1894, 'p01': -0.1282,
                  'p20': -6.5933e-04, 'p11': -6.2744e-05, 'p02': 1.6335e-04},
    source='RH_bias_correction - 20 to 95.txt',
    rh_range=(20., 95.),
    temp_range=(10., 35.),
    notes='Apply when the built-in temperature sensor has been corrected.')

#: Higher-order general fit (MATLAB poly41).
RH_POLY41_GENERAL = SurfaceCorrection(
    name='rh_poly41_general',
    coefficients={'p00': 13.68, 'p10': 0.1244, 'p01': -0.03791,
                  'p20': 0.0334, 'p11': 0.001814, 'p30': -0.0003154,
                  'p21': -6.892e-05, 'p40': 4.864e-07, 'p31': 5.873e-07},
    source='Poly41 Correction.png (sf_H_general)',
    notes='Fitted ranges were not recorded with the coefficients; '
          'extrapolation cannot be detected for this one.')

#: Corrections available by name.
CORRECTIONS = {c.name: c for c in (RH_POLY22_UNCORRECTED_T,
                                   RH_POLY22_CORRECTED_T,
                                   RH_POLY41_GENERAL)}


def get(name):
    """ A registered correction.

    :param str name: its registered name
    :rtype: SurfaceCorrection
    """
    if name not in CORRECTIONS:
        raise KeyError(f'unknown bias correction {name!r}; '
                       f'registered: {sorted(CORRECTIONS)}')
    return CORRECTIONS[name]


def apply_rh_correction(rh, temp, correction, warn_outside=0.05):
    """ Apply a humidity bias correction, warning on extrapolation.

    :param rh: relative humidity in percent
    :param temp: temperature in degrees Celsius
    :param correction: a SurfaceCorrection, or the name of a registered one
    :param float warn_outside: warn when this fraction of samples falls
       outside the fitted range
    :rtype: np.ndarray
    """
    if isinstance(correction, str):
        correction = get(correction)

    fraction = correction.applied_outside_range(rh, temp)
    if fraction > warn_outside:
        warnings.warn(
            f'{correction.name}: {fraction:.0%} of samples fall outside the '
            f'range it was fitted over (rh {correction.rh_range}, temp '
            f'{correction.temp_range} C). The correction is extrapolating.',
            stacklevel=2)

    return correction(rh, temp)
