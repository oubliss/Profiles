"""
Retrieval functions, tested without constructing anything.

Before extraction these lived inside constructors that needed a gridded
parent, a unit registry, a file path and a coefficient table, so none of
this could be checked against a hand-computed value.
"""
import numpy as np
import pytest
from metpy.units import units

from profiles.retrievals import thermo, wind

LINEAR = {'A': '32.1', 'B': '-4.2'}
QUADRATIC = {'A': '37.1', 'B': '3.8'}


def attitude(roll_deg, pitch_deg, yaw_deg):
    to_q = lambda v: np.atleast_1d(np.asarray(v, dtype=float)) * units.deg
    return to_q(roll_deg), to_q(pitch_deg), to_q(yaw_deg)


class TestTiltAndAzimuth:
    def test_level_flight_has_zero_tilt(self):
        psi, _ = wind.tilt_and_azimuth(*attitude(0, 0, 0))
        assert psi.magnitude[0] == pytest.approx(0.0, abs=1e-12)

    def test_tilt_matches_closed_form(self):
        """psi = arccos(cos(pitch) cos(roll)), independent of yaw."""
        rng = np.random.default_rng(0)
        roll = rng.uniform(-45, 45, 50)
        pitch = rng.uniform(-45, 45, 50)
        yaw = rng.uniform(-180, 180, 50)

        psi, _ = wind.tilt_and_azimuth(*attitude(roll, pitch, yaw))
        expected = np.arccos(np.cos(np.radians(pitch)) * np.cos(np.radians(roll)))
        np.testing.assert_allclose(psi.magnitude, expected, atol=1e-12)

    def test_tilt_is_independent_of_heading(self):
        psi_north, _ = wind.tilt_and_azimuth(*attitude(10, 5, 0))
        psi_east, _ = wind.tilt_and_azimuth(*attitude(10, 5, 90))
        assert psi_north.magnitude[0] == pytest.approx(psi_east.magnitude[0],
                                                       abs=1e-12)

    def test_azimuth_rotates_with_heading(self):
        """Leaning the same way in body frame points elsewhere in earth frame."""
        _, az_a = wind.tilt_and_azimuth(*attitude(0, 10, 0))
        _, az_b = wind.tilt_and_azimuth(*attitude(0, 10, 90))
        separation = abs(wind.to_compass(az_a).magnitude[0]
                         - wind.to_compass(az_b).magnitude[0])
        assert separation == pytest.approx(90.0, abs=1e-9)


class TestSpeedFromTilt:
    def test_linear_equation(self):
        psi = np.array([0.3]) * units.rad
        speed = wind.speed_from_tilt(psi, LINEAR, 'E1')
        expected = 32.1 * np.sqrt(np.tan(0.3)) - 4.2
        assert speed.magnitude[0] == pytest.approx(expected)

    def test_quadratic_equation(self):
        psi = np.array([0.3]) * units.rad
        speed = wind.speed_from_tilt(psi, QUADRATIC, 'E5')
        root = np.sqrt(np.tan(0.3))
        assert speed.magnitude[0] == pytest.approx(37.1 * root ** 2 + 3.8 * root)

    def test_negative_speeds_become_nan(self):
        """Near vertical the linear fit extrapolates below zero."""
        psi = np.array([1e-6]) * units.rad
        assert np.isnan(wind.speed_from_tilt(psi, LINEAR, 'E1').magnitude[0])

    def test_unknown_equation_is_rejected(self):
        with pytest.raises(KeyError, match='unknown wind calibration'):
            wind.speed_from_tilt(np.array([0.3]) * units.rad, LINEAR, 'E99')

    def test_registered_equations(self):
        assert set(wind.EQUATIONS) == {'E1', 'E5'}


class TestCompass:
    def test_negative_angles_wrap(self):
        wrapped = wind.to_compass(np.array([-90.0, -1.0, 10.0]) * units.deg)
        np.testing.assert_allclose(wrapped.magnitude, [270.0, 359.0, 10.0])


class TestThermoDerivation:
    def test_saturated_air_has_dewpoint_equal_to_temperature(self):
        derived = thermo.derive(np.array([1000.0]) * units.hPa,
                                np.array([293.15]) * units.kelvin,
                                np.array([100.0]) * units.percent)
        dewpoint = derived['T_d'].to(units.kelvin).magnitude[0]
        assert dewpoint == pytest.approx(293.15, abs=0.05)

    def test_potential_temperature_equals_temperature_at_reference(self):
        derived = thermo.derive(np.array([1000.0]) * units.hPa,
                                np.array([300.0]) * units.kelvin,
                                np.array([50.0]) * units.percent)
        assert derived['theta'].magnitude[0] == pytest.approx(300.0, abs=1e-6)

    def test_drier_air_has_less_water(self):
        wet = thermo.derive(np.array([1000.0]) * units.hPa,
                            np.array([293.15]) * units.kelvin,
                            np.array([80.0]) * units.percent)
        dry = thermo.derive(np.array([1000.0]) * units.hPa,
                            np.array([293.15]) * units.kelvin,
                            np.array([20.0]) * units.percent)
        assert wet['mixing_ratio'].magnitude[0] > dry['mixing_ratio'].magnitude[0]

    def test_returns_every_derived_field(self):
        derived = thermo.derive(np.array([900.0]) * units.hPa,
                                np.array([285.0]) * units.kelvin,
                                np.array([60.0]) * units.percent)
        assert set(derived) == {'mixing_ratio', 'theta', 'T_d', 'q'}
