"""
Sensor-ensemble QC.

The original _bias/_s_dev put both the rejection and the return inside the
per-sensor loop while max_diff accumulated across it. The first pass therefore
compared only sensor 0 against the others, and once anything was rejected the
stale max_diff kept the branch alive for the rest of that sweep, so a single
outlier could take otherwise good sensors down with it.
"""
import numpy as np
import pytest

from profiles.utils import _bias, _s_dev, qc

GOOD, BIAS, LAG, EMPTY = 0, 2, 3, 4


def ensemble(offsets, n=500, noise=0.0, seed=0):
    """Four sensors tracking one signal, each with a constant offset."""
    rng = np.random.default_rng(seed)
    base = np.linspace(290.0, 300.0, n)
    return [base + off + (rng.normal(0, noise, n) if noise else 0.0)
            for off in offsets]


class TestBias:
    def test_agreeing_sensors_are_all_kept(self):
        flags = _bias(ensemble([0.0, 0.05, -0.05, 0.02]), max_abs_error=0.25)
        assert list(flags) == [GOOD] * 4

    def test_single_outlier_is_rejected_alone(self):
        # Sensor 2 sits 5 K off; the other three agree to within 0.05 K.
        flags = _bias(ensemble([0.0, 0.05, 5.0, -0.02]), max_abs_error=0.25)
        assert list(flags) == [GOOD, GOOD, BIAS, GOOD], (
            'exactly one sensor should be rejected')

    def test_two_outliers_both_rejected(self):
        flags = _bias(ensemble([0.0, 4.0, -4.0, 0.03]), max_abs_error=0.25)
        assert list(flags) == [GOOD, BIAS, BIAS, GOOD]

    def test_spread_exactly_at_threshold_is_accepted(self):
        flags = _bias(ensemble([0.0, 0.25]), max_abs_error=0.25)
        assert list(flags) == [GOOD, GOOD]

    def test_spread_just_over_threshold_rejects_one(self):
        flags = _bias(ensemble([0.0, 0.26]), max_abs_error=0.25)
        assert sum(f == BIAS for f in flags) == 1

    def test_never_rejects_below_two_survivors(self):
        """Wildly disagreeing sensors must not all be thrown away."""
        flags = _bias(ensemble([0.0, 50.0, -50.0, 100.0]), max_abs_error=0.25)
        assert sum(f == GOOD for f in flags) >= 1

    def test_all_nan_sensor_does_not_poison_the_ensemble(self):
        data = ensemble([0.0, 0.05, -0.05])
        data.append(np.full(500, np.nan))
        flags = _bias(data, max_abs_error=0.25)
        assert list(flags[:3]) == [GOOD, GOOD, GOOD]

    def test_outlier_position_does_not_matter(self):
        """The old loop favoured whichever index it reached first."""
        for position in range(4):
            offsets = [0.0, 0.02, -0.02, 0.01]
            offsets[position] = 6.0
            flags = _bias(ensemble(offsets), max_abs_error=0.25)
            assert flags[position] == BIAS, f'outlier at {position} not caught'
            assert sum(f == BIAS for f in flags) == 1, (
                f'outlier at {position} took other sensors with it')


class TestStandardDeviation:
    def test_matched_variability_is_kept(self):
        flags = _s_dev(ensemble([0, 0, 0, 0], noise=0.10, seed=1),
                       max_abs_error=0.10)
        assert list(flags) == [GOOD] * 4

    def test_flatlined_sensor_is_rejected(self):
        data = ensemble([0, 0, 0], noise=0.10, seed=2)
        data.append(np.full(500, 295.0))   # no variability at all
        flags = _s_dev(data, max_abs_error=0.5)
        assert flags[3] == LAG
        assert list(flags[:3]) == [GOOD, GOOD, GOOD]


class TestQcCombination:
    def test_empty_sensor_flagged_as_empty(self):
        data = ensemble([0.0, 0.05, -0.05])
        data.append(np.zeros(500))
        flags = qc(data, 0.25, 0.1)
        assert flags[3] == EMPTY

    def test_clean_ensemble_passes(self):
        flags = qc(ensemble([0.0, 0.05, -0.05, 0.02], noise=0.01, seed=3),
                   0.25, 0.1)
        assert list(flags) == [GOOD] * 4

    def test_returns_one_flag_per_sensor(self):
        for n_sensors in (2, 3, 4, 5):
            data = ensemble([0.0] * n_sensors, noise=0.01, seed=4)
            assert len(qc(data, 0.25, 0.1)) == n_sensors


def test_terminates_on_pathological_input():
    """Previously the loop could cascade; make sure it always returns."""
    rng = np.random.default_rng(5)
    for _ in range(200):
        n_sensors = int(rng.integers(2, 6))
        offsets = rng.normal(0, 20, n_sensors)
        flags = _bias(ensemble(list(offsets), n=50), max_abs_error=0.25)
        assert len(flags) == n_sensors
        assert set(np.unique(flags)) <= {GOOD, BIAS}
