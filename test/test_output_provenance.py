"""
What a written file says about how it was made.

A c1 file used to record its values and nothing else: not the coefficients
applied, not the QC thresholds, not which sensor had been dropped from the
ensemble mean. Reprocessing an old flight and getting a different answer
was undiagnosable.
"""
import netCDF4
import numpy as np
import pytest

from profiles import qc
from profiles.processing import ProcessingConfig, all_profiles, process_flights
from test.harness import (PROFILE_START_HEIGHT, RESOLUTION, RES_UNITS,
                          staged_bin)


@pytest.fixture(scope='module')
def written(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('provenance')
    config = ProcessingConfig(
        resolution=RESOLUTION, res_units=RES_UNITS, dev=True, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT, min_levels=50,
        bias_correction='rh_poly22_corrected_t')

    profile = all_profiles(process_flights([str(staged_bin(tmp))], config))[0]
    path = tmp / 'c1.cdf'
    profile.save_netcdf(str(path))
    return profile, path


def test_records_the_coefficients_actually_applied(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        attributes = {name: handle.getncattr(name) for name in handle.ncattrs()}

    # The fixture's imet1 row is A=9.89E-04 B=2.65E-04 C=1.38E-07, E2.
    assert attributes['coef_imet1_a'] == '9.89E-04'
    assert attributes['coef_imet1_equation'] == 'E2'
    assert attributes['coef_wind_equation'] == 'E1'
    assert 'Steinhart-Hart' in attributes['coef_temperature_source']


def test_records_qc_thresholds(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        assert handle.getncattr('qc_temp_max_bias') == pytest.approx(0.25)
        assert handle.getncattr('qc_rh_max_variability') == pytest.approx(0.2)


def test_records_processing_identity(written):
    import profiles
    _, path = written
    with netCDF4.Dataset(path) as handle:
        assert handle.getncattr('processing_version') == profiles.__version__
        assert handle.getncattr('processing_level') == 'c1'
        assert handle.getncattr('tail_number') == 'N934UA'
        assert 'coefficient_directory' in handle.ncattrs()
        assert 'coefficient_revision' in handle.ncattrs()


def test_records_the_bias_correction_even_though_it_is_not_applied(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        assert handle.getncattr('rh_bias_correction') == 'rh_poly22_corrected_t'
        assert handle.getncattr('rh_bias_p10') == pytest.approx(1.1894)
        np.testing.assert_allclose(
            handle.getncattr('rh_bias_fitted_rh_range'), [20., 95.])


def test_qc_flags_reach_the_published_file(written):
    """Previously flags existed only on the thermo_ intermediate."""
    profile, path = written
    with netCDF4.Dataset(path) as handle:
        assert 'temp_qc' in handle.variables
        assert 'rh_qc' in handle.variables
        np.testing.assert_array_equal(handle.variables['temp_qc'][:],
                                      np.asarray(profile.temp_flags))


def test_qc_flags_are_cf_conventional(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        variable = handle.variables['temp_qc']
        values = list(variable.flag_values)
        meanings = variable.flag_meanings.split()

        assert len(values) == len(meanings)
        assert sorted(values) == sorted(qc.FLAG_MEANINGS)
        for value, meaning in zip(values, meanings):
            assert qc.FLAG_MEANINGS[value] == meaning


def test_flags_are_per_sensor_not_per_level(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        assert handle.variables['temp_qc'].dimensions == ('sensor',)
        assert len(handle.dimensions['sensor']) == 4


def test_thresholds_are_configurable(tmp_path):
    config = ProcessingConfig(
        resolution=RESOLUTION, res_units=RES_UNITS, dev=True, nc_level=None,
        profile_start_height=PROFILE_START_HEIGHT, min_levels=50,
        qc_thresholds={'temp': (5.0, 5.0)})

    profile = all_profiles(process_flights(
        [str(staged_bin(tmp_path))], config))[0]

    assert profile.qc_thresholds['temp'] == (5.0, 5.0)
    # A threshold that loose should reject nothing but the empty sensor.
    assert set(profile.temp_flags) <= {qc.GOOD, qc.EMPTY}
