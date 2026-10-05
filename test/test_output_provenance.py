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


def test_records_the_bias_correction_as_requested_not_applied(written):
    _, path = written
    with netCDF4.Dataset(path) as handle:
        assert 'rh_bias_correction' not in handle.ncattrs()
        assert handle.getncattr('rh_bias_correction_requested') == \
            'rh_poly22_corrected_t'
        assert handle.getncattr('rh_bias_correction_applied') == 'no'
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


class TestNaming:
    """One function builds every output name; five writers used to."""

    def test_explicit_path_always_wins(self):
        from profiles.io import naming
        assert naming.resolve('/tmp/x.cdf', None, 'c1', '') == '/tmp/x.cdf'
        assert naming.resolve('/tmp/x.nc', None, 'c1', '') == '/tmp/x.nc'

    def test_name_is_built_from_metadata(self):
        from profiles.io import naming

        class FakeMeta:
            def get(self, key):
                return {'location': 'Lake Thunderbird',
                        'platform_id': 'N934UA'}[key]

        name = naming.output_name(FakeMeta(), 'c1', '20210609.080023',
                                  resolution=10, tag='Ascending')
        assert name == 'LakeThunderbird10N934UACMTAscending.c1.20210609.080023.cdf'

    def test_a0_omits_resolution_because_it_is_not_gridded(self):
        from profiles.io import naming

        class FakeMeta:
            def get(self, key):
                return {'location': 'KAEFS', 'platform_id': 'N934UA'}[key]

        assert naming.output_name(FakeMeta(), 'a0', '20210609.080023') == \
            'KAEFSN934UACMT.a0.20210609.080023.cdf'

    def test_no_metadata_and_no_fallback_is_an_error(self):
        from profiles.io import naming
        with pytest.raises(IOError, match='specify a file name'):
            naming.resolve('/tmp/flight.BIN', None, 'c1', '')

    def test_fallback_is_used_when_offered(self):
        from profiles.io import naming
        assert naming.resolve('/tmp/f.BIN', None, 'a0', '',
                              fallback='/tmp/f.nc') == '/tmp/f.nc'


def test_specific_humidity_is_scaled_consistently(written, tmp_path):
    """The thermo_ writer wrote kg/kg under a g/kg label; c1 did not."""
    import netCDF4

    profile, c1_path = written
    thermo_path = tmp_path / 'thermo.cdf'
    profile._save_thermo_netCDF(str(thermo_path))

    with netCDF4.Dataset(c1_path) as combined, \
            netCDF4.Dataset(thermo_path) as per_variable:
        combined_q = combined.variables['q'][:]
        per_variable_q = per_variable.variables['q'][:]
        assert combined.variables['q'].units == per_variable.variables['q'].units
        np.testing.assert_allclose(per_variable_q[:len(combined_q)],
                                   combined_q, rtol=1e-12)
        # ~16 g/kg, not ~0.016
        assert 1.0 < float(np.nanmean(combined_q)) < 40.0
