"""
Coefficient tables: dated selection, honest ambiguity, bias corrections.
"""
import warnings

import numpy as np
import pytest

from profiles import bias
from profiles.Coef_Manager import (OnboardCalibration, TableCalibration,
                                   source_for_flight)
from profiles.coef_table import (AmbiguousCoefficients, CoefTable,
                                 CopterRegistry, MissingCoefficients)
from test import COEF_PATH

HEADER = 'SensorType,SerialNumber,ScoopID,Equation,A,B,C,D,Offset,SensorStatus'
DATED_HEADER = HEADER + ',ValidFrom,ValidTo'


def write_table(tmp_path, lines, header=HEADER):
    path = tmp_path / 'MasterCoefList.csv'
    path.write_text(header + '\n' + '\n'.join(lines) + '\n')
    return CoefTable(path)


class TestDatedSelection:
    def test_undated_rows_match_any_date(self, tmp_path):
        table = write_table(tmp_path, ['Imet,100,A,E2,1.0,2.0,3.0,na,na,Active'])
        assert table.lookup('Imet', 100, when='2021-06-09')['A'] == '1.0'
        assert table.lookup('Imet', 100)['A'] == '1.0'

    def test_recalibration_is_selected_by_flight_date(self, tmp_path):
        table = write_table(tmp_path, [
            'Imet,100,A,E2,1.0,2.0,3.0,na,na,Active,2020-01-01,2022-06-01',
            'Imet,100,A,E2,9.0,8.0,7.0,na,na,Active,2022-06-01,',
        ], header=DATED_HEADER)

        assert table.lookup('Imet', 100, when='2021-06-09')['A'] == '1.0'
        assert table.lookup('Imet', 100, when='2023-01-01')['A'] == '9.0'

    def test_validity_window_is_half_open(self, tmp_path):
        table = write_table(tmp_path, [
            'Imet,100,A,E2,1.0,2.0,3.0,na,na,Active,2020-01-01,2022-06-01',
            'Imet,100,A,E2,9.0,8.0,7.0,na,na,Active,2022-06-01,',
        ], header=DATED_HEADER)
        # the boundary belongs to the later row, so no date matches twice
        assert table.lookup('Imet', 100, when='2022-06-01')['A'] == '9.0'

    def test_date_before_any_row_is_an_error_not_a_guess(self, tmp_path):
        table = write_table(tmp_path, [
            'Imet,100,A,E2,1.0,2.0,3.0,na,na,Active,2022-06-01,',
        ], header=DATED_HEADER)
        with pytest.raises(MissingCoefficients, match='none valid at'):
            table.lookup('Imet', 100, when='2020-01-01')

    def test_duplicate_rows_raise_with_detail(self, tmp_path):
        """The repo's own coefs/ has five identical-equation rows per sensor."""
        table = write_table(tmp_path, [
            'Imet,100,A,E2,1.0,2.0,3.0,na,na,Active',
            'Imet,100,B,E2,1.0,2.0,3.0,na,na,Active',
        ])
        with pytest.raises(AmbiguousCoefficients) as raised:
            table.lookup('Imet', 100)
        assert 'ValidFrom' in str(raised.value)

    def test_equation_separates_otherwise_identical_rows(self, tmp_path):
        table = write_table(tmp_path, [
            'Wind,N1,na,E1,32.1,-4.2,na,na,na,Active',
            'Wind,N1,na,E5,37.1,3.8,na,na,na,Active',
        ])
        assert table.lookup('Wind', 'N1', equation='E1')['A'] == '32.1'
        assert table.lookup('Wind', 'N1', equation='E5')['A'] == '37.1'

    def test_missing_serial_names_the_file(self, tmp_path):
        table = write_table(tmp_path, ['Imet,100,A,E2,1.0,2.0,3.0,na,na,Active'])
        with pytest.raises(MissingCoefficients, match='MasterCoefList'):
            table.lookup('Imet', 999)

    def test_serials_compare_regardless_of_numeric_form(self, tmp_path):
        table = write_table(tmp_path, ['Imet,100,A,E2,1.0,2.0,3.0,na,na,Active'])
        for form in (100, '100', 100.0, '100.0'):
            assert table.lookup('Imet', form)['A'] == '1.0'


class TestCopterRegistry:
    def test_single_mapping_is_silent(self, tmp_path):
        path = tmp_path / 'copterID.csv'
        path.write_text('6,N934UA\n')
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            assert CopterRegistry(path).tail_number(6) == 'N934UA'

    def test_ambiguous_mapping_warns_once_and_says_why(self, tmp_path):
        path = tmp_path / 'copterID.csv'
        path.write_text('1,FA3TANE3MF\n1,N944UA\n')
        registry = CopterRegistry(path)

        with pytest.warns(UserWarning, match='maps to 2 tail numbers'):
            assert registry.tail_number(1) == 'FA3TANE3MF'

        with warnings.catch_warnings():
            warnings.simplefilter('error')
            registry.tail_number(1)   # already warned; must not warn again

    def test_unknown_id_raises(self, tmp_path):
        path = tmp_path / 'copterID.csv'
        path.write_text('6,N934UA\n')
        with pytest.raises(MissingCoefficients, match='not in'):
            CopterRegistry(path).tail_number(99)


class TestCalibrationSource:
    def test_table_source_reads_the_fixture(self):
        source = TableCalibration(COEF_PATH)
        assert source.get_coefs('Imet', 62275)['Equation'] == 'E2'

    def test_onboard_source_refuses_thermo_lookups(self):
        source = OnboardCalibration(COEF_PATH)
        with pytest.raises(MissingCoefficients, match='calibrated onboard'):
            source.get_coefs('Imet', 62275)

    def test_onboard_source_still_resolves_wind(self):
        source = OnboardCalibration(COEF_PATH)
        assert source.get_coefs('Wind', 'N934UA', 'E1')['A'] == '3.21E+01'

    def test_source_is_chosen_by_what_the_log_reports(self):
        legacy = {'imet1': 62275, 'rh1': 16, 'copterID': 6.0}
        onboard = {'imet1': 0, 'imet2': 0, 'rh1': 0, 'copterID': 1.0}
        assert isinstance(source_for_flight(legacy, COEF_PATH),
                          TableCalibration)
        assert isinstance(source_for_flight(onboard, COEF_PATH),
                          OnboardCalibration)

    def test_missing_directory_says_where_it_looked(self, tmp_path):
        with pytest.raises(MissingCoefficients, match='passed explicitly'):
            TableCalibration(tmp_path / 'nope')


class TestBiasCorrections:
    def test_poly22_matches_a_hand_computed_value(self):
        correction = bias.get('rh_poly22_uncorrected_t')
        rh, temp = 60.0, 22.0
        expected = (19.9512 + 1.1893 * rh - 0.1269 * temp
                    - 6.5933e-04 * rh ** 2 - 6.2418e-05 * rh * temp
                    + 1.6144e-04 * temp ** 2)
        assert correction(rh, temp) == pytest.approx(expected)

    def test_poly41_matches_a_hand_computed_value(self):
        correction = bias.get('rh_poly41_general')
        rh, temp = 50.0, 15.0
        expected = (13.68 + 0.1244 * rh - 0.03791 * temp + 0.0334 * rh ** 2
                    + 0.001814 * rh * temp - 0.0003154 * rh ** 3
                    - 6.892e-05 * rh ** 2 * temp + 4.864e-07 * rh ** 4
                    + 5.873e-07 * rh ** 3 * temp)
        assert correction(rh, temp) == pytest.approx(expected)

    def test_works_on_arrays(self):
        correction = bias.get('rh_poly22_corrected_t')
        result = correction(np.array([30., 60., 90.]), np.array([10., 20., 30.]))
        assert result.shape == (3,)

    def test_extrapolation_is_reported(self):
        correction = bias.get('rh_poly22_uncorrected_t')
        inside = correction.applied_outside_range(np.array([30., 60.]),
                                                  np.array([15., 25.]))
        outside = correction.applied_outside_range(np.array([5., 99.]),
                                                   np.array([15., 25.]))
        assert inside == 0.0
        assert outside == 1.0

    def test_extrapolation_warns_when_applied(self):
        with pytest.warns(UserWarning, match='extrapolating'):
            bias.apply_rh_correction(np.array([2., 3.]), np.array([40., 40.]),
                                     'rh_poly22_uncorrected_t')

    def test_provenance_records_every_coefficient(self):
        attributes = bias.get('rh_poly22_uncorrected_t').provenance()
        assert attributes['rh_bias_correction'] == 'rh_poly22_uncorrected_t'
        assert attributes['rh_bias_p10'] == 1.1893
        assert attributes['rh_bias_fitted_rh_range'] == [20., 95.]

    def test_malformed_coefficient_name_is_rejected(self):
        with pytest.raises(ValueError, match='p<i><j>'):
            bias.SurfaceCorrection(name='bad', coefficients={'alpha': 1.0})

    def test_unknown_correction_is_rejected(self):
        with pytest.raises(KeyError, match='unknown bias correction'):
            bias.get('nope')
