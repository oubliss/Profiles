"""
Deprecations are visible, not swallowed.

pytest.ini used to ignore every DeprecationWarning, which hid the package's
own. It now turns one into an error when it is attributed to `profiles` or to
the test suite (third-party ones stay ignored). These tests pin that policy
and that the import-time shims really do warn. (Profile_Set's own warning is
asserted in test_writer_review.)
"""
import subprocess
import sys

import pytest


def _warn_from(module_name):
    exec("import warnings\nwarnings.warn('x', DeprecationWarning)",
         {'__name__': module_name})


@pytest.mark.parametrize('module', ['profiles', 'profiles.somewhere',
                                    'test.test_something'])
def test_own_deprecations_are_errors(module):
    with pytest.raises(DeprecationWarning):
        _warn_from(module)


def test_third_party_deprecations_stay_quiet():
    _warn_from('numpy.fake')     # must not raise under pytest.ini


@pytest.mark.parametrize('module,text', [
    ('profiles.conf', 'profiles.conf is deprecated'),
    ('profiles.Raw_Profile', 'profiles.flight'),
])
def test_import_time_shims_warn(module, text):
    """A fresh interpreter: the warning fires once, at first import."""
    code = ('import warnings\n'
            'warnings.simplefilter("error", DeprecationWarning)\n'
            f'import {module}\n')
    result = subprocess.run([sys.executable, '-c', code],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert 'DeprecationWarning' in result.stderr
    assert text in result.stderr
