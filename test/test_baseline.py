"""
Characterization tests: the package must keep producing the numbers it
produced when the snapshot was captured.

These are the safety net for restructuring. They do not assert the output is
*correct* - several known defects are baked into the current snapshot on
purpose. They assert it does not change by accident. When a change is
intentional, re-capture with ``python -m test.capture_baseline`` and record
the reported differences in the commit message.
"""
import numpy as np
import pytest

from test import BASELINE_PATH
from test.capture_baseline import VARIANTS
from test.harness import run_reference_pipeline, staged_bin


@pytest.fixture(scope='module')
def bin_path(tmp_path_factory):
    return staged_bin(tmp_path_factory.mktemp('flight'))


def _describe(name, expected, actual):
    if expected.shape != actual.shape:
        return f'{name}: shape {expected.shape} -> {actual.shape}'
    with np.errstate(invalid='ignore'):
        delta = np.abs(np.asarray(actual, float) - np.asarray(expected, float))
    worst = np.nanmax(delta) if delta.size else 0.0
    idx = int(np.nanargmax(delta)) if delta.size else 0
    return (f'{name}: max |delta| = {worst:.6g} at index {idx} '
            f'(was {expected.flat[idx]:.6g}, now {actual.flat[idx]:.6g})')


@pytest.mark.parametrize('variant,lowpass', sorted(VARIANTS.items()))
def test_matches_baseline(bin_path, variant, lowpass):
    snapshot_path = BASELINE_PATH / f'{variant}.npz'
    if not snapshot_path.exists():
        pytest.skip(f'no baseline captured yet: run python -m test.capture_baseline')

    expected = np.load(snapshot_path)
    actual = run_reference_pipeline(bin_path, lowpass=lowpass)

    missing = sorted(set(expected.files) - set(actual))
    added = sorted(set(actual) - set(expected.files))
    assert not missing, f'variables disappeared from the output: {missing}'
    assert not added, f'variables appeared in the output: {added}'

    differences = []
    for name in expected.files:
        want, got = expected[name], np.asarray(actual[name])
        if want.shape != got.shape or not np.allclose(
                want, got, rtol=1e-9, atol=0, equal_nan=True):
            differences.append(_describe(name, want, got))

    assert not differences, (
        f'{len(differences)} of {len(expected.files)} arrays changed:\n  '
        + '\n  '.join(differences)
        + '\n\nIf this is intentional, re-capture with '
          '`python -m test.capture_baseline` and list these in the commit.')
