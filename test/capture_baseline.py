"""
Re-capture the baseline snapshots.

    python -m test.capture_baseline

Run this only when output is *meant* to change. Commit the regenerated .npz
files together with the code change and say in the commit message which
variables moved and why - test_baseline.py prints exactly that list when it
fails, so paste it in.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

from test import BASELINE_PATH
from test.harness import run_reference_pipeline, staged_bin

VARIANTS = {'flight616_10m': False, 'flight616_10m_lowpass': True}


def main():
    BASELINE_PATH.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        bin_path = staged_bin(Path(tmp))
        for name, lowpass in VARIANTS.items():
            print(f'capturing {name} (lowpass={lowpass}) ...')
            snapshot = run_reference_pipeline(bin_path, lowpass=lowpass)
            out_path = BASELINE_PATH / f'{name}.npz'
            np.savez_compressed(out_path, **snapshot)
            size_kb = out_path.stat().st_size / 1024
            print(f'  -> {out_path} ({len(snapshot)} arrays, {size_kb:.0f} KB)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
