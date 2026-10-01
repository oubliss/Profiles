# Test coefficient fixture

A **frozen snapshot**, not a source of truth. It contains only the sensors that
`flight616_20210609_080023.BIN` reports, so the test suite can run without a
`~/.wxuas` directory and without the `SensorCoefficients` submodule.

Rows were copied from `~/.wxuas/MasterCoefList.csv` on 2026-10-01, which is the
table the package's production processing used at that time.

Do not sync this with the live table. If a sensor here is recalibrated, the
baseline snapshots in `test/data/baseline/` are expected to change, and that
change must be deliberate — see `test/test_baseline.py`.

## Why not use the repo's top-level `coefs/`?

It cannot process this flight. `coefs/MasterCoefList.csv` carries five rows for
each of IMET 62275 / 62288 / 62298, identical except for `ScoopID`, and
`Coef_Manager.get_coefs` only disambiguates duplicates by `Equation` — all five
are `E2`. Every lookup raises:

    RuntimeError: Multiple entries found for "62275" with sensor type "Imet",
    be sure to specify an equation type: ['E2' 'E2' 'E2' 'E2' 'E2']

The top-level table also disagrees with the live table on the coefficient
values themselves (e.g. 62275 A = 1.00761568E-03 here vs 9.89E-04 live), which
is the reproducibility gap that dated coefficient records are meant to close.
