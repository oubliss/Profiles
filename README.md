# Profiles

This package was built to handle processing data primarily from the CopterSonde, a thermodynamic and kinematic profiling UAS.

See https://oucass.github.io/Profiles for detailed documentation of the API.

## Quick start

```python
from profiles.processing import ProcessingConfig, process_flights, all_profiles

config = ProcessingConfig(resolution=10, res_units='m', ascent=True,
                          profile_start_height=340, min_levels=20)

results = process_flights(['flight2859.BIN'], config)
for profile in all_profiles(results):
    profile.save_cfnetcdf('N934UA', terrain_elevation=340)
```

`.BIN` files are read directly. One `Profile` holds the whole gridded
profile — temperature, humidity, wind and the derived quantities — so there
are no separate thermo and wind objects to fetch.

Files that fail are recorded against their path rather than aborting the
batch:

```python
for result in results:
    if not result.ok:
        print(result.path, result.error)
```

## Installation

##### 1. Clone and install

```
git clone git@github.com:oucass/Profiles.git
cd Profiles
pip install -e .
```

##### 2. Point it at the coefficient tables

The package needs `MasterCoefList.csv` and `copterID.csv`. It looks for them
in this order:

1. a directory passed to `TableCalibration(...)`
2. the `WXUAS_DIR` environment variable
3. `~/.wxuas`

```
export WXUAS_DIR=/path/to/coefficients
```

The `SensorCoefficients` submodule is the canonical source:

```
git submodule init && git submodule update
ln -s $PWD/SensorCoefficients/MasterCoefList.csv ~/.wxuas/MasterCoefList.csv
```

##### 3. DroneLogBook (optional)

Only needed to download flights and pull flight metadata. The package works
without it.

```
cd dronelogbook && pip install -e . && cd ..
python dronelogbook/scripts/sync_dlb.py
```

## Calibration

Which path temperature takes is a property of the *flight*, not of what the
log happens to contain:

| the log | source | temperature |
|---|---|---|
| reports sensor serial numbers | `TableCalibration` | Steinhart-Hart from logged resistance |
| reports none | `OnboardCalibration` | `IMET.T` as logged, untouched |

Current firmware calibrates the thermistors in flight and logs the result,
so there is nothing left to apply — `OnboardCalibration` deliberately does
no arithmetic, and raises rather than returning generic coefficients if
asked for a thermodynamic lookup. Wind is a table lookup either way: the
airframe calibration is per tail number and is not applied onboard.

Detection is automatic. Override it when a log is wrong about itself:

```python
ProcessingConfig(calibration='table')    # or 'onboard', default 'auto'
```

The path taken is written to c1 files as `coef_temperature_source`.

Bias corrections (`profiles/bias.py`) are separate and remain opt-in; they
apply on top of either source.

## Processing levels

| level | what it is |
|---|---|
| a0 | every parsed message at native rate, no QC |
| b1 | calibrated and QC-flagged, still native rate |
| c1 | gridded profile, with QC flags and coefficient provenance |

c1 files record which coefficients produced them, the QC thresholds in
force, the coefficient table's git revision, and per-sensor QC flags as CF
`flag_values` / `flag_meanings`.

## Tests

```
pip install -e '.[test]'
pytest
```

The suite is hermetic: it uses `test/data/coefs` rather than `~/.wxuas` and
makes no network calls.

### References:

https://www.atmos-meas-tech.net/13/2833/2020/ \
https://www.atmos-meas-tech.net/13/3855/2020/ \
https://www.mdpi.com/1424-8220/19/12/2720


