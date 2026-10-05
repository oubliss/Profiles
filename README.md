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

With no metadata and no explicit path, the file is written beside the input
as `<input>.c1.<resolution>.cf_ascent.nc` (`save_netcdf` writes
`<input>.c1.<resolution>.ascent.nc`, so calling both keeps both). Pass
`file_path='out.nc'` to choose the name.

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

The source also decides *which coefficients*: every temperature, wind and
tail-number lookup goes through `flight.calibration_source`. Point it at a
directory with `FlightLog(..., coefficient_dir=...)`,
`ProcessingConfig(coefficient_dir=...)`, or by assigning a
`TableCalibration(directory)` to `flight.calibration_source`. Rows with
`ValidFrom`/`ValidTo` are selected by the flight's start time
(`flight.start_time`). An explicit `tail_number` always overrides the
copterID registry. The directory, source and lookup time are recorded in c1
files (`coefficient_directory`, `calibration_source`,
`coefficient_lookup_time`).

## Which barometer

Current firmware logs every barometer as `BARO` with an instance field `I`;
which one is in the scoop depends on the airframe. `ProcessingConfig(baro_instance=1)`
(default 1, also `FlightLog(..., baro_instance=)`) picks it, and
`ekf_core=0` picks the EKF core (`XKF1` field `C`). If the log has
instance-numbered barometers but not the requested one, processing raises
an error listing the instances present; it never falls back silently.

Old logs with separate `BARO` and `BAR2` messages and no instance field
always use `BAR2` (the external barometer), and the setting does not apply
to them. The instance used is written to the a0 file and to c1 files as
`baro_instance` (with `baro_message_type`); both are absent when the log
numbers no barometers. Passing an already-parsed `flight=` to
`profiles_from_flight` uses that log's own choice.

Bias corrections (`profiles/bias.py`) are **not applied** by the pipeline.
Naming one in `ProcessingConfig(bias_correction=...)` records it in the c1
file as `rh_bias_correction_requested` with `rh_bias_correction_applied = no`,
and processing warns.

## Processing levels

| level | what it is |
|---|---|
| a0 | every parsed message at native rate, no QC |
| b1 | calibrated and QC-flagged, still native rate |
| c1 | gridded profile, with QC flags and coefficient provenance |

c1 files record which coefficients produced them, the QC thresholds in
force, the coefficient table's git revision (the last commit that changed
`MasterCoefList.csv`, suffixed `-dirty` if it has uncommitted edits) and its
SHA-256 (`coefficient_sha256`), and per-sensor QC flags as CF
`flag_values` / `flag_meanings`.

`Profile.save_netcdf(lookup_place=True)` additionally records a place name
for the first fix; that asks OpenStreetMap Nominatim, so it is off by default.

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


