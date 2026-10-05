# Changelog

## 1.4.0 (unreleased)

Stages 1 and 2 of the v2 restructuring plan: run on a current stack, and fix
the defects that reach published files. No API redesign yet.

### Output changes — affects previously published files

**Altitudes in `thermo_*` and `wind_*` c1 files were half a bin too low.**
`Thermo_Profile.alt` and `Wind_Profile.alt` were assigned the bin *edges*
(`gridded_base`) while every other variable in those files was a bin *mean*.
At the 10 m resolution used for the reference flight this is a systematic
**−5 m**; in general it is −resolution/2. They now carry bin centres, matching
`Profile.save_netcdf`, which was already correct and therefore disagreed with
the per-variable files for the same flight.

**Every `thermo_*` and `wind_*` file carried a trailing fill row.** The `time`
coordinate was written from the N+1 bin edges while all variables had N
points, so the unlimited dimension grew to N+1 and NetCDF padded every
variable with one masked value. Child profiles now carry a bin-centre `time`
alongside the edges they bin with.

Anything processed with 1.3.x and submitted to an archive is affected by both.
`processing_version` moves to 1.4.0 so the two can be told apart.

**Sensor QC removed the wrong sensor.** `utils._bias` and `_s_dev` put both
the rejection and the return inside the per-sensor loop, while `max_diff`
accumulated across it. The practical effect was that the *first* sensor was
flagged almost regardless of which one was actually anomalous: with a single
large outlier placed at each of four positions in turn, the old code flagged
sensor 1 every time. The genuinely biased sensor stayed in the ensemble mean
and a good sensor was discarded.

On the reference flight this is not hypothetical:

| | old flags | new flags |
|---|---|---|
| profile 0, temp | `[3 0 0 4]` | `[0 0 2 4]` |
| profile 0, rh   | `[3 0 0 4]` | `[0 0 0 4]` |
| profile 1, temp | `[2 0 0 4]` | `[0 0 2 4]` |

Per-sensor means for profile 0 were imet1 295.606 K, imet2 295.808 K,
imet3 295.075 K. imet3 is 0.42 K from the ensemble mean and imet1 is 0.11 K
from it, so imet3 is the outlier — but imet1 was the sensor being thrown away.
Correcting this shifts gridded temperature by **+0.28 K on average and up to
+0.36 K**, dewpoint by up to +0.34 K, and RH by up to +0.30 %.

Rejection is now a single shared `_reject_outliers` helper: while the spread
across accepted sensors exceeds the threshold, drop the one furthest from
their mean and re-test, never going below two survivors.

**Empty sensor slots poisoned the ensemble.** The CopterSonde logs three
thermistors in four slots, so every flight carried an all-zero sensor into
the comparison above. `_bias` saw a ~283 K spread and `_s_dev` a ~0.9 K one,
and because rejection stops at two survivors it discarded the two *real*
sensors and kept the two empty ones. `qc()` now classifies each position
first — a series whose mean is zero or non-finite is `EMPTY` — and compares
only the populated subset.

Over the 25 OK3DM flights, 22 are bit-identical and 3 changed. Those three
(2026-05-19) have only two populated thermistors and were losing temperature
entirely:

| | old flags | new flags | temp |
|---|---|---|---|
| flight2859 p0 | `[4 3 3 4]` | `[4 0 0 4]` | NaN 36/36 → 0/36 |
| flight2860 p0 | `[4 3 3 4]` | `[4 0 0 4]` | NaN 39/39 → 0/39 |
| flight2861 p0 | `[4 3 3 4]` | `[4 0 0 4]` | NaN 39/39 → 0/39 |

RH recovered on 2860 and 2861 too. The reference flight is unchanged — it has
three populated sensors, which is why the characterization tests did not
catch this.

**Onboard-calibrated temperature was being recomputed with generic
coefficients.** `calibrate_temperature` chose its path by looking for
resistances in the log: present, recompute via Steinhart-Hart. That was right
while only legacy logs carried resistances. Current firmware logs resistances
*and* the temperature it already calibrated onboard, so the recompute branch
always won — and those logs no longer carry `USER_SENSORS` parameters, so
every serial resolved to `0` and the lookup fell to the catch-all `Imet,0`
row. One shared transfer function was applied to all four thermistors.

On flight2859 that moved imet2 by **+0.100 K** and imet3 by **−0.109 K**.
More damaging than the offsets, it compressed the inter-sensor spread from
**0.273 K to 0.064 K** — manufacturing agreement between sensors and feeding
the ensemble QC, whose bias threshold is 0.25 K, a flattened signal.

The choice now belongs to the flight's `CalibrationSource`
(`temperature_from`), resolved from whether the log reports serial numbers
and overridable with `ProcessingConfig(calibration='table'|'onboard')`.
`OnboardCalibration` performs no arithmetic on temperature — that is the
point — and refuses thermodynamic table lookups rather than returning generic
coefficients. Wind remains a table lookup; the airframe calibration is per
tail number and is not applied onboard. The path taken is recorded in the
output as `coef_temperature_source`.

**Current-firmware logs merged redundant sensors into one series.**
ArduPilot 4.x logs every barometer as `BARO` and every EKF core as `XKF1`,
distinguished only by an instance field (`I`, `C`); the schema assumed one of
each (and a separate `BAR2` for the external barometer), so the parser
interleaved all copies. On OK3DM flight2862 that meant **7688** pressure
samples (two barometers, 3844 each, no `BAR2`) differing by **94.5 Pa** on
average, and **11532** rotation samples (three EKF cores) with **7688**
duplicate timestamps. Groups now declare an `Instance(field, keep)`: pressure
keeps `BARO` instance 1, rotation keeps core 0, and each group is back to
3844 samples with strictly increasing times. The IMU's ad-hoc `I == 0` filter
uses the same mechanism. Logs without the field (flight616: `BARO` + `BAR2`,
`NKF1`) are untouched. Override with `parse(..., instances={'pres': 0})`; the
chosen instance is recorded as `source_instance` beside `source_message_type`.
**Which barometer is the scoop sensor is not recorded in the log and the
default needs confirming per airframe.** On flight2862 instance 1 is the
cooler (33.2 C vs 47.1 C die temperature) and more variable (std 1.00 vs
0.24 C) one, the same signature `BAR2` has against `BARO` in flight616
(40.1 vs 44.6 C; 0.54 vs 0.13 C).

**ESC RPM merged motors on current firmware.** Motors were `Instance % 4`;
flight2862 logs instances 8, 9, 11, 12, which became 0, 1, 3, 0 — two motors
averaged together and one empty — and the reshape assumed exactly four
messages per step. Motors are now the sorted set of observed instances,
numbered 1..N, and a timestep ends when an instance repeats, so a dropped
message leaves a NaN instead of shifting every later sample.

### Fixed

- `import profiles` no longer reads the filesystem. `utils.coef_manager` was
  constructed at import time, so the package failed to import without
  `~/.wxuas`, and repointing `conf.coef_info` after import had no effect.
- `np.NaN` → `np.nan`; removed in NumPy 2 and broke all thermo processing.
- Guard `units.define`; modern pint/MetPy already define `percent`.
- Replace `np.matrix` throughout (pending removal from NumPy, and it overloads
  `*` as matmul). Verified bit-exact over 2000 random cases.
- `Profile.trucate_to` → `truncate_to`; the method did not exist, so the
  length-reconciliation branches raised `AttributeError`.
- `Profile` passed `scoop_id` positionally into `Raw_Profile`'s `nc_level`
  slot, silently suppressing the a0 file and coefficient application.
- `Raw_Profile._read_netCDF` assigned `rot_list[6]` three times: PN held PD's
  values, and PE/PD were left as raw NetCDF variables rather than quantities.
- `Raw_Profile.wind_data` returned `rotation[6]` for all of `pos_n`, `pos_e`
  and `pos_d`.
- A file with no `WIND` messages cleared `self.events` where `self.wind` was
  meant.
- `Meta` read the flight file only when a *header* path was supplied.
- Placeholder `lat`/`lon` (when no position is available) were built at the
  edge count rather than the output length.
- `nc_level` now means one thing. `Raw_Profile` tested `== 'low'` while
  `Thermo_Profile` and `Wind_Profile` tested `is not None`, so the documented
  "write nothing" value `'none'` (a truthy string) suppressed the raw file
  but still wrote `thermo_` and `wind_` ones. All three now share
  `utils.writes_netcdf`, which accepts `None`/`'none'`/`'low'`
  case-insensitively and warns on anything else.
  **`scripts/process_data_from_dlb.py` passes `nc_level='None'`** and will
  stop writing per-variable files as a result - that was always the
  documented intent.
- `utils.rh_calib` is now an explicit pass-through. It looked up the sensor's
  `A` coefficient and then overwrote the result with `0` on the next line, so
  no offset has ever reached the data. Behaviour is unchanged; the docstring
  now records what has to be settled before a correction is reinstated.
- **a0 files corrupted temperature on read-back (about +100 K).**
  `FlightLog._save_netCDF` wrote temperatures in K as `volt<n>` labelled "mV"
  and never wrote resistances; `_read_netCDF` read them as millivolts into a
  tuple with no resistance slots, so `thermo_data()` paired the wrong
  arrays. It now writes `temp<n>` and `resi<n>` with the units the values
  carry and rebuilds the same tuples a `.BIN` gives, so reprocessing from an
  a0 file reproduces the gridded profile. RH-sensor and barometer
  temperatures were labelled "F", which pint reads as farad; they now carry
  `kelvin` and `degree_Fahrenheit`. **a0 files written by earlier 1.4.0-dev
  builds have mislabelled temperatures**: they are still read (`volt<n>` is
  taken as K, resistances come back NaN, "F" falls back to the right unit)
  but should be regenerated from the log.
- a0 writer and reader now agree on what is optional. The reader no longer
  raises `KeyError` on logs with no `wind` or `events` group (pre-2021), the
  writer no longer raises on a log with no `copterID`, every serial number
  (including `wind`), rpm and the calibrated `calib_*` groups are read back,
  and message text is stored as variable-length strings instead of `'<U13'`,
  which cut anything longer.
- The a0 fallback name could be the input log itself: `FLIGHT.JSON` or
  `x.Bin` matched none of the `.replace()` patterns, so the log was opened
  for writing and truncated. The name now comes from `os.path.splitext`, and
  `_save_netCDF` raises `ValueError` rather than write over its input.
- `Profile` matched file extensions by substring, had no `.cdf` branch
  although `FlightLog` reads one, and called `sys.exit(0)` on an unknown
  type. It now compares the extension, accepts `.cdf`, and raises
  `ValueError`. `utils.regrid_data` likewise raises instead of exiting the
  interpreter.
- The flight's calibration source now decides the *coefficients*, not only
  the path. Temperature, wind and tail-number lookups all went through a
  process-global `Coef_Manager` built from the first directory configured, so
  `TableCalibration(directory=...)` and `flight.calibration_source = ...` had
  no effect on the numbers, `OnboardCalibration`'s refusal was never reached,
  and provenance reported the global directory. Every lookup now goes through
  `FlightLog.calibration_source`; `FlightLog` and `ProcessingConfig` accept
  `coefficient_dir=`. `utils.coef_manager`/`get_coef_manager` remain but warn
  `DeprecationWarning`. c1 files record the source's directory and class
  (`coefficient_directory`, `calibration_source`).
  Numbers are unchanged for a table without validity dates.
- Dated coefficient rows are now selected. No caller passed `when=`, so a
  sensor with two dated rows raised `AmbiguousCoefficients`. The flight's
  first valid timestamp (`FlightLog.start_time`) is threaded into every
  lookup and written as `coefficient_lookup_time`.
- An explicit `tail_number` (`FlightLog(tail_number=...)`,
  `ProcessingConfig.tail_number`) now wins. `Profile` always looked up the
  log's copterID and used the explicit value only if that lookup failed, so
  a flight forced to `N944UA` was processed and labelled as `N934UA`.
- c1 files no longer claim an RH bias correction that was never applied.
  `ProcessingConfig.bias_correction` has only ever been recorded; nothing
  applies it (that is a scientific decision not yet made), yet the file
  carried `rh_bias_correction` and its coefficients as if it had been. The
  file now says `rh_bias_correction_requested` and
  `rh_bias_correction_applied = 'no'`, and processing warns. The old
  `rh_bias_correction` attribute is gone.
- **Descent processing (`ascent=False`) works.** It raised `IndexError`
  (`regrid_base` read a third index the 2-tuple leg does not have), and would
  have gridded to nothing even past that because the grid and the level walk
  assumed altitude rises along the leg. A descent is now gridded in flight
  order: the first level is the top, altitude falls down the arrays, and the
  edge times increase along the leg so `regrid_data`'s (start, end] bins are
  unchanged. Levels sit on the same `base_start + n*res` lattice as an
  ascent, so an ascent and a descent can share a grid. `Profile._base_start`
  is now the lowest edge for either direction (it was `gridded_base[0]`).
  Ascents are bit-identical to before.
- **Pressure grids (`res_units='hPa'`/`'Pa'`) work.** `regrid_base` looked
  the GPS leg times up in the barometer's clock with `list.index`, which
  raised `ValueError`. They are now mapped to the nearest sample on the
  base's own clock (`utils.nearest_index`). The resolution is also converted
  to the base's units: 5 hPa against a barometer in Pa used to step 5 Pa.
  `profile_start_height` is a height in metres and is ignored, with a
  warning, on a pressure grid; pass `base_start` as a pressure instead.
- **Phantom legs are not emitted.** `FlightLog.find_legs` takes `min_extent`
  (metres, default `FlightLog.MIN_LEG_EXTENT` = 50; `ProcessingConfig
  .min_leg_extent`, 0 disables) and `ascent`, and drops peak-detected legs
  that do not climb (descend, for `ascent=False`) that far. It is measured
  in the direction processed because the legs are valley-peak-valley
  triples: flight616's 1.1 m wiggle is a good *descent* (its "end" is the
  bottom of the real descent) and the real leg's own descent is a 27 m dip.
  Evidence for 50 m: flight616 and the 2026 OK3DM logs that parse have real
  profiles of 76 m to 1420 m; the only other legs are 1-27 m. The legacy
  finder is not filtered. flight616 now yields one ascending profile, not
  two; profile 0 is unchanged. (The descent of that flight is the *second*
  leg, and it now yields one profile too.)
- `ProcessingConfig.min_levels` counts levels that hold data
  (`Profile.n_populated_levels`: the leg reached the level and a
  temperature sample falls in its time bin), not the length of the grid.
  On a common grid a short leg used to keep a full-length vertical
  coordinate and pass.
- `profiles_from_flight` and `process_flights` grid identically.
  `profile_start_height` was honoured only by `process_flights`, so
  flight616 gridded to 141 levels through one and 142 through the other;
  `Profile` now applies it itself.
- A leg that grids to no levels is never returned as a `Profile`
  (`Profile` raises `EmptyProfileError`, a `ValueError`;
  `profiles_from_flight` skips it with a warning). Such a `Profile` used to
  crash in `__str__` and the writers on `time[0]`.
- `Profile.__init__` no longer has its own leg finder. It calls
  `FlightLog.find_legs`, so a directly constructed `Profile` and
  `process_flights` pick the same legs (the private copy lacked the 5 s
  sensor settle). A `profile_num` past the legs found raises `IndexError`
  naming the count, instead of re-running the constructor without
  `index_list`/`raw_profile`/`metadata` and recursing until
  `RecursionError`; `profile_num < 1` raises `ValueError` rather than
  indexing from the end. `Profile` still does not plot detected peaks when
  `confirm_bounds=True` (its default), as before.
- `utils.identify_profile` passed `to_return` positionally into
  `confirm_bounds` when the user rejected a selection; it now forwards
  `confirm_bounds` and re-asks for the start height.

- `save_cfnetcdf` no longer raises without metadata (the README quick
  start); with neither metadata nor a path the name is derived from the input.
  `save_netcdf` and `save_cfnetcdf` now resolve to distinct names
  (`...ascent.nc` / `...cf_ascent.nc`), so calling both no longer overwrites.
- `save_cfnetcdf` CF/WMO attributes: both references are kept; `pressure` is
  `air_pressure` with no `positive`/`axis`; `altitude` uses the CF name
  `altitude` and is the only Z axis; `featureType = trajectory` is backed by a
  `trajectory_id` variable with `cf_role`; `time` is `f8` so sub-second bin
  centres are not truncated; an unresolved `copter_id`/`tail_number` no
  longer raises.
- `Profile.q` is a true kg/kg quantity (it was kg/kg magnitudes labelled
  g/kg). The writers convert explicitly and write `g/kg` instead of `gPerKg`.
- `coefficient_revision` is the last commit that changed `MasterCoefList.csv`
  (symlinks resolved, `unknown` if untracked, `-dirty` for uncommitted edits)
  instead of the HEAD of whatever repository encloses the directory. A new
  `coefficient_sha256` attribute records the table's hash.
- Writers use `platform.node()` instead of `os.uname()` (works on Windows)
  and no longer call the deprecated `datetime.utcnow()`/`utcfromtimestamp()`.
- `Profile.__init__` no longer overwrites its `file_path` argument;
  `Profile_Set` forwards `coefficient_dir`.
- BARO `Temp`/`GndTemp` are declared and written in degC (they were declared
  degF; ArduPilot logs degC). a0 files labelled `F` read back as degC.
- Samples logged before the GPS clock was set are dropped at parse time (the
  whole row) instead of carried as NaN times, which broke the a0 writer and
  datetime comparisons.
- Per-sensor calibrated series are slot-aligned: an absent middle sensor is an
  all-NaN slot instead of shifting later sensors down, which mis-paired
  thermistors with serials and misnumbered a0 `calib_temp<n>`/`calib_rh<n>`.
- CSV reader: all six IMU channels are read, down-velocity comes from `vz`,
  `FlightLog.data` is always defined, and the removed pandas option
  `infer_datetime_format` is no longer passed.
- Importing `profiles` no longer changes the host process: no global
  RuntimeWarning suppression, no UnitStrippedWarning-as-error, no pandas
  matplotlib converter registration, and `matplotlib.pyplot` is imported
  lazily where plots are drawn.
- lat/lon/alt_MSL are binned with the same (start, end] bins as every other
  gridded variable, so they always have the same length (the old
  `[start, end)` helper could drop the final bin).
- `process_flights` no longer parses the first flight twice (re-reading the
  BIN and, with `nc_level='low'`, writing the a0 file twice).
- `Meta`: `dronelogbook` is imported only when fetching from DroneLogbook;
  `read_file` applies every field from a flight file instead of stopping at
  the first replaced one, and no longer raises on non-string values or a
  missing timestamp.

- ESC RPM samples logged before the GPS clock was set are dropped, as for
  the other groups, instead of being stamped 1970.
- `FlightLog.is_equal` no longer raises `TypeError` on pint Quantities, and
  treats matching NaNs as equal.

### Changed

- **`.BIN` files are read directly.** `profiles/readers/` streams an
  ArduPilot log into memory instead of shelling it through a vendored
  `mavlogdump` fork that wrote a newline-delimited `.json` beside the
  original and re-read it. The 12 MB reference flight produced 83 MB of
  intermediate; an 18 MB OK3DM flight produced 104 MB. Old `.json` files
  still load. `profiles/mavlogdump_Profiles.py` is deleted (459 lines).
- **Parsing is schema-driven.** `profiles/schema.py` declares each message
  type's named fields and units; `profiles/parsing.py` has one loop that
  builds an **xarray Dataset** per group on a shared time coordinate,
  available as `FlightLog.data`. This replaces ~500 lines of per-type blocks
  that addressed every measurement by slot number. The positional tuples are
  derived from the Datasets and unchanged for now.
- **`Raw_Profile` is now `FlightLog`**, in `profiles/flight.py`. It holds a
  whole flight at native rate, not a profile, and the old name made it read
  as a variant of `Profile`. `profiles.Raw_Profile` still imports, with a
  `DeprecationWarning`.
- `xarray` is a new dependency.

- **`Thermo_Profile` and `Wind_Profile` are merged into `Profile`.** Both
  were handed the same grid and produced variables on it: `pres`, `alt` and
  `time` were computed three times from the same inputs, `lat`/`lon` twice
  by two different implementations, and 78 lines of `Profile` existed only
  to reconcile the two children's lengths. Use `compute_thermo()` and
  `compute_wind()`; `get_thermo_profile()` / `get_wind_profile()` remain as
  deprecated shims returning `self`.
- **Physics and QC are plain functions**, in `profiles/retrievals/`,
  `profiles/calibration.py` and `profiles/qc.py`, testable without
  constructing a Profile. The rotation-matrix tilt loop existed in three
  identical copies; the per-sensor calibration block in two.
- **`Profile_Set` is deprecated** in favour of
  `profiles.processing.process_flights()` with a `ProcessingConfig`. It
  still works and delegates. `Profile_Set.add_profile` and `read_netCDF`
  are deleted - both raised `TypeError` unconditionally.
- Leg detection moved onto `FlightLog.find_legs()`; it previously existed
  in both `Profile_Set.add_all_profiles` and `Profile.__init__`, and the
  copies had drifted (only one passed the fan-aspiration window).
- `Profile.__eq__` defined a nested `__lt__` and fell off the end returning
  `None`, so every Profile compared unequal to every other including
  itself. Full comparison set plus `__hash__` now.

- **Coefficient rows can carry `ValidFrom`/`ValidTo`** and are selected by
  flight date. Undated rows stay valid at all times, so existing tables work
  untouched. Recalibrating a sensor no longer forces a choice between
  destroying history and breaking every lookup.
- **Ambiguity is reported.** Duplicate coefficient rows raise
  `AmbiguousCoefficients` listing each candidate and its validity window.
  `get_tail_n` still returns the first matching row — so nothing changes
  numerically — but warns once per ID naming every candidate.
- **`CalibrationSource` has two real implementations.** `TableCalibration`
  for flights that logged sensor serials; `OnboardCalibration` for current
  firmware, which refuses `Imet`/`RH` table lookups rather than silently
  applying generic coefficients over the aircraft's own calibration. Wind
  stays a table lookup. `source_for_flight()` chooses from the log.
- **Bias corrections are a separate versioned layer** (`profiles/bias.py`).
  Both lab corrections are MATLAB surface fits, so one evaluator covers
  them: `pij` multiplies `rh**i * temp**j`. Each records its source and
  fitted range and reports extrapolation. **None is applied by default.**
- **Output records what produced it**: per-sensor coefficients and equation,
  how temperature was derived, wind coefficients, tail number, QC
  thresholds, the coefficient directory and its git revision, processing
  version and level, and any bias correction's full coefficient set.
- **QC flags reach the published c1 file** as CF-conventional `temp_qc` /
  `rh_qc` variables on a `sensor` dimension, with `flag_values` and
  `flag_meanings`. They previously existed only as group attributes on the
  `thermo_` intermediate and were dropped from c1 entirely.
- **QC thresholds are configurable** via `ProcessingConfig.qc_thresholds`
  instead of being literals with a `# TODO read these from file`.
- Configuration moved to `profiles/config.py`: explicit path → `$WXUAS_DIR`
  → `~/.wxuas`. `profiles.conf` still works, still takes precedence when
  set, and now warns.
- The Azure coefficient backend is deleted. Every line was commented out and
  the branch selecting it fell through to `pass`.
- Coefficient tables are indexed once at construction; `get_coefs` was doing
  a full DataFrame copy per sensor per profile.
- **One function builds every output file name** (`profiles/io/naming.py`).
  Five writers each built their own, with unintended differences: the
  per-variable c1 names put the ascent tag in a different position from the
  combined ones, and three of the five raised `IOError` without metadata
  while `FlightLog` quietly derived a name from the input.
- **`q` in `thermo_*` files was wrong by a factor of 1000.** MetPy returns
  kg/kg and the code attached a `gPerKg` label without rescaling.
  `Profile.save_netcdf` compensated with an explicit `* 1e3`; the
  per-variable writer did not, so that file claimed g/kg while holding
  kg/kg. Both now scale.

- `Profile.save_netcdf` no longer looks up a place name over the network by
  default; pass `lookup_place=True` to record `flight_location`. The lookup
  has a 5 s timeout and an identifying User-Agent, and a failure warns
  instead of aborting the save.
- `regrid_data_group` yields one entry per bin, including empty bins.
  `Profile.lat`/`lon`/`alt_MSL` shift by up to ~5 mm (lat) and ~0.4 m
  (alt_MSL) from the previous binning.
- `FlightLog._as_datetimes` raises `ValueError` on NaT instead of returning
  NaN.
- Docs workflow runs in `python:3.12-slim-bookworm` (was end-of-life
  `debian:buster-slim`), with Sphinx tooling installed from pip.

- A requested barometer instance or EKF core that is absent from the log
  raises `ValueError` listing the instances present, instead of a generic
  "no pres messages" failure. A log with only BARO `I=0` processed with the
  default `baro_instance=1` now errors rather than silently using `I=0`.

### Added

- Tests: a trimmed current-firmware fixture
  (`test/data/flight2862_ascent_trim.BIN`, 1.8 MB: BARO `I`, XKF1 `C`, ESC
  instances, no USER_SENSORS) with end-to-end `process_flights` coverage and
  onboard-calibration selection; `process_flights` is characterized against
  the flight616 p0 snapshot, so the v2 entry point is no longer held only by
  the deprecated `Profile_Set` harness; a0 reprocessing is checked for theta,
  T_d, mixing ratio and q. `pytest.ini` turns DeprecationWarnings attributed
  to `profiles` or the test suite into errors (third-party ones stay ignored).
- `ProcessingConfig.baro_instance` (default 1, provisional) and `ekf_core`
  (default 0), also accepted by `FlightLog`, choose the BARO `I` instance and
  XKF1 core `C` on current firmware. Logs with separate BARO/BAR2 messages
  ignore them (BAR2 is still preferred). The instance used is recorded as
  `baro_instance`, with `baro_message_type`, in a0 and c1 files.

- `test/` is a real suite: baseline characterization snapshots, grid
  invariants, and an a0 NetCDF round trip. Hermetic — uses `test/data/coefs`,
  never `~/.wxuas`, and makes no network calls.
- `.github/workflows/tests.yml` runs it on Python 3.10 and 3.12.
- A log carrying none of the required message types raises a `ValueError`
  naming what is missing, rather than failing on `len(None)` deeper in.

### Known issues not addressed here

- `Thermo_Profile`/`Wind_Profile._save_netCDF` raise `IOError` under
  `nc_level='low'` unless metadata or an explicit `.cdf` path is supplied.
  Pre-existing; `Raw_Profile` falls back to a filename derived from the input
  and does not.
- RH per-sensor corrections are not applied at all (see `rh_calib`).
- **Peak detection still reports phantom legs; they are now filtered.**
  `identify_profile_peaks` itself is unchanged (a one-metre prominence), so
  on flight616 it still finds a 2-second, 1.1 m "profile" at 08:07:51 on top
  of the real one. `FlightLog.find_legs` now discards legs that climb (or
  descend) less than `ProcessingConfig.min_leg_extent` - see Fixed. Calling
  `identify_profile_peaks` directly, or `find_legs(min_extent=0)`, still
  returns them.
- `_read_csv` is still the original hand-written parser. It has no test
  coverage and no sample data in the repo, so it was left alone rather than
  migrated to the schema on faith.
- `copterID` 1 resolves to four different tail numbers in the live
  `copterID.csv` (`FA3TANE3MF`, `FA3TANFWPA`, `FA3XEX7RKR`, `N944UA`) and
  `get_tail_n` returns whichever is listed first. Their wind coefficients
  differ (A 37.6 vs 32.8, B +6.8 vs -4.5), so retrieved wind speed depends
  on row order in a CSV. Every current OK3DM flight logs `SYSID_THISMAV = 1`.
- Current firmware no longer logs `USER_SENSORS` parameters, so sensor
  serials fall back to 0 and `temp_calib` recomputes temperature from
  resistance with *generic* coefficients - discarding the already-calibrated
  `T1..T4` the autopilot logs. This is the onboard-calibration question and
  belongs to Stage 5.

- The repo's top-level `coefs/MasterCoefList.csv` cannot process the repo's
  own test flight: five rows per IMET sensor differing only by `ScoopID`,
  which `get_coefs` does not key on. It also disagrees with the live
  `~/.wxuas` table on the coefficient values themselves.
- `plotting.contour_height_time` uses `Profile.gridded_base` (N+1 edges) as
  the vertical coordinate against N-point variables.
