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

### Added

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
- **Peak detection reports phantom profiles.** `identify_profile_peaks`
  calls `find_peaks(alts, prominence=1)` - a one-metre prominence. On
  flight616 the second "profile" it finds runs 08:07:51 to 08:07:53, two
  seconds, over one metre of altitude (1746.8 to 1747.8 m): a wiggle at the
  top of the real profile. With `profile_start_height` forced to 350 m the
  grid is built from 350 to 1740 m regardless, so 119 of its 140 levels
  fall past the end of that two-second leg and every variable is NaN there.
  The `if len(p.gridded_times) > 3` guard in the processing scripts passes,
  so a mostly-empty second c1 file gets written. `ProcessingConfig
  .min_levels` is a stopgap; a real fix needs a minimum profile depth in
  the detector, which changes which profiles get emitted.
- `Thermo_Profile.q` carried kg/kg magnitudes labelled `gPerKg`.
  `Profile.save_netcdf` compensated with `* 1e3` but the per-variable
  writer did not, so `thermo_*` files are mislabelled by 1000. Unchanged
  pending Stage 5's single writer.
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
