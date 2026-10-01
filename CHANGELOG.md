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

### Added

- `test/` is a real suite: baseline characterization snapshots, grid
  invariants, and an a0 NetCDF round trip. Hermetic — uses `test/data/coefs`,
  never `~/.wxuas`, and makes no network calls.
- `.github/workflows/tests.yml` runs it on Python 3.10 and 3.12.

### Known issues not addressed here

- The repo's top-level `coefs/MasterCoefList.csv` cannot process the repo's
  own test flight: five rows per IMET sensor differing only by `ScoopID`,
  which `get_coefs` does not key on. It also disagrees with the live
  `~/.wxuas` table on the coefficient values themselves.
- `plotting.contour_height_time` uses `Profile.gridded_base` (N+1 edges) as
  the vertical coordinate against N-point variables.
- `utils.rh_calib` overwrites its coefficient lookup with `offset = 0`.
- `utils._bias` / `_s_dev` can over-flag; the return and the flagging both sit
  inside the sensor loop while `max_diff` accumulates across it.
- `nc_level` means three different things across `Raw_Profile`,
  `Thermo_Profile` and `Wind_Profile`.
