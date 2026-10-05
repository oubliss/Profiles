"""
Batch processing: settings in one object, a function over a list of files.

This replaces Profile_Set. What that class actually provided was a bag of
settings applied uniformly plus one piece of real logic - standardising the
vertical grid's starting height across every file so profiles from one
mission share levels. Both are here, without a container class whose
contents are the only thing anyone ever wanted.
"""
import warnings
from dataclasses import dataclass, field
from typing import Optional

from profiles.Profile import Profile, EmptyProfileError
from profiles.flight import FlightLog


@dataclass
class ProcessingConfig:
    """ Settings applied uniformly to every flight in a batch.

    :var resolution: vertical grid spacing
    :var res_units: units of resolution, anything pint parses. Altitude
       units grid on altitude; pressure units grid on pressure.
    :var ascent: process the ascending leg; False for the descent
    :var dev: flag the data as from a developmental flight
    :var confirm_bounds: plot detected legs for a visual check
    :var profile_start_height: force the grid to start here, in metres. Left
       as None, the first flight's own starting height is adopted and then
       applied to the rest of the batch. Both profiles_from_flight and
       process_flights honour it. It is a height, so on a pressure grid
       (res_units 'hPa' or 'Pa') it is ignored with a warning and the first
       flight's own starting pressure is adopted instead.
    :var nc_level: 'low' writes per-object NetCDF files, None or 'none'
       writes nothing
    :var legacy_peak_id: use the pre-2021 altitude-threshold leg finder
    :var tail_number: override the airframe identity, needed when the log
       does not carry a usable vehicle ID
    :var wind_algorithm: 'linear' or 'quadratic'
    :var min_levels: discard profiles with fewer than this many levels that
       actually hold data. Levels the leg never reached (NaN pressure or
       altitude) do not count, so a short leg gridded onto a long common
       grid is dropped rather than passed for its NaN padding.
    :var min_leg_extent: metres. Legs that climb (descend, for
       ``ascent=False``) less than this are not treated as profiles. Peak
       detection reports small wiggles at the top of real profiles; the
       default of 50 m removes them (see FlightLog.find_legs). 0 disables
       the check. Not applied to the legacy finder.
    :var baro_instance: which barometer supplies pressure: the ``I`` field
       of current firmware's BARO messages. Default 1, the scoop
       barometer on OK3DM/CopterSonde airframes (confirm for others). A log
       without that instance raises an error naming those it has. Old logs
       with separate BARO and BAR2 messages and no instance field always
       use BAR2, and this setting does not apply to them.
    :var ekf_core: which EKF core supplies attitude and velocity (the ``C``
       field of XKF1); default 0, the primary core
    :var calibration: 'auto' decides per flight from whether the log
       reports sensor serial numbers; 'table' forces Steinhart-Hart from
       logged resistance; 'onboard' takes IMET.T as logged.
    """
    resolution: float = 10
    res_units: str = 'm'
    ascent: bool = True
    dev: bool = False
    confirm_bounds: bool = False
    profile_start_height: Optional[float] = None
    nc_level: Optional[str] = None
    legacy_peak_id: bool = False
    tail_number: Optional[str] = None
    wind_algorithm: str = 'linear'
    min_levels: int = 0
    min_leg_extent: float = 50.0
    calibration: str = 'auto'
    baro_instance: int = 1
    ekf_core: int = 0
    #: Directory of coefficient tables for every flight's calibration
    #: source, or None to resolve it through profiles.config.
    coefficient_dir: Optional[str] = None
    #: Per-variable (max spread of sensor means, max spread of sensor
    #: standard deviations). None keeps Profile.DEFAULT_QC_THRESHOLDS.
    qc_thresholds: Optional[dict] = None
    #: Name of a registered bias correction, or None. It is recorded in the
    #: output as requested-but-NOT-applied: nothing applies it yet, and
    #: whether to is a scientific decision - see profiles/bias.py.
    bias_correction: Optional[str] = None


@dataclass
class FlightResult:
    """What came of one input file."""
    path: str
    profiles: list = field(default_factory=list)
    error: Optional[Exception] = None

    @property
    def ok(self):
        return self.error is None


def profiles_from_flight(path, config, metadata=None, base_start=None,
                         flight=None):
    """ Every profile in one flight log.

    :param str path: the log
    :param ProcessingConfig config: settings
    :param metadata: a Meta object, or None
    :param base_start: vertical coordinate of the first grid edge, for
       matching an earlier flight's levels
    :param FlightLog flight: an already-parsed log, to avoid re-reading
    :rtype: list[Profile]
    """
    if flight is None:
        flight = FlightLog(path, config.dev, nc_level=config.nc_level,
                           metadata=metadata,
                           tail_number=config.tail_number,
                           calibration=config.calibration,
                           coefficient_dir=config.coefficient_dir,
                           baro_instance=config.baro_instance,
                           ekf_core=config.ekf_core)
    elif config.calibration != 'auto':
        # The caller parsed the log themselves; the config still decides.
        from profiles import Coef_Manager
        flight.calibration_source = Coef_Manager.source_for_flight(
            flight.serial_numbers, directory=config.coefficient_dir,
            mode=config.calibration)

    legs = flight.find_legs(
        legacy=config.legacy_peak_id,
        confirm_bounds=config.confirm_bounds,
        profile_start_height=config.profile_start_height,
        min_extent=config.min_leg_extent, ascent=config.ascent)

    profiles = []
    for number in range(1, len(legs) + 1):
        try:
            # Profile applies profile_start_height itself, so this entry
            # point and process_flights grid a flight identically.
            profile = Profile(
                path, config.resolution, config.res_units, number,
                ascent=config.ascent, dev=config.dev,
                confirm_bounds=config.confirm_bounds, index_list=legs,
                raw_profile=flight,
                profile_start_height=config.profile_start_height,
                nc_level=config.nc_level, base_start=base_start)
        except EmptyProfileError as exc:
            warnings.warn(f'{path}: skipping profile {number}: {exc}',
                          UserWarning, stacklevel=2)
            continue

        if config.qc_thresholds:
            profile.qc_thresholds.update(config.qc_thresholds)
        if config.bias_correction:
            from profiles import bias
            warnings.warn(
                f'bias_correction={config.bias_correction!r} is recorded in '
                f'the output but NOT applied to the data; no bias '
                f'correction is applied by this pipeline.', UserWarning,
                stacklevel=2)
            profile.bias_correction = bias.get(config.bias_correction)

        profiles.append(profile)

    return profiles


def process_flights(paths, config=None, metadata=None, compute=True,
                    on_error='collect'):
    """ Process a batch of flights onto one common vertical grid.

    :param paths: log file paths
    :param ProcessingConfig config: settings; defaults are used if omitted
    :param metadata: a Meta object, or a {path: Meta} mapping, or None
    :param bool compute: also run compute_thermo() and compute_wind()
    :param str on_error: 'collect' records the exception against the file
       and carries on; 'raise' stops at the first failure
    :rtype: list[FlightResult]
    """
    config = config or ProcessingConfig()
    results = []

    base_start = None

    for path in paths:
        meta = metadata.get(path) if isinstance(metadata, dict) else metadata
        try:
            profiles = profiles_from_flight(path, config, metadata=meta,
                                            base_start=base_start)

            # Adopt the first flight's grid so later ones share its levels.
            # A profile_start_height already fixed it, so there is nothing
            # to adopt (and nothing to re-run) in that case.
            if base_start is None and profiles:
                base_start = profiles[0]._base_start
                if not profiles[0]._base_start_given:
                    # The log is already parsed; only the grid changes.
                    profiles = profiles_from_flight(
                        path, config, metadata=meta, base_start=base_start,
                        flight=profiles[0]._raw_profile)

            profiles = [p for p in profiles
                        if p.n_populated_levels >= config.min_levels]

            if compute:
                for profile in profiles:
                    profile.compute_thermo()
                    profile.compute_wind(algorithm=config.wind_algorithm)

            results.append(FlightResult(path, sorted(profiles)))

        except Exception as exc:
            if on_error == 'raise':
                raise
            results.append(FlightResult(path, [], error=exc))

    return results


def all_profiles(results):
    """ Flatten successful results into one time-ordered list.

    :param results: what process_flights returned
    :rtype: list[Profile]
    """
    return sorted(p for result in results if result.ok
                  for p in result.profiles)
