"""
Batch processing: settings in one object, a function over a list of files.

This replaces Profile_Set. What that class actually provided was a bag of
settings applied uniformly plus one piece of real logic - standardising the
vertical grid's starting height across every file so profiles from one
mission share levels. Both are here, without a container class whose
contents are the only thing anyone ever wanted.
"""
from dataclasses import dataclass, field
from typing import Optional

from profiles.Profile import Profile
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
       applied to the rest of the batch.
    :var nc_level: 'low' writes per-object NetCDF files, None or 'none'
       writes nothing
    :var legacy_peak_id: use the pre-2021 altitude-threshold leg finder
    :var tail_number: override the airframe identity, needed when the log
       does not carry a usable vehicle ID
    :var wind_algorithm: 'linear' or 'quadratic'
    :var min_levels: discard profiles gridding to fewer levels than this.
       Peak detection reports spurious profiles from small wiggles at the
       top of real ones - see CHANGELOG.
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
                           tail_number=config.tail_number)

    legs = flight.find_legs(
        legacy=config.legacy_peak_id,
        confirm_bounds=config.confirm_bounds,
        profile_start_height=config.profile_start_height)

    profiles = []
    for number in range(1, len(legs) + 1):
        profiles.append(Profile(
            path, config.resolution, config.res_units, number,
            ascent=config.ascent, dev=config.dev,
            confirm_bounds=config.confirm_bounds, index_list=legs,
            raw_profile=flight,
            profile_start_height=config.profile_start_height,
            nc_level=config.nc_level, base_start=base_start))

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
    if config.profile_start_height is not None:
        from profiles.unit_registry import units
        base_start = config.profile_start_height * units.m

    for path in paths:
        meta = metadata.get(path) if isinstance(metadata, dict) else metadata
        try:
            profiles = profiles_from_flight(path, config, metadata=meta,
                                            base_start=base_start)

            # Adopt the first flight's grid so later ones share its levels.
            if base_start is None and profiles:
                base_start = profiles[0]._base_start
                profiles = profiles_from_flight(path, config, metadata=meta,
                                                base_start=base_start)

            profiles = [p for p in profiles
                        if len(p.gridded_centers) >= config.min_levels]

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
