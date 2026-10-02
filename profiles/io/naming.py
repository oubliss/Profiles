"""
Output file names, built in one place.

Five writers built their own, with small differences nobody intended:
FlightLog's a0 name omits the resolution (correctly - a0 is not gridded),
the per-variable c1 names insert "thermo_"/"wind_" before the level, and
the two combined c1 writers append the ascent tag in a different position
from the per-variable ones. Three of the five then raise IOError when
there is no metadata, while FlightLog quietly derives a name from the
input file. All of that is now described by one function.

Shape:

    <Location><resolution><platform>CMT<tag>.<level>.<YYYYMMDD.HHMMSS>.cdf
"""
import os


def output_name(meta, level, timestamp, resolution=None, tag='',
                extension='cdf'):
    """ Build an output file name from metadata.

    :param meta: a Meta object, or None
    :param str level: processing level, 'a0' / 'b1' / 'c1'
    :param str timestamp: already formatted as YYYYMMDD.HHMMSS
    :param resolution: vertical resolution magnitude, omitted for a0
    :param str tag: distinguishes products at the same level, e.g.
       'thermo_Ascending'
    :param str extension: file extension without the dot
    :rtype: str
    :raises ValueError: if meta is None
    """
    if meta is None:
        raise ValueError(
            'cannot build an output name without metadata; pass an explicit '
            'file path ending in .nc or .cdf instead')

    location = str(meta.get('location')).replace(' ', '')
    platform = str(meta.get('platform_id'))
    scale = '' if resolution is None else str(resolution)

    return f'{location}{scale}{platform}CMT{tag}.{level}.{timestamp}.{extension}'


def resolve(file_path, meta, level, timestamp, resolution=None, tag='',
            fallback=None):
    """ The path to write to.

    An explicit .nc or .cdf path always wins. Otherwise a name is built
    from metadata and placed beside file_path. If there is no metadata and
    a fallback is given, that is used.

    :param str file_path: the caller's path, or the input file's path
    :param meta: a Meta object, or None
    :param str level: processing level
    :param str timestamp: formatted YYYYMMDD.HHMMSS
    :param resolution: vertical resolution magnitude, or None
    :param str tag: product tag
    :param str fallback: used when there is no metadata
    :rtype: str
    """
    if file_path and (file_path.endswith('.nc') or file_path.endswith('.cdf')):
        return file_path

    if meta is None:
        if fallback is not None:
            return fallback
        raise IOError(
            'Please specify a file name ending in .nc or .cdf, or include '
            'metadata, when saving NetCDF output')

    name = output_name(meta, level, timestamp, resolution=resolution, tag=tag)
    return os.path.join(os.path.dirname(file_path), name)
