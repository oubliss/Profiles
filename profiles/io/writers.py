"""
Processing levels and the variables that mark them.

a0  raw, native rate, no QC
b1  calibrated and QC-flagged, still native rate
c1  gridded profile, flags and provenance carried through
"""
import numpy as np

from profiles import qc

#: Recognised processing levels, in order.
LEVELS = ('a0', 'b1', 'c1')


def write_qc_variables(handle, profile, dimension='time'):
    """ Add per-sensor QC flags as CF-conventional variables.

    The flags were previously written as group attributes on the
    intermediate thermo_ file only, and dropped entirely from the c1 file
    that actually gets published - so a downstream user could not tell that
    a sensor had been removed from the ensemble mean, let alone which one.

    Flags are per sensor, not per level, so they go on their own dimension.

    :param handle: an open netCDF4.Dataset
    :param profile: the Profile being written
    :param str dimension: unused; kept for call-site symmetry
    """
    flags = {'temp': getattr(profile, 'temp_flags', None),
             'rh': getattr(profile, 'rh_flags', None)}
    flags = {name: values for name, values in flags.items()
             if values is not None}
    if not flags:
        return

    n_sensors = max(len(values) for values in flags.values())
    if 'sensor' not in handle.dimensions:
        handle.createDimension('sensor', n_sensors)

    values_list = sorted(qc.FLAG_MEANINGS)
    meanings = ' '.join(qc.FLAG_MEANINGS[value] for value in values_list)

    for name, values in flags.items():
        variable = handle.createVariable(f'{name}_qc', 'i1', ('sensor',))
        variable[:] = np.asarray(values, dtype='i1')
        variable.long_name = f'per-sensor QC flag for {name}'
        variable.flag_values = np.array(values_list, dtype='i1')
        variable.flag_meanings = meanings
        variable.comment = ('Sensors flagged non-zero were excluded from '
                            'the ensemble mean reported here.')
