"""
What produced this file.

A c1 file used to record its own values and nothing about how they were
arrived at: not which coefficients were applied, not which QC thresholds
rejected a sensor, not which version of the coefficient tables was current.
Reprocessing an old flight and getting a different answer was undiagnosable.
"""
import os
import subprocess
from datetime import datetime

import profiles


def _submodule_revision(path):
    """ git revision of the SensorCoefficients checkout, if there is one.

    :rtype: str
    """
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'], cwd=str(path),
            capture_output=True, text=True, timeout=5, check=False)
        if result.returncode == 0:
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return 'unknown'


def coefficient_attributes(record, prefix='coef'):
    """ Flatten a calibration record into NetCDF attributes.

    :param dict record: as filled in by profiles.calibration
    :param str prefix: attribute name prefix
    :rtype: dict
    """
    attributes = {}
    for sensor, row in sorted(record.items()):
        if not isinstance(row, dict):
            attributes[f'{prefix}_{sensor}'] = str(row)
            continue
        for field in ('A', 'B', 'C', 'D', 'Equation', 'ValidFrom', 'ValidTo'):
            value = row.get(field)
            if value not in (None, '', 'na'):
                attributes[f'{prefix}_{sensor}_{field.lower()}'] = str(value)
    return attributes


def provenance_attributes(profile, coefficient_dir=None):
    """ Everything needed to reproduce this file.

    :param profile: the Profile being written
    :param coefficient_dir: directory the coefficient tables came from
    :rtype: dict
    """
    from profiles import config

    directory = config.coefficient_dir(coefficient_dir)

    attributes = {
        'processing_version': profiles.__version__,
        'processing_datetime': datetime.utcnow().isoformat() + 'Z',
        'processing_machine': os.uname().nodename,
        'copter_id': getattr(profile, 'copter_id', -999),
        'tail_number': str(getattr(profile, 'tail_number', 'unknown')),
        'coefficient_directory': str(directory),
        'coefficient_revision': _submodule_revision(directory),
        'vertical_resolution': str(profile.resolution),
        'leg': 'ascent' if profile.ascent else 'descent',
    }

    thresholds = getattr(profile, 'qc_thresholds', None)
    if thresholds:
        for variable, (max_bias, max_spread) in sorted(thresholds.items()):
            attributes[f'qc_{variable}_max_bias'] = float(max_bias)
            attributes[f'qc_{variable}_max_variability'] = float(max_spread)

    attributes.update(
        coefficient_attributes(getattr(profile, 'calibration_record', {})))

    correction = getattr(profile, 'bias_correction', None)
    if correction is not None:
        attributes.update(correction.provenance())

    return attributes
