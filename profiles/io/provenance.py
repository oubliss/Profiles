"""
What produced this file.

A c1 file used to record its own values and nothing about how they were
arrived at: not which coefficients were applied, not which QC thresholds
rejected a sensor, not which version of the coefficient tables was current.
Reprocessing an old flight and getting a different answer was undiagnosable.
"""
import hashlib
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import profiles


def _git(directory, *args):
    """ Run git in directory; stdout stripped, or None on any failure.

    :rtype: str or None
    """
    try:
        result = subprocess.run(
            ['git', '-C', str(directory), *args],
            capture_output=True, text=True, timeout=5, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def coefficient_revision(csv_path):
    """ The revision of a coefficient table, '-dirty' if edited since.

    The revision is the last commit that changed the file itself, not the
    HEAD of whatever repository happens to enclose it: the coefficient
    tables usually live inside some other checkout (this one, for the test
    data), and that checkout's HEAD says nothing about them. Symlinks are
    followed, so a table linked in from a SensorCoefficients checkout is
    asked about in that checkout.

    :param csv_path: the table, possibly a symlink
    :return: a commit hash, optionally suffixed '-dirty', or 'unknown' if
       the file is not tracked by any repository
    :rtype: str
    """
    try:
        real = Path(csv_path).resolve(strict=True)
    except (OSError, RuntimeError):
        return 'unknown'

    if _git(real.parent, 'rev-parse', '--show-toplevel') is None:
        return 'unknown'
    # Owned means tracked there. A file merely sitting under a repository
    # (untracked, ignored) has no revision of its own.
    if _git(real.parent, 'ls-files', '--error-unmatch', '--',
            real.name) is None:
        return 'unknown'
    commit = _git(real.parent, 'log', '-1', '--format=%H', '--', real.name)
    if not commit:
        return 'unknown'        # staged but never committed
    if _git(real.parent, 'status', '--porcelain', '--', real.name):
        commit += '-dirty'
    return commit


def coefficient_sha256(csv_path):
    """ SHA-256 of the table itself, which identifies it with no repository.

    :rtype: str
    """
    try:
        digest = hashlib.sha256()
        with open(csv_path, 'rb') as handle:
            for block in iter(lambda: handle.read(1 << 16), b''):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
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
    :param coefficient_dir: directory the coefficient tables came from.
       Taken from the flight's calibration source when omitted - that is
       the directory the numbers were actually looked up in, which the
       process-wide configuration can disagree with.
    :rtype: dict
    """
    from profiles import config

    flight = getattr(profile, '_raw_profile', None)
    source = getattr(flight, 'calibration_source', None)
    if coefficient_dir is None and source is not None:
        try:
            coefficient_dir = source.directory
        except Exception:                 # provenance must never break
            coefficient_dir = None
    directory = config.coefficient_dir(coefficient_dir)
    table = Path(directory) / config.COEF_FILE

    attributes = {
        'processing_version': profiles.__version__,
        'processing_datetime': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S.%f') + 'Z',
        'processing_machine': platform.node(),
        'copter_id': getattr(profile, 'copter_id', None) or -999,
        'tail_number': str(getattr(profile, 'tail_number', None)
                           or 'unknown'),
        'calibration_source': type(source).__name__ if source else 'unknown',
        'coefficient_directory': str(directory),
        'coefficient_revision': coefficient_revision(table),
        'coefficient_sha256': coefficient_sha256(table),
        'vertical_resolution': str(profile.resolution),
        'leg': 'ascent' if profile.ascent else 'descent',
    }

    instance = getattr(flight, 'baro_instance', None)
    # Absent for logs that number no barometers (BARO/BAR2 firmware).
    attributes['baro_message_type'] = str(getattr(flight, 'baro', 'unknown'))
    if instance is not None:
        attributes['baro_instance'] = int(instance)

    when = getattr(flight, 'start_time', None)
    if when is not None:
        # The date dated coefficient rows were selected by.
        attributes['coefficient_lookup_time'] = when.isoformat()

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
