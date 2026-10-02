"""
Read ArduPilot logs straight into memory.

The package used to shell a .BIN through a vendored fork of pymavlink's
mavlogdump CLI, which wrote a newline-delimited .json file next to the
original and then re-read it line by line. For the 12 MB reference flight
that meant writing 83 MB to disk to parse 12; for the 18-58 MB OK3DM flights
it was proportionally worse. Nothing downstream wanted the .json.

iter_bin yields exactly the dict shape that .json held, so the two paths
share one parser and old .json files still load.
"""
import array
import json
import os

# Message types the parser understands. Everything else is read and
# discarded, so filtering here avoids decoding ~95% of a typical log.
# PARM carries sensor serial numbers and the vehicle ID; FMT is consumed by
# the DF reader itself and never surfaces to us.
WANTED_TYPES = ['PARM', 'EV', 'MSG', 'IMET', 'RHUM', 'POS', 'BARO', 'BAR2',
                'NKF1', 'XKF1', 'WIND', 'ESC', 'IMU']


def _to_text(value):
    """Decode a byte string without ever raising."""
    if isinstance(value, bytes):
        return value.decode(errors='backslashreplace')
    return value


def _normalise(message):
    """Convert a pymavlink message to the plain dict the parser expects."""
    data = message.to_dict()
    data.pop('mavpackettype', None)

    for key, value in data.items():
        if isinstance(value, array.array):
            data[key] = list(value)
        elif isinstance(value, bytes):
            data[key] = _to_text(value)

    return {'meta': {'type': message.get_type(),
                     'timestamp': getattr(message, '_timestamp', 0.0)},
            'data': data}


def iter_bin(file_path, types=WANTED_TYPES):
    """ Yield normalised messages from an ArduPilot .BIN.

    Filtering happens here rather than via ``recv_match(type=...)``. That
    looks like the obvious optimisation and is wrong: DFReader implements it
    by skipping messages at the binary level, and skipped messages no longer
    contribute to its GPS time-base estimation. On the 2021 reference flight
    it shifted 17,073 timestamps by up to 17 ms, which is enough to move
    samples between averaging bins. Decoding everything and discarding what
    we do not want is both faithful and, measurably, barely slower.

    :param str file_path: path to the .BIN
    :param list types: message types to keep, or None for all
    :rtype: Iterator[dict]
    """
    from pymavlink import mavutil

    wanted = set(types) if types else None

    log = mavutil.mavlink_connection(file_path)
    try:
        while True:
            message = log.recv_match()
            if message is None:
                return
            if wanted is None or message.get_type() in wanted:
                yield _normalise(message)
    finally:
        close = getattr(log, 'close', None)
        if close is not None:
            close()


def iter_json(file_path, types=None):
    """ Yield messages from a .json dump written by an older version.

    :param str file_path: path to the newline-delimited .json
    :param list types: message types to keep, or None for all
    :rtype: Iterator[dict]
    """
    wanted = set(types) if types else None
    with open(file_path, 'r') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            message = json.loads(line)
            if wanted is None or message['meta']['type'] in wanted:
                yield message


def iter_messages(file_path, types=WANTED_TYPES):
    """ Yield normalised messages from a .BIN or .json, chosen by extension.

    :param str file_path: path to the log
    :param list types: message types to keep, or None for all
    :rtype: Iterator[dict]
    :raises ValueError: if the extension is not recognised
    """
    extension = os.path.splitext(file_path)[1].lower()

    if extension == '.bin':
        return iter_bin(file_path, types=types)
    if extension == '.json':
        return iter_json(file_path, types=types)

    raise ValueError(
        f'{file_path!r} is not a log this reader handles (expected .BIN or '
        f'.json, got {extension!r})')
