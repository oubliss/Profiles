"""
Turn a stream of log messages into xarray Datasets, driven by schema.GROUPS.

One loop replaces the per-message-type blocks Raw_Profile used to carry. Each
group becomes a Dataset whose variables are named, carry a ``units``
attribute, and share one time coordinate - so nothing downstream addresses a
measurement by its position in a tuple.

Timestamps come from the message envelope rather than the TimeUS field, which
is what the hand-written parser did; TimeUS is microseconds since boot, while
the envelope has already been resolved to UTC against the GPS fixes.
"""
from datetime import datetime

import numpy as np
import xarray as xr

from profiles import schema

#: Timestamps earlier than this mean the clock had not yet been set from GPS.
#: The sample is kept so the series stays aligned, but its time is unusable.
MIN_VALID_YEAR = 2000


class GroupAccumulator:
    """Collects one schema Group's fields across a message stream."""

    def __init__(self, group):
        self.group = group
        self.active_type = None
        self.columns = {field.name: [] for field in group.fields}
        self.times = []

    def _preference(self, message_type):
        try:
            return self.group.types.index(message_type)
        except ValueError:
            return len(self.group.types)

    def accepts(self, message_type):
        """Should this message feed the group, and does it supersede earlier ones?

        Returns (accept, reset). ``reset`` is True when a more-preferred
        message type has appeared and everything collected so far belongs to
        a source we are abandoning.
        """
        if message_type not in self.group.types:
            return False, False

        if self.group.select == 'union' or self.active_type is None:
            return True, False

        if message_type == self.active_type:
            return True, False

        if self._preference(message_type) < self._preference(self.active_type):
            return True, True

        return False, False

    def add(self, message):
        message_type = message['meta']['type']
        accept, reset = self.accepts(message_type)
        if not accept:
            return

        if reset:
            self.columns = {field.name: [] for field in self.group.fields}
            self.times = []

        self.active_type = message_type
        data = message['data']

        for field in self.group.fields:
            # A field the firmware of the day did not log shows as NaN for
            # the whole flight rather than shortening the series.
            self.columns[field.name].append(data.get(field.source, np.nan))

        self.times.append(_timestamp(message))

    def to_dataset(self):
        """ Build the Dataset, or None if nothing was collected.

        :rtype: xarray.Dataset or None
        """
        if not self.times:
            return None

        dimension = f'{self.group.name}_time'
        variables = {}
        for field in self.group.fields:
            values = np.asarray(self.columns[field.name], dtype=float)
            attrs = {'units': field.units} if field.units else {}
            variables[field.name] = (dimension, values, attrs)

        return xr.Dataset(
            variables,
            coords={dimension: np.asarray(self.times, dtype='datetime64[ns]')},
            attrs={'source_message_type': self.active_type,
                   'group': self.group.name})


def _timestamp(message):
    """UTC datetime for a message, or NaT if the clock was not yet set."""
    moment = datetime.utcfromtimestamp(message['meta']['timestamp'])
    if moment.year < MIN_VALID_YEAR:
        return np.datetime64('NaT')
    return np.datetime64(moment, 'ns')


def _motor_index(message):
    """Which motor an ESC message belongs to."""
    return int(message['data'].get('Instance', 0)) % schema.N_MOTORS


def parse(messages):
    """ Build every Dataset the schema describes from one pass over a log.

    :param iterable messages: normalised messages from profiles.readers
    :rtype: dict
    :return: {
        'groups': {name: xarray.Dataset},
        'serial_numbers': dict,
        'events': (ids, times) or None,
        'messages': (text, times) or None,
        'rpm': (per-motor lists, times) or None,
       }
    """
    accumulators = {group.name: GroupAccumulator(group)
                    for group in schema.GROUPS}

    serial_numbers = {}
    event_ids, event_times = [], []
    message_text, message_times = [], []
    rpm = {i: [] for i in range(schema.N_MOTORS)}
    rpm_times = []

    for message in messages:
        message_type = message['meta']['type']

        if message_type == 'PARM':
            _read_parameter(message, serial_numbers)
            continue

        if message_type == 'EV':
            event_ids.append(message['data']['Id'])
            event_times.append(_timestamp(message))
            continue

        if message_type == 'MSG':
            message_text.append(message['data']['Message'])
            message_times.append(_timestamp(message))
            continue

        if message_type == 'ESC':
            # One message per motor per timestep, so the time series is
            # N_MOTORS times longer than any single motor's and gets
            # averaged down once the pass is over.
            rpm[_motor_index(message)].append(message['data'].get('RPM', np.nan))
            rpm_times.append(message['meta']['timestamp'])
            continue

        if message_type == 'IMU':
            # Only the first IMU. Older logs name them IMU/IMU2/IMU3
            # instead of carrying an instance field.
            if message['data'].get('I', 0) != 0:
                continue

        for accumulator in accumulators.values():
            accumulator.add(message)

    groups = {name: accumulator.to_dataset()
              for name, accumulator in accumulators.items()}

    return {'groups': {k: v for k, v in groups.items() if v is not None},
            'serial_numbers': serial_numbers,
            'events': (event_ids, event_times) if event_ids else None,
            'messages': (message_text, message_times) if message_text else None,
            'rpm': _finish_rpm(rpm, rpm_times)}


def _read_parameter(message, serial_numbers):
    """Pull vehicle identity and sensor serials out of a PARM message."""
    name = str(message['data'].get('Name', ''))
    value = message['data'].get('Value')

    if schema.VEHICLE_ID_PARAM in name:
        serial_numbers['copterID'] = value

    elif schema.SENSOR_SERIAL_PARAM in name:
        # USER_SENSORS1..4 are the temperature sensors, 5..8 the humidity
        # ones. Current firmware calibrates onboard and no longer logs these.
        index = int(name[-1])
        if index <= schema.N_SENSORS:
            serial_numbers[f'imet{index}'] = int(value)
        elif index <= 2 * schema.N_SENSORS:
            serial_numbers[f'rh{index - schema.N_SENSORS}'] = int(value)


def _finish_rpm(rpm, rpm_times):
    """Collapse per-motor ESC records onto one averaged time series."""
    if not rpm_times:
        return None

    lengths = {len(values) for values in rpm.values() if values}
    if not lengths:
        return None

    # Motors can report an unequal number of times if the log is cut
    # mid-sequence; trim to the shortest so the reshape below is square.
    n_messages = min(lengths)
    columns = [rpm[i][:n_messages] for i in range(schema.N_MOTORS)]

    usable = n_messages * schema.N_MOTORS
    stamps = np.reshape(rpm_times[:usable], (n_messages, schema.N_MOTORS))
    averaged = [np.datetime64(datetime.utcfromtimestamp(float(t)), 'ns')
                for t in np.nanmean(stamps, axis=1)]

    return columns + [averaged]
