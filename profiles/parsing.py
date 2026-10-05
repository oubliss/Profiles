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
from datetime import datetime, timezone

import numpy as np
import xarray as xr

from profiles import schema

#: Timestamps earlier than this mean the clock had not yet been set from GPS.
#: Such a sample has no usable time, so it is dropped (whole row, so a group's
#: variables stay aligned with its time coordinate) rather than carried as NaT:
#: NaT cannot be written to an a0 file and breaks every time comparison.
MIN_VALID_YEAR = 2000


class GroupAccumulator:
    """Collects one schema Group's fields across a message stream."""

    def __init__(self, group, keep=None):
        self.group = group
        # Which instance to retain, if the group selects one at all.
        self.keep = keep if keep is not None else (
            group.instance.keep if group.instance else None)
        self.active_type = None
        self.active_instance = None
        #: Every instance value seen on an accepted message type, so a
        #: request for one the log lacks can name those it has.
        self.seen_instances = set()
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

        data = message['data']
        instance = None
        if self.group.instance is not None:
            instance = data.get(self.group.instance.field)
            if instance is not None:
                self.seen_instances.add(int(instance))
                if instance != self.keep:
                    return

        if reset:
            self.active_instance = None
            self.columns = {field.name: [] for field in self.group.fields}
            self.times = []

        self.active_type = message_type
        if instance is not None:
            self.active_instance = int(instance)

        for field in self.group.fields:
            # A field the firmware of the day did not log shows as NaN for
            # the whole flight rather than shortening the series.
            self.columns[field.name].append(data.get(field.source, np.nan))

        self.times.append(_timestamp(message))

    def to_dataset(self):
        """ Build the Dataset, or None if nothing was collected.

        :rtype: xarray.Dataset or None
        """
        times = np.asarray(self.times, dtype='datetime64[ns]')
        valid = ~np.isnat(times)
        if not valid.any():
            return None
        times = times[valid]

        dimension = f'{self.group.name}_time'
        variables = {}
        for field in self.group.fields:
            values = np.asarray(self.columns[field.name], dtype=float)[valid]
            attrs = {'units': field.units} if field.units else {}
            variables[field.name] = (dimension, values, attrs)

        attrs = {'source_message_type': self.active_type,
                 'group': self.group.name}
        if self.active_instance is not None:
            # Lets provenance say which physical sensor or core this is.
            attrs['source_instance'] = self.active_instance

        return xr.Dataset(
            variables,
            coords={dimension: times},
            attrs=attrs)


def _utc_from_timestamp(seconds):
    """ Naive UTC datetime from epoch seconds (utcfromtimestamp is
    deprecated since Python 3.12)."""
    return datetime.fromtimestamp(seconds, timezone.utc).replace(tzinfo=None)


def _timestamp(message):
    """UTC datetime for a message, or NaT if the clock was not yet set."""
    moment = _utc_from_timestamp(message['meta']['timestamp'])
    if moment.year < MIN_VALID_YEAR:
        return np.datetime64('NaT')
    return np.datetime64(moment, 'ns')


class EscAccumulator:
    """Collects ESC RPM, one message per motor per timestep.

    Current firmware numbers ESC instances by output channel, so a
    CopterSonde logs 8, 9, 11, 12 rather than 0..3. The motors are therefore
    taken to be the sorted set of Instance values seen, numbered 1..N, rather
    than inferred from the value itself.

    A timestep is a run of messages with no instance repeated: the next
    message for an instance already in the current step opens a new one.
    That tolerates a dropped message (the step is padded with NaN) where
    assuming a fixed number of messages per step would shift every later
    sample onto the wrong motor.
    """

    def __init__(self):
        self.steps = []       # [{instance: rpm}]
        self.stamps = []      # [[timestamp, ...]] parallel to steps

    def add(self, message):
        instance = int(message['data'].get('Instance', 0))
        if not self.steps or instance in self.steps[-1]:
            self.steps.append({})
            self.stamps.append([])
        self.steps[-1][instance] = message['data'].get('RPM', np.nan)
        self.stamps[-1].append(_timestamp(message))

    def finish(self):
        """[rpm per motor..., times], or None if no ESC was logged."""
        if not self.steps:
            return None

        instances = sorted({i for step in self.steps for i in step})
        columns = [[step.get(i, np.nan) for step in self.steps]
                   for i in instances]
        # A step is stamped with the mean of its messages' times. One whose
        # messages all predate the GPS clock has no time and is dropped, as
        # in every other group; a mixed step averages the valid ones.
        averaged = []
        for stamps in self.stamps:
            valid = [t for t in stamps if not np.isnat(t)]
            if valid:
                ticks = np.array(valid, dtype='datetime64[ns]').astype('int64')
                averaged.append(np.datetime64(int(ticks.mean()), 'ns'))
            else:
                averaged.append(np.datetime64('NaT'))
        keep = np.array([not np.isnat(t) for t in averaged])
        if not keep.any():
            return None
        columns = [np.asarray(column, dtype=float)[keep].tolist()
                   for column in columns]
        return columns + [[t for t, k in zip(averaged, keep) if k]]


def _drop_unset_clock(values, times):
    """Remove entries stamped before the GPS clock was set (NaT)."""
    kept = [(value, time) for value, time in zip(values, times)
            if not np.isnat(time)]
    return [value for value, _ in kept], [time for _, time in kept]


def parse(messages, instances=None):
    """ Build every Dataset the schema describes from one pass over a log.

    :param iterable messages: normalised messages from profiles.readers
    :param dict instances: {group name: instance value to keep}, overriding
        the default in each group's ``schema.Instance`` - for example
        ``{'pres': 0}`` if another airframe has its scoop barometer on
        BARO instance 0. Raises ValueError if a required group is logged
        with instance numbers but not the requested one. Logs whose
        messages carry no instance field (old firmware, BAR2) ignore it.
    :rtype: dict
    :return: {
        'groups': {name: xarray.Dataset},
        'serial_numbers': dict,
        'events': (ids, times) or None,
        'messages': (text, times) or None,
        'rpm': (per-motor lists, times) or None,
       }
    """
    instances = instances or {}
    accumulators = {group.name: GroupAccumulator(group, instances.get(group.name))
                    for group in schema.GROUPS}

    serial_numbers = {}
    event_ids, event_times = [], []
    message_text, message_times = [], []
    esc = EscAccumulator()

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
            esc.add(message)
            continue

        for accumulator in accumulators.values():
            accumulator.add(message)

    event_ids, event_times = _drop_unset_clock(event_ids, event_times)
    message_text, message_times = _drop_unset_clock(message_text,
                                                    message_times)

    for name in schema.REQUIRED_GROUPS:
        _check_instance(accumulators.get(name))

    groups = {name: accumulator.to_dataset()
              for name, accumulator in accumulators.items()}

    return {'groups': {k: v for k, v in groups.items() if v is not None},
            'serial_numbers': serial_numbers,
            'events': (event_ids, event_times) if event_ids else None,
            'messages': (message_text, message_times) if message_text else None,
            'rpm': esc.finish()}


def _check_instance(accumulator):
    """Refuse to fall back when the requested instance is not in the log.

    Raised only if the log numbers its copies at all: older firmware logs
    one unnumbered message (or BAR2 for the second barometer), and those are
    accepted whatever instance was asked for.
    """
    if accumulator is None or accumulator.group.instance is None:
        return
    seen = accumulator.seen_instances
    if seen and accumulator.keep not in seen:
        raise ValueError(
            f"requested {accumulator.group.name} instance "
            f"{accumulator.keep} (message field "
            f"{accumulator.group.instance.field!r}) is not in the log; "
            f"instances present: {sorted(seen)}")


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
