"""
What each log message type contributes, declared once.

Raw_Profile used to carry a hand-written block per message type - roughly 800
lines of "make a fixed-length list, build a {field: index} map, append or
append NaN, then attach units by position". Every variable was addressed by
its slot number, which is how rotation[6] came to be read for all three
position components and assigned three times on read-back.

Here a message type maps to named variables with units. One generic loop in
Raw_Profile drives every type, and nothing downstream has to know an index.

Sensor-numbered fields are written with a ``{n}`` placeholder and expanded
over ``count``: ``'T{n}' -> 'temp{n}'`` with count=4 gives T1->temp1 .. T4->temp4.
"""
from collections import namedtuple

#: One logged field. ``source`` is the key in the message's data dict,
#: ``name`` the variable it becomes, ``units`` anything pint can parse
#: (None for a flag or a bare ratio).
Field = namedtuple('Field', 'source name units')

#: A group of variables sharing one time coordinate.
#:
#: ``types`` lists the log message types that can feed the group and
#: ``select`` says what to do when a log carries more than one of them:
#:
#: ``'preferred'``
#:     The types are alternatives, best first. Collection starts with
#:     whichever appears first in the file and restarts from empty if a
#:     better one shows up later. flight616 logs both BARO and BAR2, and only
#:     the BAR2 record is kept.
#: ``'union'``
#:     The types are interchangeable and are appended to one series. No file
#:     seen so far carries both NKF1 and XKF1, so this is untested in
#:     practice; it preserves what the hand-written parser did, which was to
#:     merge them.
Group = namedtuple('Group', 'name types fields select')


def _numbered(source, name, units, count):
    """Expand a {n}-templated field over sensor numbers 1..count."""
    return [Field(source.format(n=i), name.format(n=i), units)
            for i in range(1, count + 1)]


#: Sensors of each type the CopterSonde can carry. The scoop has three
#: populated positions today; the fourth is logged as zeros.
N_SENSORS = 4

#: Motors on the airframe. ESC messages arrive one per motor per timestep,
#: tagged with Instance.
N_MOTORS = 4

TEMPERATURE = Group(
    name='temp',
    types=('IMET',),
    select='union',
    fields=(_numbered('T{n}', 'temp{n}', 'kelvin', N_SENSORS)
            + _numbered('R{n}', 'resi{n}', 'ohm', N_SENSORS)
            + [Field('Fan', 'fan_flag', None)]))

HUMIDITY = Group(
    name='rh',
    types=('RHUM',),
    select='union',
    fields=tuple(_numbered('H{n}', 'rh{n}', 'percent', N_SENSORS)
                 + _numbered('T{n}', 'temp_rh{n}', 'kelvin', N_SENSORS)))

POSITION = Group(
    name='pos',
    types=('POS',),
    select='union',
    fields=(Field('Lat', 'lat', 'degree'),
            Field('Lng', 'lon', 'degree'),
            Field('Alt', 'alt_MSL', 'meter'),
            Field('RelHomeAlt', 'alt_rel_home', 'meter'),
            Field('RelOriginAlt', 'alt_rel_orig', 'meter')))

#: BAR2 wins where both are present - it is the external barometer on the
#: CopterSonde, BARO the autopilot's internal one.
PRESSURE = Group(
    name='pres',
    types=('BAR2', 'BARO'),
    select='preferred',
    fields=(Field('Press', 'pres', 'pascal'),
            Field('Temp', 'temp', 'degF'),
            Field('GndTemp', 'ground_temp', 'degF'),
            Field('Alt', 'alt', 'meter')))

#: Roll, pitch and yaw are logged in degrees already; do not convert.
ROTATION = Group(
    name='rotation',
    types=('NKF1', 'XKF1'),
    select='union',
    fields=(Field('VE', 'speed_east', 'meter / second'),
            Field('VN', 'speed_north', 'meter / second'),
            Field('VD', 'speed_down', 'meter / second'),
            Field('Roll', 'roll', 'degree'),
            Field('Pitch', 'pitch', 'degree'),
            Field('Yaw', 'yaw', 'degree'),
            Field('PN', 'pos_n', 'meter'),
            Field('PE', 'pos_e', 'meter'),
            Field('PD', 'pos_d', 'meter')))

#: The autopilot's own wind estimate plus the third column of its rotation
#: matrix. Logged by firmware from 2021 onward; absent from older files.
WIND = Group(
    name='wind',
    types=('WIND',),
    select='union',
    fields=(Field('wdir', 'wdir', 'degree'),
            Field('wspeed', 'wspeed', 'meter / second'),
            Field('R13', 'R13', None),
            Field('R23', 'R23', None),
            Field('R33', 'R33', None)))

INERTIAL = Group(
    name='imu',
    types=('IMU',),
    select='union',
    fields=(Field('GyrX', 'gyr_x', None),
            Field('GyrY', 'gyr_y', None),
            Field('GyrZ', 'gyr_z', None),
            Field('AccX', 'acc_x', None),
            Field('AccY', 'acc_y', None),
            Field('AccZ', 'acc_z', None)))

#: Every group parsed from a message stream, in the order they are built.
#: PRESSURE must come after POSITION: its altitude is stored relative to the
#: first MSL fix, so the position group has to exist first.
GROUPS = (TEMPERATURE, HUMIDITY, POSITION, PRESSURE, ROTATION, WIND, INERTIAL)

#: Groups whose absence makes a file unusable.
REQUIRED_GROUPS = ('temp', 'rh', 'pos', 'pres', 'rotation')

#: PARM entries that carry identity rather than configuration.
VEHICLE_ID_PARAM = 'SYSID_THISMAV'
SENSOR_SERIAL_PARAM = 'USER_SENSORS'

#: Message types that are handled individually rather than through GROUPS.
SPECIAL_TYPES = ('PARM', 'EV', 'MSG', 'ESC')


def all_message_types():
    """ Every log message type the parser looks at.

    :rtype: list[str]
    """
    types = set(SPECIAL_TYPES)
    for group in GROUPS:
        types.update(group.types)
    return sorted(types)
