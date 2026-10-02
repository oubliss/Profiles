"""
The sensor coefficient table, indexed and date-aware.

Two problems this addresses.

**Recalibration destroyed history.** MasterCoefList had no validity dates,
so a recalibrated sensor either overwrote its old row - making every file
processed before that point irreproducible - or added a second row, which
made every lookup raise "Multiple entries found". Rows may now carry
ValidFrom and ValidTo and be selected by the flight's date. Rows without
them are treated as always valid, so existing tables keep working.

**One short ID, several aircraft.** copterID.csv maps SYSID_THISMAV to a
tail number, and in the live table short ID 1 maps to four of them with
materially different wind coefficients (A 37.6 vs 32.8, B +6.8 vs -4.5).
get_tail_n returned whichever was listed first, silently. It still resolves
the same way when dates cannot separate them, but it says so.

The table is also indexed once at construction. get_coefs previously did
`self.coefs.copy().to_dict('list')` - a full DataFrame copy - on every
lookup, which is once per sensor per profile.
"""
import warnings
from collections import defaultdict

import pandas as pd

#: Optional validity columns. Absent or blank means "always valid".
VALID_FROM = 'ValidFrom'
VALID_TO = 'ValidTo'

#: Columns returned for a matched row.
COEF_FIELDS = ('A', 'B', 'C', 'D', 'Equation', 'Offset')


class AmbiguousCoefficients(LookupError):
    """More than one row matched and nothing could separate them."""


class MissingCoefficients(LookupError):
    """No row matched."""


def _parse_date(value):
    """A date cell to a Timestamp, or None when blank/'na'."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in ('na', 'nan', 'none'):
        return None
    try:
        return pd.Timestamp(text)
    except ValueError:
        return None


def _covers(row, when):
    """Is this row valid at `when`? Undated rows are valid at all times."""
    if when is None:
        return True

    when = pd.Timestamp(when)
    starts = _parse_date(row.get(VALID_FROM))
    ends = _parse_date(row.get(VALID_TO))

    if starts is not None and when < starts:
        return False
    if ends is not None and when >= ends:
        return False
    return True


def _describe(row):
    """One-line summary of a row, for error messages."""
    window = f"{row.get(VALID_FROM) or '-'} to {row.get(VALID_TO) or '-'}"
    return (f"Equation={row.get('Equation')} A={row.get('A')} "
            f"B={row.get('B')} valid {window}")


class CoefTable:
    """ A MasterCoefList, indexed by (sensor type, serial number)."""

    def __init__(self, path):
        """
        :param path: the MasterCoefList.csv to read
        """
        self.path = str(path)
        frame = pd.read_csv(self.path, dtype=str).fillna('')

        self._rows = defaultdict(list)
        for record in frame.to_dict('records'):
            key = (str(record.get('SensorType', '')).strip(),
                   str(record.get('SerialNumber', '')).strip())
            self._rows[key].append(record)

        self.has_validity = VALID_FROM in frame.columns

    def __len__(self):
        return sum(len(rows) for rows in self._rows.values())

    def lookup(self, sensor_type, serial_number, equation=None, when=None):
        """ The coefficients for one sensor.

        :param str sensor_type: 'Imet', 'RH' or 'Wind'
        :param serial_number: the sensor's serial, or an airframe tail number
        :param str equation: required only when several rows remain
        :param when: flight time, used to pick among dated rows
        :rtype: dict
        :raises MissingCoefficients: nothing matched
        :raises AmbiguousCoefficients: several matched and none was chosen
        """
        serial = _normalise_serial(serial_number)
        candidates = self._rows.get((sensor_type, serial), [])

        if not candidates:
            raise MissingCoefficients(
                f'no coefficients for sensor type {sensor_type!r} serial '
                f'{serial!r} in {self.path}')

        dated = [row for row in candidates if _covers(row, when)]
        if not dated:
            windows = '; '.join(_describe(row) for row in candidates)
            raise MissingCoefficients(
                f'{sensor_type} {serial} has {len(candidates)} row(s) but '
                f'none valid at {when}: {windows}')

        if len(dated) > 1 and equation is not None:
            dated = [row for row in dated
                     if str(row.get('Equation', '')).strip() == str(equation)]
            if not dated:
                raise MissingCoefficients(
                    f'{sensor_type} {serial} has no row with equation '
                    f'{equation!r} valid at {when}')

        if len(dated) > 1:
            detail = '; '.join(_describe(row) for row in dated)
            raise AmbiguousCoefficients(
                f'{len(dated)} rows match {sensor_type} {serial}'
                + (f' at {when}' if when is not None else ' (no date given)')
                + f'. Pass an equation, or add {VALID_FROM}/{VALID_TO} to '
                f'separate them: {detail}')

        row = dated[0]
        result = {name: row.get(name, '') for name in COEF_FIELDS}
        result[VALID_FROM] = row.get(VALID_FROM, '')
        result[VALID_TO] = row.get(VALID_TO, '')
        result['SerialNumber'] = serial
        result['SensorType'] = sensor_type
        return result


class CopterRegistry:
    """ copterID.csv: short vehicle ID to tail number."""

    def __init__(self, path):
        self.path = str(path)
        frame = pd.read_csv(self.path, names=['id', 'tail'], dtype=str,
                            header=None).fillna('')

        self._rows = defaultdict(list)
        for record in frame.to_dict('records'):
            self._rows[str(record['id']).strip()].append(
                str(record['tail']).strip())

        self._warned = set()

    def tail_number(self, copter_id, when=None):
        """ Tail number for a short vehicle ID.

        :param copter_id: SYSID_THISMAV from the log
        :param when: flight time, currently unused - copterID.csv has no
           date columns. Accepted so callers can start passing it.
        :rtype: str
        :raises MissingCoefficients: the ID is not in the registry
        """
        key = _normalise_serial(copter_id)
        tails = self._rows.get(key, [])

        if not tails:
            raise MissingCoefficients(
                f'vehicle ID {key!r} is not in {self.path}')

        if len(tails) > 1 and key not in self._warned:
            self._warned.add(key)
            warnings.warn(
                f'vehicle ID {key} maps to {len(tails)} tail numbers in '
                f'{self.path}: {", ".join(tails)}. Using {tails[0]!r}, which '
                f'is simply the first row - their coefficients may differ. '
                f'Pass tail_number= explicitly, or give copterID.csv '
                f'validity dates.', stacklevel=3)

        return tails[0]


def _normalise_serial(value):
    """Serials arrive as 6.0, '6', 6 or 'N934UA'; compare them as text."""
    if value is None:
        return ''
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).strip()
