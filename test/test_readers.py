"""
Log readers.

iter_bin replaced a round trip through a vendored mavlogdump fork that wrote
a newline-delimited .json beside the original and re-read it. These tests pin
the two properties that made that swap safe: the message stream is identical,
and .json files written by the old path still load.
"""
import math

import pytest

from profiles.readers import iter_bin, iter_json, iter_messages
from profiles.readers.mavlink import WANTED_TYPES
from test.harness import staged_bin


@pytest.fixture(scope='module')
def messages(tmp_path_factory):
    return list(iter_bin(str(staged_bin(tmp_path_factory.mktemp('reader')))))


def test_reads_the_expected_message_types(messages):
    found = {m['meta']['type'] for m in messages}
    # flight616 is a 2021 log: NKF1 rather than XKF1, no ESC.
    assert {'IMET', 'RHUM', 'POS', 'BARO', 'BAR2', 'NKF1', 'PARM'} <= found
    assert found <= set(WANTED_TYPES), f'unfiltered types leaked: {found}'


def test_message_shape(messages):
    for message in messages[:200]:
        assert set(message) == {'meta', 'data'}
        assert set(message['meta']) == {'type', 'timestamp'}
        assert isinstance(message['meta']['timestamp'], float)
        assert 'mavpackettype' not in message['data']


def test_timestamps_are_non_decreasing(messages):
    stamps = [m['meta']['timestamp'] for m in messages]
    assert all(b >= a for a, b in zip(stamps, stamps[1:]))


def test_type_filter_is_applied_after_decoding(tmp_path):
    """Filtering must not be pushed into recv_match.

    DFReader implements recv_match(type=...) by skipping messages at the
    binary level, and skipped messages stop contributing to its GPS
    time-base estimation. Doing that shifted 17,073 timestamps by up to
    17 ms on this flight. Decoding everything and filtering afterwards must
    give the same timestamps whatever subset is requested.
    """
    path = str(staged_bin(tmp_path))

    everything = {m['meta']['timestamp']
                  for m in iter_bin(path, types=None)
                  if m['meta']['type'] == 'IMET'}
    just_imet = {m['meta']['timestamp']
                 for m in iter_bin(path, types=['IMET'])}

    assert just_imet == everything, (
        'IMET timestamps depend on which other types were requested')


def test_json_round_trip_matches_bin(tmp_path):
    """A .json written from this .BIN must reload identically."""
    import json

    path = str(staged_bin(tmp_path))
    from_bin = list(iter_bin(path))

    json_path = tmp_path / 'dump.json'
    with open(json_path, 'w') as handle:
        for message in from_bin:
            handle.write(json.dumps(message) + '\n')

    from_json = list(iter_json(str(json_path)))
    assert len(from_json) == len(from_bin)

    def same(left, right):
        if isinstance(left, float) and isinstance(right, float):
            return left == right or (math.isnan(left) and math.isnan(right))
        return left == right

    for a, b in zip(from_bin, from_json):
        assert a['meta'] == b['meta']
        assert set(a['data']) == set(b['data'])
        assert all(same(a['data'][k], b['data'][k]) for k in a['data'])


def test_unknown_extension_is_rejected():
    with pytest.raises(ValueError, match='not a log this reader handles'):
        list(iter_messages('/tmp/nope.txt'))


def test_file_with_no_usable_messages_reports_clearly(tmp_path):
    """A failed download saved as .BIN should say what is missing."""
    from profiles.flight import FlightLog

    decoy = tmp_path / 'notalog.BIN'
    decoy.write_text('{"status":3,"description":"Flight Id does not exist."}')

    with pytest.raises(ValueError, match='contains no usable data'):
        FlightLog(str(decoy), dev=True, nc_level=None)
