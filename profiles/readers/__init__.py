"""
Readers turn a data file into a stream of log messages.

Every supported input - ArduPilot .BIN, the .json dumps the package used to
write, and anything else added later - is normalised to the same shape so
that exactly one parser consumes them:

    {"meta": {"type": str, "timestamp": float}, "data": {field: value, ...}}

``timestamp`` is UTC seconds since the Unix epoch.
"""
from profiles.readers.mavlink import iter_messages, iter_bin, iter_json

__all__ = ['iter_messages', 'iter_bin', 'iter_json']
