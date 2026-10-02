"""
Writing processed data out.

One place that knows about NetCDF attributes, processing levels and
provenance, rather than five hand-written save methods with four filename
conventions between them.
"""
from profiles.io.provenance import provenance_attributes
from profiles.io.naming import output_name, resolve
from profiles.io.writers import LEVELS, write_qc_variables

__all__ = ['provenance_attributes', 'write_qc_variables', 'LEVELS',
           'output_name', 'resolve']
