"""metrics/nmd_fragility_nearest.py
NMD fragility under the nearest-EJC model.

A CDS position is NMD-competent iff the *next* downstream exon-exon
junction is more than 50 nt past it. Equivalent per-exon formulation:
within each non-terminal exon, positions in [exon_start, exon_end - 50)
are competent. The last 50 nt of each non-terminal exon and the entire
terminal exon are excluded.

Under this model, a close downstream junction masks the NMD signal
from more distal ones — the assumption is that a junction within 50 nt
of a stop either has its EJC stripped by the terminating ribosome or
sterically inhibits SURF-EJC assembly, and more distal junctions can't
rescue the signal.

See `metrics/_nmd_fragility_core.py` for shared implementation and
METRICS.md for full column-level spec, model-choice rationale, and
edge cases.
"""
import logging
from pathlib import Path

from metrics import _nmd_fragility_core as _core

log = logging.getLogger('metrics.nmd_fragility_nearest')

OUTPUT_FILENAME = 'nmd_fragility_nearest.tsv'


def get_input_paths(paths, metric_config):
    return _core.get_input_paths(paths, metric_config)


def compute(paths, metric_config, output_path: Path) -> None:
    _core.run(paths, metric_config, output_path, model='nearest', log=log)
