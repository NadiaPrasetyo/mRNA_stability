"""metrics/nmd_fragility_any.py
NMD fragility under the any-EJC model.

A CDS position is NMD-competent iff the *last* downstream exon-exon
junction is more than 50 nt past it. Equivalent: any sufficiently
distant downstream junction fires NMD even when closer junctions
exist. The competent zone is a single contiguous CDS prefix ending at
(last_junction - 50).

Under this model, close downstream junctions don't mask more distal
ones — the assumption is that closer EJCs get stripped or occluded by
the terminating ribosome anyway, leaving more distal ones bound and
NMD-competent. Matches the semantics of junctions.tsv's
stop_dist_last_downstream > 50.

See `metrics/_nmd_fragility_core.py` for shared implementation and
METRICS.md for full column-level spec, model-choice rationale, and
edge cases.
"""
import logging
from pathlib import Path

from metrics import _nmd_fragility_core as _core

log = logging.getLogger('metrics.nmd_fragility_any')

OUTPUT_FILENAME = 'nmd_fragility_any.tsv'


def get_input_paths(paths, metric_config):
    return _core.get_input_paths(paths, metric_config)


def compute(paths, metric_config, output_path: Path) -> None:
    _core.run(paths, metric_config, output_path, model='any', log=log)
