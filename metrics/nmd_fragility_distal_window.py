"""metrics/nmd_fragility_distal_window.py
NMD fragility under the distal-window model.

This model hypothesises that the NMD machinery is most sensitive to
premature stops located in a specific "proximal" window just prior
to the last exon-exon junction.

By default, it evaluates a 120-nt window terminating 50 nt upstream
of the last exon-exon junction. Three parameters are configurable
via the dataset YAML:

  - window_size    (default 120): length of the scanning window in nt.
  - nmd_threshold  (default 50):  size of the NMD-immune safe zone;
                                  literature gives 50-55 nt.
  - apply_nmd_rule (default True): if False, the window runs right
                                   up to (but not including) the
                                   junction, ignoring the safe zone.

If the available upstream CDS is shorter than the configured window,
it scans the available length and adjusts the density denominator
downward to reflect the true length scanned.
"""
import logging
from pathlib import Path

from metrics import _nmd_fragility_core as _core

log = logging.getLogger('metrics.nmd_fragility_distal_window')

OUTPUT_FILENAME = 'nmd_fragility_distal_window.tsv'


def get_input_paths(paths, metric_config):
    return _core.get_input_paths(paths, metric_config)


def compute(paths, metric_config, output_path: Path) -> None:
    _core.run(paths, metric_config, output_path, model='distal_window', log=log)
