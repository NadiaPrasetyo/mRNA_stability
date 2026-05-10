"""metrics/
Plugins consumed by bin/01b_metrics.py.

Plugin contract
---------------
Each module under this package must expose:

    OUTPUT_FILENAME : str
        Basename of the primary output written under <paths.metrics_dir>/.
        Convention: '<plugin_name>.tsv'.

    def get_input_paths(paths, metric_config) -> Iterable[Path]:
        Return the file paths this plugin reads. The orchestrator uses
        these for staleness checking and pre-flight existence checks.

    def compute(paths, metric_config, output_path) -> None:
        Read inputs and write the output TSV. Always runs when called;
        the orchestrator decides whether to skip via mtime comparison.

    OPTIONAL (a future plugin may add):
        REQUIRES_REFERENCES : tuple[str, ...]
            Names of reference-data files that must exist under
            paths.references_root. The orchestrator can validate these
            up-front. Not used yet by junctions but reserved for CSC,
            CAI, etc.

`paths` is a frozen lib.paths.PathContext. `metric_config` is the dict
under metrics.<plugin_name> in the dataset YAML, minus the 'enabled' flag.
"""
