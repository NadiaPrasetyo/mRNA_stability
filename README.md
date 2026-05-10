# mRNA Stability Pipeline

Staged pipeline for extracting transcript regions, running structure-prediction
tools (RNAfold and friends), and preparing per-region tables for downstream
half-life analysis.

## Pipeline shape

```
01_extract.py     -> runs/<dataset>/extracted_regions/{extracted_<region>.fa, manifest.tsv}
02_stratify.sh    -> runs/<dataset>/lists/{tier_<n>/, tier_<n>.txt, .lengths.tsv}
03_calibrate.sh   -> runs/<dataset>/<tool>/calibration/<timestamp>/recommendations.tsv
04_submit.sh      -> SLURM array per tier, calibrated resources
05_collate.sh     -> runs/<dataset>/<tool>/combined.tsv
```

Three independent axes:

- **Dataset** = (genome, GFF, gene list, region selection). Defined in
  `configs/datasets/<name>.yaml`.
- **Tool** = (worker, parameters, output schema). Defined in
  `configs/tools/<name>.sh` plus `tools/<name>/`.
- **Cluster** = SLURM partition, concurrency, etc. Defined in
  `configs/cluster.sh` (machine-level).

Each pipeline step takes `--dataset` and (where applicable) `--tool`. Outputs
land under `runs/<dataset>/[<tool>/]` so multiple datasets and tools coexist
without collision.


```
01_extract.py     -> runs/<dataset>/extracted_regions/{extracted_<region>.fa,
                                                       manifest.tsv,
                                                       canonical.gff}
01b_metrics.py    -> runs/<dataset>/metrics/<metric>.tsv          (NEW)
02_stratify.sh    -> runs/<dataset>/lists/{tier_<n>/, tier_<n>.txt, .lengths.tsv}
03_calibrate.sh   -> runs/<dataset>/<tool>/calibration/<timestamp>/recommendations.tsv
04_submit.sh      -> SLURM array per tier, calibrated resources
05_collate.sh     -> runs/<dataset>/<tool>/combined.tsv
```

`canonical.gff` is now a kept artefact of `01_extract.py`: a filtered GFF
containing one transcript per gene (the same transcripts whose sequences
end up in the FASTAs). Metric plugins consume it directly.

## Metrics

Lightweight per-transcript features that don't need SLURM. Computed by
`bin/01b_metrics.py`, which dispatches to plugin modules under `metrics/`.

```bash
./bin/01b_metrics.py -d my_dataset                # run all enabled metrics
./bin/01b_metrics.py -d my_dataset -m junctions   # one metric only
./bin/01b_metrics.py -d my_dataset --force        # recompute from scratch
./bin/01b_metrics.py --list-plugins               # show discovered plugins
```

Outputs: `runs/<dataset>/metrics/<plugin>.tsv`, joinable to `manifest.tsv`
on `transcript_id` and to other metrics on `gene_id` / `transcript_id`.

### Configuring metrics

In the dataset YAML, enable plugins under a `metrics:` block:

```yaml
species: human            # selects data/references/<species>/

metrics:
  junctions:
    enabled: true
  # codon_indices:
  #   enabled: true       # (when added; needs references.csc_table etc.)

# Optional: per-plugin reference tables. Plugin docstrings list the keys
# they expect.
references:
  # csc_table: data/references/human/csc_presnyak2015.tsv
  # cai_weights: data/references/human/cai_weights.tsv
```

### Adding a new metric

Three things, mirroring the tools framework:

**1. `metrics/<name>.py`** — the plugin module:

```python
"""metrics/foo.py — short description of what's measured."""
from pathlib import Path
from typing import Iterable

OUTPUT_FILENAME = 'foo.tsv'


def get_input_paths(paths, metric_config) -> Iterable[Path]:
    """Files this metric reads (used for staleness checking)."""
    return [paths.manifest_tsv,
            paths.extract_dir / 'extracted_CDS.fa']


def compute(paths, metric_config, output_path: Path) -> None:
    """Read inputs, write output_path."""
    ...
```

**2. `configs/datasets/<dataset>.yaml`** — enable it:

```yaml
metrics:
  foo:
    enabled: true
```

**3. (Optional) `data/references/<species>/<file>`** — if the metric needs
species-specific tables, document the YAML keys it expects in its docstring,
populate `references:` in the dataset YAML, and read them in `compute`.

That's it. The orchestrator handles staleness checking, dispatch, error
reporting, and idempotency.

## `metrics/junctions.py`

Exon-junction features per transcript. Output columns:

| column | description |
|---|---|
| `transcript_id` | from `manifest.tsv` |
| `gene_id` | |
| `strand` | |
| `n_exons` | |
| `n_5UTR_junctions` | junctions falling within 5′UTR features |
| `n_CDS_junctions` | junctions falling within CDS features |
| `n_3UTR_junctions` | junctions falling within 3′UTR features |
| `stop_dist_closest_upstream` | spliced nt; `NA` if none |
| `stop_dist_closest_downstream` | spliced nt; `NA` if stop is in last exon |
| `stop_dist_last_downstream` | spliced nt to *last* downstream junction (canonical NMD metric); `NA` if stop is in last exon |
| `start_dist_closest_upstream` | spliced nt |
| `start_dist_closest_downstream` | spliced nt |

All distances are in **spliced (mature mRNA) coordinates**. NMD analyses
typically threshold `stop_dist_last_downstream` at 50 nt — kept as a raw
distance so threshold sweeps don't require recomputation.

Region junction counts assume the GFF splits UTRs and CDS per-exon (MANE,
GENCODE, modern Ensembl). Annotations that collapse a UTR into a single
feature spanning introns will report 0 UTR junctions.

## Quick start

```bash
# 1. Create configs/datasets/my_dataset.yaml from configs/datasets/example.yaml.

# 2. Extract regions
./bin/01_extract.py -d my_dataset

# 3. Stratify by length (shared across tools)
./bin/02_stratify.sh -d my_dataset

# 4. Calibrate (run on a compute node — fast now, MFE-only sampling for RNAfold)
srun --pty -p aoraki --time=00:10:00 -c 1 --mem=4G bash -c \
    './bin/03_calibrate.sh -d my_dataset -t rnafold'

# Optional: verify extrapolation accuracy on one full sample per tier
./bin/03_calibrate.sh -d my_dataset -t rnafold --verify

# 5. Submit (uses calibration recommendations)
./bin/04_submit.sh -d my_dataset -t rnafold

# 6. Resume any incomplete sequences
./bin/find_missing.sh -d my_dataset -t rnafold > /tmp/redo.txt
./bin/04_submit.sh -d my_dataset -t rnafold --list /tmp/redo.txt

# 7. Collate
./bin/05_collate.sh -d my_dataset -t rnafold
```

The `-d / -t` flags can be replaced with `DATASET=` and `TOOL=` env vars.

## Directory layout

```
configs/
  cluster.sh                     # machine-level: SLURM partition, concurrency, tier bounds
  datasets/
    example.yaml                 # dataset config template
    my_dataset.yaml              # genome + GFF + gene list + region selection
  tools/
    rnafold.sh                   # tool config: binaries, per-tier policy, calib hooks
    rnalfold.sh                  # (when added)

lib/
  paths.sh                       # shared arg parser + path resolver

bin/
  01_extract.py                  # dataset-level
  02_stratify.sh                 # dataset-level
  03_calibrate.sh                # dataset + tool
  04_submit.sh                   # dataset + tool
  05_collate.sh                  # dataset + tool
  find_missing.sh                # dataset + tool
  migrate_legacy.sh              # one-off: pre-refactor -> new layout

tools/
  rnafold/
    worker.sh                    # per-sequence work
    collate.py                   # per-tool combiner
  rnalfold/                      # (when added)
    worker.sh
    collate.py

slurm/
  array.sbatch                   # generic launcher

runs/                            # output, gitignored
  <dataset>/
    extracted_regions/           # shared across tools
    lists/                       # shared across tools
    <tool>/
      results/  errors/  tmp/  slurm_logs/  raw_shuffles/
      calibration/<timestamp>/   per-tool, no clash with other tools
      combined.tsv               from 05_collate.sh
```

## Adding a new tool

The whole point of the refactor. Three files:

**1. `configs/tools/foo.sh`** — config + hooks:

```bash
TOOL_NAME=foo
WORKER_SCRIPT="$PROJECT_ROOT/tools/foo/worker.sh"
COLLATE_SCRIPT="$PROJECT_ROOT/tools/foo/collate.py"

: "${FOO_BIN:=/path/to/foo}"
: "${FOO_WINDOW:=150}"
PATH="$(dirname "$FOO_BIN"):$PATH"
export FOO_BIN FOO_WINDOW PATH WORKER_SCRIPT COLLATE_SCRIPT TOOL_NAME

# OPTIONAL: only if this tool has a per-tier work-scaling parameter.
# tier_params <tier> -> comma-joined K=V pairs to export to sbatch
# tier_params() { ... }

# OPTIONAL: only if calibration can run a cheaper proxy (like RNAfold's
# MFE-only mode). If absent, calibration runs the worker as-is.
# calib_params() { echo "FOO_FAST=1"; }
# predict_wall_s() { local m=$1 t=$2; awk -v m=$m 'BEGIN{print int(m*N)}'; }
# predict_rss_mb() { echo "$1"; }
```

**2. `tools/foo/worker.sh`** — process one FASTA. Read `RESULTS_DIR`,
`ERRORS_DIR`, `TMP_DIR`, `TOOL_DIR` from env (set by sourcing `lib/paths.sh`).
Write `$RESULTS_DIR/results_<seqname>.csv` with a header line + data line.

**3. `tools/foo/collate.py`** — combine per-sequence CSVs into one TSV.
Read `RESULTS_DIR`, `MANIFEST_TSV`, `TOOL_DIR` from env. Write to
`$TOOL_DIR/combined.tsv`.

That's it. `02_stratify.sh` is shared, `03_calibrate.sh` and `04_submit.sh`
work for any tool, and `05_collate.sh` dispatches to your collate script.

## Calibration speedup (RNAfold)

RNAfold calibration used to run the full 1000-shuffle pipeline on each sample.
That work is dominated by 1001 RNAfold MFE folds; the shuffles themselves are
rounding error. The refactored calibration runs MFE-only (`N_SHUFFLES=0`) and
extrapolates `t_full ≈ t_mfe × (N + 1)`.

Memory is unchanged because the shuffle stream is processed one sequence at a
time, so peak RSS is set by sequence length, not by N. `predict_rss_mb` is
the identity function.

Net effect: tier 1–7 calibration runs are roughly 1000× faster. The 2× safety
factor on time and memory is unchanged.

The recommendations TSV records both measured and predicted values:

```
tier  n_samples  max_len  measured_wall_s  measured_rss_mb  predicted_wall_s  predicted_rss_mb  rec_time  rec_mem
```

so post-hoc comparison against `seff` output is straightforward.

`--verify` runs ONE full-work sample per tier and writes
`calibration/<ts>/verify.tsv` with predicted-vs-measured divergence. Useful
after binary upgrades or for new tools.

## Configuration

| File | What lives there | When you edit it |
|---|---|---|
| `configs/cluster.sh` | partition, concurrency, tier bounds, safety factors | once per machine; rarely thereafter |
| `configs/datasets/<name>.yaml` | genome, GFF, gene list, regions | once per dataset |
| `configs/tools/<name>.sh` | binary paths, tier policy, calib hooks | once per tool, plus when ViennaRNA etc. updates |

## Notes

- **Idempotency**: every step skips work that's already complete. Re-running
  any step is safe.
- **Multiple datasets in flight**: outputs are namespaced under
  `runs/<dataset>/`, so two extractions can run concurrently without collision.
- **Multiple tools per dataset**: each tool gets its own subdir under
  `runs/<dataset>/<tool>/`; tiering (`lists/`) is shared.
- **Re-tiering**: edit `TIER_BOUNDS` in `configs/cluster.sh` and re-run
  `02_stratify.sh -d <dataset>`. Recalibration is needed if breakpoints moved
  meaningfully.
- **Resume**: `find_missing.sh` produces a list of unfinished sequences;
  feed it back through `04_submit.sh --list`.
