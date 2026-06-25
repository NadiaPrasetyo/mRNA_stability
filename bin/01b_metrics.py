#!/usr/bin/env python3
"""
01b_metrics.py
Run lightweight per-transcript metrics defined as plugins under metrics/.

Each metric is a Python module under metrics/<name>.py exposing the plugin
contract documented in metrics/__init__.py.

─────────────────────────────────────────────────────────────────────────────
MODES
─────────────────────────────────────────────────────────────────────────────
Dataset mode  (original behaviour, unchanged):
  ./bin/01b_metrics.py --dataset human_test
  ./bin/01b_metrics.py -d human_test --metric junctions
  ./bin/01b_metrics.py -d human_test --force
  ./bin/01b_metrics.py -d human_test --list-plugins

Split-FASTA mode:
  ./bin/01b_metrics.py \\
      --fasta-5utr  sequences_5utr.fa \\
      --fasta-cds   baseline_cds.fa   \\
      --fasta-3utr  sequences_3utr.fa \\
      --species "Homo sapiens"

  Rules:
  • --fasta-cds must contain exactly ONE sequence.  That single CDS is
    used as the shared baseline for every sample.  Junction and NMD
    metrics therefore depend only on CDS structure and are identical
    across all samples.
  • --fasta-5utr / --fasta-3utr must have matching accession IDs (the
    first whitespace-delimited token after '>').  Use --allow-missing-utrs
    to treat absent partners as zero-length rather than erroring.
  • No GFF download is required.  CDS/UTR coordinates are derived from
    the assembled sequence lengths and written as a canonical.gff so
    GFF-dependent plugins (junctions, architecture, NMD, uORF) work
    without modification.

Shared flags:
  --metric NAME    run only this plugin (repeatable)
  --force          ignore staleness, recompute everything
  --list-plugins   enumerate available metrics/ modules
─────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import csv
import importlib
import logging
import os
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

# ── Path bootstrap ─────────────────────────────────────────────────────────
_THIS         = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from lib.paths import PathContext, resolve_paths, add_project_root_to_syspath  # noqa: E402

# Ensure lib.* is importable in every mode from the start.
add_project_root_to_syspath(_PROJECT_ROOT)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('metrics')

# Plugins that read GFF / extract_dir annotation.  In split-FASTA mode
# these still run: annotation is synthesised from sequence lengths so no
# external GFF download is needed.
_GFF_DEPENDENT_PLUGINS: frozenset[str] = frozenset({
    'architecture',
    'junctions',
    'nmd_fragility_core',
    'nmd_fragility_full',
    'nmd_fragility_window',
    'uorf',
})


# ══════════════════════════════════════════════════════════════════════════════
# Shared plugin helpers  (unchanged from original)
# ══════════════════════════════════════════════════════════════════════════════

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        log.error("PyYAML not installed (pip install pyyaml).")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f) or {}


def _load_plugin(name: str):
    """Dynamically import metrics.<name>.  Raises ImportError if missing."""
    return importlib.import_module(f'metrics.{name}')


def _is_output_current(output_path: Path, input_paths: Iterable[Path]) -> bool:
    """True if output exists and is newer than every existing input."""
    if not output_path.exists():
        return False
    out_mtime = output_path.stat().st_mtime
    for ip in input_paths:
        ip = Path(ip)
        if ip.exists() and ip.stat().st_mtime > out_mtime:
            return False
    return True


def _enabled_metrics(config: dict) -> dict:
    """Return {name: metric_config} for metrics flagged enabled: true."""
    metrics_block = config.get('metrics') or {}
    return {
        name: (mcfg or {})
        for name, mcfg in metrics_block.items()
        if (mcfg or {}).get('enabled', False)
    }


def _list_available_plugins() -> list[str]:
    """Discover plugin modules by listing metrics/*.py."""
    metrics_dir = _PROJECT_ROOT / 'metrics'
    if not metrics_dir.is_dir():
        return []
    return sorted(
        p.stem for p in metrics_dir.glob('*.py')
        if p.stem != '__init__' and not p.stem.startswith('_')
    )


def _run_one_plugin(name: str, mcfg: dict, paths, force: bool) -> bool:
    """
    Load and run one plugin.
    Returns True on success (including 'skipped as current'), False on failure.
    *paths* is a real PathContext (dataset mode) or the split-FASTA shim.
    """
    plugin_log = logging.getLogger(f'metrics.{name}')

    try:
        plugin = _load_plugin(name)
    except ImportError as e:
        plugin_log.error(f"Plugin not found: {e}")
        return False

    required = ('OUTPUT_FILENAME', 'get_input_paths', 'compute')
    missing  = [r for r in required if not hasattr(plugin, r)]
    if missing:
        plugin_log.error(f"Plugin missing attributes: {missing}")
        return False

    output_path = paths.metrics_dir / plugin.OUTPUT_FILENAME
    inputs      = list(plugin.get_input_paths(paths, mcfg))

    if not force and _is_output_current(output_path, inputs):
        plugin_log.info(f"Output current ({output_path.name}) — skipping.")
        return True

    missing_inputs = [str(p) for p in inputs if not Path(p).exists()]
    if missing_inputs:
        plugin_log.error(
            "Required input(s) missing:\n  " + "\n  ".join(missing_inputs))
        return False

    paths.metrics_dir.mkdir(parents=True, exist_ok=True)
    plugin_log.info(f"Computing → {output_path}")
    try:
        plugin.compute(paths, mcfg, output_path)
    except Exception as e:
        plugin_log.exception(f"Plugin raised: {e}")
        return False

    if not output_path.exists():
        plugin_log.error(f"Plugin completed but did not produce {output_path}")
        return False

    plugin_log.info(f"Done ({output_path.stat().st_size} bytes)")
    return True


def _run_plugins(plugins: dict, paths, force: bool) -> None:
    """Run every plugin in *plugins* {name: mcfg}; exit(1) on any failure."""
    log.info(f"Running {len(plugins)} metric(s): {', '.join(plugins)}")
    paths.metrics_dir.mkdir(parents=True, exist_ok=True)
    n_ok = sum(_run_one_plugin(n, m, paths, force) for n, m in plugins.items())
    log.info(f"Metrics complete: {n_ok}/{len(plugins)} succeeded.")
    if n_ok < len(plugins):
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# Dataset mode  (original, unchanged)
# ══════════════════════════════════════════════════════════════════════════════

def _run_dataset_mode(args: argparse.Namespace) -> None:
    try:
        paths = resolve_paths(args.dataset, _PROJECT_ROOT)
    except FileNotFoundError as e:
        log.error(str(e))
        sys.exit(1)

    config  = _load_yaml(paths.dataset_yaml)
    species = config.get('species')
    if species:
        paths = resolve_paths(args.dataset, _PROJECT_ROOT, species=species)

    log.info(f"Dataset: {args.dataset}"
             + (f" (species: {species})" if species else ""))

    enabled = _enabled_metrics(config)
    if args.metric:
        requested = {}
        for m in args.metric:
            if m in enabled:
                requested[m] = enabled[m]
            elif m in (config.get('metrics') or {}):
                log.warning(f"Metric '{m}' is in YAML but not enabled — running anyway.")
                requested[m] = config['metrics'][m] or {}
            else:
                log.warning(f"Metric '{m}' has no config block — running with empty config.")
                requested[m] = {}
        enabled = requested

    if not enabled:
        log.warning("No metrics enabled.  Add a 'metrics:' block to your dataset YAML.")
        return

    _run_plugins(enabled, paths, args.force)


# ══════════════════════════════════════════════════════════════════════════════
# Split-FASTA helpers
# ══════════════════════════════════════════════════════════════════════════════

# ── FASTA I/O ─────────────────────────────────────────────────────────────────

def _read_raw_fasta(path: Path) -> dict[str, str]:
    """Return {seq_id: uppercase_sequence} for every record in *path*."""
    seqs: dict[str, str] = {}
    cur_id: Optional[str] = None
    parts:  list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(parts).upper().replace('U', 'T')
                cur_id = line[1:].split()[0]
                parts  = []
            elif line:
                parts.append(line)

    if cur_id is not None:
        seqs[cur_id] = ''.join(parts).upper().replace('U', 'T')

    return seqs


def _read_single_sequence(path: Path) -> Tuple[str, str]:
    """
    Read a FASTA that must contain exactly one record.
    Returns (seq_id, sequence).  Raises ValueError otherwise.
    """
    seqs = _read_raw_fasta(path)
    if not seqs:
        raise ValueError(f"--fasta-cds is empty: {path}")
    if len(seqs) > 1:
        ids = ', '.join(list(seqs)[:5]) + ('…' if len(seqs) > 5 else '')
        raise ValueError(
            f"--fasta-cds must contain exactly one sequence, "
            f"found {len(seqs)}: {path}\n  IDs: {ids}"
        )
    return next(iter(seqs.items()))


def _write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    """Write (seq_id, sequence) pairs to *path*, 60-char wrapped."""
    with open(path, 'w') as fh:
        for seq_id, seq in records:
            fh.write(f'>{seq_id}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')


# ── Assembly ──────────────────────────────────────────────────────────────────

def _assemble_transcripts(
    path_5utr: Path,
    path_cds:  Path,
    path_3utr: Path,
    allow_missing: bool = False,
) -> List[Tuple[str, str, int, int, int]]:
    """
    Concatenate per-sample 5UTR + shared CDS + per-sample 3UTR.

    Returns list of (sample_id, full_seq, utr5_len, cds_len, utr3_len).
    Accession IDs in the 5UTR and 3UTR FASTAs must match unless
    *allow_missing* is True, in which case absent UTRs are zero-length.
    """
    cds_id, cds_seq = _read_single_sequence(path_cds)
    cds_len = len(cds_seq)
    log.info(f"Shared CDS: '{cds_id}'  {cds_len} nt")
    if cds_len % 3 != 0:
        log.warning(f"CDS length {cds_len} nt is not divisible by 3.")

    seqs_5 = _read_raw_fasta(path_5utr)
    seqs_3 = _read_raw_fasta(path_3utr)
    ids_5, ids_3 = set(seqs_5), set(seqs_3)

    if ids_5 != ids_3:
        only_5 = sorted(ids_5 - ids_3)
        only_3 = sorted(ids_3 - ids_5)
        parts  = []
        if only_5: parts.append(f"only in 5UTR: {only_5[:5]}")
        if only_3: parts.append(f"only in 3UTR: {only_3[:5]}")
        msg = "; ".join(parts)
        if not allow_missing:
            raise ValueError(
                f"5UTR / 3UTR accession mismatch — {msg}. "
                "Use --allow-missing-utrs to treat absent regions as zero-length."
            )
        log.warning(f"5UTR / 3UTR mismatch ({msg}); absent UTRs → zero-length.")

    records: List[Tuple[str, str, int, int, int]] = []
    for sample_id in sorted(ids_5 | ids_3):
        u5 = seqs_5.get(sample_id, '')
        u3 = seqs_3.get(sample_id, '')
        records.append((sample_id, u5 + cds_seq + u3, len(u5), cds_len, len(u3)))

    log.info(
        f"Assembled {len(records)} transcript(s)  "
        f"({len(seqs_5)} × 5UTR  +  1 shared CDS  +  {len(seqs_3)} × 3UTR)."
    )
    return records


# ── extract_dir builder ───────────────────────────────────────────────────────
#
# Reproduces the file tree that metric plugins expect, without any
# external library (no lib.fasta_to_extract, no lib.fasta_header_parser).
# All annotation is synthesised from the assembled sequence lengths.

_MANIFEST_COLS = [
    'transcript_id', 'gene_id', 'region', 'length',
    'source_file', 'seqid', 'strand',
]


def _build_extract_dir(
    records: List[Tuple[str, str, int, int, int]],
    source_fasta: Path,
    extract_dir: Path,
    force: bool = False,
) -> None:
    """
    Write the extract_dir file tree consumed by metric plugins:
      manifest.tsv          one mRNA + optional region rows per sample
      canonical.gff         GFF3 with gene/mRNA/exon/CDS/UTR features
      extracted_mRNA.fa     full assembled transcripts
      extracted_CDS.fa      CDS slices (shared sequence, one record per sample)
      extracted_5UTR.fa     5' UTR slices  (omitted if all zero-length)
      extracted_3UTR.fa     3' UTR slices  (omitted if all zero-length)

    Skips rebuild if manifest.tsv is newer than source_fasta and not *force*.
    """
    extract_dir.mkdir(parents=True, exist_ok=True)
    sentinel = extract_dir / 'manifest.tsv'

    if (not force
            and sentinel.exists()
            and sentinel.stat().st_mtime > source_fasta.stat().st_mtime):
        log.info(f"extract_dir is current ({extract_dir}) — skipping rebuild.")
        return

    log.info(f"Building extract_dir → {extract_dir}")

    manifest_rows: list[dict] = []
    gff_lines:     list[str]  = ['##gff-version 3']
    mrna_recs:     list[Tuple[str, str]] = []
    utr5_recs:     list[Tuple[str, str]] = []
    cds_recs:      list[Tuple[str, str]] = []
    utr3_recs:     list[Tuple[str, str]] = []

    src = str(source_fasta)

    for sample_id, full_seq, u5_len, cds_len, u3_len in records:
        seq_len  = len(full_seq)
        gene_id  = sample_id   # use sample ID as gene ID (no separate gene names)
        strand   = '+'

        # 0-based half-open coordinates
        u5_s,  u5_e  = 0,              u5_len
        cds_s, cds_e = u5_len,         u5_len + cds_len
        u3_s,  u3_e  = u5_len + cds_len, seq_len

        # ── manifest rows ─────────────────────────────────────────────────
        manifest_rows.append({
            'transcript_id': sample_id, 'gene_id': gene_id,
            'region': 'mRNA', 'length': seq_len,
            'source_file': src, 'seqid': sample_id, 'strand': strand,
        })
        mrna_recs.append((sample_id, full_seq))

        if u5_len > 0:
            manifest_rows.append({
                'transcript_id': sample_id, 'gene_id': gene_id,
                'region': '5UTR', 'length': u5_len,
                'source_file': src, 'seqid': sample_id, 'strand': strand,
            })
            utr5_recs.append((sample_id, full_seq[u5_s:u5_e]))

        manifest_rows.append({
            'transcript_id': sample_id, 'gene_id': gene_id,
            'region': 'CDS', 'length': cds_len,
            'source_file': src, 'seqid': sample_id, 'strand': strand,
        })
        cds_recs.append((sample_id, full_seq[cds_s:cds_e]))

        if u3_len > 0:
            manifest_rows.append({
                'transcript_id': sample_id, 'gene_id': gene_id,
                'region': '3UTR', 'length': u3_len,
                'source_file': src, 'seqid': sample_id, 'strand': strand,
            })
            utr3_recs.append((sample_id, full_seq[u3_s:u3_e]))

        # ── GFF3 features ─────────────────────────────────────────────────
        # seqid = sample_id so each transcript is its own coordinate system.
        # All positions are 1-based inclusive (GFF3 convention).
        def _gff(feat, s0, e0, attrs):
            # s0/e0 are 0-based half-open → GFF3 1-based inclusive
            return (f"{sample_id}\tsplit_fasta\t{feat}\t{s0+1}\t{e0}"
                    f"\t.\t{strand}\t.\t{attrs}")

        gff_lines.append(_gff('gene', 0, seq_len,
                              f'ID={gene_id};gene_id={gene_id}'))
        gff_lines.append(_gff('mRNA', 0, seq_len,
                              f'ID={sample_id};Parent={gene_id};'
                              f'transcript_id={sample_id};gene_id={gene_id}'))

        # Single exon spanning the whole transcript (no introns in
        # a spliced mRNA built from UTR + CDS parts).
        gff_lines.append(_gff('exon', 0, seq_len,
                              f'Parent={sample_id};transcript_id={sample_id}'))

        # CDS feature
        gff_lines.append(_gff('CDS', cds_s, cds_e,
                              f'Parent={sample_id};transcript_id={sample_id}'))

        # UTR features (only when non-zero)
        if u5_len > 0:
            gff_lines.append(_gff('five_prime_UTR', u5_s, u5_e,
                                  f'Parent={sample_id};transcript_id={sample_id}'))
        if u3_len > 0:
            gff_lines.append(_gff('three_prime_UTR', u3_s, u3_e,
                                  f'Parent={sample_id};transcript_id={sample_id}'))

        # start_codon / stop_codon (first / last 3 nt of CDS)
        gff_lines.append(_gff('start_codon', cds_s, cds_s + 3,
                              f'Parent={sample_id};transcript_id={sample_id}'))
        gff_lines.append(_gff('stop_codon',  cds_e - 3, cds_e,
                              f'Parent={sample_id};transcript_id={sample_id}'))

    # ── Write files ───────────────────────────────────────────────────────
    with open(extract_dir / 'manifest.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_MANIFEST_COLS, delimiter='\t',
                           extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(manifest_rows)

    with open(extract_dir / 'canonical.gff', 'w') as fh:
        fh.write('\n'.join(gff_lines) + '\n')

    _write_fasta(extract_dir / 'extracted_mRNA.fa', mrna_recs)
    _write_fasta(extract_dir / 'extracted_CDS.fa',  cds_recs)
    if utr5_recs:
        _write_fasta(extract_dir / 'extracted_5UTR.fa', utr5_recs)
    if utr3_recs:
        _write_fasta(extract_dir / 'extracted_3UTR.fa', utr3_recs)

    log.info(
        f"extract_dir built: {len(records)} transcripts, "
        f"{len(utr5_recs)} with 5UTR, {len(utr3_recs)} with 3UTR."
    )


# ── PathContext shim ──────────────────────────────────────────────────────────

def _build_split_paths(output_dir: Path, species: str) -> object:
    """
    Minimal PathContext shim for split-FASTA runs.

    Exposes every attribute real metric plugins read from a PathContext so
    plugins need no modification.  Attribute values mirror what the real
    PathContext would contain for an equivalent dataset-mode run.
    """
    extract_dir = output_dir / 'extracted_regions'

    class _SplitPaths:
        def __init__(self):
            # ── Paths plugins read ────────────────────────────────────────
            self.extract_dir    = extract_dir
            self.canonical_gff  = extract_dir / 'canonical.gff'
            self.manifest_tsv   = extract_dir / 'manifest.tsv'
            # sequence-only plugins read extracted_mRNA.fa via extract_dir
            self.fasta_path     = extract_dir / 'extracted_mRNA.fa'
            # ── Orchestrator ──────────────────────────────────────────────
            self.metrics_dir    = output_dir / 'metrics'
            self.output_dir     = output_dir
            self.run_dir        = output_dir
            self.species        = species
            self.engineered     = True   # no external GFF in this mode
            self.dataset_yaml   = None
            self.gff_path       = None
            self.references_root = None
            # ── Mirror real PathContext fields ────────────────────────────
            # Prevents AttributeError on any field a plugin happens to read.
            try:
                import dataclasses as _dc
                for f in _dc.fields(PathContext):
                    if not hasattr(self, f.name):
                        setattr(self, f.name, None)
            except TypeError:
                pass

    return _SplitPaths()


# ══════════════════════════════════════════════════════════════════════════════
# Split-FASTA mode orchestration
# ══════════════════════════════════════════════════════════════════════════════

def _run_split_fasta_mode(args: argparse.Namespace) -> None:
    """
    Entry point for --fasta-5utr / --fasta-cds / --fasta-3utr runs.

    1. Validate the three input paths.
    2. Assemble per-sample full transcripts (5UTR + shared CDS + 3UTR).
    3. Write assembled_transcripts.fa (for audit / downstream use).
    4. Build extract_dir with manifest.tsv, canonical.gff, and split FASTAs
       so all metric plugins can run without modification.
    5. Execute all requested plugins.
    """
    path_5utr = Path(args.fasta_5utr).expanduser().resolve()
    path_cds  = Path(args.fasta_cds).expanduser().resolve()
    path_3utr = Path(args.fasta_3utr).expanduser().resolve()

    for label, p in [('5UTR', path_5utr), ('CDS', path_cds), ('3UTR', path_3utr)]:
        if not p.exists():
            log.error(f"{label} FASTA not found: {p}")
            sys.exit(1)

    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else path_cds.parent / f"metrics_{path_cds.stem}_split"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Split-FASTA mode")
    log.info(f"  5UTR : {path_5utr}")
    log.info(f"  CDS  : {path_cds}  (shared across all samples — 1 sequence)")
    log.info(f"  3UTR : {path_3utr}")
    log.info(f"  Out  : {output_dir}")
    log.info("  Junction / NMD metrics are derived from the shared CDS only "
             "and will be identical for every sample.")

    # 1. Assemble
    try:
        records = _assemble_transcripts(
            path_5utr, path_cds, path_3utr,
            allow_missing=getattr(args, 'allow_missing_utrs', False),
        )
    except ValueError as exc:
        log.error(str(exc))
        sys.exit(1)

    # 2. Write assembled FASTA (audit copy)
    assembled_fa = output_dir / 'assembled_transcripts.fa'
    _write_fasta(assembled_fa, [(r[0], r[1]) for r in records])
    log.info(f"Assembled FASTA → {assembled_fa}")

    # 3. Build extract_dir (manifest + canonical.gff + split FASTAs)
    paths = _build_split_paths(output_dir, args.species)
    _build_extract_dir(records, assembled_fa, paths.extract_dir, force=args.force)

    # 4. Determine plugins to run
    if args.metric:
        plugins = {m: {} for m in args.metric}
    else:
        plugins = {n: {} for n in _list_available_plugins()}

    if not plugins:
        log.warning("No plugins found under metrics/.")
        return

    _run_plugins(plugins, paths, args.force)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Mode selection ────────────────────────────────────────────────────
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        '--dataset', '-d',
        default=os.environ.get('DATASET'),
        metavar='NAME',
        help="Dataset mode: name of configs/datasets/<name>.yaml",
    )
    mode.add_argument(
        '--fasta-5utr',
        metavar='FILE',
        dest='fasta_5utr',
        help="Split-FASTA mode: per-sample 5' UTR sequences.",
    )

    # ── Split-FASTA companions ────────────────────────────────────────────
    parser.add_argument(
        '--fasta-cds',
        metavar='FILE',
        dest='fasta_cds',
        help=(
            "Single shared CDS sequence (required with --fasta-5utr). "
            "Junction / NMD metrics are derived from this sequence and are "
            "identical across all samples."
        ),
    )
    parser.add_argument(
        '--fasta-3utr',
        metavar='FILE',
        dest='fasta_3utr',
        help=(
            "Per-sample 3' UTR sequences (required with --fasta-5utr). "
            "Accession IDs must match --fasta-5utr."
        ),
    )
    parser.add_argument(
        '--allow-missing-utrs',
        action='store_true',
        default=False,
        dest='allow_missing_utrs',
        help=(
            "Treat samples whose 5UTR or 3UTR partner is absent as having a "
            "zero-length UTR instead of raising an error."
        ),
    )
    parser.add_argument(
        '--species', '-s',
        metavar='BINOMIAL',
        help="Species name stored in run metadata (e.g. 'Homo sapiens'). "
             "Required in split-FASTA mode.",
    )

    # ── Shared / output options ───────────────────────────────────────────
    parser.add_argument(
        '--output-dir', '-o',
        metavar='DIR',
        help=(
            "Output directory (split-FASTA mode only). "
            "Defaults to <cds_dir>/metrics_<cds_stem>_split/"
        ),
    )
    parser.add_argument(
        '--metric', '-m',
        action='append',
        default=None,
        metavar='NAME',
        help="Run only this metric (repeatable). "
             "Default: all enabled (dataset) / all discovered (split-FASTA).",
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help="Recompute even if outputs are current.",
    )
    parser.add_argument(
        '--list-plugins',
        action='store_true',
        help="List discovered plugin modules and exit.",
    )
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── --list-plugins ────────────────────────────────────────────────────
    if args.list_plugins:
        plugins = _list_available_plugins()
        if not plugins:
            print("(no plugins found under metrics/)", file=sys.stderr)
            return
        print("Available plugins:")
        print(f"  {'Plugin':<32}  GFF-dependent?")
        print(f"  {'-'*32}  {'-'*42}")
        for p in plugins:
            note = ('yes  (annotation synthesised from lengths in split-FASTA mode)'
                    if p in _GFF_DEPENDENT_PLUGINS else 'no')
            print(f"  {p:<32}  {note}")
        return

    # ── Route ────────────────────────────────────────────────────────────
    if args.fasta_5utr:
        missing = []
        if not args.fasta_cds:  missing.append('--fasta-cds')
        if not args.fasta_3utr: missing.append('--fasta-3utr')
        if not args.species:    missing.append('--species')
        if missing:
            parser.error(f"Split-FASTA mode also requires: {', '.join(missing)}")
        _run_split_fasta_mode(args)

    elif args.dataset:
        _run_dataset_mode(args)

    else:
        parser.error(
            "Provide one of:\n"
            "  --dataset NAME                                        (dataset mode)\n"
            "  --fasta-5utr F --fasta-cds F --fasta-3utr F          (split-FASTA mode)"
        )


if __name__ == '__main__':
    main()