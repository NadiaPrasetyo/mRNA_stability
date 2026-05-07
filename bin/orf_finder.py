#!/usr/bin/env python3
"""
orf_finder - identify maximal Open Reading Frames (ORFs) in DNA sequences.

Scans all six reading frames (or only the three forward frames with
--forward-only) for ORFs bounded by user-defined start and stop codons (ATG
and TAA/TAG/TGA by default).

A "maximal" ORF is the longest in-frame interval ending at a given stop
codon: when multiple in-frame start codons share the same stop, only the
most upstream start is reported.

Coordinates are reported as 1-based inclusive positions on the FORWARD
strand regardless of which strand the ORF was found on. The strand column
('+' or '-') indicates which strand the ORF lies on, and the frame column
(1, 2, or 3) is relative to that strand.

Input may be a single FASTA file (single- or multi-FASTA) or a text file
listing FASTA paths, one per line. Sequences are processed in parallel;
the unified pipeline handles inputs from a single gene up to ~15 000 genes
in one or many files.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing
import os
import sys
from dataclasses import dataclass
from functools import partial
from typing import Iterable, Iterator


# ---------------------------------------------------------------------------
# Constants and small helpers
# ---------------------------------------------------------------------------

DEFAULT_START_CODONS = "ATG"
DEFAULT_STOP_CODONS = "TAA,TAG,TGA"
DEFAULT_MIN_LENGTH = 75

VALID_BASES = frozenset("ACGT")

# str.translate is faster than per-base dict lookup for reverse-complementing.
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")

CSV_HEADER = [
    "sequence_id",
    "strand",
    "frame",
    "start",
    "end",
    "length",
    "length_fraction",
    "total_orfs_in_sequence",
    "sequence_length",
]


@dataclass(frozen=True)
class ORFRecord:
    """
    One row of output. Coordinates are 1-based inclusive forward-strand.

    Placeholder rows (emitted under --everything for sequences with no ORFs)
    use strand="" and frame=0; to_csv_row renders both as empty CSV fields,
    which most downstream tools (pandas etc.) will read as NaN.
    """
    sequence_id: str
    strand: str       # '+', '-', or '' for placeholder rows
    frame: int        # 1, 2, or 3 (relative to strand); 0 for placeholders
    start: int
    end: int
    length: int
    length_fraction: float
    total_orfs_in_sequence: int
    sequence_length: int

    @classmethod
    def no_orf_placeholder(cls, sequence_id: str, sequence_length: int) -> "ORFRecord":
        """Build the all-zero placeholder row for a sequence with no ORFs."""
        return cls(
            sequence_id=sequence_id,
            strand="",
            frame=0,
            start=0,
            end=0,
            length=0,
            length_fraction=0.0,
            total_orfs_in_sequence=0,
            sequence_length=sequence_length,
        )

    @property
    def is_placeholder(self) -> bool:
        return self.frame == 0

    def to_csv_row(self) -> list:
        return [
            self.sequence_id,
            self.strand,
            "" if self.is_placeholder else self.frame,
            self.start,
            self.end,
            self.length,
            f"{self.length_fraction:.6f}",
            self.total_orfs_in_sequence,
            self.sequence_length,
        ]


def reverse_complement(sequence: str) -> str:
    """Reverse-complement an ACGT-only DNA string."""
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


def validate_codons(codon_csv: str, param_name: str) -> frozenset:
    """Parse and validate a comma-separated list of codons (ACGT only, len 3)."""
    codons = [c.strip().upper() for c in codon_csv.split(",") if c.strip()]
    if not codons:
        raise SystemExit(f"Error: --{param_name} cannot be empty.")
    for codon in codons:
        if len(codon) != 3:
            raise SystemExit(
                f"Error: --{param_name} value '{codon}' must be exactly "
                f"3 bases (got {len(codon)})."
            )
        invalid = set(codon) - VALID_BASES
        if invalid:
            raise SystemExit(
                f"Error: --{param_name} value '{codon}' contains invalid "
                f"base(s) {sorted(invalid)}; only A, C, G, T are allowed."
            )
    return frozenset(codons)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(file_path: str) -> Iterator:
    """
    Memory-efficient generator over a FASTA file (single or multi).
    Yields (sequence_id, sequence) tuples with bases upper-cased.

    Empty records (header with no body, or a bare '>' header) are still
    yielded, with an empty sequence_id and/or empty sequence; the caller
    decides what to do with them.
    """
    sequence_id = None
    sequence_data: list = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    yield sequence_id, "".join(sequence_data)
                tokens = line[1:].split()
                # A bare '>' becomes an empty id rather than swallowing the
                # following bases under the previous id.
                sequence_id = tokens[0] if tokens else ""
                sequence_data = []
            else:
                if sequence_id is None:
                    # Bases before any header line; skip silently.
                    continue
                sequence_data.append(line.upper())

        if sequence_id is not None:
            yield sequence_id, "".join(sequence_data)


def iter_sequences(file_paths: Iterable) -> Iterator:
    """
    Yield (source_file, sequence_id, sequence) across all input files.
    Per-file IO errors are reported on stderr and skipped, so a single
    bad file doesn't halt a batch run. Empty sequences ARE yielded; the
    consumer (process_sequence / _pool_worker) decides whether to warn
    and whether to emit a placeholder row.
    """
    for path in file_paths:
        try:
            for seq_id, seq in parse_fasta(path):
                yield path, seq_id, seq
        except OSError as exc:
            print(f"Error reading {path}: {exc}", file=sys.stderr)


# ---------------------------------------------------------------------------
# ORF finding
# ---------------------------------------------------------------------------

def find_orfs_in_frame(
    sequence: str,
    frame_offset: int,
    min_length: int,
    start_codons: frozenset,
    stop_codons: frozenset,
) -> Iterator:
    """
    Single linear pass through one reading frame of `sequence`.

    Yields (orf_start, orf_end, length) for each maximal ORF, where
    [orf_start, orf_end) is 0-based half-open and length is always a
    multiple of 3 and includes the stop codon.

    Maximal-ORF semantics: between any two adjacent in-frame stop codons
    (or sequence start and the first stop), only the first start codon
    encountered triggers an ORF; subsequent in-frame starts before the
    next stop are ignored. This matches NCBI ORFfinder's behaviour.
    """
    seq_len = len(sequence)
    current_start = None
    for i in range(frame_offset, seq_len - 2, 3):
        codon = sequence[i:i + 3]
        if current_start is None:
            if codon in start_codons:
                current_start = i
        elif codon in stop_codons:
            orf_end = i + 3
            length = orf_end - current_start
            if length >= min_length:
                yield current_start, orf_end, length
            current_start = None


def rc_to_forward_coords(rc_start: int, rc_end: int, seq_length: int):
    """
    Convert an interval [rc_start, rc_end) on the reverse-complemented
    sequence to a 1-based inclusive interval on the forward strand.
    Always returns start <= end.
    """
    fwd_start_1based = seq_length - rc_end + 1
    fwd_end_1based = seq_length - rc_start
    return fwd_start_1based, fwd_end_1based


def process_sequence(
    sequence_id: str,
    sequence: str,
    *,
    forward_only: bool,
    min_length: int,
    start_codons: frozenset,
    stop_codons: frozenset,
    longest_only: bool,
    report_no_orf: bool,
) -> list:
    """Find all maximal ORFs in a single sequence and return ORFRecord rows."""
    seq_length = len(sequence)
    if seq_length == 0:
        # Empty sequence: emit a placeholder if --everything was set,
        # otherwise return nothing. The warning about malformed input is
        # printed in _pool_worker, which has the source_file context.
        if report_no_orf:
            return [ORFRecord.no_orf_placeholder(sequence_id, 0)]
        return []

    # Each entry: (strand, frame, start_1based, end_1based, length)
    raw_hits: list = []

    # Forward strand
    for frame in range(3):
        for orf_start, orf_end, length in find_orfs_in_frame(
            sequence, frame, min_length, start_codons, stop_codons
        ):
            raw_hits.append(("+", frame + 1, orf_start + 1, orf_end, length))

    # Reverse strand
    if not forward_only:
        rc = reverse_complement(sequence)
        for frame in range(3):
            for orf_start, orf_end, length in find_orfs_in_frame(
                rc, frame, min_length, start_codons, stop_codons
            ):
                fwd_start, fwd_end = rc_to_forward_coords(
                    orf_start, orf_end, seq_length
                )
                raw_hits.append(("-", frame + 1, fwd_start, fwd_end, length))

    if not raw_hits:
        if report_no_orf:
            return [ORFRecord.no_orf_placeholder(sequence_id, seq_length)]
        return []

    total_orfs = len(raw_hits)

    if longest_only:
        # Index 4 is `length`; max() is stable on ties (keeps first).
        raw_hits = [max(raw_hits, key=lambda r: r[4])]

    return [
        ORFRecord(
            sequence_id=sequence_id,
            strand=strand,
            frame=frame,
            start=start,
            end=end,
            length=length,
            length_fraction=length / seq_length,
            total_orfs_in_sequence=total_orfs,
            sequence_length=seq_length,
        )
        for (strand, frame, start, end, length) in raw_hits
    ]


def _pool_worker(item, config):
    """Adapter so multiprocessing.Pool can call process_sequence."""
    source_file, seq_id, seq = item
    if not seq:
        # Surface malformed records to the user once per occurrence.
        print(
            f"Warning: empty sequence {seq_id!r} in {source_file}",
            file=sys.stderr,
        )
    try:
        return process_sequence(seq_id, seq, **config)
    except Exception as exc:
        print(
            f"Error processing {seq_id!r} in {source_file}: {exc}",
            file=sys.stderr,
        )
        return []


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find all maximal Open Reading Frames (ORFs) in DNA "
                    "sequences.",
        epilog=(
            "Reports the longest ORF ending at each in-frame stop codon. "
            "Coordinates are 1-based inclusive forward-strand positions; "
            "the 'strand' column ('+' / '-') indicates which strand the ORF "
            "lies on, and 'frame' (1/2/3) is relative to that strand."
        ),
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-f", "--fasta-file",
        help="Path to a single input FASTA file (may be multi-FASTA).",
    )
    input_group.add_argument(
        "-l", "--file-list",
        help="Path to a text file listing FASTA file paths, one per line.",
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output CSV file.",
    )
    parser.add_argument(
        "-m", "--min-length",
        type=int,
        default=None,
        help=f"Minimum ORF length in base pairs "
             f"(default: {DEFAULT_MIN_LENGTH}; ignored with --everything).",
    )
    parser.add_argument(
        "--longest-only",
        action="store_true",
        help="Report only the single longest ORF for each sequence.",
    )
    parser.add_argument(
        "--everything",
        action="store_true",
        help="Report a row for every input sequence, including those with no "
             "ORFs (zeros in the start/end/length/length_fraction/"
             "total_orfs_in_sequence columns; empty strand and frame). "
             "Disables --min-length filtering.",
    )
    parser.add_argument(
        "--start-codons",
        type=str,
        default=DEFAULT_START_CODONS,
        help=f"Comma-separated start codons (default: {DEFAULT_START_CODONS}).",
    )
    parser.add_argument(
        "--stop-codons",
        type=str,
        default=DEFAULT_STOP_CODONS,
        help=f"Comma-separated stop codons (default: {DEFAULT_STOP_CODONS}).",
    )
    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Skip the reverse strand (3 forward frames only).",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=None,
        help="Worker processes (default: all available CPU cores).",
    )

    return parser.parse_args()


def collect_input_paths(args: argparse.Namespace) -> list:
    """Resolve --fasta-file or --file-list to a list of validated paths."""
    if args.fasta_file:
        if not os.path.isfile(args.fasta_file):
            raise SystemExit(f"Error: FASTA file not found: {args.fasta_file}")
        return [args.fasta_file]

    if not os.path.isfile(args.file_list):
        raise SystemExit(f"Error: file list not found: {args.file_list}")

    paths: list = []
    with open(args.file_list, "r") as fh:
        for raw in fh:
            path = raw.strip()
            if not path:
                continue
            if not os.path.isfile(path):
                print(
                    f"Warning: skipping non-existent file: {path}",
                    file=sys.stderr,
                )
                continue
            paths.append(path)

    if not paths:
        raise SystemExit("Error: no valid FASTA files found in --file-list.")
    return paths


def write_records(records_iter: Iterable, output_path: str) -> int:
    """Stream batches of ORFRecord to CSV. Returns the total rows written."""
    total = 0
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(CSV_HEADER)
        for batch in records_iter:
            for record in batch:
                writer.writerow(record.to_csv_row())
                total += 1
    return total


def main() -> None:
    args = parse_args()

    start_codons = validate_codons(args.start_codons, "start-codons")
    stop_codons = validate_codons(args.stop_codons, "stop-codons")

    # Resolve the effective minimum length, respecting --everything.
    if args.everything:
        if args.min_length is not None:
            print(
                "Warning: --min-length is ignored when --everything is set.",
                file=sys.stderr,
            )
        effective_min_length = 0
    else:
        effective_min_length = (
            args.min_length if args.min_length is not None else DEFAULT_MIN_LENGTH
        )
        if effective_min_length < 3:
            raise SystemExit("Error: --min-length must be >= 3.")

    if args.threads is not None and args.threads < 1:
        raise SystemExit("Error: --threads must be >= 1.")

    print(f"Start codons:    {sorted(start_codons)}")
    print(f"Stop codons:     {sorted(stop_codons)}")
    if args.everything:
        print("Min ORF length:  none (--everything)")
    else:
        print(f"Min ORF length:  {effective_min_length} bp")
    print(f"Strands:         {'forward only' if args.forward_only else 'both'}")
    if args.longest_only:
        print("Reporting only the longest ORF per sequence.")
    if args.everything:
        print("Reporting placeholder rows for sequences with no ORFs.")
    print("Maximal-ORF semantics: longest ORF per in-frame stop codon.")

    paths = collect_input_paths(args)
    print(f"Inputs:          {len(paths)} FASTA file(s).")

    config = {
        "forward_only": args.forward_only,
        "min_length": effective_min_length,
        "start_codons": start_codons,
        "stop_codons": stop_codons,
        "longest_only": args.longest_only,
        "report_no_orf": args.everything,
    }

    num_workers = args.threads if args.threads else (os.cpu_count() or 1)
    print(f"Worker processes: {num_workers}")

    sequence_iter = iter_sequences(paths)

    if num_workers == 1:
        # Skip the pool entirely on single-core / debugging runs to avoid
        # pickling overhead and to give cleaner tracebacks.
        record_iter = (_pool_worker(item, config) for item in sequence_iter)
        total = write_records(record_iter, args.output)
    else:
        worker_func = partial(_pool_worker, config=config)
        with multiprocessing.Pool(processes=num_workers) as pool:
            # imap (ordered) keeps output deterministic per input order;
            # for 1-15 000 sequences the throughput cost vs imap_unordered
            # is negligible.
            record_iter = pool.imap(worker_func, sequence_iter, chunksize=8)
            total = write_records(record_iter, args.output)

    print(f"Wrote {total} ORF record(s) to {args.output}")


if __name__ == "__main__":
    main()
