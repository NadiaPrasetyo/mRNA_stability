# Metrics

Per-transcript metric tables produced by `bin/01b_metrics.py` and the plugin
modules under `metrics/`. This file is the reference for column names,
units, edge-case handling, and join keys. For how to add a new metric, see
the README.

## Conventions across all metric tables

- **Files**: tab-separated, header on line 1, written to
  `runs/<dataset>/metrics/<plugin>.tsv`.
- **Missing or undefined values**: written as `NA`, consistent with
  `manifest.tsv` and trivially readable by R (`read.delim(..., na.strings="NA")`).
- **Join keys**: every table has `transcript_id` (preserved verbatim from
  `manifest.tsv`, joinable to it directly) and `gene_id` (normalised:
  namespace prefixes like `gene:` / `rna-` stripped, trailing version
  suffix stripped). `gene_id` is suitable for joining across metric tables
  and across datasets if they share an annotation system.
- **Coordinates**: all distances within a transcript are in **spliced
  (mature mRNA) coordinates** unless explicitly labelled as genomic.
  Intron lengths are the one routine exception — necessarily genomic.
- **Direction**: where a metric depends on a 5'→3' reading direction
  (first / last exon, codon-relative distances), it uses transcript
  reading direction, not genomic order. For `-` strand transcripts this
  means the "first exon" is the one with the highest genomic start.
- **U vs T**: T and U are treated equivalently throughout. Output column
  names use the RNA convention (`frac_U` for the combined count).
- **Format**: each plugin chooses long or wide based on its data shape;
  the per-plugin section below documents which.

## Output format choice — long vs wide

Wide is the default (one row per transcript) and is what `junctions` and
`architecture` use. `sequence_basic` is long (one row per transcript per
region) because the same set of composition columns repeats across many
regions; a wide layout would balloon to dozens of region-prefixed columns
most of which are NA for any given analysis.

When joining a long table to a wide one in R:
```r
left_join(wide_table, long_table %>% filter(region == "3UTR"),
          by = "transcript_id")
```
or use `pivot_wider(names_from = region)` to flatten if region-prefixed
columns are wanted.

---

## `junctions.tsv`

Exon-junction features per transcript. All distances in spliced
(mature mRNA) coordinates. Wide format.

| column | description |
|---|---|
| `transcript_id` | from `manifest.tsv`, preserved verbatim |
| `gene_id` | normalised |
| `strand` | `+` or `-` |
| `n_exons` | total exon count |
| `n_5UTR_junctions` | junctions within 5'UTR features = `len(five_prime_UTR features) - 1` |
| `n_CDS_junctions` | junctions within CDS features = `len(CDS features) - 1` |
| `n_3UTR_junctions` | junctions within 3'UTR features = `len(three_prime_UTR features) - 1` |
| `stop_dist_closest_upstream` | spliced nt from stop codon to nearest upstream junction; `NA` if none |
| `stop_dist_closest_downstream` | spliced nt from stop codon to nearest downstream junction; `NA` if stop is in last exon |
| `stop_dist_last_downstream` | spliced nt to the *last* downstream junction (canonical NMD metric); `NA` if stop is in last exon |
| `start_dist_closest_upstream` | spliced nt from start codon to nearest upstream junction; `NA` if none |
| `start_dist_closest_downstream` | spliced nt from start codon to nearest downstream junction; `NA` if none |

**Reference points.** First nt of start codon, last nt of stop codon —
both in transcript reading direction. Falls back to CDS extent if explicit
`start_codon` / `stop_codon` features are absent.

**NMD threshold.** Analyses typically threshold `stop_dist_last_downstream`
at 50 nt. Kept here as a raw distance so threshold sweeps don't need
recomputation.

**Junction sign convention.** A junction located exactly at the reference
position is treated as *upstream* with distance 0. A junction one nt
downstream of the reference is *downstream* with distance 1. This matters
only at the boundary case of a codon coinciding with a splice site.

**Sharp edges:**
- Region junction counts assume the GFF splits UTRs and CDS per-exon (MANE,
  GENCODE, modern Ensembl). Annotations that collapse a UTR into one
  feature spanning introns will report 0 UTR junctions even when the UTR
  spans multiple exons.
- Falling back to CDS extents for the stop codon is off by 3 nt on legacy
  RefSeq annotations that exclude the stop codon from CDS (most modern
  annotations include it). Marginal effect on closest-junction distances.

---

## `architecture.tsv`

Exon and intron length statistics per transcript. Wide format.

| column | description |
|---|---|
| `transcript_id` | from `manifest.tsv`, preserved verbatim |
| `gene_id` | normalised |
| `strand` | `+` or `-` |
| `n_exons` | total exon count |
| `n_internal_exons` | `max(0, n_exons - 2)` |
| `first_exon_length` | nt; transcript reading direction |
| `last_exon_length` | nt; transcript reading direction; equals `first_exon_length` when `n_exons == 1` |
| `internal_exon_mean` | mean length of internal exons (nt); `NA` if `n_internal_exons == 0` |
| `internal_exon_median` | median length of internal exons (nt); `NA` if `n_internal_exons == 0` |
| `internal_exon_sd` | sample SD of internal exon lengths (nt); `NA` if `n_internal_exons < 2` |
| `intron_mean` | mean intron length (genomic nt); `NA` if `n_exons < 2` |
| `intron_median` | median intron length (genomic nt); `NA` if `n_exons < 2` |
| `intron_sd` | sample SD of intron lengths (genomic nt); `NA` if `n_exons < 3` |

**SD convention.** Sample SD (Bessel-corrected, n-1 denominator) — matches
R's `sd()`. Defined only when the underlying set has at least two values.

**Intron lengths.** Computed as `next_exon.start - this_exon.end - 1` in
genomic order. The resulting set is direction-independent, so set-level
statistics (mean / median / SD) are unaffected by strand. Reported as
genomic, not spliced, nt — they are intron lengths, not transcript distances.

**Single-exon transcripts.** `n_exons = 1`; `first_exon_length` and
`last_exon_length` both equal the exon length; all stat columns are `NA`.

---

## `sequence_basic.tsv`

Length and base composition per transcript per region. Long format — one
row per `(transcript, region)` pair.

| column | description |
|---|---|
| `transcript_id` | from `manifest.tsv`, preserved verbatim |
| `gene_id` | normalised |
| `region` | region name (`mRNA`, `CDS`, `5UTR`, `3UTR`, `start_codon_region`, `stop_codon_region`, `tail_region`, …) |
| `length` | nt |
| `gc_content` | `(G+C) / length`; `NA` if `length == 0` |
| `frac_A` | `A / length`; `NA` if `length == 0` |
| `frac_C` | `C / length`; `NA` if `length == 0` |
| `frac_G` | `G / length`; `NA` if `length == 0` |
| `frac_U` | `(T+U) / length`; `NA` if `length == 0` |
| `frac_other` | non-canonical bases (N, R, Y, …) / length; `NA` if `length == 0` |
| `gc_skew` | `(G-C) / (G+C)`; `NA` if `G+C == 0` |
| `at_skew` | `(A-U) / (A+U)`; `NA` if `A+U == 0` |
| `purine_ratio` | `(A+G) / length`; `NA` if `length == 0` |
| `amino_ratio` | `(A+C) / length`; `NA` if `length == 0` |

The five `frac_*` columns sum to 1 by construction (assuming non-zero
length). `gc_content` is computed over total length, not over canonical
bases only — so a heavily-masked sequence gets a lower `gc_content` than
a clean one with the same canonical-base ratio. This makes `frac_other` a
direct read on how much sequence was ambiguous and keeps GC interpretation
transparent.

**Skew conventions.** `gc_skew` and `at_skew` use canonical-only denominators
(standard literature definition; `(G-C)/(G+C)` and `(A-U)/(A+U)`). The two
*ratios* use total length in the denominator (matching the `frac_*` columns).
They are equal to `frac_A + frac_G` and `frac_A + frac_C` respectively,
but exposed as named columns for analysis convenience. Amino bases are
A (NH2 at C6) and C (NH2 at C4); keto bases G and U/T have C=O at the
analogous position — IUPAC code `M` = aMino.

**Region discovery.** The plugin globs `extracted_*.fa` under
`extract_dir`. Region name is the filename's middle (`extracted_<region>.fa`
→ `<region>`). `UTR_pair` is hardcoded as skipped because its 3-line-per-record
ViennaRNA-constraint format is not standard FASTA and would corrupt
composition stats if naively parsed.

**Optional config:**

```yaml
metrics:
  sequence_basic:
    enabled: true
    regions: [mRNA, CDS, 5UTR, 3UTR]    # explicit region list overrides glob
    skip_regions: [tail_region]         # additional skips beyond UTR_pair
```

**Sharp edges:**
- Length-0 records (shouldn't occur in normal extraction output but
  possible if upstream emits them) get `length = 0` and `NA` everywhere
  else.
- T (DNA letter) and U (RNA letter) both count toward `frac_U`, so the
  same plugin runs over DNA-FASTA and RNA-FASTA outputs without change.
- A FASTA header whose first whitespace-token doesn't match any manifest
  transcript_id (after tolerant ID matching) is dropped with a warning.
  Verbatim match should be the norm since extraction writes manifest IDs
  as FASTA headers.

---

## Joining tables in R

Typical downstream pattern:

```r
manifest    <- read.delim("runs/<dataset>/extracted_regions/manifest.tsv", na = "NA")
junctions   <- read.delim("runs/<dataset>/metrics/junctions.tsv", na = "NA")
architect   <- read.delim("runs/<dataset>/metrics/architecture.tsv", na = "NA")
seq_basic   <- read.delim("runs/<dataset>/metrics/sequence_basic.tsv", na = "NA")

# Wide table per transcript (long → wide on the basic-composition stats)
library(dplyr); library(tidyr)
seq_wide <- seq_basic %>%
  pivot_wider(id_cols = c(transcript_id, gene_id),
              names_from = region,
              values_from = c(length, gc_content, frac_A, frac_C, frac_G,
                              frac_U, frac_other, gc_skew, at_skew,
                              purine_ratio, amino_ratio))

merged <- manifest %>%
  left_join(junctions,  by = c("transcript_id", "gene_id")) %>%
  left_join(architect,  by = c("transcript_id", "gene_id")) %>%
  left_join(seq_wide,   by = c("transcript_id", "gene_id"))
```

`n_exons` appears in both `junctions.tsv` and `architecture.tsv`. A naive
`left_join` will produce `n_exons.x` / `n_exons.y` columns. Either pick
one source explicitly (`select(-n_exons)` from one table before joining)
or treat them as a consistency check.
