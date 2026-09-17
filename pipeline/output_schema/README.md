# Track 1 output files

The pipeline writes two tabular artifacts plus one aggregate JSON.
Coordinate columns use the VCF's bare contigs (`15`); the `chr` prefix is
added only by the submission CSV writer.

## annot_full.tsv - ranked candidate table
Produced by `prioritize.py` (tab-separated, one header line).

`rank  type  gene  chrom_1  pos_1  ref_1  alt_1  chrom_2  pos_2  ref_2  alt_2  epcr  consequence_1  consequence_2  cadd_1  cadd_2`

- `type`: `compound_het` | `hom_alt`
- `epcr`: pair score in (0, 1] (`compound_het_score` for pairs, single-variant
  score for hom_alt rows)
- `cadd_*`: empty when no CADD tile covers the site (neutral score, never
  0.0); in the shipped run the tile set was not installed, so all rows are
  empty here

The shipped run wrote 8,666 rows: 2,372 `compound_het` and 6,294 `hom_alt`,
from 5,012,204 input records. Per-stage counters for the run are in
`provenance.log` at the repository root.

## submission_track1.csv - challenge submission (root of the repo)
Produced by `make_submission_csv.py`, exactly the challenge template
header (checked by `checks/check_submission_schema.py`):

`proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2, epcr, finding_type, notes`

Contigs carry the `chr` prefix here; 10 data rows, epcr descending.

## aggregates.json - figure export (figures/)
Produced by `export_aggregates.py`; summary numbers only (histogram,
medians, top-10, consequence counts), no variant-level rows. Consumed by
`figures/make_figures.py`, so text and PNGs share the same source of truth.

## Instrumentation
Every pipeline run prints stage counters to stderr in the form
`[instrument] vcf_in=... qc_pass=...`, `coding_or_splice=...`,
`compound_het_pairs=... hom_candidates=...`; a full example is stored in
`provenance.log`.
