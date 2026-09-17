# Track 1 output files

The pipeline writes two tabular artifacts plus one aggregate JSON.
Coordinate columns use the VCF (variant call format)'s bare contigs (`15`); the `chr` prefix is added only by the submission CSV writer.

## annot_full.tsv - ranked candidate table
Produced by `prioritize.py` (tab-separated, one header line).

`rank  type  gene  chrom_1  pos_1  ref_1  alt_1  chrom_2  pos_2  ref_2  alt_2  epcr  consequence_1  consequence_2  cadd_1  cadd_2`

- `type`: `compound_het` (compound heterozygous: two heterozygous variants in the same gene on different chromosome copies) | `hom_alt` (homozygous alternate: both alleles carry the alternate variant)
- `epcr`: float — Estimated probability of causal relationship, (0, 1] (the column name is fixed by the challenge submission template)
- `cadd_*`: CADD (Combined Annotation Dependent Depletion) phred-scaled deleteriousness scores; empty when no CADD tile covers the site, which is treated as a neutral contribution (the column is never set to 0.0). In the shipped run the CADD tile set was not installed, so all rows are empty here.

The shipped run wrote 8,666 rows: 2,372 `compound_het` and 6,294 `hom_alt`,
from 5,012,204 input records. Per-stage counters for the run are in
`provenance.log` at the repository root.

## submission_track1.csv - challenge submission (root of the repo)
Produced by `make_submission_csv.py`, exactly the challenge template
header (checked by `checks/check_submission_schema.py`):

`proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2, epcr, finding_type, notes`

Contigs carry the `chr` prefix here; 10 data rows, sorted by epcr in descending order.
## aggregates.json - figure export (figures/)
Produced by `export_aggregates.py`; summary numbers only (histogram, medians, top‑10, consequence counts), no variant‑level rows. Read by `figures/make_figures.py`; the text and PNGs share the same source of truth.

## Instrumentation
Every pipeline run prints stage counters to stderr in the form
`[instrument] vcf_in=... qc_pass=...`; `qc_pass` counts records that passed
the genotype-quality control (FILTER is PASS or absent, genotype quality
at least the configured minimum, read depth at least the configured minimum,
and alternate-allele read fraction at least the configured minimum — see the
`qc` block in `pipeline/conf/track1_config.yml`); `coding_or_splice=...`,
`compound_het_pairs=... hom_candidates=...`; a full example is stored in
`provenance.log`.
