# MVA Hackathon 2026 — Track 1: variant prioritization

Team `fulopjoz` (individual). Task: find the disease-causing genotype in one
whole-genome sequencing dataset of a child with a rare disease.

**Result: BUB1B p.Leu737Ter + p.Asn1002Lys, rank 1 of 6,962 candidates
(score 0.5823)** — the clinically confirmed compound-heterozygous diagnosis,
found in one pass over 5,012,204 variants (1 h 31 min, 48-core workstation,
CPU only).

## Deliverables

| File | What it is |
|---|---|
| `submission_track1.csv` | Ranked-variants submission (challenge schema, 10 rows) |
| `methods_description_FORM_FILLED.xlsx` | Methods-description form |
| `REPORT.md` | Methods and results report (rendered: `TRACK1_REPORT.pdf`) |

## Reproduce

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install pysam cyvcf2 pyyaml
```

References (paths in `pipeline/conf/track1_config.yml`): GENCODE v47 MANE,
GRCh38 FASTA, gnomAD v4.1 exome sites, ClinVar, HPO (hp.obo +
gene-to-phenotype map). Then:

```bash
python3 pipeline/scripts/prioritize.py \
    --conf pipeline/conf/track1_config.yml --with-af --with-hpo \
    --out output/annot_full.tsv          # VCF -> ranked table (6,962 rows)

python3 pipeline/scripts/make_submission_csv.py \
    --in output/annot_full.tsv --out submission_track1.csv
```

This regenerates all typed columns (positions, refs, alts, scores, finding_type)
verbatim. The `notes` column of the shipped `submission_track1.csv` was phrased
by hand for the judges; the generator emits a compact `gene=...;consequence/...`
string in its place.

Quick checks without a full run:

```bash
python3 pipeline/scripts/annotate.py    # consequence-calling self-test (SKIP if refs absent)
python3 pipeline/scripts/_hpo.py        # HPO scoring self-test (SKIP if refs absent)
head -3 submission_track1.csv           # rank 1: chr15 40,209,701 + chr15 40,220,612
```
