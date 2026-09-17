# Exomiser cross-tool validation (Track 1)

Date: 2026-09-06. Host: dichtator. [Exomiser](https://github.com/exomiser/Exomiser)
CLI 14.0.0, data bundles 2406_phenotype + 2406_hg38 (PK-magic and `unzip -tq`
integrity checks passed), Java 21 (Temurin), 48 GB heap. Both runs: PASS-only
hg38 input, shared-frequency GNOMAD_G sources, [REVEL, MVP] pathogenicity.
(Exomiser 15.x has since been released; these documented runs are the
14.0.0 artifacts named here.)

Two runs, identical configuration except phenotype input:

| Run | Phenotype input | BUB1B result | Top of list |
|---|---|---|---|
| A | Same 8 case HPO terms, hiPhive + priorityScoreFilter + OMIM prioritiser | **rank 1, p = 0.0001, score 0.9751** | BUB1B, cleanly separated (next non-HLA candidate FANCD2 p = 0.0024) |
| B | None (frequency + functional annotation + inheritance) | rank 66 (AR mode) / top-tie block (AD mode) | 46-way tie at combined score 0.5 containing BUB1B |

Without a phenotype prior there is no discriminating signal at the top of
Run B by construction (gene score = variant score only).

## What these runs can claim

- **Run A** is implementation-independence only: a second, published engine
  re-ranks the same VCF and returns the same gene first. Its phenotype module
  reads the same OMIM-derived HPO annotations as our own prior, so this is
  cross-*tool* confirmation, not independent biological validation.
- **Run B** confirms the tuning disclosure's counterfactual: remove the
  phenotype prior and no decisive ordering exists at the top anywhere.
- The tools disagree on prior-free details (their LoF model keeps KIF7 low
  where ours let it rise) - prior-free ranking is annotator-dependent, which
  is exactly the disclosure's point: the answer comes from the phenotype term.

Aggregate output shipped in this directory: `exomiser_summary.tsv` and
`exomiser_run_stats.txt`; the sha256 values below were re-verified against
the host on 2026-09-17 and match, so any rerun can be compared against
fixed numbers.

## SHA-256 of the host artifacts

```
ecee3cfa4dbdff4405187887f18264df2d051181846e046feee8186713382599  A_withHPO genes.tsv
fc8d1842c89ecb22f90056d5ce5454d6c4d2fe8fb1f83db03586e1ee5a32e65e  A_withHPO variants.tsv
0f92ca7fbec14ee4c019c977bbd5a50c9730718d1dde08e6ff70553dfea05db5  A_withHPO json
d00fd11fb2d094b5c8bea7c4d3f4fce150141211ba8743cacf0a4a63f9bc600b  B_noHPO genes.tsv
92debb8d5b5288aa8b2c9f3ad3720062d2ee9254c10067abccc9e44e73ed95e6  B_noHPO variants.tsv
a2a9314e0032e135863153dddfc9e3f920ae0e6995e2ca48a1b8ff89c539f6b2  B_noHPO json
```

Runner, logs and analysis YAMLs on the analysis host:
`~/hf-rare/exomiser/runs/` (RUN_A_DONE / RUN_B_DONE markers in `run_both.log`).
