# Correction: minus-strand CDS assembly in the annotator (fixed 2026-09-29)

## What was wrong
`Transcript.cds_spliced` in `pipeline/scripts/annotate.py` built the coding
sequence of minus-strand transcripts in a way that cancelled itself out: it
reversed the list of CDS segments and then reverse-complemented the joined
string. The two reversals undo each other at segment level, so each segment
came out reverse-complemented but placed in ascending genomic order instead of
translation order. Variant positions were indexed in translation order, so for
variants in multi-exon minus-strand transcripts the annotator read codons from
the wrong bases. Loss-of-function changes in such genes could be labelled
`missense_variant` or rejected as a codon mismatch, so compound-heterozygous
pairs that need a loss-of-function variant could not form there. Plus-strand
transcripts and minus-strand transcripts with a single CDS segment were not
affected.

## How it was found
A synthetic known-answer test (public GRCh38 chromosome 2 sequence, GENCODE v47
MANE annotation, 10 designed variants whose effects were verified independently
at codon level) expected 2 compound-heterozygous pairs; the pipeline called 0.
Rebuilding the faulty sequence model reproduced the annotator's output on all 6
minus-strand test variants exactly, which confirmed the mechanism.

## Fix and affected versions
Fixed in [PR #2](https://github.com/fulopjoz/mva-track1/pull/2) (merge commit
`42acf36`): the plus-strand sequences of the CDS segments are joined in
ascending genomic order and the joined string is reverse-complemented once
(equivalently, each segment is reverse-complemented and the segments are joined
in translation order). Regression check: `checks/chk_cds_minus_strand.py`.
Affected versions: from commit `HUMAN: introducing commit` up to the parent of
`42acf36`. A schema checker for the submission file was added in
[PR #1](https://github.com/fulopjoz/mva-track1/pull/1).

## Influence on results
Every number below was confirmed by the authors before this note was published.

| Quantity | Before fix | After fix | Source |
|---|---|---|---|
| Rank of the true causal pair | HUMAN | HUMAN | clinical re-run |
| Submitted rows changed / total | — | HUMAN | corrected submission |
| Minus-strand transcripts reclassified (aggregate) | — | HUMAN or "not quantified: reason" | re-run |
| Rank of the true pair without the phenotype term | HUMAN: withdrawn value and where it appeared | HUMAN | fixed scorer |

This is a single case (n = 1); these numbers are sensitivity observations, not
evidence of general robustness.

## What changed in our process
- A known-answer synthetic acceptance test for the end-to-end pipeline.
- A regression check for minus-strand codon assembly.
- A schema checker for the submission file, with tests.

## Summary (for the report and methods form)
A bug in the annotator's coding-sequence assembly for multi-exon minus-strand
transcripts could mislabel loss-of-function variants in those genes. It was
found with a synthetic known-answer test and fixed in PR #2 (commit 42acf36),
with a regression check. After re-annotation the true causal pair's rank is
HUMAN; HUMAN submitted row(s) changed and the submission was corrected. The
phenotype-term ablation was re-measured (rank HUMAN without, HUMAN with). This
is a single case; the numbers are sensitivity observations.
