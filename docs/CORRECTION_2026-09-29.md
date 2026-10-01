# Correction: minus-strand CDS assembly in the annotator (fixed 2026-09-29)

## What was wrong
`Transcript.cds_spliced` in `pipeline/scripts/annotate.py` built the coding
sequence of minus-strand transcripts in a way that cancelled itself out: it
reversed the list of CDS segments and then reverse-complemented the joined
string. The two reversals undo each other at segment level, so each segment
came out reverse-complemented but placed in ascending genomic order instead of
translation order. Variant positions were indexed in translation order, so for
variants in multi-exon minus-strand transcripts the annotator read codons from
the wrong bases. Any coding variant in such transcripts could be affected, not
only loss-of-function ones (in our test, 5 of 6 were; one matched by
coincidence): a variant was either given a wrong consequence (for example a
stop-gain labelled `missense_variant`) or rejected as a codon mismatch and
dropped without a trace. Compound-heterozygous pairs that need a
loss-of-function variant could therefore not form in these genes. Plus-strand
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
Affected versions: from commit `HUMAN: introducing commit` up to `bbc00f2`
(the first parent of `42acf36`, i.e. the last `main` commit before the fix).
A schema checker for the submission file was added in
[PR #1](https://github.com/fulopjoz/mva-track1/pull/1).

## Influence on results
Every number below was confirmed by the authors before this note was published.

| Quantity | Before fix | After fix | Where the earlier value appeared | Source |
|---|---|---|---|---|
| Rank of the true causal pair | HUMAN | HUMAN | HUMAN | clinical re-run |
| Submitted rows changed / total | — | HUMAN | submission CSV | corrected submission |
| Minus-strand coding variants reclassified / dropped / newly paired | — | HUMAN or "not quantified: reason" | — | re-run |

The corrected submission passes `checks/check_submission_schema.py`, the
regression check passes, and the re-run is deterministic (HUMAN: confirm).

### Separate correction: phenotype-term ablation
Independently of the annotator bug, an earlier "rank without the phenotype term"
value is withdrawn: it came from a scorer stub that repeated the
with-phenotype rank instead of re-ranking. The ablation was re-measured with a
fixed scorer.

| Quantity | Withdrawn value | Re-measured | Where the withdrawn value appeared |
|---|---|---|---|
| Rank of the true pair without the phenotype term | HUMAN | HUMAN | HUMAN |

This is a single case (n = 1); these numbers are sensitivity observations, not
evidence of general robustness.

## What changed in our process
- A known-answer synthetic acceptance test for the end-to-end pipeline.
- A regression check for minus-strand codon assembly.
- A schema checker for the submission file, with tests.

## Summary (for the report and methods form)
A bug in the annotator's coding-sequence assembly for multi-exon minus-strand
transcripts could give coding variants in those genes wrong consequences or
drop them as codon mismatches, so loss-of-function variants there were missed. It
was found with a synthetic known-answer test and fixed in PR #2 (commit
42acf36), with a regression check. After re-annotation the true causal pair's
rank is HUMAN; HUMAN submitted row(s) changed and the corrected submission
passes the schema check. Separately, an earlier phenotype-ablation value from a
scorer stub is withdrawn and replaced by a re-measurement (rank HUMAN without
the phenotype term, HUMAN with it). This is a single case; the numbers are
sensitivity observations.
