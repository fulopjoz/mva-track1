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
only loss-of-function ones (in our test, 5 of 6 were affected; for the sixth
the misread codon happened to give the same change): a variant was either given a wrong consequence (for example a
stop-gain labelled `missense_variant`) or rejected as a codon mismatch and
dropped without a trace. Compound-heterozygous pairs that need a
loss-of-function variant could therefore not form in these genes. Plus-strand
transcripts and minus-strand transcripts with a single CDS segment were not
affected (with one segment the two operations reduce to a single, correct
reverse-complement).

## How it was found
A synthetic known-answer test (public GRCh38 chromosome 2 sequence, GENCODE v47
MANE annotation, 10 designed variants whose effects were verified independently
at codon level) expected 2 compound-heterozygous pairs; the pipeline called 0.
Rebuilding the faulty sequence model reproduced the annotator's output on all 6
minus-strand test variants exactly, and an independent full-CDS translation of
each test gene agreed with the designed effects. This is consistent with the
mechanism on the variants that detected it; no held-out validation set was
used. The acceptance test is kept in the authors' workspace outside this
repository and is available on request.

## Fix and affected versions
Fixed in [PR #2](https://github.com/fulopjoz/mva-track1/pull/2) (merge commit
`42acf36`): the plus-strand sequences of the CDS segments are joined in
ascending genomic order and the joined string is reverse-complemented once
(equivalently, each segment is reverse-complemented and the segments are joined
in translation order). Regression check: `checks/chk_cds_minus_strand.py`.
Affected versions: from commit `5b7eed9` (where `annotate.py` and `cds_spliced`
were added) up to `bbc00f2`
(the first parent of `42acf36`, i.e. the last `main` commit before the fix);
no commit in that range changed lines naming `cds_spliced`.
A schema checker for the submission file was added in
[PR #1](https://github.com/fulopjoz/mva-track1/pull/1).

## Influence on results
Values below come from the authors' records; items still pending are marked as such.

| Quantity | Before fix | After fix | Where the earlier value appeared | Source |
|---|---|---|---|---|
| Rank of the true causal pair | 1 | 1 | original submission (unchanged) | original run (`5b7eed9`); clinical re-run on fixed `main` (2026-09-29) |
| Submitted rows changed / total | — | pending: corrected submission not yet regenerated | submission CSV | corrected submission |
| Minus-strand coding variants reclassified / dropped / newly paired | — | not quantified: no aggregate count was computed in the re-run | — | re-run |

The regression check passes on `main`. The schema check and a determinism check
of the corrected submission are pending its regeneration. The rank of the true
pair is a single observation: because reclassified and dropped variants were not
counted, the extent of the bug's effect on this case's other candidates is
unknown.

### Separate correction: phenotype-term ablation
Independently of the annotator bug, an earlier "rank without the phenotype term"
value is withdrawn: it came from a scorer stub that repeated the
with-phenotype rank instead of re-ranking. The ablation was re-measured with a
fixed scorer. The report (last revised 2026-09-18) contains one statement about
the phenotype weight that predates this re-measurement; it is pending review
by the authors.

| Quantity | Withdrawn value | Re-measured | Where the withdrawn value appeared |
|---|---|---|---|
| Rank of the true pair without the phenotype term | 1 | 21 (with the phenotype term: 1) | internal scorer output only; not in the report, technical documentation or methods form |

This is a single case (n = 1); these numbers are sensitivity observations, not
evidence of general robustness.

## What changed in our process
- A known-answer synthetic acceptance test for the end-to-end pipeline.
- A regression check for minus-strand codon assembly.
- A schema checker for the submission file, with tests.

## Summary (for the report and methods form)
A bug in the annotator's coding-sequence assembly for multi-exon minus-strand
transcripts could give coding variants in those genes wrong consequences or
drop them as codon mismatches, so loss-of-function variants there could be missed. It
was found with a synthetic known-answer test and fixed in PR #2 (commit
42acf36), with a regression check. After re-annotation the true causal pair's
rank is unchanged at 1, although the bug's effect on other candidates was not
measured; the submission is being regenerated and checked. Separately, an earlier phenotype-ablation value from a
scorer stub is withdrawn and replaced by a re-measurement (rank 21 without
the phenotype term, 1 with it). This is a single case; the numbers are
sensitivity observations.

---

# Amendment 2026-10-19: independent audit, splice-classification and multiallelic corrections

## Additional findings and fixes
An independent focused audit (private package received 2026-10-19; findings
reproduced on synthetic cases, blob SHAs matching repository state at the
time) reported three further defects in the same annotate/prioritize core,
plus ordering questions. All are fixed and regression-tested on branch
`fix/audit-2026-10-19-t1` (51 new regression/integration tests, suite 154
green including the September minus-strand check):

- **Canonical intronic splice donor/acceptor calls were overwritten with
  `splice_region`** (audit T1-01). Genuine loss-of-function splice alleles
  could not seed compound-het pairs. Now classified strand-aware (donor =
  the first two intronic bases in transcription order; acceptor = the last
  two), on both strands; the region label applies only to non-canonical
  near-edge positions (3-8 nt).
- **Multiallelic records were evaluated from AD[1] regardless of the called
  allele, and every ALT was processed even when absent from the genotype**
  (audit T1-02). Quality control and downstream processing are now
  allele-faithful: the genotype is decomposed, per-allele fraction gates
  on the called allele's own AD entry, malformed/missing AD/DP or
  partial/fully missing calls skip the record, haploid calls count as
  hom_alt (hemizygote), and out-of-spec allele indices (AD longer than the
  ALT list) are skipped without crashing.
- **Candidate construction inconsistencies** (audit T1-03): reversed
  duplicate compound-het rows deduplicated with a deterministic,
  position-sorted orientation; the phantom `hom_candidates` >0.3 subset
  counter that agreed with nothing written was removed (the emitted
  `hom_singles` counter equals the written homozygous rows); tie-breaking
  made input-order independent in both strata (round-2 review caught a
  missing tie-break within the singles stratum, also fixed); non-base
  alleles are counted via `excluded_other` instead of vanishing silently.

## Ordering specification incident (recorded, not hidden)
The audit also questioned the ordering of the two output strata (candidate
pairs vs homozygous singles). An interim single-list global sort was
implemented and produced a valid but off-design interim table on the
authorized host (`annot_full_audit19.tsv`, 6,962 rows): under it, the
highest-scoring homozygous singles (~0.61) outranked the true pair (0.5823)
and the true pair's rank was 8. For this biallelic-recessive retrieval
task the declared ordering — candidate pairs first, then homozygous singles
— is the right retrieval structure (it is what this report had always
declared in prose); the fix was to put that specification explicitly into
the emission code and pin it with a test, not to keep the global ordering.
The definitive rerun (`annot_full_audit19b.tsv`) implements it.

## Updated全程 numbers from the definitive rerun (2026-10-19)

| Quantity | September 17 submission | September 29 (minus-strand fix) | Definitive 2026-10-19 |
|---|---|---|---|
| ranked rows | 8,666 | 6,921 | **6,962** |
| compound-het pairs | 2,372 | 2,170 | **2,211** (dedupe + legitimate LoF-splice seeding) |
| homozygous candidates (rows) | 6,294 | 4,751 | **4,751** |
| coding-or-splice candidates | 16,963 | 12,950 | **12,952** |
| QC pass (records) | 4,592,140 | 4,592,140 | **4,592,493** (+353: GT-faithful multiallelic rescue) |
| ClinVar hits | — | 4,067 | **4,067** |
| rank of the true causal pair | 1 | 1 | **1** |
| pair score (epcr) | 0.5823 | 0.5823 | **0.5823** |

Regenerated submission CSV from the definitive table: schema check
PASSED (run host-side on the authorized host). Rank of the true pair
without the phenotype term, on the definitive table: **14** (with the
phenotype term at weight 2.0 it remains 1) -- re-measured from
`desensitization_audit19b.json`, production-reproduction gate
QUALIFIED_REPRO_OK exact rows 6,962/6,962. (Stratification: 25 of 8,666 on
2026-09-17's submission table; 21 of 6,921 on the September-corrected table;
14 of 6,962 on the definitive splice-classification + multiallelic-corrected
table. The disclosing prose in REPORT.md carries the same drift.)

The phenotypic-ablation caution from the September correction stands:
single case (n = 1), sensitivity observation, not robustness evidence.
