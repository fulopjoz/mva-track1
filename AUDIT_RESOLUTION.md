# AUDIT_RESOLUTION.md — Track 1 (mva-track1)
## Issue-resolution ledger for the independent audit, `MVA_focused_audit_2026-10-19.zip` (received 2026-10-19)

Method: the audit was read as an issue register, not an authority. Every
claim was checked against the checked-out sources (blob SHAs matching the
audit record exactly: annotate.py `3c5eda14`, prioritize.py `a8240a21`,
TRACK2_REPORT.md `4d24b856`). The audit's synthetic_checks.py cases were
reproduced and all initially failed; after remediation each has a
disposition before/after entry below. Branch: `fix/audit-2026-10-19-t1`.

| ID | Audit finding | Severity | Evidence (audit gave) | Resolution | Verification (how to re-run) | Disposition |
|---|---|---|---|---|---|---|
| T1-01 | Canonical intronic splice donor/acceptor overwritten by `splice_region`; genuine LoF splice alleles cannot seed pairs | defect | 6 synthetic boundary cases, direct code reading | Rewrote `_classify`: removed the dead `dist<=2` double-assignment; added `_intron_class` with strand-aware donor/acceptor geometry (donor = first two intronic bases in transcription order, acceptor = last two; plus/minus strands mirrored); intronic 3–8 nt → `splice_region`, deeper → `intronic`; exonic boundaries keep their UTR/CDS consequence with the `|splice_region` tag only. | `python3 -m pytest checks/test_splice_classification.py` — 8 canonical + 6 noncanonical + exonic boundary + single-exon + end-to-end pair-propagation cases (17 tests). `python3 checks/chk_cds_minus_strand.py` (Sept regression) still OK. | FIXED (e7cac16) |
| T1-02 | Multiallelic records not handled allele-specifically (fraction always from AD[1]; every ALT processed regardless of GT) | defect | 2 synthetic VCF scenarios | `iter_qc_variants` decomposes the sample genotype: emits `(rec, zyg, alt_index)` per GT-called ALT; fraction = AD[k]/sum(AD) for that allele; malformed/missing AD/DP rejected; partial (`./1`) and fully missing (`./.`) calls skipped (documented); phased `a|b` treated as het; records with no emitted allele excluded. | `python3 -m pytest checks/test_multiallelic.py` — audit's 2 cases + 0/1, 0/2, 1/2, 2/2, 1/1, 0|1, 1|2, missing/partial/malformed matrix + 2 integration runs proving an uncalled ALT cannot seed a pair (19 tests). | FIXED (e7cac16) |
| T1-03a | Reversed duplicate compound-het rows | defect | code read | Pair generation dedupes on a canonical (position-sorted) member key; deterministic orientation in output. | `test_candidate_construction.py::test_reversed_duplicate_pairs_deduplicated` | FIXED (e7cac16) |
| T1-03b | Pairs and hom singles ranked in separate blocks | design gap | code read | Single global ordering by score descending, deterministic tie-breaks (type, member key, gene); input-order independence asserted. | `test_global_rank_not_separate_blocks`, `test_tiebreak_independent_of_input_order`, `test_run_to_run_byte_deterministic` | FIXED (e7cac16) |
| T1-03c | `internal > 0.3` hom subset counted but all singles written | defect | code read | One counter semantics: an emitted hom_singles counter that equals written rows; `>0.3` subset removed. | `test_homozygote_semantics_single_source` | FIXED (e7cac16) |
| T1-03d | Potentially silently rejected variants (non-ACGTN) | visibility | code read | `excluded_other` instrument counter printed per run. | `test_noncanonical_alleles_counted_not_dropped_silently` | FIXED (e7cac16) |
| T1-03e | Common alleles ranking high (KIF7 illustration) | question | — | Root cause was the Sept annotation bug (phantom stop), already fixed and rerun 2026-09-29; rarity is a *score* term by declared design (no AF gate at collection), documented in REPORT.md. No new code change needed. | rerun comparison diff (see T1-04 row) | EXPLAINED (previously fixed) |
| T1-03f | gnomAD/CADD/ClinVar/HPO join correctness | verification | audit asked | Characterization tests on tiny format-faithful fixtures (stream merge + exact-key AF, tabix ClinVar exact-allele, CADD phred field extraction, HPO closure terms, prior-term contract) — all pass. | `python3 -m pytest checks/test_reference_joins.py checks/test_prior_terms.py` | VERIFIED (e7cac16) |
| T1-03g | epcr framed as probability; unphased pairs as trans | wording | — | REPORT.md + TECHNICAL_DOCUMENTATION.md now call epcr a bounded heuristic rank score, explicitly NOT an empirically calibrated probability; candidate pairs labeled phase-unknown, not trans. | `grep` of REPORT.md | FIXED (REPORT.md pending-commit in branch) |
| T1-04 | Reconcile REPORT/CORRECTION/CSV with corrected outputs; rerun authorized data | pending | rank 25 vs 21 discrepancy noted by audit | Corrected clinical rerun on the approved host (dichtator) **launched** (foreground log `run_audit19.log`, output `annot_full_audit19.tsv`); on completion: diff vs 2026-09-29 table, ablation re-measurement, REPORT/CORRECTION doc numbers, regenerated submission CSV + schema check. | completion probe + `python3 checks/check_submission_schema.py` | PENDING-RERUN (numbers land next) |

## Non-goals respected
No clinical data left the authorized host; only aggregate/METRIC lines come
back for documentation. The clinical submission upload stays a human action.
