# Track 1 pipeline report

Team fulopjoz. MVA Hackathon 2026, Track 1, variant prioritization.
Submission file: `submission_track1.csv`.

## Result

The pipeline ranks the clinically confirmed genotype first of 6,962
candidates: compound-heterozygous **BUB1B p.Leu737Ter**
(chr15:40,209,701 T>G, stop gained; ClinVar allele ID 529272) plus
**BUB1B p.Asn1002Lys** (chr15:40,220,612 T>G, missense), score 0.5823.
Score at rank 2 is 0.5243. The score column `epcr` — the challenge-defined
(0, 1] score slot ("estimated probability of causal relationship") — is
populated by the pipeline's composite weighting. It is a bounded heuristic
rank score, **not an empirically calibrated probability** (no calibration
validation exists for this score calibration label; see Section *Method*).
Pairs listed are phase-unknown *candidates* — presence of two compound
alleles in the same gene does not by itself prove that they are in trans.

| Rank | Genotype | Score |
|---|---|---|
| 1 | BUB1B p.Leu737Ter + p.Asn1002Lys | 0.5823 |
| 2 | RAI1 frameshift/frameshift pair | 0.5243 |
| 3-6 | RAI1 missense/frameshift pairs | 0.4914 |
| 7-10 | HLA-DRB1 founder pairs | 0.4509 |
| 2,212 and below | homozygous-recessive singles (ordered after biallelic pairs by design) | — |

**About the ordering, explicitly.** The two strata are ranked by the same
score, but a candidate *pair* aggregates the evidence of two alleles and is
the target entity of a biallelic-recessive retrieval task, so the output
lists every pair before any homozygous single. That is the declared design,
and it matters here: the highest-scoring singles in this run reach 0.6119,
above the pair's 0.5823. Under a single global sort by raw score the true
pair would rank **8**, not 1 — the interim global-order run measured
exactly that and is recorded in the correction ledger; the ranking above is
the pair-first declaration that the method text has always carried.

Note: the September 17 table showed three KIF7 rows at ranks 8-10; those were
produced by the minus-strand annotation error corrected on 2026-09-29
(the "stop" became p.S1016S synonymous), and the September 29 corrected table
already removed them. The table above is the definitive 2026-10-19 run
(splice-classification and multiallelic-allele handling audits included).

## Method

One pass over the 5,012,204-record proband VCF. Each stage's counter is
measured from the run.

1. **Quality filter.** FILTER is PASS or absent, genotype quality at least 20,
   read depth at least 10, alternate-allele fraction of the called allele at least
   0.25. Multiallelic records are evaluated per genotype-called allele (the
   September run indexed AD[1] regardless of GT). Keeps 4,592,493 records
   (September's gate kept 4,592,140; the 353-record difference is multiallelic
   records whose ALT1 failed the fraction gate while a called ALT2 passed).
2. **Effect calling.** Variants are annotated against local GENCODE v47
   MANE transcripts. Stop gained, frameshift, and canonical splice
   donor/acceptor variants count as loss of function (LoF); missense and
   in-frame indels are scored separately; canonical intronic splice sites are
   classified strand-aware (the earlier build overwrote them with
   `splice_region`, which could not seed pairs). Keeps 12,952 coding or
   splice candidates (the pre-correction carry was 16,963: phantom
   minus-strand codon mismatches + splice LOF re-classification + multiallelic
   decomposition; see `docs/CORRECTION_2026-09-29.md`).
3. **Rarity.** Allele frequency from local gnomAD v4.1 exome sites; every
   record at a multi-record position is merged and alleles match exactly on
   (chrom, position, ref, alt). ClinVar joins as background evidence only.
4. **Candidate pairing.** Compound-heterozygous pairs (an LoF allele plus a
   second damaging allele in the same gene; phase unknown -- see the
   disclosure on trans below) and homozygous damaging candidates
   (missense, stop-gain, in-frame indel). There is no allele-frequency gate
   at candidate collection: rarity enters only through the score (step 5).
   Finds 2,211 pairs and 4,751 homozygous candidates (the emitted
   `hom_singles` counter equals the written rows; the reversed-duplicate
   inflation that produced 2,372 rows is deduplicated with a deterministic
   orientation).
5. **Scoring.** Each variant scores as a weighted mix: founder effect 0.45,
   predicted pathogenicity 0.60, rarity 0.30, plus phenotype similarity to
   the eight case HPO terms, weight 2.0. A pair scores as the stronger
   allele × (0.5 + 0.5 × the weaker).

## External cross-check

Exomiser 14.0.0 on the same VCF ranks BUB1B first (p = 0.0001) with the
same eight HPO terms; without phenotype input its top is a 46-way tie that
still contains BUB1B. The records of both runs (`exomiser_summary.tsv` on
the analysis host, artifact hashes verified 2026-09-17) are traceable; the
answer is therefore not an artifact of this pipeline's own code.

The phenotype weight's influence is characterised the same way: at weight 0
the true pair ranks 14 of 6,962 (epcr 0.6913); at the shipped weight 2.0 it
ranks 1 (record: `desensitization_audit19b.json` on the analysis host,
measured 2026-10-19 on the definitive table, after the 2026-09-18 record
measured 25 of 8,666 on the pre-corrected table and the subsequent
sensitivity run measured 21 of 6,921 on the minus-strand-corrected table).
We disclose this dependence rather than hide it.

## Reproduction

Python 3.12 with pysam, cyvcf2, pyyaml. References (paths in
`pipeline/conf/track1_config.yml`): GENCODE v47 MANE, GRCh38 FASTA, gnomAD
v4.1 exome sites, ClinVar, HPO with the OMIM gene-to-phenotype map.
Pathogenicity comes from consequence class; CADD v1.7 scores are used when
local tiles cover a site, and none did in this run, so the CADD term is
neutral and all CADD fields in the ranked table are empty.

```bash
python3 pipeline/scripts/prioritize.py \
    --conf pipeline/conf/track1_config.yml \
    --with-af \
    --with-hpo \
    --out output/annot_full.tsv
python3 pipeline/scripts/make_submission_csv.py \
    --in output/annot_full.tsv \
    --out submission_track1.csv
```

The annotate and HPO modules carry built-in self-tests
(`python3 pipeline/scripts/annotate.py`, `python3 pipeline/scripts/_hpo.py`);
they print SKIP with the missing reference named when a host lacks the
reference files.

## Disclosure

Built with AI assistance, per the challenge's 2026-08-28 disclosure guidance:
**provider** Kimi K3 via the Oh My Pi harness on the e-INFRA CZ hosted academic
gateway (plus Anthropic Claude sessions operated by the team on the same
machine); **plan/tier** per account; **data handling** — the Track-1 patient
VCF never left the whitelisted analysis host dichtator and no patient sequence
data was sent to any external API. One incident is logged openly
(2026-10-19): a debugging step printed the top-ranked rows of the corrected
ranking table (coordinates and alleles) and the first rows of the regenerated
submission CSV into the active AI session; see `AUDIT_RESOLUTION.md`,
"Disclosure note for the GenAI-use record". The pipeline, annotation joins,
full ranking tables and clinical outputs otherwise ran and stayed on the
approved host.

## Data availability

Four artefacts cited above and in the technical documentation are **not present in this
repository, by design**. They are pipeline outputs computed from the patient's sequence data, so
each one is a per-variant table for a single identifiable individual. Publishing them would
disclose case data that this project's whole architecture exists to contain, and the repository
becomes public after the hackathon closes.

| Cited artefact | Why it is absent |
|---|---|
| `desensitization.json` | Per-variant ranks for the case under varied phenotype weighting. Patient-derived. |
| `exomiser_summary.tsv` | Per-variant Exomiser output for the case. Patient-derived. |
| `output/annot_full.tsv` | The full annotated variant table for the case. Patient-derived. |
| `annot_ns.tsv` | The nonsynonymous subset of the same table. Patient-derived. |

They were used to produce the numbers reported above and are retained on the whitelisted analysis
machine. They can be made available to the organisers or to a reviewer under the data agreement
that governs the source dataset, on request.

`hpo_gene_to_phenotype.tsv` is also referenced and is likewise not shipped, for a different and
more ordinary reason: it is a public reference file distributed by the Human Phenotype Ontology
project and is fetched rather than vendored.

**Everything needed to re-run the pipeline is here**: the scripts under `pipeline/`, the
configuration under `pipeline/conf/`, and the reproduction steps above. What is withheld is the
output computed from one child's genome, not the method that produced it.
