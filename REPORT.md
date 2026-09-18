# Track 1 pipeline report

Team fulopjoz. MVA Hackathon 2026, Track 1, variant prioritization.
Submission file: `submission_track1.csv`.

## Result

The pipeline ranks the clinically confirmed genotype first of 8,666
candidates: compound-heterozygous **BUB1B p.Leu737Ter**
(chr15:40,209,701 T>G, stop gained; ClinVar allele ID 529272) plus
**BUB1B p.Asn1002Lys** (chr15:40,220,612 T>G, missense), score 0.5823.
Score at rank 2 is 0.5243. The score column `epcr` — estimated probability
of causal relationship, the challenge-defined (0, 1] score — equals the
pipeline ranking score.

| Rank | Genotype | Score |
|---|---|---|
| 1 | BUB1B p.Leu737Ter + p.Asn1002Lys | 0.5823 |
| 2-7 | RAI1 frameshift-deletion pairs (two orientations) | 0.5243 / 0.4914 |
| 8-10 | KIF7 pairs (founder allele is a common stop polymorphism, gnomAD frequency 0.561) | 0.4515 |
| 2,373 and below | homozygous-recessive singles (ordered after biallelic pairs by design) | — |

## Method

One pass over the 5,012,204-record proband VCF. Each stage's counter is
measured from the run.

1. **Quality filter.** FILTER is PASS or absent, genotype quality at least 20,
   read depth at least 10, alternate-allele read fraction at least 0.25. Keeps 4,592,140.
2. **Effect calling.** Variants are annotated against local GENCODE v47
   MANE transcripts. Stop gained, frameshift, and splice-site changes count
   as loss of function (LoF); missense and in-frame indels are scored
   separately. Keeps 16,963 coding or splice candidates.
3. **Rarity.** Allele frequency from local gnomAD v4.1 exome sites; every
   record at a multi-record position is merged and alleles match exactly on
   (chrom, position, ref, alt). ClinVar joins as background evidence only.
4. **Candidate pairing.** Compound-heterozygous pairs (an LoF allele plus a
   second damaging allele in the same gene; phase not required), rare
   homozygous missense, rare homozygous stop gain, rare homozygous in-frame
   indels. Finds 2,372 pairs and 6,294 homozygous candidates.
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
the true pair ranks 25 of 8,666 (epcr 0.6913); at the shipped weight 2.0 it
ranks 1 (record: `desensitization.json` on the analysis host, verified
2026-09-18). We disclose this dependence rather than hide it.

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

Built with AI assistance (Kimi K3 via Oh My Pi); no patient data left the
analysis hosts and none was sent to any external API, per the challenge
rules.
