# Phase probe - haplotype phase of the two pathogenic BUB1B alleles

**Genome:** WGS_EX2312012 (short-read paired-end whole-genome sequencing (WGS), 2 x 150 bp reads, mean depth about 53x)

| Variant | chr15 pos | Ref>Alt | Consequence |
|---|---|---|---|
| BUB1B p.Leu737Ter | 40,209,701 | T>G | nonsense |
| BUB1B p.Asn1002Lys | 40,220,612 | T>G | missense |

The two sites are 10,911 bp apart.

## VERDICT: measured-unresolvable

Whether the two alleles are in cis (both on the same chromosome copy) or in trans (one on each copy) cannot be determined from this sequencing library because no read or read pair connects either site to a second heterozygous site.

| Target site | nearest other heterozygous site | distance | can a fragment from this 150 bp-read library span it? |
|---|---|---|---|
| 40,209,701 | 40,216,470 | 6,769 bp | no |
| 40,220,612 | 40,221,421 | 809 bp | only 2.55% of fragments are this long; none observed |

Measured insert sizes (the fragment lengths spanned by each read pair; n = 216,196 proper pairs, i.e. pairs mapped in the expected orientation): median 383 bp, 95th percentile 726 bp.
Zero reads span either gap. [whatshap](https://whatshap.readthedocs.io/)
v2.8 read-based phasing over the BUB1B region leaves both sites as singletons (variants not grouped into any phase block) and no phase block is formed; the nearest phase block ends 7,568 bp from the downstream variant.
Allele balance in the independently extracted reads matches the VCF (variant call format; both ~50% reference / 50% alternate), so the read set represents the sites correctly.

This is a physical measurement limit, not a processing failure. Statistical panel phasing (SHAPEIT5/Eagle) – a method that infers haplotype phase by comparing variants to a reference panel of known haplotypes – does not change the verdict for variants absent from the panel. Published SHAPEIT5 validation reports a posterior probability (PP) of 0.5 (uninformative) and an approximate 35% switch error rate for singletons (variants not represented in the reference panel). Even when applying a stringent filter of PP ≥ 0.9, the singleton switch-error rate remains 12.1% (Lassen et al., Cell Genomics 2024).
Both BUB1B alleles are absent from gnomAD, so panel phasing of this pair
would be exactly that uninformative case and was not run.

## What would settle phase

- Long-read sequencing (Oxford Nanopore Technologies (ONT) or PacBio HiFi) that can span a variant and a heterozygous marker on a single DNA molecule,
- Parental genomes (trio segregation analysis),
- A barcoded or linked-read library (reads sharing a barcode from the same long DNA fragment).

## Reproduction

`phase_probe.sh` plus its helper scripts (shipped in this directory in
`phase_probe_helpers/`, all from the original host run) reproduce every
step on the analysis host: copy the helpers next to the script into a
working directory and run `bash phase_probe.sh`.
It selects the reads overlapping the BUB1B window (134,953 candidate pairs),
realigns them to the full GRCh38 reference without alternate contigs (using
minimap2; retaining only alignments with mapping quality (MAPQ) ≥ 20, which
reflects competition against the whole genome), phases the window with
whatshap, and measures insert sizes; the resulting artifacts (cand.bam,
phased_bub1b.vcf.gz, block tables) stay on the host.
