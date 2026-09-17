# Phase probe - haplotype phase of the two pathogenic BUB1B alleles

**Genome:** WGS_EX2312012 (short-read paired-end WGS, 2 x 150 bp, ~53x)

| Variant | chr15 pos | Ref>Alt | Consequence |
|---|---|---|---|
| BUB1B p.Leu737Ter | 40,209,701 | T>G | nonsense |
| BUB1B p.Asn1002Lys | 40,220,612 | T>G | missense |

The two sites are 10,911 bp apart.

## VERDICT: measured-unresolvable

Cis/trans cannot be determined from this library, because the evidence chain
breaks at both sites:

| Target site | nearest other het site | distance | can a 150 bp fragment bridge it? |
|---|---|---|---|
| 40,209,701 | 40,216,470 | 6,769 bp | no |
| 40,220,612 | 40,222,421 | 809 bp | only 2.55% of fragments are this long; none observed |

Measured insert sizes (n = 216,196 proper pairs): median 383 bp, p95 726 bp.
Zero reads span either gap. [whatshap](https://whatshap.readthedocs.io/)
v2.8 phasing over the window leaves both sites as singletons in no phase
block; the nearest block ends 7,568 bp from the downstream variant. Allele
balance in the independently extracted reads matches the VCF (both ~50/50),
so the read set represents the sites correctly.

This is a physical measurement limit, not a processing failure. Statistical
panel phasing (SHAPEIT5/Eagle) does not change the verdict for variants
absent from the panel - published SHAPEIT5 validation reports posterior
probability 0.5 (uninformative) and ~35% switch error for singletons
([Hofmeister et al., Nature Genetics 2023](https://www.nature.com/articles/s41588-023-01415-w));
even at the PP >= 0.9 filter the singleton switch-error rate is 12.1%
([Lassen et al., Cell Genomics 2024](https://doi.org/10.1016/j.xgen.2024.100602)).
Both BUB1B alleles are absent from gnomAD, so panel phasing of this pair
would be exactly that uninformative case and was not run.

## What would settle phase

- Long-read sequencing (ONT or PacBio HiFi) spanning a variant and a
  heterozygous marker on one molecule, or
- parental genomes (trio segregation), or
- a barcoded/linked-read library.

## Reproduction

`phase_probe.sh` plus its helper scripts (shipped in this directory in
`phase_probe_helpers/`, all from the original host run) reproduce every
step on the analysis host: copy the helpers next to the script into a
working directory on dichtator, adjust the `cd` line, and `bash
phase_probe.sh`. It selects the reads overlapping the BUB1B window
(134,953 candidate pairs), realigns them honestly to full GRCh38-no-alt
(minimap2, MAPQ >= 20), phases the window with whatshap, and measures
insert sizes; the resulting artifacts (cand.bam, phased_bub1b.vcf.gz,
block tables) stay on the host.
