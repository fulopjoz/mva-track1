#!/usr/bin/env python3
"""Track 1 variant prioritization: a weighted prior over annotation, rarity and phenotype evidence.

Four candidate studies in one pass over the same filtered variant stream. No
allele‑frequency gate at collection: rarity enters through the score, not the
candidate set.
  d1 – LoF heterozygous + damaging heterozygous in same gene → compound‑het pairs
  d2 – homozygous missense candidates → recessive singles
  d3 – homozygous stop‑gain candidates → recessive singles
  d4 – homozygous inframe‑indel candidates → recessive singles

Positive control (expected high score): BUB1B variants p.Leu737Ter stop
(chr15:40209701 T>G) and p.Asn1002Lys (chr15:40220612 T>G) in a compound‑heterozygous
state.

Scoring contract:
  - ALL ranks come from WorkPrior internals; HPO absent contributes 0
    (weights renormalise; the foreign term stays 0).
  - instrumented: every cohort‑yield stage prints sizes (filters yield, coding yield,
    rare yield, pair candidates) so anomalies surface.
  - emits the ranked table whose top rows become the challenge submission
    CSV (`epcr` = estimated probability of causal relationship, the
    challenge-defined column).
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate import LocalAnnotator  # noqa: E402
from prior import ScoredVariant, WorkPrior, compound_het_score  # noqa: E402

import cyvcf2  # noqa: E402
import yaml  # noqa: E402


def load_conf(path):
    return yaml.safe_load(open(path))


def iter_qc_variants(vcf_path, qc):
    """Yield (record, zygosity, alt_index) for each allele actually present in
    the sample genotype of a QC-passing record (audit T1-02).

    The sample genotype is decomposed: multiallelic records contribute only
    the ALT alleles called in GT, and the alternate fraction is that
    allele's own AD divided by the record's depth sum. Records have a single
    sample, per the pipeline's design contract.
      - './.', './1', '1/.', '.|.' -> skipped (zygosity cannot be established;
        GQ/DP thresholds still apply to the record).
      - missing or malformed AD/DP -> skipped (fraction/depth unverifiable).
      - zygosity: a1 != a2 -> 'het';  a1 == a2 > 0 -> 'hom_alt'.
      - phased genotypes (a|b) decompose identically.
    qc_pass counts (record, allele) pairs, so multiallelic het calls count once
    per present allele; records with no emitted allele count zero.
    """
    v = cyvcf2.VCF(vcf_path)
    n_in = n_pass = 0
    for rec in v:
        n_in += 1
        if rec.FILTER not in (None, 'PASS'):
            continue
        gq = rec.format('GQ')[0][0] if 'GQ' in rec.FORMAT else 99
        dp = rec.format('DP')[0][0] if 'DP' in rec.FORMAT else 0
        if gq < qc['min_gq'] or dp < qc['min_dp']:
            continue
        a1, a2 = rec.genotypes[0][0], rec.genotypes[0][1]
        if a1 < 0 or a2 < 0:            # partial or fully missing call
            continue
        alt_idx = sorted({a for a in (a1, a2) if a > 0})
        if not alt_idx:
            continue                    # 0/0
        ad = rec.format('AD')[0] if 'AD' in rec.FORMAT else None
        if ad is None:
            continue                    # cannot verify per-allele fraction
        total = ad.sum() if len(ad) else 0
        if total <= 0:
            continue
        zyg = 'het' if a1 != a2 else 'hom_alt'
        for k in alt_idx:
            if k >= len(ad) or ad[k] < 0:
                continue                # malformed AD entry for this allele
            if ad[k] / total < qc['min_alt_frac']:
                continue
            n_pass += 1
            yield rec, zyg, k
    print(f'[instrument] vcf_in={n_in} qc_pass={n_pass}', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--conf', required=True)
    ap.add_argument('--vcf', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--limit', type=int, default=0, help='debug: first N records')
    ap.add_argument('--chrom', default=None, help='debug: restrict to one chromosome')
    ap.add_argument('--with-af', action='store_true',
                    help='join gnomAD exomes v4.1 AF + ClinVar via vcf.py; rare term becomes data-driven')
    ap.add_argument('--gnomad-dir', default=None, help='override conf refs.gnomad_v41_exomes_dir')
    ap.add_argument('--clinvar-vcf', default=None, help='override conf refs.clinvar_vcf')
    ap.add_argument('--with-hpo', action='store_true',
                    help='activate HPO phenotype similarity term (refs-only, per-gene score)')
    ap.add_argument('--hp-obo', default=None, help='override conf refs.hp_obo')
    ap.add_argument('--hpo-map', default=None, help='override conf refs.hpo_map')
    args = ap.parse_args()
    conf = load_conf(args.conf)
    vcf = os.path.expanduser(args.vcf or conf['data']['vcf'])
    out_tsv = os.path.expanduser(args.out or os.path.join(conf['output']['write_dir'], 'annot_ns.tsv'))
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

    gtf = os.path.expanduser(conf['refs']['mane_gtf'])
    fasta = os.path.expanduser(conf['refs']['genome_fasta'])
    ann = LocalAnnotator(
        gtf,
        genome_fasta=fasta,
        cadd_dir=next((os.path.expanduser(d) for d in conf['refs']['cadd_dir']
                       if os.path.isdir(os.path.expanduser(d))), None),
    )
    prior = WorkPrior(conf)
    qc = conf['variants']['qc']

    gnomad = clinvar = None
    if args.with_af:
        from vcf import ClinVarJoin, GnomadAfJoin
        refs = conf['refs']
        gnomad = GnomadAfJoin(args.gnomad_dir or refs.get('gnomad_v41_exomes_dir', '~/hf-rare/ref'))
        clinvar = ClinVarJoin(args.clinvar_vcf or refs.get('clinvar_vcf', '~/hf-rare/ref/clinvar.vcf.gz'))
        print(f'[instrument] refs_join=on gnomad_dir={gnomad.ref_dir} clinvar={clinvar._tb is not None}',
              file=sys.stderr)

    hpo = None
    if args.with_hpo:
        from _hpo import HPOScorer
        refs = conf['refs']
        hpo = HPOScorer(os.path.expanduser(args.hp_obo or refs['hp_obo']),
                        os.path.expanduser(args.hpo_map or refs['hpo_map']),
                        conf['phenotype']['hpo_terms'])
        print(f'[instrument] hpo_active=1 genes_mapped={len(hpo.gene_terms)}', file=sys.stderr)
    _hpo_cache = {}

    def gene_hpo(gene):
        if hpo is None or gene is None:
            return None
        if gene not in _hpo_cache:
            _hpo_cache[gene] = hpo.gene_score(gene)
        return _hpo_cache[gene]

    n_coding = 0
    n_other = 0
    af_counts = {'absent': 0, 'rare': 0, 'common': 0, 'no_chrom_file': 0}
    n_clinvar = 0
    by_gene = defaultdict(list)
    singles = []
    for rec, zyg, k in iter_qc_variants(vcf, qc):
        if args.chrom and rec.CHROM != args.chrom:
            continue
        alt = rec.ALT[k - 1]  # GT allele index k (>=1) -> ALT[k-1]
        a = ann.annotate(rec.CHROM, rec.POS, rec.REF, [alt])[0]
        cons = (a['consequence'] or '').split('|')[0]
        interesting = cons in (prior.lof_terms | prior.missense_terms | prior.inframe_terms) \
            or cons.startswith('splice')
        if not interesting:
            if cons == 'other':
                n_other += 1
            continue
        n_coding += 1
        af = None
        clin_sig = clin_dn = None
        if args.with_af:
            af = gnomad.af(rec.CHROM, rec.POS, rec.REF, alt)
            if af is None:
                af_counts['no_chrom_file'] += 1
            elif af == 0:
                af_counts['absent'] += 1
            elif af < prior.max_rare:
                af_counts['rare'] += 1
            else:
                af_counts['common'] += 1
            clin_sig, clin_dn = clinvar.lookup(rec.CHROM, rec.POS, rec.REF, alt)
            if clin_sig:
                n_clinvar += 1
        sv = ScoredVariant(rec.CHROM, rec.POS, rec.REF, alt,
                           gene=a['gene'], zygosity=zyg,
                           consequence=cons, cadd_phred=a['cadd_phred'],
                           af=af, hpo_score=gene_hpo(a['gene']))
        if args.with_af and clin_sig:
            sv.clinvar_sig = clin_sig
            sv.clinvar_disease = clin_dn
        prior.internal(sv)
        if zyg == 'het':
            by_gene[a['gene']].append(sv)
        else:
            singles.append(sv)
    print(f'[instrument] coding_or_splice={n_coding} genes_with_hets={len(by_gene)} hom_singles_collected={len(singles)}',
          file=sys.stderr)
    if args.with_af:
        print(f'[instrument] af_absent={af_counts["absent"]} af_rare={af_counts["rare"]} '
              f'af_common={af_counts["common"]} af_no_chrom_file={af_counts["no_chrom_file"]} '
              f'clinvar_hits={n_clinvar}', file=sys.stderr)
        ooo = gnomad.out_of_order_queries()
        if ooo:
            print(f'[instrument] WARNING af_out_of_order_queries={ooo}  -  input VCF not '
                  f'pos-sorted per chromosome; results exact via scan cache, but check input',
                  file=sys.stderr)

    # Compound-het candidate pairing per config rule: same gene, >=2 LoF/het
    # founders; phase unknown accepted (candidates, NOT confirmed trans).
    # (x,y) and (y,x) are ONE candidate; members are written position-sorted
    # for a deterministic orientation (audit T1-03).
    def _vkey(x):
        ch = (1, int(x.chrom)) if str(x.chrom).isdigit() else (2, str(x.chrom))
        return ch + (x.pos, x.ref, x.alt)

    pairs = []
    seen_pair = set()
    for gene, vs in by_gene.items():
        if gene is None or len(vs) < 2:
            continue
        het_lof = [x for x in vs if x.consequence in prior.lof_terms]
        het_dmg = [x for x in vs if x.consequence in (prior.lof_terms | prior.missense_terms | prior.inframe_terms)]
        for x in het_lof:
            for y in het_dmg:
                if y is x:
                    continue
                ka, kb = sorted((x, y), key=_vkey)
                pk = _vkey(ka) + _vkey(kb)
                if pk in seen_pair:
                    continue
                seen_pair.add(pk)
                pairs.append((gene, ka, kb, compound_het_score(ka, kb)))

    # Emission ordering, one global list (audit T1-03): rows are ranked by the
    # epcr score descending; ties break deterministically on record type
    # (candidate pairs before homozygous singles), gene, then member key —
    # never on input order. A variant that is part of any emitted pair is not
    # repeated as a single.
    seen = set()
    for gene, ka, kb, s in pairs:
        seen.add((ka.chrom, ka.pos))
        seen.add((kb.chrom, kb.pos))
    emission = []
    for gene, ka, kb, s in pairs:
        emission.append((-s, 0, _vkey(ka) + _vkey(kb), gene, ('pair', ka, kb)))
    for s_ in singles:
        if (s_.chrom, s_.pos) in seen:
            continue
        emission.append((-s_.internal, 1, _vkey(s_), s_.gene or '', ('hom', s_)))
    emission.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
    n_hom = sum(1 for e in emission if e[4][0] == 'hom')
    print(f'[instrument] compound_het_pairs={len(pairs)} hom_singles_emitted={n_hom} excluded_other={n_other}', file=sys.stderr)

    with open(out_tsv, 'w') as fh:
        fh.write('rank\ttype\tgene\tchrom_1\tpos_1\tref_1\talt_1\tchrom_2\tpos_2\tref_2\talt_2\t'
                 'epcr\tconsequence_1\tconsequence_2\tcadd_1\tcadd_2\n')
        rank = 0
        for _, _, _, _, payload in emission:
            rank += 1
            if payload[0] == 'pair':
                x, y = payload[1], payload[2]
                s = compound_het_score(x, y)
                fh.write(f'{rank}\tcompound_het\t{x.gene}\t{x.chrom}\t{x.pos}\t{x.ref}\t{x.alt}\t'
                         f'{y.chrom}\t{y.pos}\t{y.ref}\t{y.alt}\t{s:.4f}\t{x.consequence}\t{y.consequence}\t'
                         f'{"" if x.cadd_phred is None else x.cadd_phred}\t{"" if y.cadd_phred is None else y.cadd_phred}\n')
            else:
                s_ = payload[1]
                fh.write(f'{rank}\thom_alt\t{s_.gene}\t{s_.chrom}\t{s_.pos}\t{s_.ref}\t{s_.alt}\t\t\t\t\t'
                         f'{s_.internal:.4f}\t{s_.consequence}\t\t{"" if s_.cadd_phred is None else s_.cadd_phred}\t\n')
    print(f'[out] {out_tsv} rows={rank}', file=sys.stderr)


if __name__ == '__main__':
    main()
