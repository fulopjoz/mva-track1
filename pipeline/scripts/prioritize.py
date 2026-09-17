#!/usr/bin/env python3
"""Track 1 stage‑0 collector: genome‑wide honest prior prioritization.

Study collection (conf.variants.input_candidate_defs) runs in parallel against the
same filtered variant stream. No allele‑frequency gate at collection: rarity enters
through the score, not the candidate set.
  d1 – LoF heterozygous + damaging heterozygous in same gene → compound‑het pairs
  d2 – homozygous missense candidates → recessive singles
  d3 – homozygous stop‑gain candidates → recessive singles
  d4 – homozygous inframe‑indel candidates → recessive singles

Positive control (expected high score): BUB1B variants p.Leu737Ter stop
(chr15:40209701 T>G) and p.Asn1002Lys (chr15:40220612 T>G) in a compound‑heterozygous
state.

V0 contract:
  - honest first: ALL ranks come from WorkPrior internals; HPO absent = 0
    contribution (weights renormalise; foreign term stays 0).
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
    v = cyvcf2.VCF(vcf_path)
    n_in = n_pass = 0
    for rec in v:
        n_in += 1
        if rec.FILTER not in (None, 'PASS'):
            continue
        gq = rec.format('GQ')[0][0] if 'GQ' in rec.FORMAT else 99
        dp = rec.format('DP')[0][0] if 'DP' in rec.FORMAT else 0
        ad = rec.format('AD')[0] if 'AD' in rec.FORMAT else None
        if gq < qc['min_gq'] or dp < qc['min_dp']:
            continue
        alt_frac = (ad[1] / ad.sum()) if ad is not None and ad.sum() else 0.0
        if alt_frac < qc['min_alt_frac']:
            continue
        gt = rec.gt_types[0]  # 0=HOM_REF 1=HET 2=UNKNOWN 3=HOM_ALT
        if gt not in (1, 3):
            continue
        n_pass += 1
        yield rec, ('het' if gt == 1 else 'hom_alt')
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
    af_counts = {'absent': 0, 'rare': 0, 'common': 0, 'no_chrom_file': 0}
    n_clinvar = 0
    by_gene = defaultdict(list)
    singles = []
    for rec, zyg in iter_qc_variants(vcf, qc):
        if args.chrom and rec.CHROM != args.chrom:
            continue
        for ai, alt in enumerate(rec.ALT):
            a = ann.annotate(rec.CHROM, rec.POS, rec.REF, [alt])[0]
            cons = (a['consequence'] or '').split('|')[0]
            interesting = cons in (prior.lof_terms | prior.missense_terms | prior.inframe_terms) \
                or cons.startswith('splice')
            if not interesting:
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
    print(f'[instrument] coding_or_splice={n_coding} genes_with_hets={len(by_gene)} hom_singles={len(singles)}',
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

    # compound-het pairing per config rule: same gene, >=2 founders(het), phase-unknown ok
    pairs = []
    for gene, vs in by_gene.items():
        if gene is None or len(vs) < 2:
            continue
        het_lof = [x for x in vs if x.consequence in prior.lof_terms]
        het_dmg = [x for x in vs if x.consequence in (prior.lof_terms | prior.missense_terms | prior.inframe_terms)]
        for x in het_lof:
            for y in het_dmg:
                if y is x:
                    continue
                pairs.append((gene, x, y, compound_het_score(x, y)))
    pairs.sort(key=lambda p: -p[3])
    n_hom = [s for s in singles if s.internal > 0.3]
    print(f'[instrument] compound_het_pairs={len(pairs)} hom_candidates={len(n_hom)}', file=sys.stderr)

    with open(out_tsv, 'w') as fh:
        fh.write('rank\ttype\tgene\tchrom_1\tpos_1\tref_1\talt_1\tchrom_2\tpos_2\tref_2\talt_2\t'
                 'epcr\tconsequence_1\tconsequence_2\tcadd_1\tcadd_2\n')
        rank = 0
        seen = set()
        for gene, x, y, s in pairs:
            rank += 1
            fh.write(f'{rank}\tcompound_het\t{gene}\t{x.chrom}\t{x.pos}\t{x.ref}\t{x.alt}\t'
                     f'{y.chrom}\t{y.pos}\t{y.ref}\t{y.alt}\t{s:.4f}\t{x.consequence}\t{y.consequence}\t'
                     f'{"" if x.cadd_phred is None else x.cadd_phred}\t{"" if y.cadd_phred is None else y.cadd_phred}\n')
            seen.add((x.chrom, x.pos))
            seen.add((y.chrom, y.pos))
        for s_ in sorted(singles, key=lambda z: -z.internal):
            if (s_.chrom, s_.pos) in seen:
                continue
            rank += 1
            fh.write(f'{rank}\thom_alt\t{s_.gene}\t{s_.chrom}\t{s_.pos}\t{s_.ref}\t{s_.alt}\t\t\t\t\t'
                     f'{s_.internal:.4f}\t{s_.consequence}\t\t{"" if s_.cadd_phred is None else s_.cadd_phred}\t\n')
    print(f'[out] {out_tsv} rows={rank}', file=sys.stderr)


if __name__ == '__main__':
    main()
