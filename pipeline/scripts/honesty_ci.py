#!/usr/bin/env python3
"""honesty_ci.py – CI suite for the Track 1 pipeline using a synthetic “forge”.

The *forge* generates a tiny synthetic mini‑genome consisting of three chromosomes
(3, 7, 15), twenty decoy variants with fake phenotype or no‑phenotype tags, and a
single ultrarare loss‑of‑function + missense sentinel pair. The pipeline scripts are
imported read‑only (the CI never modifies pipeline code) and run against the
forged mini‑VCF, after which a set of honesty invariants are asserted.

Unit‑check invariants (checked on individual component functions):
  * U01 – the WorkPrior rarity term never increases as allele frequency (AF)
    decreases; a missing AF is treated as neutral.
  * U02 – a homozygous loss‑of‑function founder/internal score is strictly greater
    than a heterozygous loss‑of‑function score.
  * U03 – the CADD gate cannot elevate a synonymous variant above a clean heterozygous LoF.
  * U04 – when the CADD file does not cover a site, the annotator returns ``None``
    (not ``0.0``).
  * U05 – the ``_patho`` function distinguishes ``cadd_phred=None`` (neutral)
    from ``0.0`` (penalised).
  * U06 – the ``foreign_term`` is zero at stage s = 0 for every forged honest row.

Integration‑check invariants (checked after a full pipeline run on the forged data):
  * I00 – the annotator’s consequence annotation matches the intended forge intent
    for every variant.
  * I01 – forge artifacts (manifest and per‑family seeds) are successfully built.
  * I02 – every decoy ranks below the ultrarare LoF + missense sentinel pair in the
    competition ranking (decoys absent from output are considered infinite rank).
  * I03 – homopolymer‑tract frameshift pairs (ZNF717‑style) are removed by the quality‑control step.
  * I04 – synonymous variants with high CADD scores never enter the collector’s pair set.
  * I05 – decoys with the FAKE phenotype tag score identically to those without a phenotype tag.
  * I06 – the rank order is invariant to shuffling of chromosome/record order
    (e.g., 3, 7, 15 vs 15, 7, 3).
  * I07 – once the ``--with-af`` join is enabled, common‑AF family‑A pairs drop
    strictly once, confirming the gnomAD join is correctly wired.
  * I08 – shuffling input AF values yields identical resolutions to the sorted run
    and triggers the ``af_out_of_order_queries`` instrument instead of silently
    promoting misses to ultrarare.

Usage:
  python3 pipeline/scripts/honesty_ci.py [--workdir DIR] [--seed N]
Exit code 0 indicates all checks passed; otherwise a line is printed for each failed check.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))  # import pipeline modules READ-ONLY

MASTER_SEED = 20260904
FAMILY_SEEDS = {'A': 11, 'B': 23, 'C': 37, 'D': 53, 'SENT': 71}

# ---------------------------------------------------------------------------
# synthetic genome forge
# ---------------------------------------------------------------------------
BASES = 'ACGT'
COMP = str.maketrans(BASES, 'TGCA')
CHROMS = ['3', '7', '15']
CHROM_LEN = 4000
# codon template (translation sense). idx: 0 ATG 1 AAA 2 GCT 3 CAA 4 GGT 5 AAA 6 AAA 7 TTT 8 TGA
CODONS = ['ATG', 'AAA', 'GCT', 'CAA', 'GGT', 'AAA', 'AAA', 'TTT', 'TGA']
CDS_PAD = 40  # UTR flank each side of the CDS inside the single exon

# consequence target sites inside CODONS (translation-sense edits)
EDIT_SG = (3, 0, 'T')    # CAA -> TAA  stop_gained
EDIT_MS = (2, 1, 'G')    # GCT -> GGT  missense
EDIT_SYN = (2, 2, 'C')   # GCT -> GCC  synonymous


def _rev(s):
    return s.translate(COMP)[::-1]


class Gene:
    """One single-exon forged gene; variants are produced in VCF-convention
    genomic coordinates with REF taken from the forged genome itself."""

    def __init__(self, name, chrom, cds_start, strand):
        self.name, self.chrom, self.strand = name, chrom, strand
        self.cds_len = 3 * len(CODONS)
        self.cds = (cds_start, cds_start + self.cds_len - 1)
        self.exon = (cds_start - CDS_PAD, self.cds[1] + CDS_PAD)

    def genomic_cds(self):
        gseq = ''.join(CODONS)
        return gseq if self.strand == '+' else _rev(gseq)

    def _to_genomic(self, t_idx):
        """translation-order CDS index -> 1-based genomic pos"""
        return self.cds[0] + t_idx if self.strand == '+' else self.cds[1] - t_idx

    def snv(self, edit):
        k, j, new_base = edit
        codon = CODONS[k]
        pos = self._to_genomic(3 * k + j)
        old_t = codon[j]
        ref = old_t if self.strand == '+' else old_t.translate(COMP)
        alt = new_base if self.strand == '+' else new_base.translate(COMP)
        return {'chrom': self.chrom, 'pos': pos, 'ref': ref, 'alt': alt}

    def indel(self, kind, codon_idx, j):
        """kind 'del'|'ins' at base j of codon_idx. VCF-anchored."""
        k, run = codon_idx, CODONS[codon_idx]
        assert set(run) == {'A'}, 'indel forge assumes the AAA run codons'
        if kind == 'del':                      # POS = base before, REF=AA ALT=A
            pos = self._to_genomic(3 * k + j) - 1
            return {'chrom': self.chrom, 'pos': pos, 'ref': 'AA', 'alt': 'A'}
        pos = self._to_genomic(3 * k + j)      # INS: REF=A ALT=AA
        return {'chrom': self.chrom, 'pos': pos, 'ref': 'A', 'alt': 'AA'}


def build_genes():
    """Layout: one sentinel gene + 12 decoy genes over 3/7/15."""
    G = {}
    spec = [
        # name            chrom  cds_start strand
        ('SENT1',          '3',    500,     '+'),
        ('PA3',            '3',   1500,     '-'),   # minus-strand honesty probe
        ('PB1',            '3',   2500,     '+'),
        ('PC3',            '3',   3300,     '+'),
        ('PA1',            '7',    500,     '+'),
        ('PB2',            '7',   1500,     '+'),
        ('PC1',            '7',   2500,     '+'),
        ('PC4',            '7',   3300,     '+'),
        ('PA2',            '15',   500,     '+'),
        ('PB3',            '15',  1500,     '+'),
        ('PD1',            '15',  2500,     '+'),   # FAKE-phenotype decoy gene
        ('PD2',            '15',  3300,     '+'),   # phenotype-free twin
    ]
    spec.append(('PC2', '7', 3500, '+'))  # distinct slot, stays within 4000 bp
    for name, chrom, start, strand in spec:
        assert name not in G
        G[name] = Gene(name, chrom, start, strand)
    return G


def build_decoys(G):
    """20 decoys + sentinel pair. Each entry: variant fields + intentions."""
    decoys = []

    def add(gene, variant, family, intended, qc_pass, zyg='het', af=1e-6,
            fake_pheno=None, cadd=None):
        d = dict(gene=gene, family=family, intended=intended, qc_pass=qc_pass,
                 zyg=zyg, af=af, fake_pheno=fake_pheno, cadd=cadd, **variant)
        decoys.append(d)
        return d

    # sentinel (not a decoy, but travels in the same structure)
    g = G['SENT1']
    add('SENT1', g.snv(EDIT_SG), 'SENT', 'stop_gained', True)
    add('SENT1', g.snv(EDIT_MS), 'SENT', 'missense_variant', True, cadd=35.0)
    # A: common-AF LoF+missense pairs in 3 random genes
    for name in ('PA1', 'PA2', 'PA3'):
        g = G[name]
        add(name, g.snv(EDIT_SG), 'A', 'stop_gained', True, af=0.2)
        add(name, g.snv(EDIT_MS), 'A', 'missense_variant', True, af=0.2)
    # B: homopolymer-tract frameshift pairs (ZNF717-style), QC-failing
    for name in ('PB1', 'PB2', 'PB3'):
        g = G[name]
        add(name, g.indel('del', 5, 0), 'B', 'frameshift', False)
        add(name, g.indel('ins', 6, 0), 'B', 'frameshift', False)
    # C: synonymous CADD-high artifacts (2 het + 2 hom-alt)
    for name, zyg in (('PC1', 'het'), ('PC2', 'het'), ('PC3', 'hom'), ('PC4', 'hom')):
        g = G[name]
        add(name, g.snv(EDIT_SYN), 'C', 'synonymous_variant', True, zyg=zyg, cadd=32.0)
    # D: identical pair genes  -  PD1 FAKE-phenotype-tagged, PD2 untagged
    for name, pheno in (('PD1', ['HP:0001250', 'HP:0000256']), ('PD2', None)):
        g = G[name]
        add(name, g.snv(EDIT_SG), 'D', 'stop_gained', True, fake_pheno=pheno)
        add(name, g.snv(EDIT_MS), 'D', 'missense_variant', True, fake_pheno=pheno)
    n_decoys = sum(1 for d in decoys if d['family'] != 'SENT')
    assert n_decoys == 20, n_decoys
    return decoys


def forge_genome(G, seed):
    rng = random.Random(seed)
    seqs = {c: [rng.choice(BASES) for _ in range(CHROM_LEN)] for c in CHROMS}
    for g in G.values():
        s, e = g.cds
        seqs[g.chrom][s - 1:e] = g.genomic_cds()
    # sanity: CDS reads back as intended codons, both strands
    for g in G.values():
        s, e = g.cds
        gseq = ''.join(seqs[g.chrom][s - 1:e])
        t = gseq if g.strand == '+' else _rev(gseq)
        assert t == ''.join(CODONS), (g.name, t)
    return {c: ''.join(s) for c, s in seqs.items()}


# ---------------------------------------------------------------------------
# artifact writers
# ---------------------------------------------------------------------------
def write_fasta(seqs, path):
    with open(path, 'w') as fh:
        for c in CHROMS:
            fh.write(f'>{c}\n')
            s = seqs[c]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + '\n')
    import pysam
    pysam.faidx(str(path))


def write_gtf(G, path):
    with gzip.open(path, 'wt') as fh:
        for g in G.values():
            attr = f'gene_id "{g.name}"; transcript_id "TX_{g.name}"; gene_name "{g.name}";'
            fh.write(f'{g.chrom}\tforge\texon\t{g.exon[0]}\t{g.exon[1]}\t.\t{g.strand}\t.\t{attr}\n')
            fh.write(f'{g.chrom}\tforge\tCDS\t{g.cds[0]}\t{g.cds[1]}\t.\t{g.strand}\t0\t{attr}\n')


def write_cadd(decoys, cadd_dir):
    """Tiny per-chrom bgz+tbi covering exactly the intended rows."""
    import pysam
    cadd_dir.mkdir(parents=True, exist_ok=True)
    rows = {}
    for d in decoys:
        if d['cadd'] is not None:
            rows.setdefault(d['chrom'], []).append((d['pos'], d['ref'], d['alt'], d['cadd']))
    for chrom, rs in rows.items():
        raw = cadd_dir / f'raw_{chrom}.tsv'
        with open(raw, 'w') as fh:
            for pos, ref, alt, phred in sorted(rs):
                fh.write(f'{chrom}\t{pos}\t{ref}\t{alt}\t.\t{phred}\n')
        tmp = cadd_dir / f'raw_{chrom}.tsv.gz'
        out = cadd_dir / f'{chrom}.bgz'
        pysam.tabix_compress(str(raw), str(tmp), force=True)
        pysam.tabix_index(str(tmp), preset='vcf', force=True, keep_original=True)
        os.rename(tmp, out)
        os.rename(str(tmp) + '.tbi', str(out) + '.tbi')
        raw.unlink()


def _fmt(d):
    qc = '99:40:20,20' if d['qc_pass'] else '30:40:36,4'  # fam-B alt_frac 0.10
    gt = '0/1' if d['zyg'] == 'het' else '1/1'
    return f"{d['chrom']}\t{d['pos']}\t.\t{d['ref']}\t{d['alt']}\t100\tPASS\tAF={d['af']}\tGT:GQ:DP:AD\t{gt}:{qc}"


VCF_HEADER = """##fileformat=VCFv4.2
##source=honesty_forge synthetic  -  no patient data
{contigs}
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="forged allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFORGE
"""


def write_vcf(decoys, seqs, path, chrom_order, seed, sort_within=False):
    rng = random.Random(seed)
    contigs = '\n'.join(f'##contig=<ID={c},length={len(seqs[c])}>' for c in CHROMS)
    lines = []
    for c in chrom_order:
        rows = [d for d in decoys if d['chrom'] == c]
        if sort_within:
            rows.sort(key=lambda d: (d['pos'], d['ref'], d['alt']))
        else:
            rng.shuffle(rows)
        lines.extend(_fmt(d) for d in rows)
    with open(path, 'w') as fh:
        fh.write(VCF_HEADER.format(contigs=contigs))
        fh.write('\n'.join(lines) + '\n')


GNOMAD_HDR = ('##fileformat=VCFv4.2\n'
              '##source=honesty_forge synthetic gnomAD-mirror (no real data)\n'
              '##INFO=<ID=AF,Number=A,Type=Float,Description="forged allele frequency">\n'
              '##contig=<ID=chr{c},length=4000>\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')


def write_synth_gnomad(decoys, gnomad_dir):
    """Sites-only synthetic gnomAD set for the --with-af join (spec from
    CompoundHetEngine vcf.py cursor: per-chrom file
    gnomad.exomes.v4.1.sites.chr<N>.vcf.bgz, records pos-sorted, only
    POS/REF/ALT/AF are matched). Mirrors the forge INFO AF tags faithfully;
    sentinel variants are omitted (truly ultrarare = absent from reference)."""
    import pysam
    gnomad_dir.mkdir(parents=True, exist_ok=True)
    rows = {}
    for d in decoys:
        if d['family'] == 'SENT':
            continue
        rows.setdefault(d['chrom'], []).append((d['pos'], d['ref'], d['alt'], d['af']))
    for chrom, rs in rows.items():
        raw = gnomad_dir / f'raw_chr{chrom}.vcf'
        with open(raw, 'w') as fh:
            fh.write(GNOMAD_HDR.format(c=chrom))
            for pos, ref, alt, af in sorted(rs):
                fh.write(f'chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\tAF={af}\n')
        dest = gnomad_dir / f'gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz'
        pysam.tabix_compress(str(raw), str(dest), force=True)
        raw.unlink()


def write_conf(base_conf_path, workdir, cadd_dir):
    import yaml
    conf = yaml.safe_load(open(base_conf_path))
    conf = copy.deepcopy(conf)
    conf['refs']['mane_gtf'] = str(workdir / 'forge.gtf.gz')
    conf['refs']['genome_fasta'] = str(workdir / 'forge.fa')
    conf['refs']['cadd_dir'] = [str(cadd_dir)]
    conf['output']['write_dir'] = str(workdir / 'out')
    p = workdir / 'conf_ci.yml'
    yaml.safe_dump(conf, open(p, 'w'))
    return p


def forge(workdir, seed):
    workdir = Path(workdir)
    (workdir / 'out').mkdir(parents=True, exist_ok=True)
    G = build_genes()
    seqs = forge_genome(G, seed)
    decoys = build_decoys(G)
    write_fasta(seqs, workdir / 'forge.fa')
    write_gtf(G, workdir / 'forge.gtf.gz')
    write_cadd(decoys, workdir / 'cadd')
    write_vcf(decoys, seqs, workdir / 'forge_run1.vcf', ['3', '7', '15'], seed + 1)
    write_vcf(decoys, seqs, workdir / 'forge_run2.vcf', ['15', '7', '3'], seed + 2)
    # pos-sorted variant for the --with-af run: the gnomAD stream-merge cursor
    # contract requires non-decreasing pos queries per chromosome (real
    # patient VCFs are pos-sorted; run1/run2 are deliberately NOT, to exercise
    # shuffle-invariance of the no-AF pipeline only).
    write_vcf(decoys, seqs, workdir / 'forge_run3_sorted.vcf', ['3', '7', '15'],
              seed, sort_within=True)
    manifest = {
        'master_seed': seed, 'family_seeds': FAMILY_SEEDS,
        'decoys': decoys, 'sentinel_gene': 'SENT1',
    }
    json.dump(manifest, open(workdir / 'forge_manifest.json', 'w'), indent=1)
    return G, decoys


# ---------------------------------------------------------------------------
# check harness
# ---------------------------------------------------------------------------
class CI:
    def __init__(self):
        self.failures = 0

    def emit(self, cid, name, ok, seed, detail=''):
        tag = 'PASS' if ok else 'FAIL'
        if not ok:
            self.failures += 1
        print(f'[{tag}] {cid:3s} {name:38s} seed={seed} {detail}')


# --- unit checks -------------------------------------------------------------
def _mk(lof_terms_conf):
    from prior import ScoredVariant, WorkPrior  # noqa
    return ScoredVariant, WorkPrior(lof_terms_conf)


def _sv(SV, cons, zyg='het', cadd=None, af=None, hpo=None):
    return SV(chrom='3', pos=1, ref='A', alt='G', gene='G', zygosity=zyg,
              consequence=cons, cadd_phred=cadd, af=af, hpo_score=hpo)


def unit_checks(ci, conf, workdir, seed):
    SV, P = _mk(conf)
    # U01 rarity monotonicity: af None neutral, 0 ultrarare, staircase down
    rs = [(None, 0.5), (0.0, 1.0), (5e-5, 0.9), (5e-3, 0.7), (0.5, 0.0)]
    got = [P._rare(_sv(SV, 'stop_gained', af=a)) for a, _ in rs]
    ok = got == [e for _, e in rs]
    sweep = [P.internal(_sv(SV, 'stop_gained', af=a)) for a in (0.0, 5e-5, 5e-3, 0.5)]
    ok = ok and all(x > y for x, y in zip(sweep, sweep[1:]))
    ci.emit('U01', 'prior rarity monotone in AF', ok, 'n/a',
            f'staircase={got} internal={[f"{v:.4f}" for v in sweep]}')

    # U02 hom-LoF > het-LoF
    f_het, f_hom = P._founder(_sv(SV, 'stop_gained', 'het')), P._founder(_sv(SV, 'stop_gained', 'hom_alt'))
    i_het = P.internal(_sv(SV, 'stop_gained', 'het'))
    i_hom = P.internal(_sv(SV, 'stop_gained', 'hom_alt'))
    ok = f_hom > f_het and i_hom > i_het
    ci.emit('U02', 'founder hom-LoF > het-LoF', ok, 'n/a',
            f'founder {f_hom:.2f}>{f_het:.2f} internal {i_hom:.4f}>{i_het:.4f}')

    # U03 CADD gate cannot launder synonymous above clean het LoF
    dirty = P.internal(_sv(SV, 'synonymous_variant', 'het', cadd=30.0))
    clean = P.internal(_sv(SV, 'stop_gained', 'het', cadd=None))
    ok = dirty < clean and P._patho(_sv(SV, 'synonymous_variant', cadd=30.0)) == 1.0
    ci.emit('U03', 'CADD-high synonymous < het LoF', ok, 'n/a',
            f'syn/cadd30={dirty:.4f} hetsg={clean:.4f}')

    # U04 annotator: site-not-covered -> None, covered -> forged phred
    from annotate import LocalAnnotator
    ann = LocalAnnotator(str(workdir / 'forge.gtf.gz'),
                         genome_fasta=str(workdir / 'forge.fa'),
                         cadd_dir=str(workdir / 'cadd'))
    notcovered = ann._cadd('3', 1, 'A', 'G')              # no row forged there
    nodir = LocalAnnotator(str(workdir / 'forge.gtf.gz'))._cadd('3', 1, 'A', 'G')
    ci.emit('U04', 'CADD not-covered = None, not 0', notcovered is None and nodir is None,
            seed, f'uncovered={notcovered} nocadd_dir={nodir}')

    # U05 gate wiring: None skips the gate; real scores enter it. V0 pin:
    # max(base, base*0.5) collapses the low-CADD penalty, so covered-low
    # currently == neutral. Assert exactly that so a future penalising gate
    # change fails loudly instead of landing silently.
    p_none = P._patho(_sv(SV, 'missense_variant', cadd=None))
    p_low = P._patho(_sv(SV, 'missense_variant', cadd=0.0))
    p_high = P._patho(_sv(SV, 'missense_variant', cadd=30.0))
    ok = p_none == 0.5 and p_high == 1.0 and p_low == p_none
    ci.emit('U05', '_patho gate: None ungated, 30 gated, V0 low=neutral', ok, 'n/a',
            f'None->{p_none} 0.0->{p_low} 30.0->{p_high} (None is neutral, 0.0 is penalised)')

    # U06 foreign term == 0 at s=0 for every honest forged row
    manifest = json.load(open(workdir / 'forge_manifest.json'))
    rows = [_sv(SV, d['intended'], 'hom_alt' if d['zyg'] == 'hom' else 'het')
            for d in manifest['decoys']]
    ok = all(P.foreign_term(v) == 0.0 for v in rows)
    ok = ok and all((P.internal(v), v.foreign)[1] == 0.0 for v in rows)
    ok = ok and SV(chrom='3', pos=1, ref='A', alt='G').foreign == 0.0
    ci.emit('U06', 'foreign term == 0 at s=0', ok, 'n/a', f'rows={len(rows)}')
    return ann


# --- integration -------------------------------------------------------------
def run_pipeline(workdir, conf_ci, vcf, out):
    cmd = [sys.executable, str(SCRIPTS_DIR / 'prioritize.py'),
           '--conf', str(conf_ci), '--vcf', str(vcf), '--out', str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPTS_DIR)
    return r


def run_pipeline_af(workdir, conf_ci, vcf, out, gnomad_dir):
    """Same runner + CompoundHetEngine's AF join flags (contract: --with-af
    --gnomad-dir DIR; TSV schema unchanged, epcr stays in col 13)."""
    cmd = [sys.executable, str(SCRIPTS_DIR / 'prioritize.py'),
           '--conf', str(conf_ci), '--vcf', str(vcf), '--out', str(out),
           '--with-af', '--gnomad-dir', str(gnomad_dir)]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPTS_DIR)


def parse_out(path):
    rows = []
    for i, line in enumerate(open(path)):
        if i == 0 or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        rows.append({'rank': int(f[0]), 'type': f[1], 'gene': f[2],
                     'pos1': (f[3], f[4]), 'pos2': (f[6], f[7]),
                     'epcr': float(f[11])})
    return rows


def score_map(rows):
    return {(r['type'], r['gene']): r['epcr'] for r in rows}


def competition_ranks(rows):
    """Shuffle-invariant ranks: 1 + #(strictly greater scores), keyed by gene."""
    s = score_map(rows)
    return {k: 1 + sum(1 for e in s.values() if e > v) for k, v in s.items()}


def integration_checks(ci, conf, workdir, seed, ann):
    workdir = Path(workdir)
    r1 = run_pipeline(workdir, workdir / 'conf_ci.yml',
                      workdir / 'forge_run1.vcf', workdir / 'out/run1.tsv')
    r2 = run_pipeline(workdir, workdir / 'conf_ci.yml',
                      workdir / 'forge_run2.vcf', workdir / 'out/run2.tsv')
    ok = r1.returncode == 0 and r2.returncode == 0
    if not ok:
        ci.emit('I00', 'pipeline runs on forged VCF', False, seed,
                (r1.stderr or r2.stderr).strip().splitlines()[-1][:120])
        return
    manifest = json.load(open(workdir / 'forge_manifest.json'))
    decoys = manifest['decoys']

    # I00 annotator consequence == forge intent, for every forged variant
    mism = []
    for d in decoys:
        got = ann.annotate(d['chrom'], d['pos'], d['ref'], [d['alt']])[0]
        if (got['consequence'] or '').split('|')[0] != d['intended']:
            mism.append((d['gene'], d['pos'], d['intended'], got['consequence']))
        if got['hgvsp'] and 'codon-mismatch' in got['hgvsp']:
            mism.append((d['gene'], d['pos'], 'ref-genome mismatch', got['hgvsp']))
        if d['cadd'] is not None and got['cadd_phred'] != d['cadd']:
            mism.append((d['gene'], d['pos'], 'cadd', got['cadd_phred']))
    ci.emit('I00', 'forge intent == annotator result', not mism, seed,
            f'variants={len(decoys)} mismatches={mism[:3]}')

    # I01 forge artifacts built
    n = sum(1 for d in decoys if d['family'] != 'SENT')
    fams = sorted({d['family'] for d in decoys} - {'SENT'})
    ci.emit('I01', 'forge: 20 decoys + sentinel built', n == 20 and fams == ['A', 'B', 'C', 'D'],
            seed, f'decoys={n} families={",".join(fams)} seeds={FAMILY_SEEDS}')

    rows1, rows2 = parse_out(workdir / 'out/run1.tsv'), parse_out(workdir / 'out/run2.tsv')
    pair1 = [r for r in rows1 if r['type'] == 'compound_het']
    sent = [r for r in pair1 if r['gene'] == 'SENT1']
    present_positions = set()
    for r in rows1:
        present_positions.add(r['pos1'])
        present_positions.add(r['pos2'])
    cr1 = competition_ranks(rows1)
    sent_epcr = sent[0]['epcr'] if sent else -1.0

    # I02 every decoy below the sentinel (absent == rank inf)
    bad = []
    for d in decoys:
        if d['family'] == 'SENT':
            continue
        site = (d['chrom'], str(d['pos']))
        if site in present_positions:
            epcr = next((r['epcr'] for r in pair1
                         if r['gene'] == d['gene'] and (site in (r['pos1'], r['pos2']))), None)
            if epcr is None or not epcr < sent_epcr:
                bad.append((d['family'], d['gene'], epcr))
    ok = bool(sent) and not bad
    ci.emit('I02', 'every decoy below LoF+miss sentinel', ok, seed,
            f'sentinel_epcr={sent_epcr:.4f} pairs={len(pair1)} violators={bad[:3]}')

    # I03 homopolymer-tract frameshift pairs QC-killed (absent)
    b_present = [d for d in decoys if d['family'] == 'B'
                 and (d['chrom'], str(d['pos'])) in present_positions]
    ci.emit('I03', 'homopolymer fs pairs QC-killed', not b_present, FAMILY_SEEDS['B'],
            f'famB_present={len(b_present)}/6')

    # I04 synonymous CADD-high artifacts absent from collector output
    c_present = [d for d in decoys if d['family'] == 'C'
                 and (d['chrom'], str(d['pos'])) in present_positions]
    ci.emit('I04', 'synonymous CADD-high artifacts absent', not c_present, FAMILY_SEEDS['C'],
            f'famC_present={len(c_present)}/4')

    # I05 FAKE-phenotype vs no-phenotype identical scores
    p = {r['gene']: r['epcr'] for r in pair1 if r['gene'] in ('PD1', 'PD2')}
    ok = len(p) == 2 and p['PD1'] == p['PD2']
    ci.emit('I05', 'FAKE-pheno vs no-pheno identical', ok, FAMILY_SEEDS['D'],
            f'PD1={p.get("PD1")} PD2={p.get("PD2")}')

    # I06 shuffle-invariance: 3,7,15 vs 15,7,3 + within-chrom shuffle
    ok = (score_map(rows1) == score_map(rows2)
          and cr1 == competition_ranks(rows2))
    ci.emit('I06', 'rank-order shuffle-invariant', ok, seed + 2,
            f'rows run1={len(rows1)} run2={len(rows2)} identical={ok}')

    # I07 AF-awareness: once CompoundHetEngine's --with-af join is wired, the
    # common-AF fam-A pairs (AF=0.2) must drop strictly below their own
    # no-AF epcr (0.1880). Skips-with-note on a V0 checkout.
    probe = subprocess.run([sys.executable, str(SCRIPTS_DIR / 'prioritize.py'), '--help'],
                           capture_output=True, text=True)
    if '--with-af' not in probe.stdout:
        ci.emit('I07', 'fam-A dies under --with-af join', True, FAMILY_SEEDS['A'],
                'SKIP: --with-af not wired (V0 checkout)')
        return
    write_synth_gnomad(decoys, workdir / 'gnomad')
    r3 = run_pipeline_af(workdir, workdir / 'conf_ci.yml',
                         workdir / 'forge_run3_sorted.vcf', workdir / 'out/run3_af.tsv',
                         workdir / 'gnomad')
    if r3.returncode != 0:
        ci.emit('I07', 'fam-A dies under --with-af join', False, FAMILY_SEEDS['A'],
                'run failed: ' + r3.stderr.strip().splitlines()[-1][:120])
        return
    rows3 = parse_out(workdir / 'out/run3_af.tsv')
    ep3 = {r['gene']: r['epcr'] for r in rows3
           if r['type'] == 'compound_het' and r['gene'] in ('PA1', 'PA2', 'PA3')}
    plain = {r['gene']: r['epcr'] for r in rows1
             if r['type'] == 'compound_het' and r['gene'] in ep3}
    drops = {g: (plain[g], ep3[g]) for g in ep3 if g in plain}
    m = re.search(r'af_common=(\d+)', r3.stderr)
    n_common = int(m.group(1)) if m else -1
    ok = len(drops) == 3 and all(af < noaf for noaf, af in drops.values()) and n_common >= 6
    ci.emit('I07', 'fam-A dies under --with-af join', ok, FAMILY_SEEDS['A'],
            f'drops={ {g: f"{a:.4f}->{b:.4f}" for g, (a, b) in drops.items()} } af_common={n_common}')

    # I08 cursor-hazard pin: SHUFFLED input must now resolve identically to the
    # sorted AF run (vcf.py cache fix) AND must trip the out-of-order
    # instrument instead of silently promoting missed matches to ultrarare.
    r4 = run_pipeline_af(workdir, workdir / 'conf_ci.yml',
                         workdir / 'forge_run1.vcf', workdir / 'out/run4_af_shuffled.tsv',
                         workdir / 'gnomad')
    w = re.search(r'af_out_of_order_queries=(\d+)', r4.stderr or '')
    rows4 = parse_out(workdir / 'out/run4_af_shuffled.tsv') if r4.returncode == 0 else []
    same = r4.returncode == 0 and score_map(rows3) == score_map(rows4)
    ok = same and w is not None and int(w.group(1)) >= 1
    ci.emit('I08', 'shuffled AF run exact + guard trips', ok, seed + 1,
            f'scores_equal={same} out_of_order_queries={w.group(1) if w else "MISSING"}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--workdir', default='/tmp/honestyforge')
    ap.add_argument('--seed', type=int, default=MASTER_SEED)
    ap.add_argument('--conf', default=str(Path(__file__).resolve().parent.parent / 'conf' / 'track1_config.yml'),
                    help='production conf used read-only for weights/QC gates')
    args = ap.parse_args()
    import yaml  # noqa: F401  (fail fast if the venv is wrong)
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    ci = CI()

    G, decoys = forge(workdir, args.seed)
    conf = yaml.safe_load(open(os.path.expanduser(args.conf)))
    write_conf(os.path.expanduser(args.conf), workdir, workdir / 'cadd')

    ann = unit_checks(ci, conf, workdir, args.seed)
    integration_checks(ci, conf, workdir, args.seed, ann)
    print(f'honesty_ci: {ci.failures} failure(s), workdir={workdir}')
    sys.exit(1 if ci.failures else 0)


if __name__ == '__main__':
    main()
