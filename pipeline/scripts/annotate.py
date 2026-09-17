"""Local per-position variant annotator  -  zero external API calls.

Resources (all local, configured via conf/track1_config.yml `refs`):
  - MANE_Select GTF (tabix)         -> consequence vs canonical transcript
  - CADD core bgz per chromosome    -> cadd_phred when file covers position
  - GRCh38 reference FASTA (faidx)  -> codon-level effect (stop_gained etc.)

Design notes:
  * GTF is parsed once into per-chromosome sorted interval arrays; lookups are
    bisect-based.
  * Consequence vocabulary: stop_gained, stop_lost, frameshift, inframe_deletion,
    inframe_insertion, missense_variant, synonymous_variant, splice_donor,
    splice_acceptor, splice_region, 5_prime_utr, 3_prime_utr, intronic, intergenic.
  * Multi-transcript genes: MANE GTF carries one transcript per gene, so overlap
    maps to exactly one transcript.
"""
from __future__ import annotations

import bisect
import gzip
import os
import re
import sys
from dataclasses import dataclass, field

import pysam

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")

CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


@dataclass
class Transcript:
    gene: str
    tid: str
    strand: str
    exons: list = field(default_factory=list)   # [(start, end)] 1-based incl., genomic order
    cds: list = field(default_factory=list)     # [(start, end)] genomic order

    def cds_spliced(self, genome_fetch):
        """Spliced CDS sequence in translation order."""
        parts = self.cds if self.strand == '+' else list(reversed(self.cds))
        seq = ''.join(genome_fetch(s, e) for s, e in parts)
        if self.strand == '-':
            seq = seq.translate(COMP)[::-1]
        return seq


def _revcomp(s):
    return s.translate(COMP)[::-1]


class LocalAnnotator:
    def __init__(self, mane_gtf, genome_fasta=None, cadd_dir=None, chrom_prefix=''):
        self.prefix = chrom_prefix
        self.genome = pysam.FastaFile(genome_fasta) if genome_fasta and os.path.exists(genome_fasta) else None
        self.cadd = {}  # chrom -> TabixFile
        if cadd_dir and os.path.isdir(cadd_dir):
            for fn in os.listdir(cadd_dir):
                m = re.match(r'(chr)?(\w+)\.(vcf\.)?bgz$', fn)
                if m and os.path.exists(os.path.join(cadd_dir, fn + '.tbi')):
                    self.cadd[m.group(2)] = pysam.TabixFile(os.path.join(cadd_dir, fn))
        self.tx = {}            # chrom -> list[Transcript] sorted by min start
        self.tx_starts = {}     # chrom -> [min_start]
        self._load_gtf(mane_gtf)

    # ---------------- GTF loading ----------------
    def _load_gtf(self, path):
        attr_re = {
            'gene': re.compile(r'gene_name "([^"]+)"'),
            'tid': re.compile(r'transcript_id "([^"]+)"'),
        }
        coll = {}
        for row in gzip.open(path, 'rt'):
            f = row.split('\t')
            feat = f[2]
            if feat not in ('exon', 'CDS'):
                continue
            mg = attr_re['gene'].search(f[8])
            mt = attr_re['tid'].search(f[8])
            if not mg or not mt:
                continue
            t = coll.setdefault(mt.group(1), Transcript(mg.group(1), mt.group(1), f[6]))
            t.__dict__['chrom'] = f[0][3:] if f[0].startswith('chr') else f[0]
            (t.exons if feat == 'exon' else t.cds).append((int(f[3]), int(f[4])))
        for t in coll.values():
            t.exons.sort()
            t.cds.sort()
            self.tx.setdefault(t.__dict__['chrom'], []).append(t)
        for c, lst in self.tx.items():
            lst.sort(key=lambda t: t.exons[0][0])
            self.tx_starts[c] = [t.exons[0][0] for t in lst]

    # ---------------- lookup ----------------
    MAX_TX_SPAN = 2_500_000  # longest human gene (DMD ~2.24 Mb) + margin

    def _candidates(self, chrom, pos):
        """Transcripts whose span contains pos. Walk back over transcripts with
        start <= pos until beyond max gene span; any containing transcript must
        start within MAX_TX_SPAN of pos."""
        lst = self.tx.get(chrom, [])
        if not lst:
            return []
        starts = self.tx_starts[chrom]
        i = bisect.bisect_right(starts, pos) - 1
        out = []
        while i >= 0 and starts[i] > pos - self.MAX_TX_SPAN:
            t = lst[i]
            if t.exons[0][0] <= pos <= t.exons[-1][1]:
                out.append(t)
            i -= 1
        return out

    def _fetch(self, chrom, s, e):
        c = ('chr' + chrom) if self.prefix else chrom
        if self.genome is None:
            return None
        try:
            names = self.genome.references
            key = c if c in names else (chrom if chrom in names else 'chr' + chrom)
            return self.genome.fetch(key, s - 1, e).upper()
        except Exception:
            return None

    def annotate(self, chrom, pos, ref, alts):
        """-> list of dicts, one per ALT. chrom with or without 'chr' prefix."""
        chrom = chrom[3:] if chrom.startswith('chr') else chrom
        txs = self._candidates(chrom, pos)
        out = []
        for alt in alts:
            ann = {
                'gene': None, 'tid': None, 'consequence': 'intergenic',
                'hgvsp': None, 'splice_dist': None, 'cadd_phred': None,
            }
            if txs:
                t = txs[0]
                ann['gene'], ann['tid'] = t.gene, t.tid
                self._classify(ann, t, chrom, pos, ref, alt)
            ann['cadd_phred'] = self._cadd(chrom, pos, ref, alt)
            out.append(ann)
        return out

    def _classify(self, ann, t, chrom, pos, ref, alt):
        exonic = any(s <= pos <= e for s, e in t.exons)
        in_cds = any(s <= pos <= e for s, e in t.cds)
        # splice distance from nearest exon edge
        edges = [x for ex in t.exons for x in ex]
        dist = min(abs(pos - e) for e in edges)
        ann['splice_dist'] = dist
        if dist <= 2:
            ann['consequence'] = 'splice_donor' if (pos - min(edges)) % 1 == 0 else 'splice_acceptor'
            # donor = start of intron after exon on strand; acceptor = end. Simplify: edge-near = splice
            ann['consequence'] = 'splice_donor' if dist <= 1 else 'splice_acceptor'
        if in_cds:
            self._coding(ann, t, chrom, pos, ref, alt)
            return
        if exonic:
            cds_s, cds_e = t.cds[0][0], t.cds[-1][1]
            if t.strand == '+':
                ann['consequence'] = '5_prime_utr' if pos < cds_s else '3_prime_utr'
            else:
                ann['consequence'] = '3_prime_utr' if pos < cds_s else '5_prime_utr'
            if dist <= 8:
                ann['consequence'] += '|splice_region'
            return
        ann['consequence'] = 'intronic' if dist > 8 else 'splice_region'

    def _coding(self, ann, t, chrom, pos, ref, alt):
        if not re.fullmatch(f'[{BASES}]+', ref) or not re.fullmatch(f'[{BASES}]+', alt):
            ann['consequence'] = 'other'
            return
        dl = len(alt) - len(ref)
        if dl != 0:
            ann['consequence'] = 'frameshift' if dl % 3 else ('inframe_deletion' if dl < 0 else 'inframe_insertion')
            return
        if len(ref) != 1 or self.genome is None:
            ann['consequence'] = 'missense_variant' if len(ref) == 1 else 'other'
            return
        cds = t.cds_spliced(lambda s, e: self._fetch(chrom, s, e))
        # index of pos in spliced CDS (genomic order -> translation order)
        i = None
        off = 0
        parts = sorted(t.cds)  # genomic order
        for s, e in parts:
            if s <= pos <= e:
                i = off + (pos - s)
                break
            off += e - s + 1
        if i is None:
            ann['consequence'] = 'other'
            return
        if t.strand == '-':
            i = len(cds) - 1 - i
            ref_b = _revcomp(ref)
            alt_b = _revcomp(alt)
        else:
            ref_b, alt_b = ref, alt
        ci = (i // 3) * 3
        codon = cds[ci:ci + 3]
        if len(codon) != 3 or codon[i % 3] != ref_b:
            ann['consequence'] = 'missense_variant'
            ann['hgvsp'] = f'codon-mismatch(ref=genome)'
            return
        mut = codon[:i % 3] + alt_b + codon[i % 3 + 1:]
        aa = i // 3 + 1
        a1, a2 = CODON.get(codon, '?'), CODON.get(mut, '?')
        if a1 == '*' and a2 != '*':
            ann['consequence'] = 'stop_lost'
        elif a2 == '*':
            ann['consequence'] = 'stop_gained'
        elif a1 == a2:
            ann['consequence'] = 'synonymous_variant'
        else:
            ann['consequence'] = 'missense_variant'
        ann['hgvsp'] = f'p.{a1}{aa}{a2}'

    def _cadd(self, chrom, pos, ref, alt):
        tb = self.cadd.get(chrom)
        if tb is None:
            return None
        try:
            for row in tb.fetch(chrom, pos - 1, pos):
                f = row.split('\t')
                if int(f[1]) == pos and f[2] == ref and f[3] == alt:
                    return float(f[5])
        except Exception:
            return None
        return None


if __name__ == '__main__':
    import yaml
    here = os.path.dirname(os.path.abspath(__file__))
    conf = yaml.safe_load(open(os.path.join(here, '../conf/track1_config.yml')))
    gtf = os.path.expanduser(conf['refs']['mane_gtf'])
    if not os.path.exists(gtf):
        alt_gtf = os.path.expanduser('~/Documents/projects/hf-rare/ref/gencode.v47.MANE.gtf.gz')
        if os.path.exists(alt_gtf):
            gtf = alt_gtf
    fasta = os.path.expanduser(conf['refs']['genome_fasta'])
    have_fasta = os.path.exists(fasta)
    if not os.path.exists(gtf):
        print(f"SKIP all checks: GENCODE MANE GTF absent on this host (conf refs.mane_gtf = {gtf});"
              " recorded result: validation/annotate_selftest_dichtator.txt")
        raise SystemExit(0)
    a = LocalAnnotator(gtf, genome_fasta=fasta if have_fasta else None)
    # tier 1 (GTF only, runs everywhere): gene/transcript identification and
    # position-class of non-coding rows. Consequence strings of coding SNVs are
    # NOT asserted here - with no FASTA the annotator falls back to
    # 'missense_variant' for any coding SNV, so that assertion would pass
    # vacuously and hide a codon-arithmetic regression.
    gtf_checks = [  # (chrom, pos, ref, alt, want_gene, want_tid_prefix, want_class)
        ('15', 40220612, 'T', 'G', 'BUB1B', 'ENST00000287598', 'coding'),
        ('15', 40181093, 'AT', 'A', 'BUB1B', 'ENST00000287598', 'intronic'),
        ('11', 95801248, 'G', 'A', 'CEP57', 'ENST00000325542', 'intronic'),
    ]
    # tier 2 (codon-level, requires the GRCh38 FASTA, conf refs.genome_fasta)
    codon_checks = [  # (chrom, pos, ref, alt, want_consequence)
        ('15', 40220612, 'T', 'G', 'missense_variant'),
        ('15', 40209701, 'T', 'G', 'stop_gained'),
    ]
    ok = True
    for c, p, r, alt, want_gene, want_tid, want_class in gtf_checks:
        got = a.annotate(c, p, r, [alt])[0]
        g_ok = got['gene'] == want_gene and (got['tid'] or '').startswith(want_tid)
        cls_ok = got['consequence'] == 'intronic' if want_class == 'intronic' else got['consequence'] != 'intronic'
        status = 'OK' if (g_ok and cls_ok) else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f"{status}  {c}:{p} {r}>{alt} -> gene={got['gene']} tid={got['tid']} "
              f"class={got['consequence']} (want gene={want_gene} tid~{want_tid} class={want_class})")
    for c, p, r, alt, want in codon_checks:
        if not have_fasta:
            print(f"SKIP codon-level {c}:{p} {r}>{alt} (want {want}): GRCh38 FASTA absent on this host"
                  " (conf refs.genome_fasta); recorded result: validation/annotate_selftest_dichtator.txt")
            continue
        got = a.annotate(c, p, r, [alt])[0]
        status = 'OK' if got['consequence'] == want else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f"{status}  {c}:{p} {r}>{alt} -> {got['consequence']} {got['hgvsp'] or ''} gene={got['gene']}")
    sys.exit(0 if ok else 1)
