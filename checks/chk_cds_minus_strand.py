#!/usr/bin/env python3
"""Regression check: multi-exon minus-strand CDS assembly (fix/minus-strand-cds-assembly).

Builds a tiny synthetic MANE-shaped GTF + FASTA with
  - one minus-strand gene with THREE CDS segments whose true spliced CDS is
    ATG (GCT)^n ... TAA  (start codon + readable frame + stop codon), and
  - one plus-strand single-exon control.
Asserts the annotator (pipeline/scripts/annotate.py):
  1. cds_spliced starts with ATG for both genes (translation-order assembly),
  2. the minus-strand spliced CDS ends with a stop codon,
  3. annotate() of a variant designed to produce a premature stop on the
     minus-strand transcript returns consequence == 'stop_gained' (not the
     'codon-mismatch(ref=genome)' fallback),
  4. a plus-strand synonymous-designed variant stays 'synonymous_variant'.

Exits 0 when all assertions pass, 1 with a diff otherwise.
Deterministic, self-contained, no clinical data.
"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pipeline', 'scripts'))
import pysam  # noqa: E402
from annotate import LocalAnnotator  # noqa: E402

COMP = str.maketrans('ACGT', 'TGCA')

def revcomp(s):
    return s.translate(COMP)[::-1]

# transcript CDS (translation order) fragments for the minus gene:
TFRAGS = ['ATGGCTGCC', 'TCCAAGCAT', 'CGTTAA']  # 21 nt -> 7 codons, ends TAA
assert len(''.join(TFRAGS)) % 3 == 0
# genomic layout (ascending): last transcript fragment first, each revcomp'd
GSEGS = [(2001, 2001 + len(TFRAGS[2]) - 1, revcomp(TFRAGS[2])),   # 2001-2006
         (3001, 3001 + len(TFRAGS[1]) - 1, revcomp(TFRAGS[1])),   # 3001-3009
         (4001, 4001 + len(TFRAGS[0]) - 1, revcomp(TFRAGS[0]))]   # 4001-4009
# plus single-exon control: CDS ATG AAA TAA (9 nt) at 5001-5009, genome provides
# exactly that plus-strand sequence.
PSEG = (5001, 5009, 'ATGAAATAA')

CHROM_LEN = 6000
genome = ['N'] * CHROM_LEN
for s, e, frag in GSEGS + [PSEG]:
    genome[s - 1:e] = frag
GENOME = ''.join(genome)

GTF = ('1\tHAVANA\tCDS\t{s}\t{e}\t.\t-\t0\tgene_name "GMINUS"; transcript_id "TXM1";\n'
       .format)

def write_inputs(d):
    fa = os.path.join(d, 'genome.fa')
    with open(fa, 'w') as fh:
        fh.write('>1\n')
        for i in range(0, CHROM_LEN, 60):
            fh.write(GENOME[i:i + 60] + '\n')
    pysam.faidx(fa)
    gtf = os.path.join(d, 'mane.gtf.gz')
    rows = []
    for s, e, _ in GSEGS:
        rows.append(f'1\tHAVANA\tCDS\t{s}\t{e}\t.\t-\t0\tgene_name "GMINUS"; transcript_id "TXM1";')
    rows.append(f'1\tHAVANA\texon\t1991\t2010\t.\t-\t.\tgene_name "GMINUS"; transcript_id "TXM1";')
    rows.append(f'1\tHAVANA\texon\t2991\t3010\t.\t-\t.\tgene_name "GMINUS"; transcript_id "TXM1";')
    rows.append(f'1\tHAVANA\texon\t3991\t4010\t.\t-\t.\tgene_name "GMINUS"; transcript_id "TXM1";')
    rows.append(f'1\tHAVANA\tCDS\t5001\t5009\t.\t+\t0\tgene_name "GPLUS"; transcript_id "TXP1";')
    rows.append(f'1\tHAVANA\texon\t4991\t5010\t.\t+\t.\tgene_name "GPLUS"; transcript_id "TXP1";')
    rows.sort(key=lambda r: (int(r.split('\t')[3]), int(r.split('\t')[4])))
    with __import__('gzip').open(gtf, 'wt') as fh:
        fh.write('\n'.join(rows) + '\n')
    return fa, gtf

def main():
    with tempfile.TemporaryDirectory() as d:
        fa, gtf = write_inputs(d)
        ann = LocalAnnotator(gtf, genome_fasta=fa)
        t_minus = next(t for c in ann.tx.values() for t in c if t.gene == 'GMINUS')
        t_plus = next(t for c in ann.tx.values() for t in c if t.gene == 'GPLUS')
        fetch = lambda s, e: ann._fetch('1', s, e)
        cds_m = t_minus.cds_spliced(fetch)
        cds_p = t_plus.cds_spliced(fetch)
        checks = [
            ('minus CDS starts ATG', cds_m[:3] == 'ATG'),
            ('minus CDS equals design', cds_m == ''.join(TFRAGS)),
            ('minus CDS ends stop', cds_m[-3:] in ('TAA', 'TAG', 'TGA')),
            ('plus CDS equals design', cds_p == 'ATGAAATAA'),
        ]
        # TFRAGS[1]='TCCAAGCAT' -> codon 4 'TCC'... design a premature stop:
        # transcript codon positions: ATG(1) GCT(2) GCC(3) TCC(4) AAG(5) CAT(6) CGT(7) TAA(8)
        # Mutating transcript CAG? pick codon 5 'AAG'(Lys) -> 'TAG'(stop): position 2 of codon 5, A->T on transcript
        # codon 5 starts at nt index 4*3=12 (0-based) -> within TFRAGS[1] idx 3 (12-9=3), transcript base A
        # TFRAGS[1][3] = 'A' (AAG), mutate idx 3 A->T gives TAG stop? TCC TAG CAT -> need codon 'TAG': yes A->T.
        # transcript index 12 corresponds to genome pos: reverse layout; find pos via CDS scan.
        i = 12
        remaining = i
        tx_pos = None
        for s, e, _frag in [(2001, 2006, None)] + []:
            pass
        # transcript order segments are the *reversed* genomic ascending list:
        for s, e in reversed(sorted([(g[0], g[1]) for g in GSEGS])):
            ln = e - s + 1
            if remaining < ln:
                tx_pos = e - remaining  # minus strand: transcript 5'->3' runs descending in genome
                break
            remaining -= ln
        assert tx_pos is not None
        # transcript base = A (on transcription strand) -> genomic REF = revcomp('A') = 'T'
        # transcript ALT 'T' -> genomic ALT = revcomp('T') = 'A'
        out = ann.annotate('1', tx_pos, 'T', 'A')[0]
        checks.append(('stop design -> stop_gained', out['consequence'] == 'stop_gained'))
        checks.append(('no codon-mismatch trap', out['hgvsp'] != 'codon-mismatch(ref=genome)'))
        checks.append(('minus gene name', out['gene'] == 'GMINUS'))
        # plus synonymous: ATG AAA TAA -> aa2 AAA(K)->AAG(K) synonymous via 3rd base of codon 2
        outsyn = ann.annotate('1', 5006, 'A', 'G')[0]  # genome 5006 = PSEG idx 5 = codon2 base 3
        checks.append(('plus synonymous stays synonymous', outsyn['consequence'] == 'synonymous_variant'))
        bad = [name for name, ok in checks if not ok]
        for name, ok in checks:
            print(f'{"PASS" if ok else "FAIL"}  {name}')
        if bad:
            print('CHK_CDS_MINUS: FAIL', bad, file=sys.stderr)
            return 1
        print('CHK_CDS_MINUS_OK')
        return 0

if __name__ == '__main__':
    sys.exit(main())
