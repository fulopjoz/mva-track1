"""T1-01: canonical splice donor/acceptor classification (strand-aware, no overwrite).

Regression guard for the audit finding: intronic +1/+2 (donor) and -2/-1
(acceptor, 3' side in transcription order) were labeled `splice_region` by
console overwrite, so genuine LoF splice alleles could not seed compound-het
pairs. Every case runs against the real LocalAnnotator with a synthesized GTF.
Cases and expectations follow Ensembl/VEP geometry conventions.
All coordinates are invented.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from conftest import SCRIPT_DIR, write_conf, write_gtf, write_vcf

sys.path.insert(0, str(SCRIPT_DIR))
from annotate import LocalAnnotator  # noqa: E402


def annotator(tmp_path):
    gtf = write_gtf(tmp_path)
    return LocalAnnotator(str(gtf))   # no FASTA: coding SNVs fall back, introns unaffected


# ---------- intronic canonical sites, both strands ----------

@pytest.mark.parametrize('chrom,pos,want', [
    # + strand: donor = first two intronic bases after the exon end
    ('1', 200, 'splice_donor'),
    ('1', 201, 'splice_donor'),
    # + strand: acceptor = last two intronic bases before the next exon
    ('1', 299, 'splice_acceptor'),
    ('1', 298, 'splice_acceptor'),
    # - strand: transcription runs high->low; donor sits at the HIGH end
    ('2', 299, 'splice_donor'),
    ('2', 298, 'splice_donor'),
    # - strand: acceptor = LOW end of the intron
    ('2', 200, 'splice_acceptor'),
    ('2', 201, 'splice_acceptor'),
])
def test_intronic_canonical_sites(tmp_path, chrom, pos, want):
    a = annotator(tmp_path)
    got = a.annotate(chrom, pos, 'A', ['G'])[0]['consequence']
    assert got == want, f'pos {chrom}:{pos} -> {got}, want {want}'


@pytest.mark.parametrize('chrom,pos,want', [
    ('1', 204, 'splice_region'),   # noncanonical near edge (3..8 nt)
    ('1', 296, 'splice_region'),
    ('1', 210, 'intronic'),        # deeper than the region band
    ('2', 204, 'splice_region'),
    ('2', 296, 'splice_region'),
    ('2', 290, 'intronic'),
])
def test_intronic_noncanonical_and_deep(tmp_path, chrom, pos, want):
    a = annotator(tmp_path)
    assert a.annotate(chrom, pos, 'A', ['G'])[0]['consequence'] == want


def test_exonic_boundary_variants_keep_region_tag_not_canonical(tmp_path):
    a = annotator(tmp_path)
    # exonic, UTR, within the 8nt band -> UTR with splice_region tag
    got = a.annotate('1', 105, 'A', ['G'])[0]
    assert got['consequence'] == '5_prime_utr|splice_region'
    # exonic boundary inside CDS (dist<=2): must NOT be called splice_donor;
    # coding branch owns it (no FASTA -> generic coding fallback).
    got = a.annotate('1', 111, 'A', ['G'])[0]
    assert not got['consequence'].startswith('splice')
    assert got['consequence'] in ('missense_variant', 'other')
    # exonic edge on the minus strand, UTR band
    got = a.annotate('2', 395, 'A', ['G'])[0]
    assert got['consequence'] == '5_prime_utr|splice_region'


def test_single_exon_no_intronic_calls(tmp_path):
    a = annotator(tmp_path)
    assert a.annotate('3', 101, 'A', ['G'])[0]['consequence'] == '5_prime_utr|splice_region'
    assert a.annotate('3', 99, 'A', ['G'])[0]['consequence'] == 'intergenic'
    assert a.annotate('3', 300, 'A', ['G'])[0]['consequence'] == 'intergenic'


def test_splice_donor_seeds_compound_het_pair(tmp_path):
    """End-to-end: an intronic donor + a coding missense in the same gene,
    both heterozygous, must yield one candidate pair. Pre-fix, the donor
    degraded to splice_region and no pair was emitted."""
    gtf = write_gtf(tmp_path)
    rows = [
        '1\t200\t.\tA\tG\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20',
        '1\t150\t.\tA\tG\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20',
    ]
    vcf = write_vcf(tmp_path, rows)
    out = tmp_path / 'out.tsv'
    conf = write_conf(tmp_path, gtf, vcf)
    from prioritize import main as prio_main  # noqa: E402
    sys.argv = ['prioritize.py', '--conf', str(conf), '--out', str(out)]
    prio_main()
    lines = out.read_text().strip().split('\n')
    assert len(lines) == 2, f'expected header+1 pair, got: {lines}'
    cols = lines[1].split('\t')
    assert cols[1] == 'compound_het' and cols[2] == 'GENEP'
    assert cols[4] == '200' or cols[8] == '200'
    assert 'splice_donor' in (cols[12], cols[13])
