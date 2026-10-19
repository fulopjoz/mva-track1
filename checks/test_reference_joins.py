"""Reference-join verification: gnomAD AF merge, ClinVar tabix lookup,
CADD line parsing, HPO ancestor scoring (audit T1-03 verification sweep).

Tiny synthetic reference files in the exact on-disk formats the joins expect
(plain VCF for gnomAD stream readers; bgzipped+tabix-indexed for tabix hits).
All values invented. These are characterization tests: a FAIL indicates a real
join defect, before any clinical-scale rerun is trusted.
"""
from __future__ import annotations

import sys

import pysam
import pytest

from conftest import SCRIPT_DIR, write_gtf

sys.path.insert(0, str(SCRIPT_DIR))
from _hpo import HPOScorer  # noqa: E402
from annotate import LocalAnnotator  # noqa: E402
from vcf import ClinVarJoin, GnomadAfJoin  # noqa: E402

VHC = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID={chrom}>\n'
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Alelle frequency">\n'
    '##INFO=<ID=CLNSIG,Number=A,Type=String,Description="significance">\n'
    '##INFO=<ID=CLNDN,Number=A,Type=String,Description="disease">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')


def _bgz(ltmp, name, lines, chrom='1'):
    p = ltmp / name
    p.write_text(VHC.format(chrom=chrom) + '\n'.join(lines) + '\n')
    return str(p)


def _tabixed(ltmp, name, lines, chrom='1'):
    raw = _bgz(ltmp, name, lines, chrom)
    import shutil
    traw = raw + '.raw'
    shutil.move(raw, traw)          # tabix_compress must not write onto its own input
    pysam.tabix_compress(traw, raw, force=True)
    pysam.tabix_index(raw, force=True, seq_col=0, start_col=1, end_col=1)
    return raw


# ---------------- gnomAD AF join ----------------

def _gnomad_dir(tmp_path):
    lines = [
        '1\t100\t.\tA\tC\t.\tPASS\tAF=0.01',
        '1\t100\t.\tA\tG\t.\tPASS\tAF=0.3',     # split line at the same position
        '1\t200\t.\tT\tG\t.\tPASS\tAF=0.00005',  # ultrarare
    ]
    _bgz(tmp_path, 'gnomad.exomes.v4.1.sites.chr1.vcf.bgz', lines)
    return GnomadAfJoin(str(tmp_path))


def test_gnomad_exact_allele_match_and_position_merge(tmp_path):
    g = _gnomad_dir(tmp_path)
    # both split alleles at one position are merged into one position map
    assert g.af('1', 100, 'A', 'C') == pytest.approx(0.01)
    assert g.af('1', 100, 'A', 'G') == pytest.approx(0.3)
    # exact REF/ALT keying: wrong ALT reverses to absent
    assert g.af('1', 100, 'A', 'T') == pytest.approx(0.0)
    assert g.af('1', 200, 'T', 'G') == pytest.approx(0.00005)
    # position scanned past without a match -> absent (0.0), NOT None
    assert g.af('1', 250, 'A', 'C') == pytest.approx(0.0)
    # missing chromosome file -> None (neutral), never falsely ultrarare
    assert g.af('2', 100, 'A', 'C') is None


def test_gnomad_out_of_order_detection_keeps_exact_answers(tmp_path):
    g = _gnomad_dir(tmp_path)
    assert g.af('1', 200, 'T', 'G') == pytest.approx(0.00005)
    assert g.af('1', 100, 'A', 'C') == pytest.approx(0.01)   # backward query
    assert g.out_of_order_queries() >= 1


# ---------------- ClinVar join ----------------

def test_clinvar_exact_lookup(tmp_path):
    lines = [
        '1\t150\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;CLNDN=TestDiseaseA',
        '1\t160\t.\tA\tC,G\t.\tPASS\tCLNSIG=Benign;CLNDN=TestDiseaseB',
    ]
    cb = _tabixed(tmp_path, 'clinvar.vcf.gz', lines)
    j = ClinVarJoin(cb)
    assert j.lookup('1', 150, 'A', 'G') == ('Pathogenic', 'TestDiseaseA')
    # multiallelic clinvar record: exact ALT match only
    assert j.lookup('1', 160, 'A', 'G') == ('Benign', 'TestDiseaseB')
    assert j.lookup('1', 160, 'A', 'T') == (None, None)
    # different position / missing region
    assert j.lookup('1', 500, 'A', 'G') == (None, None)


# ---------------- CADD extraction ----------------

def test_cadd_phred_field_extraction(tmp_path):
    gtf = write_gtf(tmp_path)
    lines = ['1\t500\tA\tG\t3.2\t25.5', '1\t500\tA\tC\t0.1\t5.0']
    cb = _tabixed(tmp_path, 'chr1.bgz', lines)
    a = LocalAnnotator(str(gtf), cadd_dir=str(tmp_path))
    # CADD layout: pos, ref, alt, raw, phred (indices 2,3,4,5) -> 25.5 exact key
    ann = a.annotate('1', 500, 'A', ['G'])[0]
    assert ann['cadd_phred'] == pytest.approx(25.5)
    ann = a.annotate('1', 500, 'A', ['T'])[0]
    assert ann['cadd_phred'] is None
    ann = a.annotate('1', 900, 'A', ['G'])[0]
    assert ann['cadd_phred'] is None
    # chromosome without a tile: neutral None
    ann = a.annotate('2', 105, 'A', ['G'])[0]
    assert ann['cadd_phred'] is None


# ---------------- HPO ancestor semantics ----------------

def _hpo_files(tmp_path):
    obo = tmp_path / 'hp.obo'
    obo.write_text(
        'format-version: 1.2\n\n'
        '[Term]\nid: HP:0001\nname: Root\n\n'
        '[Term]\nid: HP:0002\nname: Parent\nis_a: HP:0001\n\n'
        '[Term]\nid: HP:0003\nname: Exact\nis_a: HP:0002\n\n'
        '[Term]\nid: HP:0004\nname: Other\nis_a: HP:0005\n\n'
        '[Term]\nid: HP:0005\nname: OtherBranch\nis_a: HP:0001\n')
    gmap = tmp_path / 'map.tsv'
    gmap.write_text('GENEA\tHP:0003\tExact\nGENEB\tHP:0002\tParent\n')
    return str(obo), str(gmap)


def test_hpo_ancestor_closure_and_parent_level(tmp_path):
    obo, gmap = _hpo_files(tmp_path)
    s = HPOScorer(obo, gmap, ['HP:0003'])
    assert s.gene_score('GENEA') == pytest.approx(1.0)      # same term
    assert s.gene_score('GENEB') == pytest.approx(0.5)      # direct-parent-level ancestor
    assert s.gene_score('GENEC') == pytest.approx(0.0)      # unmapped
    # multiple proband terms average per gene
    s2 = HPOScorer(obo, gmap, ['HP:0003', 'HP:0004'])
    assert s2.gene_score('GENEA') == pytest.approx((1.0 + 0.0) / 2)
