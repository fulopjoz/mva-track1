"""T1-02: genotype-aware, allele-specific multiallelic processing.

Regression guard for the audit findings:
  * alternate fraction was always computed from AD[1] regardless of GT;
  * every ALT in a record was processed downstream, even ALTs absent from
    the sample genotype.
Design now: iterate the sample genotype's ALT allele indices only; fraction
is that allele's AD / sum(AD); zygosity comes from GT (het vs hom).
Cases: 0/1, 0/2, 1/2, 2/2, phased 0|1, partial missing, fully missing,
malformed/missing AD or DP. All coordinates invented; uses the REAL modules.
"""
from __future__ import annotations

import sys

import cyvcf2
import pytest

from conftest import SCRIPT_DIR, VCF_HEADER, write_conf, write_gtf, write_vcf

sys.path.insert(0, str(SCRIPT_DIR))
import prioritize  # noqa: E402

QC = {'min_gq': 20, 'min_dp': 10, 'min_alt_frac': 0.25}


def iter_cases(tmp_path, rows):
    """(chrom,pos,gt) pairs yielded by the real QC iterator for the given rows."""
    vcf = write_vcf(tmp_path, rows)
    return [(rec.CHROM, rec.POS, k, zyg) for rec, zyg, k in prioritize.iter_qc_variants(str(vcf), QC)]


def r(pos, ref='A', alt='G', gt='0/1', ad='20,20', dp=40, gq=99, fmt='GT:GQ:DP:AD'):
    return f'1\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\t{fmt}\t{gt}:{gq}:{dp}:{ad}'


# ---------- audit's two headline failures ----------


def test_second_alt_fraction_uses_own_ad(tmp_path):
    # GT 0/2, AD 10,0,10: called ALT2 fraction is 0.5 -> keep. Old code used AD[1]=0 -> drop.
    got = iter_cases(tmp_path, [r(500, alt='C,G', gt='0/2', ad='10,0,10')])
    assert got == [('1', 500, 2, 'het')]


def test_only_called_alt_is_processed(tmp_path):
    # GT 0/1 with ALT=C,G: downstream must see ALT1 only, never ALT2.
    got = iter_cases(tmp_path, [r(500, alt='C,G', gt='0/1', ad='10,10,0')])
    assert got == [('1', 500, 1, 'het')]


# ---------- genotype matrix ----------

@pytest.mark.parametrize('gt,ad,expected', [
    ('0/1', '20,20', [('1', 500, 1, 'het')]),
    ('0/2', '10,0,10', [('1', 500, 2, 'het')]),
    ('1/2', '10,10,10', [('1', 500, 1, 'het'), ('1', 500, 2, 'het')]),
    ('2/2', '10,0,25', [('1', 500, 2, 'hom_alt')]),
    ('1/1', '8,30', [('1', 500, 1, 'hom_alt')]),
    ('0|1', '20,20', [('1', 500, 1, 'het')]),      # phased het
    ('1|2', '10,10,10', [('1', 500, 1, 'het'), ('1', 500, 2, 'het')]),
])
def test_genotype_matrix(tmp_path, gt, ad, expected):
    got = iter_cases(tmp_path, [r(500, alt='C,G', gt=gt, ad=ad)])
    assert got == expected


@pytest.mark.parametrize('gt', ['./.', './1', '1/.', '.|.'])
def test_missing_genotypes_skipped(tmp_path, gt):
    # half- or fully-missing calls cannot establish zygosity: skip (documented).
    got = iter_cases(tmp_path, [r(500, alt='C,G', gt=gt, ad='10,10,5')])
    assert got == []


def test_allele_fraction_is_called_alleles_fraction(tmp_path):
    # Called ALT1's fraction is 5/40 = 0.125 -> drop, even though ALT2 would pass.
    assert iter_cases(tmp_path, [r(500, alt='C,G', gt='0/1', ad='20,5,15')]) == []
    # Same depths with GT 0/2: called ALT2's fraction is 15/40 = 0.375 -> keep ALT2 only.
    got = iter_cases(tmp_path, [r(500, alt='C,G', gt='0/2', ad='20,5,15')])
    assert got == [('1', 500, 2, 'het')]


@pytest.mark.parametrize('fmt,gq,dp,ad', [
    ('GT', 99, 40, ''),          # no FORMAT AD/DP/GQ columns at all
    ('GT:AD', 99, 40, '20,20'),  # DP missing from FORMAT -> dp=0 < min_dp
])
def test_malformed_or_missing_dp_rejected(tmp_path, fmt, gq, dp, ad):
    val = f'0/1' + (f':{gq}' if 'GQ' in fmt else '') + (f':{dp}' if 'DP' in fmt else '') + (f':{ad}' if 'AD' in fmt else '')
    row = f'1\t500\t.\tA\tG\t99\tPASS\t.\t{fmt}\t{val}'
    assert iter_cases(tmp_path, [row]) == []


def test_missing_ad_rejected(tmp_path):
    # DP present, AD field declared but value '.' -> cannot verify fraction -> drop.
    row = '1\t500\t.\tA\tG\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:.'
    assert iter_cases(tmp_path, [row]) == []


# ---------- integration: downstream sees only genotype-decomposed ALT ----------

def run_pipeline(tmp_path, rows):
    gtf = write_gtf(tmp_path)
    vcf = write_vcf(tmp_path, rows)
    out = tmp_path / 'out.tsv'
    conf = write_conf(tmp_path, gtf, vcf)
    sys.argv = ['prioritize.py', '--conf', str(conf), '--out', str(out)]
    prioritize.main()
    return out.read_text().strip().split('\n')


def test_downstream_only_called_alt_reaches_table(tmp_path):
    """Multiallelic record in GENEP at coding pos 150, ALT=C,G called 0/1.
    Only the C row may ever reach the output (the G allele was never in the
    patient's genotype). Pre-fix both were annotated and could seed pairs."""
    # R1: multiallelic called ALT C only -> 1 het allele in GENEP
    # R2: second variant (donor) -> pair possible only with ALT C
    lines = run_pipeline(tmp_path, [
        '1\t150\t.\tA\tC,G\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20,0',
        '1\t200\t.\tA\tG\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20',
    ])
    body = [ln for ln in lines[1:]]
    assert len(body) == 1, f'expected exactly one pair row, got {body}'
    cols = body[0].split('\t')
    members = {(cols[4], cols[6]), (cols[8], cols[10])}  # {(pos, alt)}
    assert ('150', 'G') not in members and ('150', 'C') in members


def test_uncalled_alt_cannot_seed_pair(tmp_path):
    """A single multiallelic splice-donor record, GT 0/1: only ALT C is in the
    genotype, so the gene holds exactly ONE het donor allele and no pair can
    form. Pre-fix both ALTs were processed and the record paired with itself."""
    lines = run_pipeline(tmp_path, [
        '1\t200\t.\tA\tC,G\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20,0',
    ])
    assert len(lines) == 1, f'expected header only (no pair), got {lines}'
