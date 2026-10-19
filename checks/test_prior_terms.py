"""Scoring-term contract locks (audit T1-03): founder/patho/rare/hpo component
semantics and the internal normalisation, including neutral-when-missing AF
and CADD. Locks decision-critical values before the clinical rerun.
"""
from __future__ import annotations

import sys

import pytest

from conftest import SCRIPT_DIR, write_conf

sys.path.insert(0, str(SCRIPT_DIR))
from prior import ScoredVariant, WorkPrior, compound_het_score  # noqa: E402


def _prior(tmp_path):
    gtf = tmp_path / 'mini.gtf.gz'
    # reuse conftest conf builder for an authoritative conf
    from conftest import write_conf as wc, write_gtf
    gtf = write_gtf(tmp_path)
    vcf = tmp_path / 'fixture.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    conf = wc(tmp_path, gtf, vcf)
    import yaml
    return WorkPrior(yaml.safe_load(open(conf)))


def test_founder_term(tmp_path):
    p = _prior(tmp_path)
    v = lambda c, z: ScoredVariant('1', 1, 'A', 'G', gene='G', zygosity=z, consequence=c)
    assert p._founder(v('stop_gained', 'hom_alt')) == 1.0
    assert p._founder(v('stop_gained', 'het')) == 0.65
    assert p._founder(v('missense_variant', 'het')) == 0.55
    assert p._founder(v('synonymous_variant', 'het')) == 0.2
    # splice_donor is LoF for the founder term (audit T1-01 downstream contract)
    assert p._founder(v('splice_donor', 'het')) == 0.65


def test_patho_term_and_cadd_gate(tmp_path):
    p = _prior(tmp_path)
    mk = lambda c, cad=None: ScoredVariant('1', 1, 'A', 'G', consequence=c, cadd_phred=cad)
    assert p._patho(mk('stop_gained')) == 1.0
    assert p._patho(mk('missense_variant')) == 0.5
    assert p._patho(mk('scroll')) == 0.1
    # CADD >= 25 lifts to 1.0; below the threshold it must NOT lower the base
    assert p._patho(mk('missense_variant', 30.0)) == 1.0
    assert p._patho(mk('missense_variant', 5.0)) == 0.5


def test_rare_term_neutral_when_missing(tmp_path):
    p = _prior(tmp_path)
    mk = lambda af: ScoredVariant('1', 1, 'A', 'G', af=af)
    assert p._rare(mk(None)) == 0.5      # unavailable: neutral, never ultrarare
    assert p._rare(mk(0.0)) == 1.0       # absent: ultrarare
    assert p._rare(mk(0.00005)) == 0.9   # < max_ultrarare
    assert p._rare(mk(0.005)) == 0.7     # < max_rare
    assert p._rare(mk(0.5)) == 0.0       # common


def test_internal_normalisation_and_symmetric_pair(tmp_path):
    p = _prior(tmp_path)
    a = ScoredVariant('1', 1, 'A', 'G', consequence='splice_donor', zygosity='het')
    b = ScoredVariant('1', 2, 'A', 'G', consequence='missense_variant', zygosity='het')
    ia, ib = p.internal(a), p.internal(b)
    # weights founder .45 patho .60 rare .30 hpo 2.0 with no refs -> rare 0.5, hpo 0
    assert ia == pytest.approx((0.45 * 0.65 + 0.60 * 1.0 + 0.30 * 0.5) / 3.35)
    assert ib == pytest.approx((0.45 * 0.55 + 0.60 * 0.5 + 0.30 * 0.5) / 3.35)
    assert compound_het_score(a, b) == compound_het_score(b, a)  # symmetric: orientation-free dedupe basis
