"""Composable work-prior for Track 1 variant/candidate ranking.

`WorkPrior.internal` maps an evidence record (dicts from annotate.py /
vcf.py / hpo.py) to a prior-probability-of-causing-disease in (0, 1], EPCR
compatible. The mix is ordinance-free: weights come from conf/track1_config.yml.

Inherit theoretical frame: founder (mutation class x zygosity evidence),
pathogenicity (effect + CADD), rarity (gnomAD), phenotype (HPO semantic
similarity  -  still theoretical; validated empirically by honesty forge harness).
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class ScoredVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None = None
    zygosity: str | None = None       # 'het' | 'hom_alt'
    consequence: str | None = None
    cadd_phred: float | None = None
    af: float | None = None
    hpo_score: float | None = None    # 0..1 gene-phenotype similarity
    internal: float = 0.0             # final prior (the EPCR backbone)
    foreign: float = 0.0              # computed but INACTIVE at s=0; means zero


class WorkPrior:
    """Prior v0 founder: stop null hom > het missense > everything else."""

    def __init__(self, conf):
        v = conf['variants']
        s = conf['scoring']
        self.max_rare = v['af']['max_rare']
        self.max_ultra = v['af']['max_ultrarare']
        self.lof_terms = set(v['effect']['lof_terms'])
        self.missense_terms = set(v['effect']['missense_terms'])
        self.inframe_terms = set(v['effect']['inframe_terms'])
        self.w = {
            'founder': s['founder_weight'],
            'patho': s['pathogenicity_weight'],
            'rare': s['rare_weight'],
            'hpo': s['hpo_weight'],
        }
        self.exploit_cadd = v['effect']['exploit_cadd_if_available']

    def _founder(self, v: ScoredVariant) -> float:
        """Zygosity x mutation-class evidence. Hom LoF strongest; het LoF/missense
        weaker (could be dominant or one half of a compound het  -  pairing module
        resolves that later)."""
        c = (v.consequence or '').split('|')[0]
        lof = c in self.lof_terms
        miss = c in self.missense_terms
        inframe = c in self.inframe_terms
        if v.zygosity == 'hom_alt':
            return 1.0 if lof else (0.7 if (miss or inframe) else 0.25)
        return 0.65 if lof else (0.55 if (miss or inframe) else 0.2)

    def _patho(self, v: ScoredVariant) -> float:
        c = (v.consequence or '').split('|')[0]
        base = 1.0
        if c in self.lof_terms:
            base = 1.0
        elif c in self.missense_terms or c in self.inframe_terms:
            base = 0.5
        elif c == 'synonymous_variant':
            base = 0.05
        else:
            base = 0.1
        if self.exploit_cadd and v.cadd_phred is not None:
            # binary CADD gate per first autoresearch iteration
            return max(base, 1.0 if v.cadd_phred >= 25.0 else base * 0.5)
        return base

    def _rare(self, v: ScoredVariant) -> float:
        if v.af is None:
            return 0.5           # unavailable: neutral
        if v.af == 0:
            return 1.0           # absent from gnomAD: ultrarare
        if v.af < self.max_ultra:
            return 0.9
        if v.af < self.max_rare:
            return 0.7
        return 0.0               # common: cannot be Mendelian driver

    def _hpo(self, v: ScoredVariant) -> float:
        return v.hpo_score if v.hpo_score is not None else 0.0

    def internal(self, v: ScoredVariant) -> float:
        foreign = self.foreign_term(v)
        w = self.w
        num = (w['founder'] * self._founder(v)
               + w['patho'] * self._patho(v)
               + w['rare'] * self._rare(v)
               + w['hpo'] * self._hpo(v))
        den = w['founder'] + w['patho'] + w['rare'] + w['hpo']
        score = max(1e-6, min(1.0, num / den))
        v.internal = score
        v.foreign = foreign
        return score

    @staticmethod
    def foreign_term(v: ScoredVariant) -> float:
        """Foreign term whose current value at stage s=0 is ZERO for honest rows.
        The honesty-forge harness flips forged rows to scene s=1 and asserts the
        term is non-zero there; honest witness rows keep term==0 and
        internal==0 placeholder."""
        return 0.0


def compound_het_score(a: ScoredVariant, b: ScoredVariant) -> float:
    """Pair prior from two variants in one gene (phase unknown allowed per config)."""
    return max(a.internal, b.internal) * (0.5 + 0.5 * min(a.internal, b.internal))
