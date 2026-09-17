"""HPO phenotype similarity – reference-only, zero patient data.

Reads at runtime (files may appear after code is written):
  refs.hp_obo   : hp.obo  – Human Phenotype Ontology (HPO) is_a hierarchy for ancestor closure
  refs.hpo_map  : hpo_gene_to_phenotype.tsv  – gene \t hpo_id \t hpo_name (built by fetch script)

Version 0 scorer: for each proband term, credit 1.0 if the gene's ancestor‑closed term set
contains it, 0.5 if it contains a direct parent (one generalisation step), otherwise 0.
Gene score is the mean over proband terms. This simple overlap method is implemented now;
information‑content weighting is planned for a future calibration pass.
"""
from __future__ import annotations

import os
import re
from functools import lru_cache


class HPOScorer:
    def __init__(self, hp_obo, hpo_map, proband_terms):
        self.proband = [t for t in proband_terms]
        self.parents = {}      # term -> set(direct is_a parents)
        self.gene_terms = {}   # gene -> set(hpo ids)
        if os.path.exists(hp_obo):
            self._parse_obo(hp_obo)
        if os.path.exists(hpo_map):
            self._parse_map(hpo_map)

    def _parse_obo(self, path):
        cur_id = None
        for line in open(path):
            line = line.strip()
            if line == '[Term]':
                cur_id = None
            elif line.startswith('id: HP:'):
                cur_id = line[4:]
                self.parents.setdefault(cur_id, set())
            elif line.startswith('is_a: HP:') and cur_id:
                self.parents[cur_id].add(line.split()[1])

    def _parse_map(self, path):
        for line in open(path):
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 2:
                continue
            self.gene_terms.setdefault(f[0], set()).add(f[1])

    @lru_cache(maxsize=200_000)
    def closure(self, term):
        """Term + all ancestors (is_A transitive)."""
        seen = {term}
        stack = [term]
        while stack:
            t = stack.pop()
            for p in self.parents.get(t, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        return frozenset(seen)

    def gene_score(self, gene: str) -> float:
        if not gene or gene not in self.gene_terms or not self.proband:
            return 0.0
        gset = set()
        for t in self.gene_terms[gene]:
            gset |= self.closure(t)
        tot = 0.0
        for pt in self.proband:
            if pt in gset:
                tot += 1.0
            elif self.parents.get(pt) and self.parents[pt] & gset:
                tot += 0.5
        return tot / len(self.proband)


if __name__ == '__main__':
    import sys
    import yaml
    conf = yaml.safe_load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../conf/track1_config.yml')))
    hp_obo = os.path.expanduser(conf['refs']['hp_obo'])
    hpo_map = os.path.expanduser(conf['refs']['hpo_map'])
    if not (os.path.exists(hp_obo) and os.path.exists(hpo_map)):
        print(f"SKIP hpo scoring: refs absent on this host ({hp_obo}, {hpo_map});"
              " recorded result: validation/hpo_selftest_dichtator.txt")
        sys.exit(0)
    terms = conf['phenotype']['hpo_terms']
    s = HPOScorer(hp_obo, hpo_map, terms)
    ok = True
    # positive controls: proband-phenotype genes must score > 0
    for g in ['BUB1B', 'CEP57']:
        sc = s.gene_score(g)
        print(g, f'{sc:.3f}')
        if sc <= 0.0:
            ok = False
    # negative control: unrelated gene must score 0
    sc = s.gene_score('ZNF717')
    print('ZNF717', f'{sc:.3f}')
    if sc != 0.0:
        ok = False
    sys.exit(0 if ok else 1)
