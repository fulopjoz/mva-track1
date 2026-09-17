#!/usr/bin/env python3
"""Aggregate exporter (Track 1).

Reads the ranked table (annot TSV, patient-derived; lives on the compute
host) and writes figures/aggregates.json: summary numbers only, no
variant-level rows. The JSON is what the shipped figures are built from
(figures/make_figures.py), and it is what the README's grep/awk checks
read. Run on the compute host after a pipeline pass:

  python3 pipeline/scripts/export_aggregates.py \
      ~/hf-rare/pipeline/output/annot_full.tsv > figures/aggregates.json

Fields (all integers/floats/strings, all aggregates):
  n_rows, n_genes: output table size.
  funnel: record counts per pipeline stage (measured at run time).
  epcr_hist_edges/counts: histogram of row scores, 20 bins of 0.05.
  median_epcr, max_epcr: distribution landmarks.
  type_counts: candidate classes (compound_het pairs vs hom-alt singles).
  consequence_top: top consequence-class pairs with counts.
  top10: first 10 ranked rows as [gene, epcr, type].
"""
from __future__ import annotations

import json
import statistics
import sys
from collections import Counter


def main():
    path = sys.argv[1]
    scores = []
    genes = set()
    types = Counter()
    consequences = Counter()
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            gene = f[col['gene']]
            epcr = float(f[col['epcr']])
            ftype = f[col['type']]
            c1 = f[col['consequence_1']] or '-'
            c2 = f[col['consequence_2']] or '-'
            if ftype == 'hom_alt':
                cons = f'{c1} (hom-alt)'
            else:
                cons = c1 if c2 == '-' else f'{c1} + {c2}'
            genes.add(gene)
            types[ftype] += 1
            consequences[cons] += 1
            scores.append(epcr)
            rows.append([gene, round(epcr, 4), ftype])
    scores.sort()
    edges = [i / 20 for i in range(21)]
    counts = [0] * 20
    for s in scores:
        counts[min(int(s * 20), 19)] += 1
    out = {
        'n_rows': len(scores),
        'n_genes': len(genes),
        'epcr_hist_edges': edges,
        'epcr_hist_counts': counts,
        'median_epcr': round(statistics.median(scores), 4),
        'max_epcr': round(scores[-1], 4),
        'type_counts': dict(types),
        'consequence_top': [[k, v] for k, v in consequences.most_common(12)],
        'top10': rows[:10],
    }
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == '__main__':
    main()
