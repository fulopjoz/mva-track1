#!/usr/bin/env python3
"""Export aggregated summary statistics for a Track 1 pipeline run.

The script reads the ranked annotation TSV (generated on the compute host) and
writes a JSON file containing only summary numbers – no per‑variant rows. The
JSON output is read by ``figures/make_figures.py`` to produce the final
submission figures and is also validated by the README’s grep/awk checks.

Typical usage on the compute host after the pipeline has completed:

  python3 pipeline/scripts/export_aggregates.py \
      ~/hf-rare/pipeline/output/annot_full.tsv > figures/aggregates.json

Fields (all integers, floats, or strings):
  * ``n_rows`` – total number of rows in the input table.
  * ``n_genes`` – count of distinct genes represented.
  * ``funnel`` – record counts per pipeline stage (measured at runtime).
  * ``epcr_hist_edges`` / ``epcr_hist_counts`` – histogram of EPCR scores
    (20 bins of width 0.05).
  * ``median_epcr`` and ``max_epcr`` – distribution landmarks.
  * ``type_counts`` – counts of candidate classes (compound‑het pairs vs
    hom‑alt singles).
  * ``consequence_top`` – the most frequent consequence‑class pairs with
    their counts (up to twelve entries).
  * ``top10`` – the first ten ranked rows as ``[gene, epcr, type]``.
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
