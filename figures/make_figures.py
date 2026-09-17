#!/usr/bin/env python3
"""Rebuild the submission figures from the shipped aggregate JSON file.

The script loads ``aggregates.json`` (produced by ``export_aggregates.py``) and
generates four PNG images:
  * ``fig1_funnel.png`` – pipeline funnel counts.
  * ``fig2_score_distribution.png`` – histogram of EPCR scores.
  * ``fig3_top10.png`` – the top‑10 ranked genes with scores.
  * ``fig4_consequences.png`` – most frequent consequence‑class pairs.

All plots use a non‑interactive Agg backend suitable for headless environments.
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

BLUE = '#1f77b4'
ORANGE = '#DD7A00'
GREY = '#999999'


def load_aggregates():
    with open(os.path.join(HERE, 'aggregates.json')) as fh:
        return json.load(fh)


def funnel_counts():
    """First (final-run) set of instrumented counts in provenance.log."""
    want = ('vcf_in', 'qc_pass', 'coding_or_splice', 'af_absent', 'af_common',
            'compound_het_pairs', 'hom_candidates')
    got = {}
    with open(os.path.join(ROOT, 'provenance.log')) as fh:
        for line in fh:
            m = re.match(r'\[instrument\] (.*)', line)
            if not m:
                continue
            for k, v in re.findall(r'(\w+)=(\d+)', m.group(1)):
                if k in want and k not in got:
                    got[k] = int(v)
            if len(got) == len(want):
                break
    if len(got) != len(want):
        raise SystemExit('provenance.log: instrument counts incomplete')
    return got


def load_weights():
    import yaml
    with open(os.path.join(ROOT, 'pipeline', 'conf', 'track1_config.yml')) as fh:
        conf = yaml.safe_load(fh)
    s = conf['scoring']
    return s['founder_weight'], s['pathogenicity_weight'], s['rare_weight'], s['hpo_weight']


def fm(n):
    return f'{n:,}'


def pipeline_diagram(agg):
    c = funnel_counts()
    fw, pw, rw, hw = load_weights()
    bub = [row for row in agg['top10'] if row[0] == 'BUB1B'][0][1]
    boxes = [
        ('INPUT', f'WGS VCF (GRCh38)\n{fm(c["vcf_in"])} records, ~80 GB', '#5b7f95'),
        ('1. Quality filter', f'GQ>=20, depth>=10, alt fraction>=0.25\n{fm(c["qc_pass"])} records kept', BLUE),
        ('2. Effect calling', f'GENCODE v47 MANE transcripts, local\nCADD v1.7 term neutral (tiles not installed)\n{fm(c["coding_or_splice"])} coding or splice candidates', BLUE),
        ('3. Frequency and background', f'gnomAD v4.1 exome sites, exact allele key\nClinVar background only, no rescue\nAF absent {fm(c["af_absent"])}, common {fm(c["af_common"])}', BLUE),
        ('4. Candidate studies (4 parallel)', f'compound-het pairs: LoF + second allele\n3 homozygous recessive models\n{fm(c["compound_het_pairs"])} pairs + {fm(agg["type_counts"].get("hom_alt", c["hom_candidates"]))} hom singles', BLUE),
        ('5. Weighted scoring', f'founder {fw}, pathogenicity {pw}, rarity {rw}\n+ HPO phenotype term, weight {hw}', '#7b5aa6'),
        ('RANKED OUTPUT', f'{fm(agg["n_rows"])} rows across {fm(agg["n_genes"])} genes\nrank 1: BUB1B stop + missense, score {bub:.4f}', ORANGE),
        ('SUBMISSION', 'submission_track1.csv, 10 rows\nschema-checked', '#4e8a5a'),
    ]
    side = [
        ('GUARDS', 'honesty CI: 15 invariants\n(repro in reruns at equal weight)', '#4e8a5a', -1.05),
        ('CROSS-CHECK', 'Exomiser 14.0.0 A/B\nBUB1B rank 1', '#4e8a5a', 1.05),
    ]
    fig, ax = plt.subplots(figsize=(9, 13))
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-0.6, len(boxes) + 0.4)
    ax.axis('off')
    n = len(boxes)
    for i, (title, body, color) in enumerate(boxes):
        y = n - 1 - i
        round_bbox = dict(boxstyle='round,pad=0.55', facecolor='white', edgecolor=color, linewidth=2.5)
        lines = body.split('\n')
        text = '\n'.join([title, *lines])
        ax.text(0, y, text, ha='center', va='center', fontsize=11,
                fontweight='bold', bbox=round_bbox,
                color=color, linespacing=1.5)
        if i < n - 1:
            ax.annotate('', xy=(0, y - 0.62), xytext=(0, y - 0.28),
                        arrowprops=dict(arrowstyle='-|>', color='#444444', lw=2))
    # side boxes hang level with the scoring stage (index 5)
    y_side = n - 1 - 5
    for title, body, color, x in side:
        ax.text(x * 1.32, y_side, f'{title}\n{body}', ha='center', va='center',
                fontsize=10, fontweight='bold', color=color,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                          edgecolor=color, linewidth=2), linespacing=1.5)
    ax.text(0, -0.45, 'Single pass, CPU only. Numbers from provenance.log and figures/aggregates.json.',
            ha='center', fontsize=10, color='#555555')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'pipeline_diagram.png'), dpi=110)
    plt.close(fig)


def fig1_funnel(agg):
    c = funnel_counts()
    vcf_in, qc_pass, coding = c['vcf_in'], c['qc_pass'], c['coding_or_splice']
    ranked = agg['n_rows']
    pairs = agg['type_counts'].get('compound_het', 0)
    homes = agg['type_counts'].get('hom_alt', 0)
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = ['VCF records\n(input WGS)', 'QC pass\n(GQ>=20, DP>=10,\nAF>=0.25)',
              'Coding / splice\ncandidates',
              f'Ranked rows\n({pairs:,} pairs +\n{homes:,} hom-alt)']
    values = [vcf_in, qc_pass, coding, ranked]
    colors = [GREY, GREY, BLUE, ORANGE]
    bars = ax.bar(labels, values, color=colors)
    ax.set_yscale('log')
    ax.set_ylim(100, 5e7)
    ax.set_ylabel('records (log scale)')
    ax.set_title(f'Track 1 analysis funnel: {vcf_in / 1e6:.1f} M WGS records '
                 f'to {ranked:,} ranked rows')
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.15, f'{v:,}',
                ha='center', fontsize=12, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fig1_funnel.png'), dpi=100)
    plt.close(fig)


def fig2_distribution(agg):
    edges = agg['epcr_hist_edges']
    counts = agg['epcr_hist_counts']
    med = agg['median_epcr']
    bub = [row for row in agg['top10'] if row[0] == 'BUB1B'][0][1]
    centers = [(a + b) / 2 for a, b in zip(edges, edges[1:])]
    width = edges[1] - edges[0]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(centers, counts, width=width * 0.9, color=BLUE)
    ax.axvline(med, ls=':', color=GREY, lw=2)
    ax.text(med + 0.02, max(counts) * 0.9, f'median {med}', color=GREY)
    ax.axvline(bub, ls='--', color=ORANGE, lw=2)
    ax.text(bub - 0.02, max(counts) * 0.75,
            f'BUB1B compound-het pair\nrank 1, epcr {bub:.4f}',
            ha='right', color=ORANGE)
    ax.set_xlabel('pathogenicity rank score (epcr)')
    ax.set_ylabel('ranked rows')
    ax.set_title(f"Distribution of {agg['n_rows']:,} ranked-row scores")
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fig2_score_distribution.png'), dpi=100)
    plt.close(fig)


def fig3_top10(agg):
    labels, scores, genes = [], [], []
    for rank, (gene, epcr, ftype) in enumerate(reversed(agg['top10'])):
        r = len(agg['top10']) - rank
        labels.append(f'{gene} (rank {r})')
        genes.append(gene)
        scores.append(epcr)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    colors = [ORANGE if g == 'BUB1B' else BLUE for g in genes]
    bars = ax.barh(labels, scores, color=colors)
    for b, s in zip(bars, scores):
        ax.text(s + 0.005, b.get_y() + b.get_height() / 2, f'{s:.4f}',
                va='center', fontsize=11)
    ax.set_xlim(0, 0.66)
    ax.set_xlabel('compound-het pair score (epcr)')
    ax.set_title('Top 10 compound-het candidates')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fig3_top10.png'), dpi=100)
    plt.close(fig)


def fig4_consequences(agg):
    rows = agg['consequence_top'][:8]
    total = agg['n_rows']
    shown = sum(v for _, v in rows)
    labels = [k for k, _ in reversed(rows)]
    values = [v for _, v in reversed(rows)]
    colors = [ORANGE if 'stop_gained' in k else BLUE for k in labels]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(labels, values, color=colors)
    for b, v in zip(bars, values):
        ax.text(v + 40, b.get_y() + b.get_height() / 2, f'{v:,}',
                va='center', fontsize=11)
    ax.text(0.97, 0.05, f'largest 8 categories shown ({shown:,} of {total:,} rows)',
            transform=ax.transAxes, ha='right', color=GREY)
    ax.set_xlabel('ranked rows')
    ax.set_title('Functional composition of the ranked candidate set')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fig4_consequences.png'), dpi=100)
    plt.close(fig)


def main():
    agg = load_aggregates()
    pipeline_diagram(agg)
    fig1_funnel(agg)
    fig2_distribution(agg)
    fig3_top10(agg)
    fig4_consequences(agg)
    print('wrote pipeline_diagram.png, fig1_funnel.png, '
          'fig2_score_distribution.png, fig3_top10.png, fig4_consequences.png')


if __name__ == '__main__':
    main()
