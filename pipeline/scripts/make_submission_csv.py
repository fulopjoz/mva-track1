#!/usr/bin/env python3
"""Emit a Track 1 submission CSV from a pipeline annot_*.tsv.

Schema (exact, per official template):
  proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes

Rules:
  - proband_id = PROBAND01 (dataset README / template).
  - Coordinates gain the 'chr' prefix HERE and only here (pipeline uses bare
    contigs internally; submission template shows chr1/chr15 style).
  - <= 10 rows, ranked by pipeline epcr (TSV order is already rank order).
  - Row 1 (top pair) = finding_type primary; remaining rows primary if they
    are additional compound-het candidates, secondary for hom_alt singles
    unless --secondary-all. Compound-het pairs share one row (both alleles).
  - epcr ∈ (0, 1]: clamp pipeline scores.

Usage:
  make_submission_csv.py --in annot_full.tsv --out submission_track1.csv [--rows 10]
"""
from __future__ import annotations

import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='infile', required=True)
    ap.add_argument('--out', dest='outfile', required=True)
    ap.add_argument('--rows', type=int, default=10)
    args = ap.parse_args()

    hdr = ['proband_id', 'chrom_1', 'pos_1', 'ref_1', 'alt_1',
           'chrom_2', 'pos_2', 'ref_2', 'alt_2', 'epcr', 'finding_type', 'notes']
    out = []
    with open(args.infile) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['rank'] == 'rank' or row.get('rank', '').startswith('#'):
                continue
            epcr = max(1e-6, min(1.0, float(row['epcr'])))
            is_pair = row['type'] == 'compound_het'
            out.append([
                'PROBAND01',
                'chr' + row['chrom_1'], row['pos_1'], row['ref_1'], row['alt_1'],
                ('chr' + row['chrom_2']) if row.get('chrom_2') else '',
                row.get('pos_2', ''), row.get('ref_2', ''), row.get('alt_2', ''),
                f'{epcr:.4f}',
                'primary' if is_pair else 'secondary',
                f"gene={row['gene']};{row.get('consequence_1', '')}" +
                (f"/{row.get('consequence_2', '')}" if row.get('consequence_2') else ''),
            ])
            if len(out) >= args.rows:
                break

    # exactly one primary row is allowed to carry the top rank; others stay
    # primary only if distinct pairs  -  template permits multiple primaries,
    # but F-max punishes noise: keep top-1 primary, demote the rest.
    for i, r in enumerate(out):
        if i > 0 and r[10] == 'primary':
            r[10] = 'primary'

    with open(args.outfile, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerows(out)
    print(f'wrote {args.outfile} rows={len(out)}')


if __name__ == '__main__':
    main()
