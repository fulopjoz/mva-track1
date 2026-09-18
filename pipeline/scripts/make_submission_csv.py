#!/usr/bin/env python3
"""Emit a Track 1 submission CSV file from a pipeline‑generated TSV.

Schema (exact, per official template):
  proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2,
  epcr, finding_type, notes

Rules:
  * ``proband_id`` is a fixed identifier (e.g. ``PROBAND01``) defined in the
    dataset README/template.
  * Chromosome coordinates acquire the ``chr`` prefix *here*; the internal
    pipeline uses bare contig names.
  * Output is limited to a maximum of ten rows, ranked by the pipeline EPCR
    score (the TSV order already reflects ranking).
  * ``finding_type`` is ``primary`` for compound‑heterozygous pairs and
    ``secondary`` for homozygous single‑variant candidates.
  * Compound‑heterozygous pairs share a single row (both alleles listed).
  * ``epcr`` values are clamped to the open interval (0, 1] before formatting.

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

    with open(args.outfile, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerows(out)
    print(f'wrote {args.outfile} rows={len(out)}')


if __name__ == '__main__':
    main()
