#!/usr/bin/env python3
"""Validate a Track-1 predictions CSV against the challenge's published schema.

Usage: python3 checks/check_submission_schema.py [CSV]
Default CSV: submission_track1.csv at the repo root.
Prints one line per problem, then "ALL CHECKS PASSED" (exit 0) or exits 1.

Error messages name the row, column and rule only, never a cell value, so the
check can run on the real submission without echoing its contents.

Schema (public "Submit - Track 1" tab): 12 exact fields; up to 10 rows; GRCh38
coordinates on `chrN` contigs; second-variant fields for compound-het pairs
only, blank otherwise; epcr in (0, 1]; finding_type primary|secondary.
Additional house rules: one proband_id, no duplicate rows, rows ordered by
non-increasing epcr (the file's order matches its ranking).
"""

import csv
import itertools
import math
import re
import sys
from pathlib import Path

HEADER = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]
MAX_ROWS = 10
CHROM = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")
ALLELE = re.compile(r"^[ACGTN]+$")
POS = re.compile(r"^[1-9][0-9]*$")
FINDING_TYPES = {"primary", "secondary"}
DEFAULT_CSV = Path(__file__).resolve().parents[1] / "submission_track1.csv"


def _variant_errors(row, n, k):
    errs = []
    chrom, pos, ref, alt = (row[f"{c}_{k}"] for c in ("chrom", "pos", "ref", "alt"))
    if not CHROM.match(chrom):
        errs.append(f"row {n}: chrom_{k}: not a GRCh38 chrN contig")
    if not POS.match(pos):
        errs.append(f"row {n}: pos_{k}: not a positive integer")
    for col, allele in ((f"ref_{k}", ref), (f"alt_{k}", alt)):
        if not ALLELE.match(allele):
            errs.append(f"row {n}: {col}: not an uppercase ACGTN allele")
    if ref and ref == alt:
        errs.append(f"row {n}: ref_{k}: equals alt_{k}")
    return errs


def validate(path):
    """Return a list of problems (empty = valid). Messages never contain cell values."""
    with open(path, newline="") as fh:
        records = list(csv.reader(fh))
    if not records or records[0] != HEADER:
        return ["header: not the 12 exact schema fields in order"]
    body = records[1:]
    errs = []
    if not 1 <= len(body) <= MAX_ROWS:
        errs.append(f"row count: {len(body)} (must be 1..{MAX_ROWS})")
    rows, epcrs = [], []
    for n, rec in enumerate(body, 1):
        if len(rec) != len(HEADER):
            errs.append(f"row {n}: field count {len(rec)} (must be {len(HEADER)})")
            continue
        row = dict(zip(HEADER, rec))
        rows.append((n, row))
        if not row["proband_id"]:
            errs.append(f"row {n}: proband_id: empty")
        errs += _variant_errors(row, n, 1)
        second = [row[c] for c in ("chrom_2", "pos_2", "ref_2", "alt_2")]
        if all(second):
            errs += _variant_errors(row, n, 2)
        elif any(second):
            errs.append(f"row {n}: second variant: partially filled (all four fields or none)")
        try:
            e = float(row["epcr"])
            ok = math.isfinite(e) and 0 < e <= 1
        except ValueError:
            e, ok = None, False
        if ok:
            epcrs.append((n, e))
        else:
            errs.append(f"row {n}: epcr: not a number in (0, 1]")
        if row["finding_type"] not in FINDING_TYPES:
            errs.append(f"row {n}: finding_type: not primary|secondary")
    if len({r["proband_id"] for _, r in rows}) > 1:
        errs.append("proband_id: more than one value across rows")
    seen = {}
    for n, r in rows:
        key = tuple(r[c] for c in HEADER[1:9])
        if key in seen:
            errs.append(f"row {n}: duplicate of row {seen[key]} (same variant fields)")
        seen.setdefault(key, n)
    for (n0, e0), (n1, e1) in itertools.pairwise(epcrs):
        if e1 > e0:
            errs.append(f"row {n1}: epcr order: higher than row {n0} (rows must be non-increasing)")
    return errs


def main(argv):
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_CSV
    if not path.is_file():
        print("input: file not found")
        return 1
    errs = validate(path)
    for e in errs:
        print(e)
    if errs:
        print(f"FAILED: {len(errs)} problem(s)")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
