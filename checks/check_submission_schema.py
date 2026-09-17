#!/usr/bin/env python3
"""Verify submission_track1.csv matches the challenge schema and invariants.

Checks (all must pass):
  1. Header is exactly the challenge template header.
  2. Exactly 10 data rows, all for PROBAND01.
  3. Chromosomes carry the `chr` prefix; positions are positive integers.
     (Allele 1 / allele 2 order within a row is not constrained by the template.)
  4. epcr is a probability in [0, 1] and rows are sorted by epcr descending.
  5. finding_type is in {primary, secondary}.
  6. Rank 1 is the expected BUB1B compound-heterozygous finding
     (chr15:40209701 T>G + chr15:40220612 T>G).

Exit 0 and print OK lines; exit 1 on the first failed check.
"""
import csv
import re
import sys

CSV_PATH = "submission_track1.csv"
EXPECTED_HEADER = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes",
]
EXPECTED_RANK1 = ("chr15", "40209701", "T", "G", "chr15", "40220612", "T", "G")


def fail(msg):
    print(f"FAIL: {msg}")
    sys.exit(1)


with open(CSV_PATH, newline="") as fh:
    rows = list(csv.reader(fh))

if not rows or rows[0] != EXPECTED_HEADER:
    fail(f"header mismatch: {rows[0] if rows else 'empty file'}")
print("OK header")

data = rows[1:]
if len(data) != 10:
    fail(f"expected 10 data rows, found {len(data)}")
print("OK row count (10)")

prev_epcr = 1.0
for i, r in enumerate(data, start=1):
    if r[0] != "PROBAND01":
        fail(f"row {i}: proband_id={r[0]}")
    for col in ("chrom_1", "chrom_2"):
        idx = EXPECTED_HEADER.index(col)
        if not re.fullmatch(r"chr[0-9XYMT]+", r[idx]):
            fail(f"row {i}: {col}={r[idx]} lacks chr prefix")
    for col in ("pos_1", "pos_2"):
        idx = EXPECTED_HEADER.index(col)
        if not r[idx].isdigit() or int(r[idx]) <= 0:
            fail(f"row {i}: {col}={r[idx]} not a positive integer")
    try:
        epcr = float(r[9])
    except ValueError:
        fail(f"row {i}: epcr={r[9]} not numeric")
    if not 0.0 <= epcr <= 1.0:
        fail(f"row {i}: epcr={epcr} outside [0,1]")
    if epcr > prev_epcr + 1e-12:
        fail(f"row {i}: epcr not sorted descending ({epcr} > {prev_epcr})")
    prev_epcr = epcr
    if r[10] not in {"primary", "secondary"}:
        fail(f"row {i}: finding_type={r[10]}")
print("OK row-level invariants (ids, chr prefix, positions, epcr, type)")

if tuple(data[0][1:9]) != EXPECTED_RANK1:
    fail(
        "rank 1 is not the expected BUB1B compound-heterozygote "
        f"(chr15:40209701 T>G + chr15:40220612 T>G); got {data[0][1:9]}"
    )
print("OK rank 1 = BUB1B p.Leu737Ter + p.Asn1002Lys")

print("ALL CHECKS PASSED")
