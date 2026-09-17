#!/usr/bin/env python3
"""Interleave index-aligned R1 and R2 FASTQs into one interleaved FASTQ."""
import sys

def read(fn):
    out = []
    with open(fn) as f:
        lines = f.read().split("\n")
    i = 0
    while i + 3 < len(lines) and lines[i].startswith("@"):
        out.append(lines[i:i+4]); i += 4
    return out

r1, r2 = read(sys.argv[1]), read(sys.argv[2])
assert len(r1) == len(r2), "R1/R2 length mismatch"
with open(sys.argv[3], "w") as f:
    for a, b in zip(r1, r2):
        f.write("\n".join(a) + "\n")
        f.write("\n".join(b) + "\n")
print(f"interleaved pairs: {len(r1)}")
