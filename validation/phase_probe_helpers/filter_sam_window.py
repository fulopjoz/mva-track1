#!/usr/bin/env python3
"""Triage filter: read SAM on stdin, emit (stdout) unique read names (qnames)
whose chr15 alignment overlaps the BUB1B window with MAPQ>=MINQ.

Used ONLY for candidate-read selection; selected reads are later realigned to
the whole genome for honest MAPQ/phasing.
"""
import sys

WIN_START = 40_000_000
WIN_END = 40_400_000
MINQ = 20

names = set()
for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 11 or f[5] == "*":
        continue
    flag, rname, pos, mapq, seq = int(f[1]), f[2], int(f[3]), int(f[4]), f[9]
    if rname != "chr15":
        continue
    if pos > WIN_END or pos + len(seq) < WIN_START:
        continue
    if mapq < MINQ:
        continue
    names.add(f[0])  # qname (full first token)

for n in sorted(names):
    sys.stdout.write(n + "\n")
sys.stderr.write(f"triage: unique candidate reads={len(names)}\n")
