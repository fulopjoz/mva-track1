#!/usr/bin/env python3
"""Produce the phase-block table from a whatshap phased VCF (rows file)."""
import sys
blocks = {}
singletons = []
targets = {40209701, 40220612, 40216470, 40221421}
with open(sys.argv[1]) as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 10:
            continue
        pos = int(c[1]); fmt = c[8].split(":"); samp = c[9].split(":")
        d = dict(zip(fmt, samp))
        ps = d.get("PS")
        if ps:
            blocks.setdefault(ps, []).append(pos)
        else:
            singletons.append(pos)

print("PHASE BLOCKS (variants sharing a PS value):")
print("block_anchor | n_sites | range(chr15)")
for ps in sorted(blocks, key=lambda x: int(x)):
    p = blocks[ps]
    print(f"{ps} | {len(p)} | {min(p)}-{max(p)}")
print(f"\nUnphased singletons (no PS): {len(singletons)}")
print("Target variant sites:")
for pos in sorted(targets):
    inblk = "singleton (unphased)"
    for ps, p in blocks.items():
        if pos in p:
            inblk = f"in block {ps}"
    print(f"  {pos}: {inblk}")
