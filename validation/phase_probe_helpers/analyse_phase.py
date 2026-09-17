#!/usr/bin/env python3
"""Analyse whatshap phased VCF + BAM to produce the phase-block evidence table."""
import sys, subprocess, re, statistics

VCF = sys.argv[1]

# Parse phased VCF: position, ref>alt, GT, PS(if present)
variants = []  # (pos, reffull, gt, ps)
with open(VCF) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10:
            continue
        pos = int(c[1]); ref = c[3]; alt = c[4]; fmt = c[8].split(":"); samp = c[9].split(":")
        d = dict(zip(fmt, samp))
        gt = d.get("GT", "./.")
        ps = d.get("PS", None)
        variants.append((pos, f"{ref}>{alt}", gt, ps))

variants.sort()
print(f"total variants: {len(variants)}")

# --- block structure (variants that share a PS are in the same phase block) ---
by_ps = {}
phased = 0
for pos, change, gt, ps in variants:
    if ps:  # has phase set -> phased
        by_ps.setdefault(ps, []).append((pos, change, gt))
        phased += 1
print(f"variants with a PS block: {phased}; distinct blocks: {len(by_ps)}")

key_pos = {40209701: "variant1_Leu737Ter", 40220612: "variant2_Asn1002Lys"}
for pos, change, gt, ps in variants:
    if pos in key_pos or pos in (40216470, 40221421):
        tag = key_pos.get(pos, f"het_neighbor@{pos}")
        ps_str = ps if ps else "NONE(singleton)"
        print(f"  {tag}: chrs15:{pos} {change} GT={gt} PS={ps_str}")

# nearest het neighbor distance for each pathogenic variant (exclude self)
hets = [v for v in variants if v[2].startswith("0/1") or (v[2] != "./." and "|" in v[2])]
print("\nnearest heterozygous site to each pathogenic variant:")
for tpos in (40209701, 40220612):
    dmin = min(abs(tpos - h[0]) for h in hets if h[0] != tpos)
    print(f"  {tpos}: nearest OTHER het at {dmin} bp")

# fragment (insert) length: minimap2 does not set TLEN, so estimate from mate
# coordinates: span = |pos - mate_pos| + read_length for pairs with a mapped mate.
bam = sys.argv[2]
print("\nfragment (insert) length stats (pairs in BUB1B window, from mate coords):")
p = subprocess.run(["samtools","view",bam,"chr15:40100000-40270000"],
                   capture_output=True, text=True).stdout
spans = []
for line in p.splitlines():
    f = line.split("\t")
    start = int(f[3]); mate = int(f[7]); rlen = len(f[9])
    if mate > 0:
        d = abs(start - mate) + rlen
        if d > 0:
            spans.append(d)
spans.sort()
n = len(spans)
print(f"  n pairs: {n}")
if n:
    import math
    def q(frac): return spans[int((n-1)*frac)]
    print(f"  p50={q(.5)} p90={q(.9)} p99={q(.99)} max={spans[-1]} mean={statistics.mean(spans):.0f}")

# confirm no read spans both variant sites (start<=40209701 and end>=40220612)
p2 = subprocess.run(["samtools","view",bam,"chr15:40209701-40220612"],
                    capture_output=True, text=True).stdout
span=0
for line in p2.splitlines():
    f = line.split("\t")
    start=int(f[3]); cigar=f[5]; seq=f[9]
    # approximate end via M ops
    end = start + sum(int(x) for x in re.findall(r"(\d+)M", cigar))
    if end is not None and start <= 40209701 and end >= 40220612:
        span += 1
print(f"\nreads spanning BOTH variant positions (start<=40209701 & end>=40220612): {span}")
