#!/usr/bin/env python3
"""Measure fragment (insert) span from a minimap2 SAM on stdin (paired reads)."""
import sys, statistics
spans=[]; tot=0
for line in sys.stdin:
    if line.startswith("@"): continue
    f=line.rstrip("\n").split("\t")
    if len(f)<11: continue
    flag=int(f[1]); pos=int(f[3]); mate=int(f[7]); rlen=len(f[9])
    if flag & 0x4: continue   # unmapped
    tot+=1
    if mate>0:
        spans.append(abs(pos-mate)+rlen)
spans.sort()
n=len(spans)
print(f"total mapped alignments: {tot}, with mate info: {n}")
if n:
    def q(x): return spans[int((n-1)*x)]
    print(f"insert(bp) all: p50={q(.5)} p90={q(.9)} max={spans[-1]}")
    # proper pairs = span within a sane cap (<2000 bp)
    proper=[s for s in spans if s<2000]
    proper.sort()
    m=len(proper)
    if m:
        def qq(x): return proper[int((m-1)*x)]
        print(f"proper pairs (span<2000bp): n={m} median={qq(.5)} p90={qq(.9)} p95={qq(.95)} max={qq(1)}")
        frac809=sum(1 for s in proper if s>=809)/m*100
        frac700=sum(1 for s in proper if s>=700)/m*100
        print(f"proper pairs spanning >=809 bp (the 40220612-40221421 gap): {frac809:.2f}%")
        print(f"proper pairs spanning >=700 bp: {frac700:.2f}%")
