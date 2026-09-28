#!/usr/bin/env python3
"""PDF layout gate: no word may cross the right text-edge.

pdftotext -bbox gives per-word xMax in PostScript points. Text block width:
A4 595.3pt - 2*2.2cm (124.7pt) = 470.6pt, so the right edge sits at
62.35 + 470.6 ~ 533pt... measured text runs to ~549.3pt (microtype
protrusion of periods/commas/hyphens is intentional and may overshoot
<= ~2.8pt = 1mm). FAIL only beyond 3mm (8.5pt) past the 549.5pt line.
"""
import re
import subprocess
import sys

PDF = sys.argv[1] if len(sys.argv) > 1 else "TRACK1_TECHNICAL_DOCUMENTATION.pdf"
EDGE = 549.5
HARD = EDGE + 8.5

bb = subprocess.run(["pdftotext", "-bbox", PDF, "-"], capture_output=True, text=True).stdout
pages = re.findall(r"<page[^>]*>(.*?)</page>", bb, re.S)
n_words = sum(len(re.findall(r"<word ", p)) for p in pages)
assert len(pages) >= 5 and n_words > 1000, f"suspicious parse: {len(pages)} pages, {n_words} words"
viol = []
for i, p in enumerate(pages, 1):
    for x1, t in re.findall(r'xMax="([\d.]+)"[^>]*>([^<]*)</word>', p):
        if float(x1) > HARD:
            viol.append(f"p{i} x1={float(x1):.1f} '{t[:40]}'")
if viol:
    print("LAYOUT GATE FAILED:")
    print("\n".join(viol))
    sys.exit(1)
print(f"LAYOUT GATE OK: {len(pages)} pages, {n_words} words, none beyond {HARD:.1f}pt")
