#!/usr/bin/env python3
"""PDF layout gate: no word may cross the right text edge.

The right text edge is derived from the geometry settings in docs/pdf_style.tex
(A4 width minus the right margin). `pdftotext -bbox` gives each word's xMax in
PostScript points. microtype deliberately lets punctuation and hyphens protrude
slightly into the margin, so words may overshoot the edge by at most
TOLERANCE_PT (1 mm); anything further fails the gate.

Usage: check_pdf_layout.py [PDF] [--min-pages N] [--min-words N]
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STYLE = HERE / "pdf_style.tex"
A4_WIDTH_PT = 595.28
TOLERANCE_PT = 2.8  # 1 mm: microtype protrusion of punctuation/hyphens
UNIT_PT = {"mm": 72 / 25.4, "cm": 72 / 2.54, "in": 72.0, "pt": 72 / 72.27}


def _length_pt(value):
    m = re.fullmatch(r"\s*([\d.]+)\s*(mm|cm|in|pt)\s*", value)
    if not m:
        raise ValueError(f"unsupported geometry length: {value!r}")
    return float(m.group(1)) * UNIT_PT[m.group(2)]


def text_edge_pt(style_tex, page_width_pt=A4_WIDTH_PT):
    """x position of the right text edge from a \\usepackage[...]{geometry} line."""
    m = re.search(r"\\usepackage\[([^\]]*)\]\{geometry\}", style_tex)
    if not m:
        raise ValueError("no \\usepackage[...]{geometry} line in the style file")
    opts = dict(o.split("=", 1) for o in (p.strip() for p in m.group(1).split(",")) if "=" in o)
    right = opts.get("right", opts.get("margin"))
    if right is None:
        raise ValueError("geometry sets neither right= nor margin=")
    return page_width_pt - _length_pt(right)


def _pages(bbox_xml):
    return re.findall(r"<page[^>]*>(.*?)</page>", bbox_xml, re.S)


def parse_ok(bbox_xml, min_pages, min_words):
    pages = _pages(bbox_xml)
    return len(pages) >= min_pages and sum(len(re.findall(r"<word ", p)) for p in pages) >= min_words


def violations(bbox_xml, edge, tolerance):
    """Words whose xMax lies more than `tolerance` points beyond `edge`."""
    out = []
    for i, page in enumerate(_pages(bbox_xml), 1):
        for x1, text in re.findall(r'xMax="([\d.]+)"[^>]*>([^<]*)</word>', page):
            if float(x1) > edge + tolerance:
                out.append(f"p{i} x1={float(x1):.1f} '{text[:40]}'")
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf", nargs="?", default="TRACK1_TECHNICAL_DOCUMENTATION.pdf")
    ap.add_argument("--min-pages", type=int, default=5)
    ap.add_argument("--min-words", type=int, default=1000)
    a = ap.parse_args(argv)
    if shutil.which("pdftotext") is None:
        print("check_pdf_layout: pdftotext not found (install poppler-utils)", file=sys.stderr)
        return 2
    edge = text_edge_pt(STYLE.read_text())
    bbox = subprocess.run(["pdftotext", "-bbox", a.pdf, "-"], capture_output=True, text=True, check=True).stdout
    if not parse_ok(bbox, a.min_pages, a.min_words):
        print(f"check_pdf_layout: suspicious parse of {a.pdf} (below --min-pages/--min-words)", file=sys.stderr)
        return 1
    viol = violations(bbox, edge, TOLERANCE_PT)
    if viol:
        print(f"LAYOUT GATE FAILED (text edge {edge:.1f}pt + {TOLERANCE_PT}pt):")
        print("\n".join(viol))
        return 1
    print(f"LAYOUT GATE OK: none beyond {edge + TOLERANCE_PT:.1f}pt (text edge {edge:.1f}pt + {TOLERANCE_PT}pt)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
