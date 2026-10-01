#!/usr/bin/env bash
# Rebuild TRACK1_TECHNICAL_DOCUMENTATION.pdf (run from repo root).
set -euo pipefail

# Pipeline figure v2 (docs/figures/pipeline_figure.py): tools, freshness of the committed
# SVG, exactly one figure reference in the document, then render the PDF pandoc embeds.
command -v dot >/dev/null || { echo "build_pdf.sh: graphviz 'dot' not found" >&2; exit 1; }
fc-match "DejaVu Sans" | grep -q DejaVuSans || { echo "build_pdf.sh: DejaVu Sans font missing" >&2; exit 1; }
python3 docs/figures/pipeline_figure.py --check
python3 docs/figures/pipeline_figure.py --count-refs docs/TECHNICAL_DOCUMENTATION.md --expect 1 >/dev/null \
  || { echo "build_pdf.sh: expected exactly 1 pipeline figure reference in the document" >&2; exit 1; }
python3 docs/figures/pipeline_figure.py --out docs/figures

pandoc docs/TECHNICAL_DOCUMENTATION.md \
  -o TRACK1_TECHNICAL_DOCUMENTATION.pdf \
  --toc --toc-depth=2 --pdf-engine=xelatex \
  --resource-path=docs:docs/figures \
  -H docs/pdf_style.tex \
  --metadata title="MVA Hackathon 2026 — Track 1 Pipeline · Technical Documentation" \
  --metadata mainfont="DejaVu Serif" \
  --metadata fontsize=10pt

# layout gate: fail the build if any word crosses the text edge
python3 docs/check_pdf_layout.py TRACK1_TECHNICAL_DOCUMENTATION.pdf
echo "built and layout-gated TRACK1_TECHNICAL_DOCUMENTATION.pdf"
