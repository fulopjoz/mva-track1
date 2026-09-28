#!/usr/bin/env bash
# Rebuild TRACK1_TECHNICAL_DOCUMENTATION.pdf (run from repo root).
set -euo pipefail
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
