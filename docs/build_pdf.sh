#!/usr/bin/env bash
# Rebuild TRACK1_TECHNICAL_DOCUMENTATION.pdf (run from repo root).
set -euo pipefail
pandoc docs/TECHNICAL_DOCUMENTATION.md \
  -o TRACK1_TECHNICAL_DOCUMENTATION.pdf \
  --toc --pdf-engine=xelatex \
  --resource-path=docs:docs/figures \
  --metadata title="MVA Track 1 Pipeline - Technical Documentation" \
  --metadata mainfont="DejaVu Serif"
echo "built TRACK1_TECHNICAL_DOCUMENTATION.pdf"
