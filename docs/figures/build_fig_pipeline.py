#!/usr/bin/env python3
"""Generate F_pipeline_overview.png — Track 1 pipeline overview diagram.

Layout rules (print legibility is the contract):
- single centered column of stage boxes; no side panels (their content was
  already carried inside the stage bodies)
- canvas 8.6 in wide, so a label at F points prints at F x 176/(25.4x8.6) pt
  on A4 portrait: body 9.3 pt -> ~7.5 pt, smallest text 8.8 pt -> ~7.1 pt
- every drawn text is measured against its container after draw; any escape
  aborts before savefig (self-gate)
All numbers come from REPORT.md / track1_config.yml (see figures appendix of
docs/TECHNICAL_DOCUMENTATION.md).
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# palette: white stages, navy ruling, neutral gray artifact — structure
# (solid flow arrows, dashed artifact outline, chip bands), not color, encodes roles
STAGE_FILL, STAGE_EDGE, STAGE_HEAD = "#FFFFFF", "#2F4B6E", "#1F3B5C"
ART_FILL, ART_EDGE = "#FFFFFF", "#9A9A9A"
CHIP_FILL, CHIP_EDGE = "#EFF3F7", "#2F4B6E"
BODY = "#262626"
MUTED = "#4A4A4A"

W, H = 100.0, 100.0
FIG_W_IN, FIG_H_IN = 8.6, 8.6
X0, X1 = 5.0, 95.0
BOXH = 8.9
ROW = [91.5, 82.2, 72.9, 63.6, 54.3, 45.0, 35.7, 19.2]  # y-top per stage
CHIP_Y0, CHIP_H = 0.55, 1.30
CHIP_BODY_FLOOR = 2.25
ART_X0, ART_X1 = 30.0, 70.0
ART_TOP, ART_H = 26.2, 3.4
FOOT_Y0, FOOT_Y1 = 0.8, 8.6

_GATE = []
GATE_MARGIN = 0.2


def _gate(text, box, name):
    _GATE.append((text, box, name))


def run_gate(fig, ax):
    """Measure every registered text; SystemExit(1) on any overflow."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = []
    for text, (x0, y0, x1, y1), name in _GATE:
        bb = text.get_window_extent(renderer=r)
        (u0, v0), (u1, v1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        if (u0 < x0 - GATE_MARGIN or u1 > x1 + GATE_MARGIN
                or v0 < y0 - GATE_MARGIN or v1 > y1 + GATE_MARGIN):
            bad.append((name, round(u0, 2), round(u1, 2), x0, x1))
    if bad:
        raise SystemExit(f"GATE FAIL: {bad}")
    print(f"GATE OK: {len(_GATE)} text artists within their boxes "
          f"(margin {GATE_MARGIN})")


def stage_box(ax, i, title, lines, chip=None):
    t, b = ROW[i], ROW[i] - BOXH
    ax.add_patch(FancyBboxPatch(
        (X0, b), X1 - X0, BOXH, boxstyle="round,pad=0,rounding_size=0.7",
        facecolor=STAGE_FILL, edgecolor=STAGE_EDGE, lw=1.4, zorder=2))
    tt = ax.text(X0 + 1.6, t - 1.55, title, ha="left", va="center",
                 fontsize=11.0, fontweight="bold", color=STAGE_HEAD, zorder=3)
    _gate(tt, (X0, b, X1, t), f"stage {i} title")
    by = t - 3.55
    for k, ln in enumerate(lines):
        bt = ax.text(X0 + 1.6, by - 1.35 * k, ln, ha="left", va="center",
                     fontsize=9.2, color=BODY, zorder=3)
        _gate(bt, (X0, b + CHIP_BODY_FLOOR, X1, t), f"stage {i} body {k}")
    if chip:
        cy = b + CHIP_Y0
        ct = ax.text(X0 + 1.6, cy + CHIP_H / 2, chip, ha="left", va="center",
                     fontsize=8.6, style="italic", color=STAGE_EDGE, zorder=3)
        _gate(ct, (X0, b, X1, b + CHIP_BODY_FLOOR + 0.2), f"stage {i} chip")
        ax.plot([X0 + 1.6, X1 - 1.6], [b + CHIP_BODY_FLOOR - 0.25] * 2,
                color=CHIP_EDGE, lw=0.6, zorder=2)


def flow_arrow(ax, y_from, y_to):
    ax.add_patch(FancyArrowPatch((50, y_from), (50, y_to),
                                 arrowstyle="-|>", mutation_scale=13, lw=1.6,
                                 color=STAGE_EDGE, shrinkA=0, shrinkB=0,
                                 zorder=1))


def main() -> None:
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "figure.dpi": 150,
        "savefig.dpi": 150,
    })
    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    canvas_box = (0.0, 0.0, W, H)

    tt = ax.text(W / 2, 98.1,
                 "MVA Hackathon 2026 — Track 1 variant-prioritization pipeline",
                 ha="center", va="center", fontsize=12.5, fontweight="bold",
                 color=STAGE_HEAD)
    _gate(tt, canvas_box, "figure title")
    sub = ax.text(W / 2, 94.9,
                  "patient WGS VCF  →  ranked variant table  →  submission CSV   ·   "
                  "single pass, all joins against local reference files",
                  ha="center", va="center", fontsize=9.0, color=MUTED)
    _gate(sub, canvas_box, "figure subtitle")

    stages = [
        ("Proband WGS VCF",
         ["5,012,204 records · GRCh38 · bgzipped VCF · bare contig names (data.vcf)"],
         None),
        ("Quality-control filter",
         ["FILTER PASS or absent · GQ ≥ 20 · DP ≥ 10 · alt fraction ≥ 0.25",
          "GT het or hom-alt · keeps 4,592,140 records"],
         "prioritize.py :: iter_qc_variants    [cyvcf2.VCF]"),
        ("Consequence annotation",
         ["MANE transcript overlap (bisect) · splice classes · codon-level effect",
          "from FASTA · CADD PHRED lookup · keeps 16,963 coding/splice candidates"],
         "annotate.py :: LocalAnnotator    [pysam.FastaFile · pysam.TabixFile]"),
        ("Rarity and ClinVar join",
         ["gnomAD v4.1 exomes AF, exact (pos, ref, alt), sorted per-chromosome merge;",
          "AF: hom ≤ 0.01 · het ≤ 0.0001 · None → neutral · ClinVar: background only"],
         "vcf.py :: GnomadAfJoin · ClinVarJoin    [pysam.TabixFile]"),
        ("HPO phenotype term",
         ["gene score vs 8 proband terms: 1.0 ancestor-closed exact match,",
          "0.5 parent-only match, mean over proband terms"],
         "_hpo.py :: HPOScorer.gene_score"),
        ("Scoring (WorkPrior)",
         ["epcr = Σ wᵢ·termᵢ / Σ wᵢ, clamped to (0, 1] · pair = stronger × (0.5 + 0.5 × weaker)",
          "weights: founder 0.45 · pathogenicity 0.60 · rarity 0.30 · HPO 2.0"],
         "prior.py :: WorkPrior.internal"),
        ("Candidate pairing and ranking",
         ["compound-het pairs: LoF het × damaging het, same gene, phase unknown",
          "2,372 pairs · hom-alt singles, internal > 0.3 → 6,294 · 8,666 rows"],
         "prioritize.py :: main"),
        ("Submission CSV",
         ["challenge schema, 12 columns · top 10 rows · chr prefix added",
          "epcr formatted to 4 decimals, clamped to (0, 1]"],
         "make_submission_csv.py"),
    ]
    for i, (title, lines, chip) in enumerate(stages):
        stage_box(ax, i, title, lines, chip)
    for i in range(6):
        flow_arrow(ax, ROW[i] - BOXH - 0.25, ROW[i + 1] + 0.25)

    # intermediate artifact on the flow line between ranking and CSV
    b6 = ROW[6] - BOXH
    ax.add_patch(FancyBboxPatch(
        (ART_X0, ART_TOP - ART_H), ART_X1 - ART_X0, ART_H,
        boxstyle="round,pad=0,rounding_size=0.6", facecolor=ART_FILL,
        edgecolor=ART_EDGE, lw=1.2, linestyle=(0, (3, 2)), zorder=2))
    at1 = ax.text((ART_X0 + ART_X1) / 2, ART_TOP - 1.15, "output/annot_full.tsv",
                  ha="center", va="center", fontsize=9.0, fontweight="bold",
                  color=BODY, family="monospace", zorder=3)
    _gate(at1, (ART_X0, ART_TOP - ART_H, ART_X1, ART_TOP), "artifact filename")
    at2 = ax.text((ART_X0 + ART_X1) / 2, ART_TOP - 2.35,
                  "ranked table · 16 columns · 8,666 rows",
                  ha="center", va="center", fontsize=8.4, color=MUTED, zorder=3)
    _gate(at2, (ART_X0, ART_TOP - ART_H, ART_X1, ART_TOP), "artifact subtitle")
    flow_arrow(ax, b6 - 0.25, ART_TOP + 0.25)
    flow_arrow(ax, ART_TOP - ART_H - 0.25, ROW[7] + 0.25)

    # footer: configuration + reference inputs (was the side panels)
    ax.add_patch(FancyBboxPatch(
        (X0, FOOT_Y0), X1 - X0, FOOT_Y1 - FOOT_Y0,
        boxstyle="round,pad=0,rounding_size=0.7", facecolor="#F5F5F5",
        edgecolor="#7A7A7A", lw=1.0, zorder=2))
    f0 = ax.text(X0 + 1.6, FOOT_Y1 - 1.35,
                 "Reference inputs via refs.* (pipeline/conf/track1_config.yml):",
                 ha="left", va="center", fontsize=8.6, fontweight="bold",
                 color=BODY, zorder=3)
    _gate(f0, (X0, FOOT_Y0, X1, FOOT_Y1), "footer head")
    f1 = ax.text(X0 + 1.6, FOOT_Y1 - 2.85,
                 "GRCh38 FASTA + FAI · GENCODE v47 MANE GTF · CADD v1.7 tiles · "
                 "gnomAD v4.1 exomes · ClinVar VCF · hp.obo + HPO map",
                 ha="left", va="center", fontsize=8.4, color=BODY, zorder=3)
    _gate(f1, (X0, FOOT_Y0, X1, FOOT_Y1), "footer refs")
    f2 = ax.text(X0 + 1.6, FOOT_Y1 - 4.35,
                 "dashed roundrect = intermediate artifact · chips = script :: function · "
                 "tool libraries in brackets",
                 ha="left", va="center", fontsize=8.4, color=MUTED, zorder=3)
    _gate(f2, (X0, FOOT_Y0, X1, FOOT_Y1), "footer legend")
    f3 = ax.text(X0 + 1.6, FOOT_Y1 - 5.85,
                 "rarity AF hom ≤ 0.01 · AF het ≤ 0.0001 (config rarity block)",
                 ha="left", va="center", fontsize=8.4, color=MUTED, zorder=3)
    _gate(f3, (X0, FOOT_Y0, X1, FOOT_Y1), "footer af")

    run_gate(fig, ax)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "F_pipeline_overview.png")
    fig.savefig(out, dpi=220, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
