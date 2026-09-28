#!/usr/bin/env python3
"""Generate F_pipeline_overview.png - Track 1 pipeline overview diagram.

Deterministic, argv-free matplotlib generator. All numbers shown are quoted
from repo artifacts (REPORT.md, README.md, pipeline/conf/track1_config.yml,
pipeline/scripts/*.py). Output is written next to this script.

Self-gate: after drawing, every registered text artist is measured
(renderer window extent -> data coordinates via ax.transData.inverted())
against its container box. Any overflow raises SystemExit(1) BEFORE the
file is written, so a layout regression can never ship silently.

Run:  python3 docs/figures/build_fig_pipeline.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# ---------------------------------------------------------------- palette ---
STAGE_FILL, STAGE_EDGE, STAGE_HEAD = "#DFEAF5", "#3F6FA8", "#1F3B5C"
REF_FILL, REF_EDGE, REF_TEXT = "#FBF2DC", "#B58A2A", "#5C4A14"
CONF_FILL, CONF_EDGE, CONF_TEXT = "#E4F0DC", "#6E9A57", "#2F4A22"
ART_FILL, ART_EDGE = "#F4F4F4", "#8A8A8A"
CHIP_FILL, CHIP_EDGE = "#FFFFFF", "#8FB1D9"
MUTED = "#5A6B7D"

# canvas is 130 x 116 coordinate units (13 x 11.6 in at 150 dpi = 1950 x 1740 px)
W, H = 130.0, 116.0

# main pipeline column
X0, X1 = 26.0, 82.0
BOXH = 10.8
# explicit y-top of each stage box, 0-based top row first; the wide gap
# before the last row hosts the intermediate-artifact node
ROWS = [104.0, 92.0, 80.0, 68.0, 56.0, 44.0, 32.0, 14.2]

# chip band geometry: chip occupies b+0.8 .. b+2.3; body text must stay
# above b + CHIP_BODY_FLOOR so it can never strike through the chip
CHIP_Y0, CHIP_H = 0.8, 1.5
CHIP_BODY_FLOOR = 2.6

# intermediate artifact node between ranking (row 6) and CSV (row 7)
ART_X0, ART_X1 = 30.0, 78.0
ART_TOP, ART_H = 19.8, 4.4

# ------------------------------------------------------------------- gate ---
_GATE: list[tuple] = []
GATE_MARGIN = 0.2  # data units of slack allowed outside a container


def _gate(text, box, name):
    """Register a text artist for the post-draw overflow gate.

    box = (x0, y0, x1, y1) in data coordinates; the text's measured extent
    must stay inside the box expanded by GATE_MARGIN on every side.
    """
    _GATE.append((text, box, name))


def run_gate(fig, ax):
    """Measure every registered text; SystemExit(1) on any overflow.

    Window extents come from the Agg renderer and are converted to data
    coordinates with ax.transData.inverted(), so the comparison is exact
    in the coordinate system the boxes are laid out in. Runs BEFORE
    savefig: on failure no PNG is written.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    failures = []
    for text, (bx0, by0, bx1, by1), name in _GATE:
        bb = inv.transform_bbox(text.get_window_extent(renderer=renderer))
        if (bb.x0 < bx0 - GATE_MARGIN or bb.x1 > bx1 + GATE_MARGIN
                or bb.y0 < by0 - GATE_MARGIN or bb.y1 > by1 + GATE_MARGIN):
            failures.append(
                f"{name}: text extent ({bb.x0:.2f}, {bb.y0:.2f})-"
                f"({bb.x1:.2f}, {bb.y1:.2f}) escapes box "
                f"({bx0:.2f}, {by0:.2f})-({bx1:.2f}, {by1:.2f})")
    if failures:
        for f in failures:
            print(f"GATE FAIL: {f}", file=sys.stderr)
        raise SystemExit(1)
    print(f"GATE OK: {len(_GATE)} text artists within their boxes "
          f"(margin {GATE_MARGIN})")


# ------------------------------------------------------------------ layout ---
def row(i: int) -> tuple[float, float]:
    """(y_top, y_bottom) of pipeline stage box i (0-based, top row first)."""
    t = ROWS[i]
    return t, t - BOXH


def stage_box(ax, i, title, lines, chip=None, tool=None):
    """Stage box with a title line, body lines, and a dedicated chip band.

    Layout inside the box (height 10.8 units, 1 unit = 0.1 in):
      title + tool annotation  top ~3.7 units
      body text                next ~3.5 units (two 9.3 pt lines)
      chip band                bottom 2.3 units, reserved - never hosts body
    The overflow gate enforces both containment and body/chip separation:
    body text is registered against the box EXCLUDING the chip strip.
    """
    t, b = row(i)
    ax.add_patch(FancyBboxPatch(
        (X0, b), X1 - X0, t - b,
        boxstyle="round,pad=0,rounding_size=0.9",
        facecolor=STAGE_FILL, edgecolor=STAGE_EDGE, lw=1.4, zorder=2))
    ttl = ax.text(X0 + 1.6, t - 1.6, title, ha="left", va="top",
                  fontsize=11, fontweight="bold", color=STAGE_HEAD, zorder=3)
    _gate(ttl, (X0, b, X1, t), f"stage {i} title")
    if tool:
        tol = ax.text(X1 - 1.6, t - 1.6, tool, ha="right", va="top",
                      fontsize=9, style="italic", color=MUTED, zorder=3)
        _gate(tol, (X0, b, X1, t), f"stage {i} tool")
    body_bottom = b + CHIP_BODY_FLOOR if chip else b
    bod = ax.text(X0 + 1.6, t - 4.2, "\n".join(lines), ha="left", va="top",
                  fontsize=9.3, color="#24303C", linespacing=1.35, zorder=3)
    _gate(bod, (X0, body_bottom, X1, t), f"stage {i} body")
    if chip:
        # dedicated bottom strip: chip rect and its gated text
        cx0, cy0 = X0 + 1.4, b + CHIP_Y0
        chip_box = (cx0, cy0, cx0 + 27.0, cy0 + CHIP_H)
        ax.add_patch(FancyBboxPatch(
            (cx0, cy0), 27.0, CHIP_H,
            boxstyle="round,pad=0,rounding_size=0.5",
            facecolor=CHIP_FILL, edgecolor=CHIP_EDGE, lw=1.0, zorder=3))
        ct = ax.text(cx0 + 0.9, cy0 + CHIP_H / 2, chip, ha="left", va="center",
                     fontsize=8.6, color=STAGE_EDGE, zorder=4)
        _gate(ct, chip_box, f"stage {i} chip")


def flow_arrow(ax, y_from, y_to, x=None):
    """Vertical data-flow arrow between two y positions on the main column."""
    x = (X0 + X1) / 2 if x is None else x
    ax.add_patch(FancyArrowPatch(
        (x, y_from), (x, y_to),
        arrowstyle="-|>", mutation_scale=13, lw=1.6,
        color=STAGE_EDGE, shrinkA=0, shrinkB=0, zorder=1))


def ref_box(ax, y_top, h, lines, target_i, dy=0.0):
    x0, x1 = 88.0, 128.0
    ax.add_patch(FancyBboxPatch(
        (x0, y_top - h), x1 - x0, h,
        boxstyle="round,pad=0,rounding_size=0.7",
        facecolor=REF_FILL, edgecolor=REF_EDGE, lw=1.2, zorder=2))
    rt = ax.text(x0 + 1.4, y_top - 1.3, "\n".join(lines), ha="left", va="top",
                 fontsize=8.8, color=REF_TEXT, linespacing=1.35, zorder=3)
    _gate(rt, (x0, y_top - h, x1, y_top), f"ref box {lines[0][:18]!r}")
    # dashed feed arrow from reference box to the pipeline stage
    t, b = row(target_i)
    yc = (t + b) / 2 + dy
    ax.add_patch(FancyArrowPatch(
        (x0 - 0.2, (y_top - h / 2)), (X1 + 0.35, yc),
        arrowstyle="-|>", mutation_scale=11, lw=1.2, linestyle=(0, (4, 2.5)),
        color=REF_EDGE, shrinkA=1, shrinkB=0, zorder=1))


def main() -> None:
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "figure.dpi": 150,
        "savefig.dpi": 150,
    })
    fig = plt.figure(figsize=(13, 11.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    canvas_box = (0.0, 0.0, W, H)

    # ---- title -------------------------------------------------------------
    ttl = ax.text(W / 2, 112.5, "MVA Hackathon 2026 - Track 1 variant-prioritization pipeline",
                  ha="center", va="center", fontsize=14, fontweight="bold",
                  color="#16222E")
    _gate(ttl, canvas_box, "figure title")
    sub = ax.text(W / 2, 108.8, "patient WGS VCF  \u2192  ranked variant table  \u2192  submission CSV  "
                 "\u00b7  single pass, all joins against local reference files",
                 ha="center", va="center", fontsize=10, color=MUTED)
    _gate(sub, canvas_box, "figure subtitle")

    # ---- pipeline stages (top to bottom) -----------------------------------
    stage_box(
        ax, 0, "Proband WGS VCF",
        ["5,012,204 records \u00b7 GRCh38 \u00b7 bgzipped VCF \u00b7 bare contig names (data.vcf)"],
        chip=None, tool=None)
    stage_box(
        ax, 1, "Quality-control filter",
        ["FILTER PASS or absent \u00b7 GQ \u2265 20 \u00b7 DP \u2265 10 \u00b7 alt fraction \u2265 0.25",
         "GT het or hom-alt \u00b7 keeps 4,592,140 records"],
        chip="prioritize.py :: iter_qc_variants", tool="cyvcf2.VCF")
    stage_box(
        ax, 2, "Consequence annotation",
        ["MANE transcript overlap (bisect) \u00b7 splice classes \u00b7 codon-level effect",
         "from FASTA \u00b7 CADD PHRED lookup \u00b7 keeps 16,963 coding/splice candidates"],
        chip="annotate.py :: LocalAnnotator", tool="pysam.FastaFile \u00b7 pysam.TabixFile")
    stage_box(
        ax, 3, "Rarity and ClinVar join",
        ["gnomAD v4.1 exomes AF, exact (pos, ref, alt), sorted per-chromosome merge;",
         "af=None \u21d2 neutral, 0.0 \u21d2 ultrarare \u00b7 ClinVar: background evidence only"],
        chip="vcf.py :: GnomadAfJoin \u00b7 ClinVarJoin", tool="pysam.TabixFile")
    stage_box(
        ax, 4, "HPO phenotype term",
        ["gene score vs 8 proband terms: 1.0 ancestor-closed exact match,",
         "0.5 parent-only match, mean over proband terms"],
        chip="_hpo.py :: HPOScorer.gene_score", tool=None)
    stage_box(
        ax, 5, "Scoring (WorkPrior)",
        ["epcr = \u03a3 w\u1d62\u00b7term\u1d62 / \u03a3 w\u1d62, clamped to (0, 1] \u00b7 pair = stronger",
         "\u00d7 (0.5 + 0.5 \u00d7 weaker) \u00b7 weights: 0.45 / 0.60 / 0.30 / 2.0"],
        chip="prior.py :: WorkPrior.internal", tool=None)
    stage_box(
        ax, 6, "Candidate pairing and ranking",
        ["compound-het pairs: LoF het \u00d7 damaging het, same gene, phase unknown",
         "2,372 pairs \u00b7 hom-alt singles, internal > 0.3 \u2192 6,294 \u00b7 8,666 rows"],
        chip="prioritize.py :: main", tool=None)
    stage_box(
        ax, 7, "Submission CSV",
        ["challenge schema, 12 columns \u00b7 top 10 rows \u00b7 chr prefix added",
         "epcr formatted to 4 decimals, clamped to (0, 1]"],
        chip="make_submission_csv.py", tool=None)

    # arrows between consecutive stage boxes (rows 0-5 -> 1-6)
    for i in range(6):
        _, b = row(i)
        t2, _ = row(i + 1)
        flow_arrow(ax, b - 0.25, t2 + 0.25)

    # ---- intermediate artifact node (connected, on the flow line) -----------
    # ranking (row 6) writes the ranked TSV; make_submission_csv.py reads it
    _, b6 = row(6)
    ax.add_patch(FancyBboxPatch(
        (ART_X0, ART_TOP - ART_H), ART_X1 - ART_X0, ART_H,
        boxstyle="round,pad=0,rounding_size=0.6",
        facecolor=ART_FILL, edgecolor=ART_EDGE, lw=1.2,
        linestyle=(0, (3, 2)), zorder=2))
    at1 = ax.text((ART_X0 + ART_X1) / 2, ART_TOP - 1.5, "output/annot_full.tsv",
                  ha="center", va="center", fontsize=9.4, fontweight="bold",
                  color="#3A3A3A", family="monospace", zorder=3)
    _gate(at1, (ART_X0, ART_TOP - ART_H, ART_X1, ART_TOP), "artifact filename")
    at2 = ax.text((ART_X0 + ART_X1) / 2, ART_TOP - 3.0,
                  "ranked table \u00b7 16 columns \u00b7 8,666 rows",
                  ha="center", va="center", fontsize=8.8, color=MUTED, zorder=3)
    _gate(at2, (ART_X0, ART_TOP - ART_H, ART_X1, ART_TOP), "artifact subtitle")
    flow_arrow(ax, b6 - 0.25, ART_TOP + 0.25)          # ranking -> artifact
    flow_arrow(ax, ART_TOP - ART_H - 0.25, row(7)[0] + 0.25)  # artifact -> CSV

    # ---- reference inputs (right column) -------------------------------------
    ref_box(ax, 84.0, 8.0,
            ["GRCh38 no-alt-analysis-set FASTA (FAI-indexed)",
             "+ GENCODE v47 MANE GTF (tabix-indexed)"], 2)
    ref_box(ax, 74.0, 8.0,
            ["CADD v1.7 per-chromosome tiles (tabix)",
             "no site covered in shipped run \u2192 term neutral"], 2, dy=-2.2)
    ref_box(ax, 66.0, 7.5,
            ["gnomAD v4.1 exome sites,", "per-chromosome .vcf.bgz"], 3)
    ref_box(ax, 57.0, 7.0,
            ["ClinVar VCF (tabix-indexed)"], 3, dy=2.2)
    ref_box(ax, 52.0, 7.5,
            ["hp.obo is_a hierarchy +", "gene-to-phenotype map"], 4)

    # ---- configuration (left column) -----------------------------------------
    cx0, cx1, cy0, cy1 = 2.0, 22.0, 30.0, 88.0
    conf_box = (cx0, cy0, cx1, cy1)
    ax.add_patch(FancyBboxPatch(
        (cx0, cy0), cx1 - cx0, cy1 - cy0,
        boxstyle="round,pad=0,rounding_size=0.9",
        facecolor=CONF_FILL, edgecolor=CONF_EDGE, lw=1.4, zorder=2))
    ct1 = ax.text((cx0 + cx1) / 2, cy1 - 1.8, "pipeline/conf/", ha="center",
                  va="top", fontsize=9.6, fontweight="bold", color=CONF_TEXT,
                  family="monospace")
    _gate(ct1, conf_box, "config title line 1")
    ct2 = ax.text((cx0 + cx1) / 2, cy1 - 3.4, "track1_config.yml", ha="center",
                  va="top", fontsize=9.6, fontweight="bold", color=CONF_TEXT,
                  family="monospace")
    _gate(ct2, conf_box, "config title line 2")
    conf_lines = [
        "QC thresholds",
        "(GQ 20 / DP 10 / alt 0.25)",
        "",
        "rarity cutoffs",
        "0.01 / 0.0001",
        "",
        "LoF / missense /",
        "inframe term lists",
        "",
        "scoring weights",
        "0.45 0.60 0.30 2.0",
        "",
        "8 proband HPO terms",
        "",
        "reference paths",
        "(refs.*)",
    ]
    ctb = ax.text(cx0 + 1.5, cy1 - 6.2, "\n".join(conf_lines), ha="left",
                  va="top", fontsize=8.4, color=CONF_TEXT, linespacing=1.45,
                  zorder=3)
    _gate(ctb, conf_box, "config body")
    for target_i in (1, 4, 5, 6):
        t, b = row(target_i)
        yc = (t + b) / 2
        ax.add_patch(FancyArrowPatch(
            (cx1 + 0.15, min(max(yc, cy0 + 1), cy1 - 1)), (X0 - 0.35, yc),
            arrowstyle="-|>", mutation_scale=11, lw=1.2, linestyle=(0, (4, 2.5)),
            color=CONF_EDGE, shrinkA=1, shrinkB=0, zorder=1))

    # ---- legend ---------------------------------------------------------------
    lx, ly = X0 + 1.0, 0.9
    for dx, fc, ec, lbl in (
            (0.0, STAGE_FILL, STAGE_EDGE, "pipeline stage"),
            (14.5, REF_FILL, REF_EDGE, "local reference input"),
            (33.0, CONF_FILL, CONF_EDGE, "configuration"),
            (46.5, ART_FILL, ART_EDGE, "output artifact")):
        ax.add_patch(FancyBboxPatch(
            (lx + dx, ly), 2.6, 1.8, boxstyle="round,pad=0,rounding_size=0.3",
            facecolor=fc, edgecolor=ec, lw=1.1))
        lt = ax.text(lx + dx + 3.4, ly + 0.9, lbl, ha="left", va="center",
                     fontsize=8.0, color="#24303C")
        _gate(lt, canvas_box, f"legend label {lbl!r}")
    ax.add_patch(FancyArrowPatch(
        (lx + 60.5, ly + 0.9), (lx + 65.0, ly + 0.9), arrowstyle="-|>",
        mutation_scale=11, lw=1.6, color=STAGE_EDGE))
    lt = ax.text(lx + 66.0, ly + 0.9, "data flow", ha="left", va="center",
                 fontsize=8.0, color="#24303C")
    _gate(lt, canvas_box, "legend label 'data flow'")
    ax.add_patch(FancyArrowPatch(
        (lx + 74.5, ly + 0.9), (lx + 79.0, ly + 0.9), arrowstyle="-|>",
        mutation_scale=11, lw=1.2, linestyle=(0, (4, 2.5)), color=REF_EDGE))
    lt = ax.text(lx + 80.0, ly + 0.9, "ref / config feed", ha="left",
                 va="center", fontsize=8.0, color="#24303C")
    _gate(lt, canvas_box, "legend label 'ref / config feed'")

    # ---- self-gate, then write -------------------------------------------------
    run_gate(fig, ax)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "F_pipeline_overview.png")
    fig.savefig(out, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
