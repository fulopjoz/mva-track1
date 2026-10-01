#!/usr/bin/env python3
"""Pipeline overview figure (Graphviz) for the Track-1 technical documentation.

Labels come from a fixed stage model plus parameter values read from
pipeline/conf/track1_config.yml at build time, so the figure always matches the
shipped configuration. Only the SVG is committed; docs/build_pdf.sh renders the
PDF used by pandoc.

Usage:
  pipeline_figure.py --out DIR                 render F_pipeline_overview_v2.{svg,pdf}
  pipeline_figure.py --check                   exit 1 if the committed SVG != a fresh render
  pipeline_figure.py --count-refs MD [--expect N]
                                               print image references to the figure;
                                               with --expect, exit 1 unless the count is N
"""

import argparse
import itertools
import re
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CONFIG = REPO / "pipeline/conf/track1_config.yml"
NAME = "F_pipeline_overview_v2"

KEYS = (
    "variants.qc.min_gq",
    "variants.qc.min_dp",
    "variants.qc.min_alt_frac",
    "variants.af.max_rare",
    "variants.af.max_ultrarare",
    "scoring.hpo_weight",
    "scoring.pathogenicity_weight",
    "scoring.rare_weight",
    "scoring.founder_weight",
)

# Okabe-Ito palette: (fill, border) per role
ROLE = {
    "input": ("#E8F4FA", "#0072B2"),
    "process": ("#F2F2F2", "#555555"),
    "decision": ("#FFF4E0", "#E69F00"),
    "output": ("#E6F5EF", "#009E73"),
    "verify": ("#F7E8F0", "#CC79A7"),
}
STAGE_ROLE = {
    "inputs": "input",
    "qc": "process",
    "annotation": "process",
    "rarity": "decision",
    "inheritance": "decision",
    "scoring": "decision",
    "output": "output",
}
# check file -> (guarded stage, stage the check node is levelled with, label)
LANE = {
    "checks/chk_cds_minus_strand.py": "annotation",
    "checks/check_submission_schema.py": "output",
}
LANE_TEXT = {
    "checks/chk_cds_minus_strand.py": "chk_cds_minus_strand.py\nminus-strand codon\nregression",
    "checks/check_submission_schema.py": "check_submission_schema.py\nsubmission CSV schema\n(+ tests)",
}

VOCAB = frozenset(
    [
        "Inputs",
        "Proband",
        "VCF",
        "GRCh38",
        "HPO",
        "phenotype",
        "terms",
        "References",
        "GENCODE",
        "v47",
        "MANE",
        "gnomAD",
        "v4.1",
        "exomes",
        "CADD",
        "v1.7",
        "ClinVar",
        "QC",
        "filter",
        "GQ",
        "DP",
        "alt",
        "fraction",
        "Annotation",
        "consequence",
        "on",
        "transcripts",
        "strand-aware",
        "CDS",
        "assembly",
        "when",
        "available",
        "Rarity",
        "population",
        "AF",
        "rare",
        "ultra-rare",
        "Inheritance",
        "compound-heterozygous",
        "pairs",
        "same",
        "gene",
        "homozygous",
        "candidates",
        "Scoring",
        "weights",
        "phenotype",
        "pathogenicity",
        "rarity",
        "founder",
        "Output",
        "ranked",
        "submission",
        "CSV",
        "rows",
        "Verification",
        "checks",
        "chk_cds_minus_strand.py",
        "check_submission_schema.py",
        "tests",
        "minus-strand",
        "codon",
        "regression",
        "schema",
    ]
)
FIXED_NUMBERS = frozenset({"10"})


def get_value(config, path):
    node = config
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"missing config key: {path}")
        node = node[part]
    return node


def fmt(v):
    """Plain decimal text for a config number (YAML may give int, float or str)."""
    if isinstance(v, bool):
        raise TypeError("boolean is not a number")
    v = float(v)  # ValueError for non-numeric strings
    if v.is_integer():
        return str(int(v))
    return format(Decimal(repr(v)), "f")


def labels(config):
    """Stage node id -> label text (config values substituted)."""
    v = {}
    for k in KEYS:
        try:
            v[k] = fmt(get_value(config, k))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"config key {k}: not a number ({exc})") from None
    return {
        "inputs": "Inputs\nProband VCF (GRCh38)\nHPO phenotype terms\n"
        "References: GENCODE v47 MANE\ngnomAD v4.1 exomes · CADD v1.7\nClinVar · HPO",
        "qc": f"QC filter\nGQ ≥ {v['variants.qc.min_gq']} · DP ≥ {v['variants.qc.min_dp']}\n"
        f"alt fraction ≥ {v['variants.qc.min_alt_frac']}",
        "annotation": "Annotation\nconsequence on MANE transcripts\nstrand-aware CDS assembly\nCADD when available",
        "rarity": f"Rarity (population AF)\nrare ≤ {v['variants.af.max_rare']}\n"
        f"ultra-rare ≤ {v['variants.af.max_ultrarare']}",
        "inheritance": "Inheritance\ncompound-heterozygous pairs\n(same gene)\nhomozygous candidates",
        "scoring": f"Scoring weights\nphenotype {v['scoring.hpo_weight']}"
        f" · pathogenicity {v['scoring.pathogenicity_weight']}\n"
        f"rarity {v['scoring.rare_weight']} · founder {v['scoring.founder_weight']}",
        "output": "Output\nranked candidates →\nsubmission CSV (≤ 10 rows)",
    }


def _q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'


def build_dot(config):
    lab = labels(config)
    order = list(STAGE_ROLE)
    out = [
        "digraph pipeline {",
        (
            '  graph [rankdir=TB, size="5.4,20", pad=0, nodesep=0.3, ranksep=0.28,'
            ' fontname="DejaVu Sans", fontsize=11];'
        ),
        (
            '  node [shape=box, style="rounded,filled", fontname="DejaVu Sans", fontsize=11,'
            ' margin="0.14,0.06", penwidth=1.3];'
        ),
        '  edge [color="#555555", penwidth=1.1, arrowsize=0.7];',
    ]
    for node, role in STAGE_ROLE.items():
        fill, border = ROLE[role]
        out.append(f'  {node} [label={_q(lab[node])}, fillcolor="{fill}", color="{border}"];')
    out += [f"  {a} -> {b};" for a, b in itertools.pairwise(order)]
    fill, border = ROLE["verify"]
    out.append(
        f'  subgraph cluster_verify {{ label="Verification checks"; style="dashed,rounded";'
        f' color="{border}"; fontcolor="{border}"; fontsize=10;'
    )
    for i, check in enumerate(LANE):
        out.append(
            f'    check{i} [label={_q(LANE_TEXT[check])}, fontsize=10, fillcolor="{fill}", color="{border}"];'
        )
    out.append("  }")
    for i, target in enumerate(LANE.values()):
        if order.index(target) == 0:
            raise ValueError(
                f"lane check cannot guard the first stage ({target!r}): nothing to level it with"
            )
        anchor = order[order.index(target) - 1]  # level the check with the stage it guards
        out.append(f"  {anchor} -> check{i} [style=invis];")
        out.append(f'  check{i} -> {target} [style=dashed, color="{border}", constraint=false];')
    out.append("}")
    return "\n".join(out) + "\n"


def render(dot, out, fmt):
    """Render DOT text with Graphviz to out (fmt: 'svg' or 'pdf')."""
    out = Path(out)
    subprocess.run(["dot", f"-T{fmt}", "-o", str(out)], input=dot.encode(), check=True)
    return out


FENCE = re.compile(r"^(```|~~~).*?^\1[^\n]*$", re.MULTILINE | re.DOTALL)
IMAGE_REF = re.compile(r"!\[[^\]]*\]\([^)]*F_pipeline_overview_v2\.pdf[^)]*\)")
XML_COMMENT = re.compile(rb"<!--.*?-->", re.DOTALL)


def _canonical(svg_bytes):
    """SVG bytes without XML comments (Graphviz writes its version and the graph title there)."""
    return XML_COMMENT.sub(b"", svg_bytes)


def count_refs(md_text):
    """Number of Markdown image references to the pipeline figure, outside fenced code."""
    return len(IMAGE_REF.findall(FENCE.sub("", md_text)))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(CONFIG))
    ap.add_argument("--svg", default=str(HERE / f"{NAME}.svg"))
    ap.add_argument("--expect", type=int)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out")
    g.add_argument("--check", action="store_true")
    g.add_argument("--count-refs")
    a = ap.parse_args(argv)
    if a.count_refs:
        n = count_refs(Path(a.count_refs).read_text())
        print(n)
        return 0 if a.expect is None or n == a.expect else 1
    if not Path(a.config).is_file():
        print(f"pipeline_figure: config not found: {a.config}", file=sys.stderr)
        return 1
    dot = build_dot(yaml.safe_load(Path(a.config).read_text()))
    if a.check:
        svg = Path(a.svg)
        if not svg.is_file():
            print("figure check: committed SVG missing", file=sys.stderr)
            return 1
        fresh = render(dot, svg.with_name(svg.stem + ".check.svg"), "svg")
        same = _canonical(fresh.read_bytes()) == _canonical(svg.read_bytes())
        fresh.unlink()
        print("figure check: " + ("ok" if same else "DRIFT (re-render and commit)"), file=sys.stderr)
        return 0 if same else 1
    out = Path(a.out)
    render(dot, out / f"{NAME}.svg", "svg")
    render(dot, out / f"{NAME}.pdf", "pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
