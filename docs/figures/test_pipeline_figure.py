"""Tests for docs/figures/pipeline_figure.py (synthetic config only)."""

import copy
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import pipeline_figure as pf

CFG = {
    "variants": {
        "qc": {"min_gq": 20, "min_dp": 10, "min_alt_frac": 0.25},
        "af": {"max_rare": 0.01, "max_ultrarare": 0.0001},
    },
    "scoring": {"hpo_weight": 2.0, "pathogenicity_weight": 0.6, "rare_weight": 0.3, "founder_weight": 0.45},
}
STAGES = ["inputs", "qc", "annotation", "rarity", "inheritance", "scoring", "output"]
LABEL = re.compile(r'label="((?:[^"\\]|\\.)*)"')


def dot_labels(dot):
    return [m.replace("\\n", "\n").replace('\\"', '"') for m in LABEL.findall(dot)]


def test_dot_has_every_stage_and_lane_edge():
    dot = pf.build_dot(CFG)
    for s in STAGES:
        assert re.search(rf"^\s*{s}\s*\[", dot, re.M), s
    for i, (check, target) in enumerate(pf.LANE.items()):
        assert re.search(rf"check{i}\s*->\s*{target}\b[^\n]*style=dashed", dot), check


def test_config_value_drives_label():
    assert "GQ ≥ 20" in pf.labels(CFG)["qc"]
    cfg = copy.deepcopy(CFG)
    cfg["variants"]["qc"]["min_gq"] = 30
    assert "GQ ≥ 30" in pf.labels(cfg)["qc"]


@pytest.mark.parametrize("path", pf.KEYS)
def test_missing_key_names_path(path):
    cfg = copy.deepcopy(CFG)
    *parents, leaf = path.split(".")
    node = cfg
    for p in parents:
        node = node[p]
    del node[leaf]
    with pytest.raises(KeyError, match=re.escape(path)):
        pf.build_dot(cfg)


@pytest.mark.parametrize(
    "value,expected",
    [
        (20, "20"),
        (2.0, "2"),
        (0.25, "0.25"),
        (0.0001, "0.0001"),
        (1e-05, "0.00001"),
        ("1e-04", "0.0001"),  # PyYAML loads exponent scalars without a dot as strings
        (1e-08, "0.00000001"),
    ],
)
def test_number_formatting(value, expected):
    assert pf.fmt(value) == expected


def test_build_dot_is_deterministic():
    assert pf.build_dot(CFG) == pf.build_dot(CFG)


def test_every_dot_label_uses_only_allowlisted_tokens():
    allowed_numbers = {pf.fmt(pf.get_value(CFG, k)) for k in pf.KEYS} | pf.FIXED_NUMBERS
    labels = dot_labels(pf.build_dot(CFG))
    assert len(labels) >= len(STAGES) + len(pf.LANE) + 1  # stages, checks, cluster label
    for text in labels:
        for tok in re.findall(r"[A-Za-z][A-Za-z0-9.+_-]*|\d+(?:\.\d+)?", text):
            if re.fullmatch(r"\d+(?:\.\d+)?", tok):
                assert tok in allowed_numbers, tok
            else:
                assert tok in pf.VOCAB, tok


def test_lane_check_files_exist():
    for check in pf.LANE:
        assert (REPO / check).is_file(), check


# ---- Task 2: rendering ---------------------------------------------------------------------
import subprocess  # noqa: E402
import xml.etree.ElementTree as ET  # noqa: E402

import yaml  # noqa: E402

SIZE_ATTR = re.compile(r'size="[^"]*",\s*')


def _render(tmp_path, fmt, cfg=CFG, name="fig"):
    return pf.render(pf.build_dot(cfg), tmp_path / f"{name}.{fmt}", fmt)


def _pdf_width_pt(pdf):
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
    return float(re.search(r"Page size:\s+([\d.]+) x", info).group(1))


def test_svg_parses_and_has_no_raster(tmp_path):
    root = ET.parse(_render(tmp_path, "svg")).getroot()
    assert not [e for e in root.iter() if e.tag.endswith("image")]


def test_pdf_is_vector_with_embedded_non_type3_fonts(tmp_path):
    pdf = _render(tmp_path, "pdf")
    imgs = subprocess.run(["pdfimages", "-list", str(pdf)], capture_output=True, text=True, check=True)
    assert len(imgs.stdout.strip().splitlines()) <= 2  # header lines only
    fonts = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True, check=True).stdout
    assert "Type 3" not in fonts and "DejaVu" in fonts


def test_width_and_effective_font_budget(tmp_path):
    capped = _pdf_width_pt(_render(tmp_path, "pdf"))
    natural_dot = SIZE_ATTR.sub("", pf.build_dot(CFG), count=1)
    assert not re.search(r"\bsize=", natural_dot)
    natural = _pdf_width_pt(pf.render(natural_dot, tmp_path / "nat.pdf", "pdf"))
    assert capped <= 470
    scale = capped / natural if natural > capped else 1.0
    assert 10 * scale >= 7, (capped, natural)  # smallest font in the figure is 10 pt


def test_svg_render_is_deterministic(tmp_path):
    a = _render(tmp_path, "svg").read_bytes()
    b = _render(tmp_path, "svg", name="again").read_bytes()
    assert a == b


def _cfgfile(tmp_path, cfg):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_check_passes_on_fresh_and_fails_on_drift(tmp_path):
    svg = tmp_path / f"{pf.NAME}.svg"
    pf.render(pf.build_dot(CFG), svg, "svg")
    assert pf.main(["--check", "--config", _cfgfile(tmp_path, CFG), "--svg", str(svg)]) == 0
    drift = copy.deepcopy(CFG)
    drift["scoring"]["hpo_weight"] = 3.0
    assert pf.main(["--check", "--config", _cfgfile(tmp_path, drift), "--svg", str(svg)]) == 1


def test_check_fails_when_svg_missing(tmp_path):
    missing = str(tmp_path / "nope.svg")
    assert pf.main(["--check", "--config", _cfgfile(tmp_path, CFG), "--svg", missing]) == 1


def test_out_writes_svg_and_pdf(tmp_path):
    assert pf.main(["--out", str(tmp_path), "--config", _cfgfile(tmp_path, CFG)]) == 0
    assert (tmp_path / f"{pf.NAME}.svg").is_file() and (tmp_path / f"{pf.NAME}.pdf").is_file()


# ---- Task 3: build integration ----------------------------------------------------------------
import os  # noqa: E402
import shutil  # noqa: E402

FENCE = "`" * 3
REF = "![Pipeline overview.](F_pipeline_overview_v2.pdf){width=100%}"
BUILD = REPO / "docs/build_pdf.sh"


def test_count_refs_image_lines_mid_line_and_wrapped_caption():
    md = f"{REF}\nSee text ![inline](F_pipeline_overview_v2.pdf) here.\n![Wrapped\ncaption](F_pipeline_overview_v2.pdf)\n"
    assert pf.count_refs(md) == 3


def test_count_refs_skips_code_blocks_and_plain_mentions():
    md = f"{FENCE}\n{REF}\n{FENCE}\nplain mention of F_pipeline_overview_v2.pdf\n"
    assert pf.count_refs(md) == 0


@pytest.mark.parametrize("n,rc", [(0, 1), (1, 0), (2, 1)])
def test_count_refs_expect_gate(tmp_path, n, rc):
    md = tmp_path / "doc.md"
    md.write_text("intro\n" + "\n".join([REF] * n) + "\n")
    assert pf.main(["--count-refs", str(md), "--expect", "1"]) == rc


def test_build_script_wires_the_figure_gate():
    s = BUILD.read_text()
    assert "command -v dot" in s and "fc-match" in s
    assert "pipeline_figure.py --check" in s
    assert "--count-refs docs/TECHNICAL_DOCUMENTATION.md --expect 1" in s
    assert "pipeline_figure.py --out docs/figures" in s
    assert s.index("pipeline_figure.py --out") < s.index("\npandoc ")


def test_build_requires_dot(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for tool in ("fc-match", "grep", "python3", "pandoc"):
        found = shutil.which(tool)
        if found:
            os.symlink(found, bindir / tool)
    env = {"PATH": str(bindir), "HOME": str(tmp_path)}
    res = subprocess.run(["/bin/bash", str(BUILD)], cwd=REPO, env=env, capture_output=True, text=True)
    assert res.returncode != 0
    assert "graphviz 'dot' not found" in res.stderr


def test_rendered_pdf_is_gitignored():
    assert "docs/figures/F_pipeline_overview_v2.pdf" in (REPO / ".gitignore").read_text().splitlines()
