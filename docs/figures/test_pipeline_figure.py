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
