"""Tests for docs/check_pdf_layout.py (synthetic inputs only)."""

import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_pdf_layout as cl

MM = 72 / 25.4


def bbox(words_per_page):
    """Minimal pdftotext -bbox XHTML: list of pages, each a list of (xMax, text)."""
    pages = "".join(
        "<page>"
        + "".join(f'<word xMin="10" yMin="1" xMax="{x}" yMax="2">{t}</word>' for x, t in words)
        + "</page>"
        for words in words_per_page
    )
    return f"<doc>{pages}</doc>"


def test_text_edge_from_uniform_margin():
    assert cl.text_edge_pt("\\usepackage[a4paper,margin=17mm,top=16mm]{geometry}") == pytest.approx(
        595.28 - 17 * MM, abs=0.01
    )


def test_text_edge_right_overrides_margin():
    tex = "\\usepackage[a4paper,margin=17mm,right=20mm]{geometry}"
    assert cl.text_edge_pt(tex) == pytest.approx(595.28 - 20 * MM, abs=0.01)


def test_text_edge_supports_cm():
    assert cl.text_edge_pt("\\usepackage[margin=2.2cm]{geometry}") == pytest.approx(
        595.28 - 22 * MM, abs=0.01
    )


def test_text_edge_without_geometry_is_an_error():
    with pytest.raises(ValueError, match="geometry"):
        cl.text_edge_pt("% no geometry here")


def test_protrusion_within_tolerance_passes_beyond_fails():
    edge = 547.09
    pages = [[(edge + 2.5, "ok."), (edge + 3.0, "over")]]
    viol = cl.violations(bbox(pages), edge, cl.TOLERANCE_PT)
    assert len(viol) == 1 and "over" in viol[0] and viol[0].startswith("p1 ")


def test_ten_point_overflow_is_caught():
    edge = 547.09
    assert cl.violations(bbox([[(edge + 10.0, "wide")]]), edge, cl.TOLERANCE_PT)


def test_parse_floor_is_configurable():
    xml = bbox([[(100, "a")]])
    assert cl.parse_ok(xml, min_pages=1, min_words=1)
    assert not cl.parse_ok(xml, min_pages=5, min_words=1000)


def test_missing_pdftotext_gives_clear_message(tmp_path):
    res = subprocess.run(
        [sys.executable, str(HERE / "check_pdf_layout.py"), str(tmp_path / "x.pdf")],
        env={"PATH": str(tmp_path)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 2
    assert "pdftotext not found" in res.stderr
