"""Tests for checks/check_submission_schema.py (synthetic rows only).

Schema source: the challenge's public "Submit - Track 1" tab (12 exact fields,
up to 10 rows, GRCh38 `chrN` contigs, epcr in (0, 1], finding_type
primary|secondary; second-variant fields blank unless compound-het).
"""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from check_submission_schema import HEADER, validate

SCRIPT = Path(__file__).parent / "check_submission_schema.py"
PAIR = ["PROBAND01", "chr1", "100000", "A", "T", "chr1", "200000", "G", "C", "0.90", "primary", ""]
SINGLE = ["PROBAND01", "chr7", "300000", "C", "T", "", "", "", "", "0.30", "secondary", "Incidental finding"]


def write(tmp_path, rows, header=HEADER):
    p = tmp_path / "sub.csv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return p


def with_field(row, col, value):
    r = list(row)
    r[HEADER.index(col)] = value
    return r


def test_valid_file_has_no_errors(tmp_path):
    assert validate(write(tmp_path, [PAIR, SINGLE])) == []


def test_cli_prints_all_checks_passed_and_exits_zero(tmp_path):
    res = subprocess.run(
        [sys.executable, SCRIPT, write(tmp_path, [PAIR, SINGLE])], capture_output=True, text=True, check=False
    )
    assert res.returncode == 0
    assert res.stdout.strip().splitlines()[-1] == "ALL CHECKS PASSED"


def test_cli_exits_nonzero_on_error(tmp_path):
    bad = write(tmp_path, [with_field(PAIR, "epcr", "1.5")])
    res = subprocess.run([sys.executable, SCRIPT, bad], capture_output=True, text=True, check=False)
    assert res.returncode == 1
    assert "ALL CHECKS PASSED" not in res.stdout


def test_wrong_header_is_rejected(tmp_path):
    hdr = list(HEADER)
    hdr[0], hdr[1] = hdr[1], hdr[0]
    assert any("header" in e for e in validate(write(tmp_path, [PAIR], header=hdr)))


@pytest.mark.parametrize("n", [0, 11])
def test_row_count_outside_1_to_10_is_rejected(tmp_path, n):
    rows = [with_field(PAIR, "pos_1", str(100000 + i)) for i in range(n)]
    assert any("row count" in e for e in validate(write(tmp_path, rows)))


def test_row_with_wrong_field_count_is_rejected(tmp_path):
    assert any("field count" in e for e in validate(write(tmp_path, [PAIR[:-1]])))


@pytest.mark.parametrize(
    "col,value",
    [
        ("chrom_1", "15"),
        ("chrom_1", "chr23"),
        ("chrom_1", ""),
        ("pos_1", "0"),
        ("pos_1", "-5"),
        ("pos_1", "12.5"),
        ("pos_1", "abc"),
        ("ref_1", "Z"),
        ("ref_1", ""),
        ("alt_1", "a"),
        ("epcr", "0"),
        ("epcr", "1.5"),
        ("epcr", "nan"),
        ("epcr", "high"),
        ("finding_type", "incidental"),
        ("finding_type", ""),
        ("proband_id", ""),
    ],
)
def test_invalid_field_value_is_rejected(tmp_path, col, value):
    errs = validate(write(tmp_path, [with_field(PAIR, col, value)]))
    assert any(col in e for e in errs), errs


def test_ref_equal_to_alt_is_rejected(tmp_path):
    assert any("ref_1" in e for e in validate(write(tmp_path, [with_field(PAIR, "alt_1", "A")])))


def test_partially_filled_second_variant_is_rejected(tmp_path):
    row = with_field(SINGLE, "pos_2", "400000")  # pos_2 set, rest of variant 2 blank
    assert any("second variant" in e for e in validate(write(tmp_path, [row])))


def test_mixed_proband_ids_are_rejected(tmp_path):
    other = with_field(SINGLE, "proband_id", "PROBAND02")
    assert any("proband_id" in e for e in validate(write(tmp_path, [PAIR, other])))


def test_duplicate_rows_are_rejected(tmp_path):
    assert any("duplicate" in e for e in validate(write(tmp_path, [PAIR, PAIR])))


def test_epcr_increasing_down_the_file_is_rejected(tmp_path):
    assert any("order" in e for e in validate(write(tmp_path, [SINGLE, PAIR])))  # 0.30 then 0.90


def test_equal_epcr_values_are_allowed(tmp_path):
    tie = with_field(with_field(SINGLE, "epcr", "0.90"), "pos_1", "300001")
    assert validate(write(tmp_path, [PAIR, tie])) == []


def test_error_messages_never_echo_cell_values(tmp_path):
    secret = "ZZSECRETZZ"
    rows = [with_field(PAIR, "chrom_1", secret), with_field(SINGLE, "finding_type", secret)]
    res = subprocess.run(
        [sys.executable, SCRIPT, write(tmp_path, rows)], capture_output=True, text=True, check=False
    )
    assert res.returncode == 1
    assert secret not in res.stdout + res.stderr
