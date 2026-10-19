"""T1-03: candidate construction, ordering, determinism.

Guards against the audit's candidate-construction findings:
  * reversed duplicate compound-het pairs (same two alleles both orders);
  * pairs and homozygous singles emitted in separate blocks with no
    specified global ordering;
  * homozygote count/write inconsistency (an internal >0.3 subset was
    counted while every single was written);
  * output must be deterministic byte-for-byte across runs and not depend
    on input row order (stable tie-breaking).
All coordinates invented; runs the real prioritize.main end to end.
"""
from __future__ import annotations

import re
import subprocess
import sys

from conftest import SCRIPT_DIR, write_conf, write_gtf, write_vcf

sys.path.insert(0, str(SCRIPT_DIR))
import prioritize  # noqa: E402

QC_ROW = 'GT:GQ:DP:AD\t0/1:99:40:20,20'


def run(tmp_path, rows, arg0='fixture.vcf'):
    gtf = write_gtf(tmp_path)
    vcf = write_vcf(tmp_path, rows)
    out = tmp_path / 'out.tsv'
    conf = write_conf(tmp_path, gtf, vcf)
    sys.argv = ['prioritize.py', '--conf', str(conf), '--out', str(out)]
    prioritize.main()
    return out.read_text().strip().split('\n')


def het(pos, chrom='1', ref='A', alt='G'):
    return f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\t{QC_ROW}'


def hom(pos, chrom='1', ref='A', alt='G'):
    return f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.\tGT:GQ:DP:AD\t1/1:99:40:10,30'


def body_pairs(lines):
    return [ln.split('\t') for ln in lines[1:] if ln.split('\t')[1] == 'compound_het']


def test_reversed_duplicate_pairs_deduplicated(tmp_path):
    """Three het donors + one het missense in GENEP: (x,y) and (y,x) are one
    candidate pair, deterministic orientation (position-sorted members)."""
    lines = run(tmp_path, [
        het(200), het(201), het(299),   # donor, donor, acceptor (all LoF)
        het(150),                       # missense (coding fallback)
    ])
    rows = body_pairs(lines)
    keys = {tuple(sorted(((c[3], c[4]), (c[7], c[8])))) for c in rows}
    # 3 donor-missense + 3 distinct donor-donor = 6 unique pairs
    assert len(rows) == len(keys) == 6, f'expected 6 unique pairs, got {len(rows)}/{len(keys)}: {rows}'
    for c in rows:
        # orientation: member 1 sorts <= member 2 by position
        assert (c[3], int(c[4]), c[5], c[6]) <= (c[7], int(c[8]), c[9], c[10])


def test_global_rank_not_separate_blocks(tmp_path):
    """A homozygous LoF single with internal > pair score must rank ABOVE the
    het pair: output is one global ordering by score, not pairs-then-singles."""
    lines = run(tmp_path, [
        het(200), het(150),             # pair scores 0.1880
        hom(111),                       # hom_alt coding LoF-ish internal 0.3582
    ])
    order = [ln.split('\t')[1] for ln in lines[1:]]
    assert order[0] == 'hom_alt', f'higher-scoring single must rank first, got {order}'


def test_tiebreak_independent_of_input_order(tmp_path):
    """Two equal-score pairs on chr1/chr2: regardless of listing order the
    chr1 pair sorts first (numeric chrom, then position)."""
    rows_fwd = [het(200), het(150), het(200, chrom='2', ref='A', alt='C'), het(150, chrom='2', ref='A', alt='C')]
    rows_rev = list(reversed(rows_fwd))
    a = run(tmp_path, rows_fwd)
    b = run(tmp_path, rows_rev)
    assert a == b, 'output depends on input VCF row order'
    first = a[1].split('\t')
    assert first[1] == 'compound_het' and int(first[3]) == 1, f'expected chr1 pair first, got {first}'


def test_run_to_run_byte_deterministic(tmp_path):
    rows = [het(200), het(201), het(150), hom(111), het(200, chrom='2')]
    assert run(tmp_path, rows) == run(tmp_path, rows)


def _run_capture_stderr(tmp_path, rows):
    gtf = write_gtf(tmp_path)
    vcf = write_vcf(tmp_path, rows)
    conf = write_conf(tmp_path, gtf, vcf)
    out = tmp_path / 'out.tsv'
    p = subprocess.run(
        [sys.executable, str(SCRIPT_DIR / 'prioritize.py'), '--conf', str(conf), '--out', str(out)],
        capture_output=True, text=True)
    return p, out.read_text()


def test_homozygote_semantics_single_source(tmp_path):
    """One counter semantics: instrument count == written hom_alt rows
    (no hidden >0.3 subset)."""
    p, text = _run_capture_stderr(tmp_path, [hom(111), hom(150)])
    m = re.search(r'hom_singles=(\d+)', p.stderr)
    written = [ln for ln in text.strip().split('\n')[1:] if ln.split('\t')[1] == 'hom_alt']
    assert m is None or int(m.group(1)) == len(written), \
        f'instrument hom_singles={m.group(1) if m else None} vs written {len(written)}'
    # and no phantom '>0.3 subset' counter survives
    assert 'hom_candidates' not in p.stderr


def test_noncanonical_alleles_counted_not_dropped_silently(tmp_path):
    """A record whose REF/ALT are not pure ACGTN is excluded from interesting
    sets by consequence 'other'; the exclusion must be instrumented."""
    p, text = _run_capture_stderr(tmp_path, [
        '1\t150\t.\tN\t<A1>\t99\tPASS\t.\tGT:GQ:DP:AD\t0/1:99:40:20,20'])
    m = re.search(r'excluded_other=(\d+)', p.stderr)
    assert m is not None and int(m.group(1)) >= 1, \
        f'non-base allele excluded without instrumentation: {p.stderr}'
