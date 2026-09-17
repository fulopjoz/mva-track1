"""vcf.py  -  local reference-VCF joins for the Track 1 variant stream.
All sources are local files on the analysis host.

Two joins:
  1. gnomAD exomes v4.1 per-chromosome sites -> per-alt allele frequency.
     INFO field 'AF' (Number=A, joint callset AF) matched on (pos, ref, alt).
     Implemented as a same-chrom sorted stream merge: the patient QC stream is
     chromosome-ordered, each gnomAD sites file is pos-sorted, so a single
     forward cursor per chromosome resolves every query with O(1) memory.
     Chromosomes whose sites file is absent (still downloading, or not
     reference-covered) yield af=None, which WorkPrior treats as neutral  - 
     never falsely ultrarare.
  2. ClinVar VCF (tabix-indexed) -> CLNSIG + CLNDN for exact
     (chrom, pos, ref, alt) matches. Diagnostic only in this iteration:
     recorded on ScoredVariant.clinvar_sig / .clinvar_disease and used by
     analyses, deliberately NOT wired into WorkPrior weights yet.

Contig naming: patient = bare ('15'), gnomAD = 'chr15', ClinVar = bare.
All keys are normalised with _bare() so joins are convention-agnostic.
"""
from __future__ import annotations

import os

import cyvcf2
import pysam

GNOMAD_PATTERN = "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"


def _bare(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


# ---------------------------------------------------------------- gnomAD ----
class _ChromCursor:
    """Forward-only sorted merge cursor over one gnomAD sites VCF.
    Contract: `af()` must be called with non-decreasing `pos`."""

    def __init__(self, path: str):
        self._it = iter(cyvcf2.VCF(path))
        self._cache: dict[int, dict[tuple[str, str], float]] = {}
        self._hi = -1     # highest gnomAD position consumed so far
        self._eof = False
        self._last_q = -1  # last query position (order-violation detection)
        self.out_of_order = 0
        self._front = -1     # highest position the iterator has consumed
        self._pending = None # reserved; always None in the single-path design

    def af(self, pos: int, ref: str, alt: str) -> float:
        """AF for exact (pos, ref, alt); 0.0 = scanned past without a match.

        Completeness invariant (2026-09-07, second fix): a cached position only
        counts as fully consumed once the iterator has moved PAST it
        (pos < self._front). The earlier drain loop cached the look-ahead
        record, leaked partial first-line maps into the early-hit path, and
        made later sibling-line alleles answer absent. Now every record flows
        through one merge path and no partial position ever answers a query."""
        if pos < self._last_q:
            self.out_of_order += 1  # contract violated; cache keeps the answer
        self._last_q = max(self._last_q, pos)
        if pos < self._front:
            m = self._cache.get(pos)
            return m.get((ref, alt), 0.0) if m is not None else 0.0
        while not self._eof and self._front <= pos:
            if self._pending is not None:
                rec = self._pending
                self._pending = None
            else:
                rec = next(self._it, None)
                if rec is None:
                    self._eof = True
                    break
            afs = rec.INFO.get("AF")
            if afs is None:
                vals = []
            elif isinstance(afs, (list, tuple)):
                vals = [float(x) for x in afs]
            else:
                vals = [float(afs)]
            m = {(rec.REF, a): f for a, f in zip(rec.ALT, vals)}
            # gnomAD splits multiallelic sites into several lines at one
            # position; merge all lines at a position keyed (REF, ALT)
            self._cache.setdefault(rec.POS, {}).update(m)
            self._front = rec.POS
            if len(self._cache) > 500_000:
                self._cache = {q: x for q, x in self._cache.items() if q >= pos}
        m = self._cache.get(pos)
        if m is not None:
            return m.get((ref, alt), 0.0)
        return 0.0


class GnomadAfJoin:
    """gnomAD exomes v4.1 AF for patient stream positions, per chromosome."""

    def __init__(self, ref_dir: str, pattern: str = GNOMAD_PATTERN):
        self.ref_dir = os.path.expanduser(ref_dir)
        self.pattern = pattern
        self._cursors: dict[str, _ChromCursor | None] = {}

    def path_for(self, chrom: str) -> str:
        return os.path.join(self.ref_dir, self.pattern.format(chrom=_bare(chrom)))

    def has_chrom(self, chrom: str) -> bool:
        return os.path.exists(self.path_for(chrom))

    def out_of_order_queries(self) -> int:
        """Order-contract violations across all opened chromosomes (>0 means the
        input VCF was not pos-sorted per chromosome; results remain exact via
        the scan cache but the caller should investigate)."""
        return sum(c.out_of_order for c in self._cursors.values() if c is not None)

    def af(self, chrom: str, pos: int, ref: str, alt: str) -> float | None:
        """AF for exact (pos, ref, alt). None only when the chromosome's sites
        file is unavailable (neutral); 0.0 = scanned and absent (ultrarare)."""
        b = _bare(chrom)
        if b not in self._cursors:
            p = self.path_for(b)
            self._cursors[b] = _ChromCursor(p) if os.path.exists(p) else None
        cur = self._cursors[b]
        if cur is None:
            return None
        return cur.af(pos, ref, alt)


# ---------------------------------------------------------------- ClinVar ----
class ClinVarJoin:
    """Exact-allele ClinVar lookup via tabix. Returns (CLNSIG, CLNDN) or
    (None, None). Both returned raw, ';' / '|' separated as stored."""

    def __init__(self, path: str | None):
        p = os.path.expanduser(path) if path else None
        # encoding="utf-8": CLNDN disease names contain non-ASCII (e.g. é);
        # pysam tuples default to ASCII and UnicodeDecodeError otherwise.
        self._tb = (pysam.TabixFile(p, parser=pysam.asTuple(), encoding="utf-8")
                    if p and os.path.exists(p) else None)

    @staticmethod
    def _info_field(info: str, key: str) -> str | None:
        for kv in info.split(";"):
            if kv.startswith(key + "="):
                return kv[len(key) + 1:]
        return None

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> tuple[str | None, str | None]:
        if self._tb is None:
            return (None, None)
        b = _bare(chrom)
        try:
            rows = self._tb.fetch(b, pos - 1, pos)
        except (ValueError, OSError):
            try:
                rows = self._tb.fetch(f"chr{b}", pos - 1, pos)
            except (ValueError, OSError):
                return (None, None)
        for f in rows:  # asTuple -> VCF columns: 0=CHROM 1=POS 2=ID 3=REF 4=ALT … 7=INFO
            if int(f[1]) != pos or f[3] != ref:
                continue
            if alt in f[4].split(","):
                return (self._info_field(f[7], "CLNSIG"), self._info_field(f[7], "CLNDN"))
        return (None, None)


# ---------------------------------------------------------------- stream ----
def annotate_stream(rows, gnomad: GnomadAfJoin | None = None,
                    clinvar: ClinVarJoin | None = None):
    """Enrich (chrom, pos, ref, alt) rows with af / CLNSIG / CLNDN.
    `rows` must be non-decreasing in pos within each chromosome for the gnomAD
    sorted merge (clinvar is index-based and order-agnostic).
    Yields dict(chrom, pos, ref, alt, af, clinvar_sig, clinvar_disease).
    """
    for chrom, pos, ref, alt in rows:
        af = gnomad.af(chrom, pos, ref, alt) if gnomad else None
        sig, dn = clinvar.lookup(chrom, pos, ref, alt) if clinvar else (None, None)
        yield {
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "af": af, "clinvar_sig": sig, "clinvar_disease": dn,
        }
