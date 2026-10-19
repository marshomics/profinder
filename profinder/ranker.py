"""
Shortlist ranking for promoter calls.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

_SPEC_PATH = Path(__file__).resolve().parent / "promoter_ranker.json"
_SPEC = None
_FEATURE_NAMES = ("neg_log_p10", "neg_log_p35_raw", "score_std", "spacer", "n_pairs",
                  "gc", "igr_len", "dist_to_cds", "m10_rank", "m35_rank",
                  "is_pathA", "is_pathB", "ext10", "is_co_r")


_WARNED = False


def spec():
    """Load and cache the shipped ranker definition (``promoter_ranker.json``
    next to this module)"""
    global _SPEC, _WARNED
    if _SPEC is None:
        if not _SPEC_PATH.is_file():
            if not _WARNED:
                print(f"  WARNING: ranker spec not found at {_SPEC_PATH}; "
                      f"shortlist_score will be blank and tables fall back to "
                      f"path-then-significance ordering.", file=sys.stderr)
                _WARNED = True
            return None
        with open(_SPEC_PATH, encoding="utf-8") as fh:
            loaded = json.load(fh)
        n = len(loaded["features"])
        if not (len(loaded["coef"]) == len(loaded["median_impute"]) == n):
            raise ValueError(f"ranker spec inconsistent: {n} features, "
                             f"{len(loaded['coef'])} coefficients, "
                             f"{len(loaded['median_impute'])} medians")
        unknown = set(loaded["features"]) - set(_FEATURE_NAMES)
        if unknown:
            raise ValueError(f"ranker spec names features this module does not "
                             f"compute: {sorted(unknown)}")
        _SPEC = loaded
    return _SPEC


def available() -> bool:
    """True when a ranker spec is loadable."""
    return spec() is not None


def _m_rank(v):
    """Architecture frequency rank from an M-id ('M007' -> 7). M-ids are assigned by
    descending promoter count, so the integer is comparable across catalogues."""
    s = str(v).strip()
    if len(s) > 1 and s[0] in "Mm" and s[1:].isdigit():
        return float(int(s[1:]))
    return float("nan")


def _p35_raw(rec):
    """The unadjusted -35 p-value.
    """
    v = rec.get("emp_p_35_uncorrected")
    try:
        f = float(v)
        if f == f and f > 0.0:
            return f
    except (TypeError, ValueError):
        pass
    try:
        adj = float(rec.get("emp_p_35"))
        k = float(rec.get("p35_multiplicity_k") or 1.0)
    except (TypeError, ValueError):
        return rec.get("emp_p_35")
    if adj != adj or k <= 1.0:
        return adj
    return -math.expm1(math.log1p(-min(max(adj, 0.0), 1.0 - 1e-15)) / k)


def _neg_log10(p):
    try:
        f = float(p)
    except (TypeError, ValueError):
        return float("nan")
    if f != f or f <= 0.0:
        return float("nan") if f != f else 300.0
    return -math.log10(max(f, 1e-300))


def features(rec, seq) -> dict:
    """Feature values for one call record.
    """
    def num(k):
        try:
            return float(rec.get(k))
        except (TypeError, ValueError):
            return float("nan")

    s = str(seq or "").upper()
    L = float(len(s)) if s else num("length")
    gc = (s.count("G") + s.count("C")) / len(s) if s else float("nan")
    p3 = num("pos_10") + 5.0
    path = str(rec.get("path", "")).strip()
    ext = rec.get("has_ext10")
    ext = 1.0 if (ext is True or str(ext).strip().lower() in ("true", "1", "yes")) else 0.0
    return {
        "neg_log_p10": _neg_log10(rec.get("emp_p_10")),
        "neg_log_p35_raw": _neg_log10(_p35_raw(rec)),
        "score_std": num("score_std"),
        "spacer": num("spacer"),
        "n_pairs": num("n_pairs"),
        "gc": gc,
        "igr_len": L,
        "dist_to_cds": L - p3,
        "m10_rank": _m_rank(rec.get("source_10")),
        "m35_rank": _m_rank(rec.get("source_35")),
        "is_pathA": 1.0 if path == "A" else 0.0,
        "is_pathB": 1.0 if path == "B" else 0.0,
        "ext10": ext,
        "is_co_r": 1.0 if str(rec.get("orientation", "")).strip() == "CO_R" else 0.0,
    }


def score(rec, seq) -> float:
    """Shortlist score in (0, 1); higher ranks first. Missing or non-finite features fall back
    to the fitted training median for that feature, so a call is never dropped from the ranking
    for want of one value."""
    sp = spec()
    if sp is None:
        return float("nan")
    f = features(rec, seq)
    z = float(sp["intercept"])
    for name, w, med in zip(sp["features"], sp["coef"], sp["median_impute"]):
        v = f.get(name, float("nan"))
        if v is None or v != v or v in (float("inf"), float("-inf")):
            v = float(med)
        z += float(w) * float(v)
    return 1.0 / (1.0 + math.exp(-max(min(z, 60.0), -60.0)))
