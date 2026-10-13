"""Calibrated positive predictive value (PPV) per call, fitted on held-out
genomes.
"""

import json
import numpy as np


# ----------------------------------------------------------------------
# Isotonic regression (pool-adjacent-violators), non-decreasing
# ----------------------------------------------------------------------
def pava(y, w=None):
    """Non-decreasing isotonic fit of y (optionally weighted). Returns fitted
    values aligned to the input order (which must be sorted by the predictor)."""
    y = np.asarray(y, float)
    w = np.ones_like(y) if w is None else np.asarray(w, float)
    val, wt, cnt = [], [], []
    for yi, wi in zip(y, w):
        cy, cw, cn = yi, wi, 1
        while val and val[-1] > cy:
            py, pw, pn = val.pop(), wt.pop(), cnt.pop()
            cy = (py * pw + cy * cw) / (pw + cw); cw += pw; cn += pn
        val.append(cy); wt.append(cw); cnt.append(cn)
    out = np.empty(len(y))
    i = 0
    for v, n in zip(val, cnt):
        out[i:i + n] = v; i += n
    return out


def fit_isotonic_ppv(scores, labels):
    """Fit monotone P(label=1 | score). Returns (score_knots, ppv_knots)
    describing a right-continuous step function usable by Calibration.
    """
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, float)
    n = max(len(labels), 1)
    eps = 1.0 / (2.0 * n)
    order = np.argsort(scores, kind="mergesort")
    s = scores[order]; y = labels[order]
    fit = np.clip(pava(y), eps, 1.0 - eps)
    # collapse to unique score knots (take the last fitted value per score)
    keep = np.concatenate([np.diff(s) > 0, [True]])   # last index of each run
    return s[keep], fit[keep]


# ----------------------------------------------------------------------
# Precision / recall summaries (reproduce the validation numbers)
# ----------------------------------------------------------------------
def precision_recall_at_threshold(scores, labels, thr):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    sel = scores >= thr
    tp = int((labels[sel] == 1).sum()); n_sel = int(sel.sum())
    total_true = int((labels == 1).sum())
    prec = tp / n_sel if n_sel else float("nan")
    rec = tp / total_true if total_true else float("nan")
    return prec, rec, n_sel, tp


def per_tier_precision(tiers, labels, tier_sets):
    """PPV and recall for named tiers. tier_sets maps name -> set of path
    labels selected (e.g. {'strict_core': {'A','C'}, 'path_B': {'A','B','C'}})."""
    tiers = np.asarray(tiers); labels = np.asarray(labels, int)
    total_true = int((labels == 1).sum())
    out = {}
    for name, paths in tier_sets.items():
        sel = np.isin(tiers, list(paths))
        n_sel = int(sel.sum()); tp = int((labels[sel] == 1).sum())
        out[name] = {"ppv": (tp / n_sel if n_sel else float("nan")),
                     "recall": (tp / total_true if total_true else float("nan")),
                     "n": n_sel, "tp": tp}
    return out


def logo_isotonic_ppv(scores, labels, genomes):
    """Leave-one-genome-out isotonic PPV: for each genome, fit on the others
    and predict its calls. Returns per-call held-out PPV (validates that the
    calibration generalises out of sample). Reliability: bin held-out PPV vs
    realised label rate."""
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    genomes = np.asarray(genomes)
    held = np.full(len(scores), np.nan)
    for g in np.unique(genomes):
        te = genomes == g; tr = ~te
        if tr.sum() < 10 or te.sum() == 0:
            continue
        sk, pk = fit_isotonic_ppv(scores[tr], labels[tr])
        idx = np.clip(np.searchsorted(sk, scores[te], "right") - 1, 0, pk.size - 1)
        held[te] = pk[idx]
    return held


# ----------------------------------------------------------------------
# Applied calibration (what the pipeline loads and calls per prediction)
# ----------------------------------------------------------------------
class Calibration:
    """Per-path calibration the pipeline loads and calls per prediction.
    """

    def __init__(self, per_path=None, ppv_iqr=None, domain=None, meta=None):
        self.per_path = per_path or {}      # path -> {score_knots, ppv_knots, pi0}
        self.ppv_iqr = ppv_iqr              # [lo, hi] PPV across calibration genomes
        self.domain = domain or {}          # {gc_min, gc_max, phyla:[...]}
        self.meta = dict(meta or {})

    @property
    def calibrated(self):
        return bool(self.per_path)

    def _entry(self, path):
        return self.per_path.get(path) or self.per_path.get("_all")

    def _ppv_true(self, path, score):
        c = self._entry(path)
        if not c or not c.get("score_knots"):
            return None
        sk = np.asarray(c["score_knots"], float); pk = np.asarray(c["ppv_knots"], float)
        idx = int(np.clip(np.searchsorted(sk, float(score), "right") - 1, 0, pk.size - 1))
        return float(pk[idx])

    def lr_pos(self, path, score):
        """Positive likelihood ratio (base-rate-free): posterior odds at this
        score / the fit set's prior odds. None if uncalibrated."""
        c = self._entry(path); p = self._ppv_true(path, score)
        if c is None or p is None:
            return None
        pi0 = c.get("pi0")
        if not pi0 or not (0.0 < pi0 < 1.0) or not (0.0 < p < 1.0):
            return None
        return (p / (1.0 - p)) / (pi0 / (1.0 - pi0))

    def ppv_from_prior(self, path, score, prior):
        """PPV against an explicitly stated prevalence prior. None if
        uncalibrated or no prior given."""
        lr = self.lr_pos(path, score)
        if lr is None or prior is None or not (0.0 < prior < 1.0):
            return None
        po = lr * (prior / (1.0 - prior))
        return po / (1.0 + po)

    def in_domain(self, gc=None):
        """Whether a genome's GC sits in the calibration's applicability domain
        (None if unknown / uncalibrated)."""
        if gc is None or "gc_min" not in self.domain:
            return None
        return bool(self.domain["gc_min"] <= gc <= self.domain["gc_max"])

    @classmethod
    def fit(cls, path_labels, scores, is_true, genomes=None, gc_values=None,
            phyla=None, **meta):
        """Fit per-path isotonic P(true|score) + prior pi0 from labeled
        calibration calls; ppv_iqr from per-genome precision; domain from
        GC/phyla. Used by calibrate_ppv.py."""
        path_labels = np.asarray(path_labels)
        scores = np.asarray(scores, float); is_true = np.asarray(is_true, int)
        per_path = {}
        for p in list(np.unique(path_labels)) + ["_all"]:
            m = np.ones(len(scores), bool) if p == "_all" else (path_labels == p)
            if m.sum() < 20:
                continue
            sk, pk = fit_isotonic_ppv(scores[m], is_true[m])
            per_path[str(p)] = {"score_knots": sk.tolist(), "ppv_knots": pk.tolist(),
                                "pi0": float(is_true[m].mean())}
        iqr = None
        if genomes is not None:
            precs = [float(is_true[np.asarray(genomes) == g].mean())
                     for g in np.unique(genomes) if (np.asarray(genomes) == g).any()]
            if precs:
                iqr = [float(np.quantile(precs, 0.25)), float(np.quantile(precs, 0.75))]
        domain = {}
        if gc_values is not None and len(gc_values):
            domain["gc_min"] = float(np.min(gc_values))
            domain["gc_max"] = float(np.max(gc_values))
        if phyla is not None:
            domain["phyla"] = sorted(set(map(str, phyla)))
        return cls(per_path, iqr, domain, meta)

    def to_json(self, path):
        with open(path, "w") as fh:
            json.dump({"per_path": self.per_path, "ppv_iqr": self.ppv_iqr,
                       "domain": self.domain, "meta": self.meta}, fh, indent=2)

    @classmethod
    def from_json(cls, path):
        with open(path) as fh:
            d = json.load(fh)
        return cls(d.get("per_path"), d.get("ppv_iqr"), d.get("domain"), d.get("meta"))
