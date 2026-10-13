"""Interval-level empirical-FDR statistical core.

One pooled construction per (length x dinucleotide-composition) stratum:
  - AE dinucleotide shuffles of the stratum's intervals, scored identically
    to the observed intervals (same motif stack, same strands);
  - pooled per-(motif, element) null score arrays -> element thresholds at
    p10 / p35 / p35_relaxed AND the within-motif empirical-p (minP) function;
  - the interval statistic (best combined -log10 minP over the spacer-linked
    search) evaluated on the shuffles -> the interval-null;
  - observed interval statistic vs that null -> add-one empirical p;
  - Benjamini-Hochberg across intervals -> q; optional scipy GPD tail.
"""

import math
import numpy as np
from .ae_shuffle import ae_shuffle
from . import fdr_core as fc

try:
    from scipy.stats import genpareto
    _HAVE_SCIPY = True
except Exception:                    # pragma: no cover
    _HAVE_SCIPY = False

_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


# ----------------------------------------------------------------------
# Stratification: length x dinucleotide composition
# ----------------------------------------------------------------------
def dinuc_freq(seq):
    v = np.zeros(16)
    n = 0
    for i in range(len(seq) - 1):
        a = _IDX.get(seq[i]); b = _IDX.get(seq[i + 1])
        if a is not None and b is not None:
            v[a * 4 + b] += 1; n += 1
    if n:
        v /= n
    return v


def _kmeans(X, k, seed=0, iters=30):
    rng = np.random.default_rng(seed)
    if len(X) <= k:
        return np.arange(len(X)) % max(1, k)
    cent = X[rng.choice(len(X), k, replace=False)].copy()
    lab = np.zeros(len(X), dtype=int)
    for _ in range(iters):
        d = ((X[:, None, :] - cent[None, :, :]) ** 2).sum(2)
        newlab = d.argmin(1)
        if (newlab == lab).all():
            break
        lab = newlab
        for j in range(k):
            if (lab == j).any():
                cent[j] = X[lab == j].mean(0)
    return lab


def stratify(seqs, n_len_bins=4, n_comp=8, min_size=300, seed=0):
    """Assign each interval a stratum label (length bin x composition cluster);
    strata smaller than min_size are merged into a shared catch-all (-1)."""
    lens = np.array([len(s) for s in seqs])
    edges = np.quantile(lens, np.linspace(0, 1, n_len_bins + 1))
    lbin = np.clip(np.digitize(lens, edges[1:-1]), 0, n_len_bins - 1)
    X = np.vstack([dinuc_freq(s) for s in seqs])
    cbin = _kmeans(X, n_comp, seed)
    strat = lbin * n_comp + cbin
    labels, counts = np.unique(strat, return_counts=True)
    small = set(labels[counts < min_size].tolist())
    return np.array([-1 if s in small else int(s) for s in strat])


# ----------------------------------------------------------------------
# Scoring / search
# ----------------------------------------------------------------------


def _emp_p_row(scores_row, null_sorted):
    """add-one empirical p for each score in scores_row vs ascending null."""
    cg = null_sorted.size - np.searchsorted(null_sorted, scores_row, side="left")
    return (cg + 1.0) / (null_sorted.size + 1.0)


def sidak(p, k):
   
    if k is None or k <= 1:
        return float(p)
    p = min(max(float(p), 0.0), 1.0)
    if p < 1e-12:                       # avoid catastrophic cancellation in 1-(1-p)^k
        return float(min(1.0, p * k))
    return float(-np.expm1(k * np.log1p(-p)))


class StratumNull:
    """Pooled per-motif null for one stratum: thresholds + emp-p + interval null."""

    def __init__(self, pwm10, pwm35, p10, p35, p35_relaxed,
                 spacer=(16, 18), cap=200_000, seed=0,
                 max_dist_to_cds_start=None, subs=None):
        self.pwm10, self.pwm35 = pwm10, pwm35
        self.p10, self.p35, self.p35r = p10, p35, p35_relaxed
        self.spacer = spacer
        self.cap = cap
        self.max_dist = max_dist_to_cds_start   # -10 right edge within this
        self.rng = np.random.default_rng(seed)  # many bp of the 3' end
        self.S10 = pwm10.shape[0]
        self.S35 = pwm35.shape[0]
        self.subs = subs                        # subgroup names, for source_10/35

    # -- pass 1: pool position scores from shuffles, derive thr + emp-p --
    def build(self, shuffle_strand_codes):
        pool10 = [[] for _ in range(self.S10)]
        pool35 = [[] for _ in range(self.S35)]
        for codes in shuffle_strand_codes:            # one shuffled interval
            for code in codes:                        # its strand(s)
                s10 = fc.score_stack(code, self.pwm10)
                s35 = fc.score_stack(code, self.pwm35)
                for m in range(self.S10):
                    r = s10[m][np.isfinite(s10[m])]
                    if r.size:
                        pool10[m].append(r)
                for m in range(self.S35):
                    r = s35[m][np.isfinite(s35[m])]
                    if r.size:
                        pool35[m].append(r)
        self.null10 = [self._finalize(pool10[m]) for m in range(self.S10)]
        self.null35 = [self._finalize(pool35[m]) for m in range(self.S35)]
        self.set_thresholds(self.p10, self.p35, self.p35r)
        return self

    def set_thresholds(self, p10, p35, p35_relaxed):
        """(Re)derive the element thresholds from the already-built null at
        new target levels. Cheap (just quantiles of the stored null score
        arrays), so a p-value sweep can reuse one built null across many
        threshold points without re-shuffling."""
        self.p10, self.p35, self.p35r = p10, p35, p35_relaxed
        self.thr10 = np.array([np.quantile(a, 1 - p10) for a in self.null10])
        self.thr35s = np.array([np.quantile(a, 1 - p35) for a in self.null35])
        self.thr35r = np.array([np.quantile(a, 1 - p35_relaxed) for a in self.null35])
        return self

    def _finalize(self, chunks):
        if not chunks:
            return np.array([-np.inf], dtype=np.float32)
        a = np.concatenate(chunks)
        if a.size > self.cap:                         # bound memory; subsample
            a = self.rng.choice(a, self.cap, replace=False)
        a.sort()
        return a

    # -- interval statistic: best combined -log10 minP over spacer-linked pairs --
    def interval_T(self, strand_codes):
        best = -np.inf
        for code in strand_codes:
            s10 = fc.score_stack(code, self.pwm10)
            s35 = fc.score_stack(code, self.pwm35)
            P = s10.shape[1]
            if P <= 0:
                continue
            minp10 = np.ones(P); minp35 = np.ones(P)
            for m in range(self.S10):
                ep = _emp_p_row(s10[m], self.null10[m])
                ep = np.where(s10[m] >= self.thr10[m], ep, 1.0)
                np.minimum(minp10, ep, out=minp10)
            for m in range(self.S35):
                ep = _emp_p_row(s35[m], self.null35[m])
                ep = np.where(s35[m] >= self.thr35s[m], ep, 1.0)
                np.minimum(minp35, ep, out=minp35)
            # distance-to-CDS filter: -10 right edge within max_dist of the
            # 3' (gene-proximal) end of the already-oriented sequence.
            min_pos_10 = 0
            if self.max_dist is not None:
                min_pos_10 = max(0, len(code) - fc._MOTIF_WIDTH - int(self.max_dist))
            hits10 = np.where(minp10 < 1.0)[0]
            hits10 = hits10[hits10 >= min_pos_10]
            for i in hits10:
                for sp in range(self.spacer[0], self.spacer[1] + 1):
                    j = i - fc._MOTIF_WIDTH - sp
                    if j < 0:
                        continue
                    if minp35[j] < 1.0:
                        comb = -math.log10(minp10[i]) - math.log10(minp35[j])
                        if comb > best:
                            best = comb
        return best

    def interval_null(self, shuffle_strand_codes):
        self.nullT = np.sort(np.array(
            [self.interval_T(c) for c in shuffle_strand_codes], dtype=float))
        # optional GPD tail fit on the upper exceedances
        self._fit_tail()
        return self.nullT

    def _fit_tail(self):
        self.gpd = None
        if not _HAVE_SCIPY or self.nullT.size < 500:
            return
        finite = self.nullT[np.isfinite(self.nullT)]
        if finite.size < 500:
            return
        u = np.quantile(finite, 0.90)             # threshold for POT fit
        exc = finite[finite > u] - u
        if exc.size < 50:
            return
        try:
            c, loc, scale = genpareto.fit(exc, floc=0)
            self.gpd = (u, c, scale, finite.size, exc.size)
        except Exception:
            self.gpd = None

    def emp_p_interval(self, T):
        n = self.nullT.size
        cg = n - np.searchsorted(self.nullT, T, side="left")
        p = (cg + 1.0) / (n + 1.0)
        # extend into the tail with the GPD fit when at/near the empirical floor
        if self.gpd is not None and cg <= 1 and np.isfinite(T):
            u, c, scale, ntot, nexc = self.gpd
            if T > u:
                tail = (nexc / ntot) * genpareto.sf(T - u, c, loc=0, scale=scale)
                p = max(min(p, tail), 1.0 / (n + 1.0) / 1e6)  # floor, avoid 0
        return p

    # -- length normalisation: continuous statistic + null-location fit --
    def continuous_stat(self, strand_codes):
        """Continuous best combined -log10 minP over ALL valid pairs (no
        element-threshold gate) plus the number of pairs tested. This is the
        length-normalisable score; path A/B/C classification still uses the
        element thresholds separately. Unlike the threshold-gated statistic it
        does not saturate, so it exposes (and lets us remove) the length bias."""
        best = -np.inf; ntest = 0
        w = fc._MOTIF_WIDTH
        for code in strand_codes:
            s10 = fc.score_stack(code, self.pwm10)
            s35 = fc.score_stack(code, self.pwm35)
            P = s10.shape[1]
            if P <= 0:
                continue
            mp10 = np.ones(P); mp35 = np.ones(P)
            for m in range(self.S10):
                np.minimum(mp10, _emp_p_row(s10[m], self.null10[m]), out=mp10)
            for m in range(self.S35):
                np.minimum(mp35, _emp_p_row(s35[m], self.null35[m]), out=mp35)
            nl10 = -np.log10(mp10); nl35 = -np.log10(mp35)
            lo = 0
            if self.max_dist is not None:
                lo = max(0, len(code) - w - int(self.max_dist))
            for sp in range(self.spacer[0], self.spacer[1] + 1):
                i0 = max(lo, w + sp)
                if i0 > P - 1:
                    continue
                i = np.arange(i0, P)
                comb = nl10[i] + nl35[i - w - sp]
                ntest += i.size
                if comb.size:
                    best = max(best, float(comb.max()))
        return best, ntest

    def fit_length_location(self, shuffle_strand_codes):
        """Fit the null statistic's LOCATION on log(n_pairs): T ~ a + b*log n
        (least squares over finite null draws). Stores (a, b) for standardize."""
        Ts, logn = [], []
        for codes in shuffle_strand_codes:
            T, n = self.continuous_stat(codes)
            if np.isfinite(T) and n > 0:
                Ts.append(T); logn.append(math.log(n))
        if len(Ts) >= 10:
            b, a = np.polyfit(np.asarray(logn), np.asarray(Ts), 1)
        else:
            a, b = 0.0, 0.0
        self._loc_a, self._loc_b = float(a), float(b)
        return self._loc_a, self._loc_b

    def standardized_stat(self, strand_codes):
        """Length-standardised score: continuous T minus the fitted null
        location a + b*log(n_pairs). Threshold and calibrate on THIS."""
        T, n = self.continuous_stat(strand_codes)
        if not np.isfinite(T) or n <= 0:
            return -np.inf, n
        a = getattr(self, "_loc_a", 0.0); b = getattr(self, "_loc_b", 0.0)
        return T - (a + b * math.log(n)), n

    def classify_and_score(self, seq, spacer=None):
        """Full path A/B/C classification on the empirical element thresholds
        (A: linked strict -35 same subgroup; B: extended -10 + linked relaxed
        -35; C: unlinked strict -35), plus the length-standardised continuous
        score. Returns a call record, or None if there is no A/B/C hit."""
        sp0, sp1 = spacer or self.spacer
        code = fc.encode(seq)
        s10 = fc.score_stack(code, self.pwm10)
        s35 = fc.score_stack(code, self.pwm35)
        P = s10.shape[1]
        if P <= 0:
            return None
        ep10 = empp_stack(s10, self.null10); ep35 = empp_stack(s35, self.null35)
        ep10m = np.where(s10 >= self.thr10[:, None], ep10, np.inf)
        best10_p = ep10m.min(0); best10_sub = ep10m.argmin(0)
        ep35sm = np.where(s35 >= self.thr35s[:, None], ep35, np.inf)
        best35s_p = ep35sm.min(0); best35s_sub = ep35sm.argmin(0)
        m35rhit = s35 >= self.thr35r[:, None]
        w = fc._MOTIF_WIDTH
        lo = 0 if self.max_dist is None else max(0, len(code) - w - int(self.max_dist))
        name = lambda k: self.subs[k] if self.subs is not None else str(k)
        best = None
        for i in np.where(np.isfinite(best10_p))[0]:
            if i < lo:
                continue
            sub10 = int(best10_sub[i]); p10 = best10_p[i]
            ext10 = i >= 3 and seq[i - 3:i - 1] == "TG"
            bA = bB = bC = None
            for spc in range(sp0, sp1 + 1):
                j = i - w - spc
                if j < 0:
                    continue
                if np.isfinite(best35s_p[j]):
                    if int(best35s_sub[j]) == sub10:
                        if bA is None or best35s_p[j] < bA[0]:
                            bA = (best35s_p[j], j, spc, int(best35s_sub[j]))
                    elif bC is None or best35s_p[j] < bC[0]:
                        bC = (best35s_p[j], j, spc, int(best35s_sub[j]))
                if m35rhit[sub10, j] and (bB is None or ep35[sub10, j] < bB[0]):
                    bB = (ep35[sub10, j], j, spc, sub10)
            if bA is not None:
                path, ch = "A", bA
            elif ext10 and bB is not None:
                path, ch = "B", bB
            elif bC is not None:
                path, ch = "C", bC
            else:
                continue
            k35 = (self.S35 - 1) if path == "C" else 1
            p35_raw = float(ch[0])
            p35_adj = sidak(p35_raw, k35)
            comb = -math.log10(p10) - math.log10(max(p35_adj, 1e-300))
            key = ({"A": 0, "B": 1, "C": 2}[path], -comb)
            if best is None or key < best[0]:
                best = (key, {"path": path, "emp_p_10": float(p10),
                              "emp_p_35": p35_adj,
                              "emp_p_35_uncorrected": p35_raw,
                              "p35_multiplicity_k": int(k35),
                              "pos_10": int(i),
                              "pos_35": int(ch[1]), "spacer": int(ch[2]),
                              "source_10": name(sub10), "source_35": name(ch[3]),
                              "has_ext10": bool(ext10)})
        if best is None:
            return None
        rec = best[1]
        sT, n = self.standardized_stat([code])
        rec["score_std"] = float(sT); rec["n_pairs"] = int(n)
        return rec


def empp_stack(scores, nulls):
    """(S, P) empirical p per motif per position (each row vs its motif null)."""
    return np.vstack([_emp_p_row(scores[m], nulls[m]) for m in range(scores.shape[0])])


def benjamini_hochberg(pvals):
    """BH q-values (monotone) for a 1-D array of p-values."""
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.clip(q_sorted, 0, 1)
    return q
