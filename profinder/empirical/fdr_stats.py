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
    """Šidák-adjust a p-value that is the minimum over k comparisons:
    P(min p <= a) = 1 - (1 - a)^k. Used for Path C, whose -35 is the best
    match among all subgroups other than the one supplying the -10; Paths A
    and B fix the subgroup from the -10 and are reported unadjusted."""
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
        """Full path A/B/C classification on the empirical element thresholds,
        plus the length-standardised continuous score. Returns a call record,
        or None if there is no A/B/C hit.

        At every -10 position i and spacer-linked -35 position j = i - 6 - spacer:

        * **Path A** -- some subgroup s passes the strict -10 threshold at i AND
          its own strict -35 threshold at j (a linked pair). Every -10-passing
          subgroup is tested, not only the single most significant one: many
          catalogue subgroups share a TATAAT-like -10 and differ in the -35, so
          testing only the sharpest -10 subgroup mis-routed canonical linked
          promoters (e.g. TTGACA-N17-TATAAT) to Path C.
        * **Path B** -- i carries an extended -10 ("TGn" 5' of the hexamer) and
          some -10-passing subgroup has a relaxed -35 at j.
        * **Path C** -- the most significant -10 subgroup at i (ties on the
          empirical quantile broken on raw log-odds) with the best strict -35
          from any OTHER subgroup at j; emp_p_35 is Sidak-adjusted for the
          S35-1 subgroups searched.

        Within a path the (subgroup, spacer) pair with the largest combined
        -log10(p10) - log10(p35) is kept; across positions the call minimising
        (path rank A<B<C, -combined) is reported."""
        sp0, sp1 = spacer or self.spacer
        code = fc.encode(seq)
        s10 = fc.score_stack(code, self.pwm10)
        s35 = fc.score_stack(code, self.pwm35)
        P = s10.shape[1]
        if P <= 0:
            return None
        w = fc._MOTIF_WIDTH
        ep10 = empp_stack(s10, self.null10); ep35 = empp_stack(s35, self.null35)
        nl10 = -np.log10(ep10); nl35 = -np.log10(ep35)
        pass10 = s10 >= self.thr10[:, None]          # (S, P) strict -10 pass
        m35s = s35 >= self.thr35s[:, None]           # strict -35 pass
        m35r = s35 >= self.thr35r[:, None]           # relaxed -35 pass
        lo = 0 if self.max_dist is None else max(0, len(code) - w - int(self.max_dist))
        name = lambda k: self.subs[k] if self.subs is not None else str(k)
        NEG = -np.inf

        # Greedy most-significant -10 subgroup per position (Path C anchor).
        ep10m = np.where(pass10, ep10, np.inf)
        best10_p = ep10m.min(0)
        best10_sub = np.where(ep10m == best10_p[None, :], s10, NEG).argmax(0)
        # Best and second-best strict -35 per position (to exclude the -10's
        # own subgroup for Path C without a per-position rescan).
        ep35sm = np.where(m35s, ep35, np.inf)
        ord35 = np.argsort(ep35sm, axis=0, kind="stable")
        top1_sub = ord35[0]; top1_p = ep35sm[top1_sub, np.arange(P)]
        if self.S35 > 1:
            top2_sub = ord35[1]; top2_p = ep35sm[top2_sub, np.arange(P)]
        else:
            top2_sub = top1_sub; top2_p = np.full(P, np.inf)

        # Extended -10 flag per position: "TG" at i-3, i-2.
        ext10 = np.zeros(P, dtype=bool)
        for i in range(3, P):
            ext10[i] = seq[i - 3:i - 1] == "TG"

        # Per-position bests for each path: combined score, subgroup(s), j, spacer.
        combA = np.full(P, NEG); subA = np.zeros(P, int); jA = np.zeros(P, int); spA = np.zeros(P, int)
        combB = np.full(P, NEG); subB = np.zeros(P, int); jB = np.zeros(P, int); spB = np.zeros(P, int)
        combC = np.full(P, NEG); subC = np.zeros(P, int); jC = np.zeros(P, int); spC = np.zeros(P, int)
        for spc in range(sp0, sp1 + 1):
            off = w + spc
            if off >= P:
                continue
            ii = np.arange(off, P); jj = ii - off
            raw = s10[:, ii] + s35[:, jj]                      # tie-break on raw log-odds
            # -- Path A: linked strict pair, any subgroup --
            mA = pass10[:, ii] & m35s[:, jj]
            cA = np.where(mA, nl10[:, ii] + nl35[:, jj], NEG)
            bestA = cA.max(0)
            kA = np.where(cA == bestA[None, :], raw, NEG).argmax(0)
            upd = bestA > combA[ii]
            combA[ii[upd]] = bestA[upd]; subA[ii[upd]] = kA[upd]; jA[ii[upd]] = jj[upd]; spA[ii[upd]] = spc
            # -- Path B: strict -10 + relaxed -35, any subgroup (ext10 gate applied later) --
            mB = pass10[:, ii] & m35r[:, jj]
            cB = np.where(mB, nl10[:, ii] + nl35[:, jj], NEG)
            bestB = cB.max(0)
            kB = np.where(cB == bestB[None, :], raw, NEG).argmax(0)
            upd = bestB > combB[ii]
            combB[ii[upd]] = bestB[upd]; subB[ii[upd]] = kB[upd]; jB[ii[upd]] = jj[upd]; spB[ii[upd]] = spc
            # -- Path C: greedy -10 subgroup with best strict -35 of another subgroup --
            g = best10_sub[ii]
            use2 = top1_sub[jj] == g
            pC = np.where(use2, top2_p[jj], top1_p[jj])
            sC = np.where(use2, top2_sub[jj], top1_sub[jj])
            valid = np.isfinite(best10_p[ii]) & np.isfinite(pC)
            cC = np.where(valid, -np.log10(best10_p[ii]) - np.log10(np.maximum(pC, 1e-300)), NEG)
            upd = cC > combC[ii]
            combC[ii[upd]] = cC[upd]; subC[ii[upd]] = sC[upd]; jC[ii[upd]] = jj[upd]; spC[ii[upd]] = spc

        best = None
        for i in range(lo, P):
            if combA[i] > NEG:
                path, k, j, spc = "A", int(subA[i]), int(jA[i]), int(spA[i])
                p10 = float(ep10[k, i]); p35_raw = float(ep35[k, j]); sub35 = k
            elif ext10[i] and combB[i] > NEG:
                path, k, j, spc = "B", int(subB[i]), int(jB[i]), int(spB[i])
                p10 = float(ep10[k, i]); p35_raw = float(ep35[k, j]); sub35 = k
            elif combC[i] > NEG:
                path, k, j, spc = "C", int(best10_sub[i]), int(jC[i]), int(spC[i])
                sub35 = int(subC[i]); p10 = float(best10_p[i]); p35_raw = float(ep35[sub35, j])
            else:
                continue
            k35 = (self.S35 - 1) if path == "C" else 1
            p35_adj = sidak(p35_raw, k35)
            comb = -math.log10(p10) - math.log10(max(p35_adj, 1e-300))
            key = ({"A": 0, "B": 1, "C": 2}[path], -comb)
            if best is None or key < best[0]:
                best = (key, {"path": path, "emp_p_10": p10,
                              "emp_p_35": p35_adj,
                              "emp_p_35_uncorrected": p35_raw,
                              "p35_multiplicity_k": int(k35),
                              "pos_10": int(i),
                              "pos_35": int(j), "spacer": int(spc),
                              "source_10": name(k), "source_35": name(sub35),
                              "has_ext10": bool(ext10[i])})
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
