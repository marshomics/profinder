"""Per-genome orchestration: build the shuffle null per (length x
dinucleotide-composition) stratum, classify each promoter-orientation IGR
(path A/B/C on the empirical element thresholds), attach the length-
standardised score and the base-rate-free
LR+ and (against a stated prior) ppv_interval. Returns a calls DataFrame.
"""
import random
import numpy as np
import pandas as pd
from . import fdr_core as fc
from . import fdr_stats as fs
from .ae_shuffle import ae_shuffle
from .. import ranker as _ranker

_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def _oriented(seq, orientation):
    seq = str(seq).upper()
    return seq.translate(_RC)[::-1] if orientation == "CO_R" else seq


def _seqs_for(co):
    if "sequence_5p_to_3p" in co.columns and co["sequence_5p_to_3p"].notna().all():
        return [str(s).upper() for s in co["sequence_5p_to_3p"]]
    return [_oriented(s, o) for s, o in zip(co["sequence"], co["orientation"])]


def build_scan_engine(igr_df, pwm10, pwm35, subs, *, max_dist, n_shuffles=10,
                      n_len_bins=2, n_comp=3, min_stratum=150, seed=0):
    """Build the per-stratum shuffle null ONCE (the expensive step: shuffle +
    score). Returns (co, seqs, strata, engines), where ``engines`` maps each
    stratum label -> (member_indices, StratumNull). The nulls are
    threshold-independent, so a p-value sweep can re-threshold them cheaply
    via reclassify() without re-shuffling."""
    if n_shuffles < 1:
        raise ValueError(f"n_shuffles must be >= 1 (got {n_shuffles}); the empirical "
                         f"null needs at least one shuffle per interval")
    co = igr_df[igr_df["orientation"].isin(["CO_F", "CO_R"])].copy().reset_index(drop=True)
    if co.empty:
        return co, [], np.array([]), {}
    seqs = _seqs_for(co)
    strata = fs.stratify(seqs, n_len_bins=n_len_bins, n_comp=n_comp,
                         min_size=min_stratum, seed=seed)
    rng = random.Random(seed)
    engines = {}
    for lab in np.unique(strata):
        idx = np.where(strata == lab)[0]
        shuf = [[fc.encode(ae_shuffle(seqs[i], rng=rng))]
                for i in idx for _ in range(n_shuffles)]
        sn = fs.StratumNull(pwm10, pwm35, 1e-3, 1e-3, 3e-2, cap=150_000,
                            seed=(int(lab) & 0x7fffffff),
                            max_dist_to_cds_start=max_dist, subs=subs)
        sn.build(shuf); sn.fit_length_location(shuf)
        engines[int(lab)] = (idx, sn)
    return co, seqs, strata, engines


def _finish_calls(co, seqs, engines, calibration, prior, genome_gc):
    """Classify every interval with each stratum's CURRENT thresholds and
    assemble the calls DataFrame (one row per interval with an A/B/C hit)."""
    dom = calibration.in_domain(genome_gc) if calibration is not None else None
    recs = []
    for lab, (idx, sn) in engines.items():
        for i in idx:
            r = sn.classify_and_score(seqs[i])
            if r is None:
                continue
            row = co.iloc[i]
            p10i, p35i = r["pos_10"], r["pos_35"]
            r.update(igr_id=row.get("igr_id"), contig=row.get("contig"),
                     start=int(row["start"]), end=int(row["end"]),
                     length=int(row["length"]), orientation=row["orientation"],
                     left_gene=row.get("left_gene"), right_gene=row.get("right_gene"),
                     seq_10=seqs[i][p10i:p10i + 6], seq_35=seqs[i][p35i:p35i + 6],
                     sequence_5p_to_3p=seqs[i], stratum=int(lab))
            if calibration is not None:
                r["lr_pos"] = calibration.lr_pos(r["path"], r["score_std"])
                r["ppv_interval"] = calibration.ppv_from_prior(r["path"], r["score_std"], prior)
            else:
                r["lr_pos"] = None; r["ppv_interval"] = None
            r["in_domain"] = dom
            r["shortlist_score"] = _ranker.score(r, seqs[i])
            recs.append(r)
    return pd.DataFrame(recs)


def reclassify(co, seqs, engines, *, p10, p35, p35_relaxed,
               calibration=None, prior=None, genome_gc=None):
    """Re-threshold a built engine at new (p10, p35, p35_relaxed) levels and
    re-classify all intervals. Cheap relative to build_scan_engine (no
    re-shuffling); the p-value sweep calls this once per sweep point."""
    for _idx, sn in engines.values():
        sn.set_thresholds(p10, p35, p35_relaxed)
    if co is None or len(co) == 0:
        return pd.DataFrame()
    return _finish_calls(co, seqs, engines, calibration, prior, genome_gc)


def scan_igrs(igr_df, pwm10, pwm35, subs, *, p10, p35, p35_relaxed, max_dist,
              n_shuffles=10, n_len_bins=2, n_comp=3, min_stratum=150, seed=0,
              calibration=None, prior=None, genome_gc=None):
    """Score all CO_F/CO_R IGRs. Returns one row per CALL (interval with a
    path A/B/C hit) with columns: igr_id, contig, start, end, length,
    orientation, path, source_10/35, has_ext10, emp_p_10, emp_p_35, pos_10/35,
    spacer, score_std, n_pairs, lr_pos, ppv_interval, in_domain, sequence_5p_to_3p,
    stratum. Uncalibrated -> lr_pos / ppv_interval are None (never fabricated)."""
    co, seqs, strata, engines = build_scan_engine(
        igr_df, pwm10, pwm35, subs, max_dist=max_dist, n_shuffles=n_shuffles,
        n_len_bins=n_len_bins, n_comp=n_comp, min_stratum=min_stratum, seed=seed)
    if co.empty:
        return pd.DataFrame()
    return reclassify(co, seqs, engines, p10=p10, p35=p35, p35_relaxed=p35_relaxed,
                      calibration=calibration, prior=prior, genome_gc=genome_gc)
