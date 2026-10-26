"""Fit the ProFinder PPV calibration on held-out genomes.

The pipeline reports a length-standardised score (``motif_score_std``) for
every call but, by default, leaves ``lr_pos`` / ``ppv_interval`` blank —
those calibrated quantities require truth data and are never fabricated.
This tool PRODUCES that calibration; the pipeline only APPLIES the resulting
JSON via ``--ppv-calibration``.

Input: a labeled predictions TSV, one row per candidate promoter-orientation
IGR across the calibration genomes, with columns:

    genome     genome id (calibration is fitted leave-one-genome-out)
    score      per-call length-standardised score (the pipeline's
               ``motif_score_std``); blank for intervals with no A/B/C hit
    path       A | B | C  (blank if no hit)
    is_true    1 if the interval is a known true promoter, else 0
    gc         (optional) genome GC fraction, for the applicability domain
    phylum     (optional) genome phylum, for the applicability domain

IMPORTANT: the predictions must come from THIS pipeline's interval
definition (the same ``--max-dist-to-cds-start`` and IGR bounds). A
calibration fitted on a different interval set does not transfer, so a
generic TSS-truth table cannot be reused as-is.

Output: ``ppv_calibration.json`` — a per-path isotonic score->P(true) curve
with each path's fit-set prior, the cross-genome PPV IQR, and the
applicability domain (GC range, phyla). Place it at
``profinder/ppv_calibration.json`` to bundle it, or pass it with
``--ppv-calibration``. The run also prints a leave-one-genome-out
reliability check and a per-tier precision@recall summary; verify these
before shipping.

Usage:
    python -m profinder.calibrate_ppv labeled_predictions.tsv -o ppv_calibration.json
"""
import argparse
import numpy as np
import pandas as pd
from profinder.empirical import Calibration, per_tier_precision, logo_isotonic_ppv

# Reporting tiers: the "strict core" (linked/unlinked strict -35, paths A+C)
# and the full set including the extended--10 path B.
TIERS = {"strict_core": {"A", "C"}, "with_path_B": {"A", "B", "C"}}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fit the ProFinder per-path PPV/LR+ calibration on held-out genomes.")
    ap.add_argument("labeled_tsv", help="labeled predictions TSV (see module docstring)")
    ap.add_argument("-o", "--out", default="ppv_calibration.json",
                    help="output calibration JSON (default: ppv_calibration.json)")
    ap.add_argument("--score-col", default="score")
    ap.add_argument("--path-col", default="path")
    ap.add_argument("--truth-col", default="is_true")
    ap.add_argument("--genome-col", default="genome")
    ap.add_argument("--gc-col", default="gc")
    ap.add_argument("--phylum-col", default="phylum")
    a = ap.parse_args(argv)

    df = pd.read_csv(a.labeled_tsv, sep="\t")
    for col in (a.score_col, a.path_col, a.truth_col, a.genome_col):
        if col not in df.columns:
            ap.error(f"required column {col!r} not found in {a.labeled_tsv} "
                     f"(columns: {list(df.columns)})")

    score = pd.to_numeric(df[a.score_col], errors="coerce").to_numpy()
    path = df[a.path_col].astype(str).str.strip().str.upper().to_numpy()
    y = pd.to_numeric(df[a.truth_col], errors="coerce").fillna(0).astype(int).to_numpy()
    genome = df[a.genome_col].astype(str).to_numpy()
    # Fit only on intervals with a qualifying A/B/C hit and a finite score.
    hit = np.isin(path, list("ABC")) & np.isfinite(score)
    if hit.sum() < 20:
        ap.error(f"only {int(hit.sum())} hit-bearing calls; need at least 20 to fit.")

    gc = None
    if a.gc_col in df.columns:
        g = pd.to_numeric(df[a.gc_col], errors="coerce").to_numpy()
        g = g[np.isfinite(g)]
        gc = g if g.size else None
    phyla = None
    if a.phylum_col in df.columns:
        p = df[a.phylum_col].astype(str)
        p = p[p.str.strip() != ""].tolist()
        phyla = p or None

    n_genomes = int(len(np.unique(genome)))
    print(f"{len(df)} candidate IGRs across {n_genomes} genomes; "
          f"{int(hit.sum())} with an A/B/C hit; {int(y.sum())} true promoters "
          f"(overall baseline {y.mean():.3f})")

    # Per-tier precision @ recall, pooled over genomes.
    print("\nper-tier precision @ recall (pooled over genomes):")
    tp_all = per_tier_precision(path, y, TIERS)
    for name in TIERS:
        tp = tp_all[name]
        print(f"  {name:12s} PPV={tp['ppv']:.3f}  recall={tp['recall']:.3f}  n={tp['n']}")

    # Leave-one-genome-out isotonic PPV on the hit-bearing calls: this is the
    # out-of-sample check that the calibration generalises across genomes.
    held = logo_isotonic_ppv(score[hit], y[hit], genome[hit])
    ok = np.isfinite(held)
    if ok.any():
        yh = y[hit]
        edges = np.unique(np.quantile(held[ok], np.linspace(0, 1, 6)))
        print("\nheld-out (leave-one-genome-out) reliability, predicted vs realised:")
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = ok & (held >= lo) & (held <= hi)
            if m.sum() > 20:
                print(f"  [{lo:.2f},{hi:.2f}] predicted {held[m].mean():.3f} "
                      f"realised {yh[m].mean():.3f}  n={int(m.sum())}")

    # Fit the bundled per-path calibration on the hit-bearing calls.
    cal = Calibration.fit(path[hit], score[hit], y[hit],
                          genomes=genome[hit], gc_values=gc, phyla=phyla,
                          source="calibrate_ppv.py", n_genomes=n_genomes,
                          n_calls=int(hit.sum()))
    cal.to_json(a.out)

    paths = sorted(k for k in cal.per_path if k != "_all")
    print(f"\nwrote {a.out}")
    print(f"  paths calibrated: {paths}")
    for p in paths:
        print(f"    path {p}: prior pi0={cal.per_path[p]['pi0']:.3f}, "
              f"{len(cal.per_path[p]['score_knots'])} isotonic knots")
    if cal.ppv_iqr:
        print(f"  cross-genome PPV IQR: {cal.ppv_iqr[0]:.3f}-{cal.ppv_iqr[1]:.3f}")
    if cal.domain.get("gc_min") is not None:
        dom = f"  applicability domain GC: {cal.domain['gc_min']:.3f}-{cal.domain['gc_max']:.3f}"
        if cal.domain.get("phyla"):
            dom += f"; phyla={cal.domain['phyla']}"
        print(dom)
    print("  NOTE: check the reliability table above before shipping, and "
          "confirm the predictions came from this pipeline's interval definition.")


if __name__ == "__main__":
    main()
