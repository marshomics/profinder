"""ProFinder empirical scoring + PPV calibration (interval-level, per-genome
shuffle null; replaces the analytic per-position p-value). Comment 2 stays
resolved: significance is empirical, not tied to the analytic null."""
from .ae_shuffle import ae_shuffle, dinucleotide_counts
from .fdr_core import (parse_meme, paired_subgroups, freq_to_logodds,
                       encode, score_stack)
from .fdr_stats import (StratumNull, stratify, benjamini_hochberg, empp_stack,
                        sidak)
from .ppv_calibration import (Calibration, fit_isotonic_ppv,
                              per_tier_precision, logo_isotonic_ppv, pava)
from .scan import scan_igrs, build_scan_engine, reclassify
