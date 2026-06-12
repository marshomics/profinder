#!/usr/bin/env python3
"""
ProFinder — bacterial and archaeal promoter identification pipeline.

Takes a single genome FASTA, annotates it with Prokka, extracts
intergenic regions, identifies operons, screens for HMM marker genes,
and outputs promoter sequences in 5'-to-3' orientation.

Promoter verification is performed by scanning for -10/-35 hexamer
motifs using position weight matrices from the bundled
``all_unique_subgroups.meme`` file. PWM log-odds, score thresholds,
and empirical p-values are all calibrated against per-genome A/C/G/T
frequencies derived from the input FASTA, so significance reflects
the composition of the genome being scanned rather than a fixed
uniform background. Before scoring, a low-complexity mask is computed
over each strand and 6-mer anchor positions inside masked tracts
(simple repeats, homopolymers, strongly biased composition) are
skipped, so the AT-rich -10 PWM does not pile up spurious hits in
poly-A runs. The motif scan is also orientation-aware: tandem-
cooriented IGRs (CO_F, CO_R) are scanned on the ``+`` strand of the
already-oriented ``sequence_5p_to_3p`` only, because a ``-`` strand
motif on such an IGR would describe a promoter pointing away from
the downstream marker gene; divergent IGRs (DP) are scanned on both
strands since each strand carries the promoter for one of the two
flanking genes. When multiple subgroups exceed threshold at the same
position, the most-significant one (lowest empirical p-value under
the per-genome background) is reported, with raw log-odds as the
tiebreaker. Paths A–C require a -35 hit with a 16–18 bp spacer to the
-10 hit (canonical σ⁷⁰ spacer is 17 ± 1); Path D requires only an
extended -10. -10 motifs whose right edge sits further than
``cfg.max_dist_to_cds_start`` bp from the downstream CDS start
(default 200) are dropped before the spacer search — distal motifs
in long IGRs are unlikely to drive the downstream gene and contribute
disproportionately to false-positive predictions. Each -10 hit is
classified as Path A (linked strict -35, same subgroup), Path B
(extended -10 with a linked relaxed -35, same subgroup), Path C
(unlinked strict -35, different subgroups), or Path D (extended -10
with no -35 in the spacer window). The extended -10 anchor accepts
"TG" at either of the two literature-reported positions immediately
upstream of the -10 hexamer (-2/-1 or -3/-2 of the motif, equivalent
to TSS positions -14/-13 and -15/-14 for a -10 anchored at the
textbook -12 to -7). Anything else is no hit. Within each path group,
rows are ranked by combined significance -log10(p10) + -log10(p35).

Usage
-----
    # Minimal (Prokka + IGR extraction + operon ID + promoter output):
    profinder -i genome.fasta -o results/

    # With HMM marker screening:
    profinder -i genome.fasta -o results/ \\
        --tigrfam /path/to/tigrfam.hmm --pfam /path/to/Pfam-A.hmm

    # Resume after interruption (skips completed steps):
    profinder -i genome.fasta -o results/

    # Force re-run from scratch:
    profinder -i genome.fasta -o results/ --force

    # Run specific steps:
    profinder -i genome.fasta -o results/ --start 3 --end 5

    # List steps:
    profinder --list

Steps
-----
    1.  Run Prokka (genome annotation)
    2.  Extract intergenic regions from GFF + FASTA
    3.  Identify operons from GFF
    4.  Run hmmsearch (TIGRfam + Pfam)       [bundled HMMs used by default]
    5.  Filter HMM output
    6.  Filter operons and add marker info
    7.  Match IGRs to marker operons
    8.  Scan promoter motifs (-10/-35 hexamer verification, A/B/C/D)
    9.  Annotate associated CDS (Prokka product names, gene, locus tag)
    10. Build final output table (all promoters, all metadata)
    11. Generate HTML report

Diagnostic mode
---------------
    --diagnose-p-sweep replaces steps 8, 10, 11 with a p-value sweep.
    The pipeline runs steps 1-7 and step 9 (Prokka, IGR extraction,
    operons, HMM, marker matching, CDS annotation — all p-value
    independent), then scans every CO_F/CO_R IGR at a range of
    (--p10, --p35) settings in three modes — matched, p10_only,
    p35_only — and writes:

      diagnostics/p_sweep.tsv         per-sweep summary counts
      diagnostics/p_sweep.png         3-panel figure (if matplotlib)
      diagnostics/sweep_tables/       one profinder_results.tsv per
                                      (mode, p10, p35), pre-filtered
                                      to motif_confirmed = True and
                                      motif_path in {A, B}

    Each per-sweep-point table has the same column layout as the
    canonical profinder_results.tsv produced by step 10, so it can be
    consumed directly by downstream benchmarking (e.g. Section 5).
    Useful for choosing defensible p-value thresholds for a given
    genome before committing to a full run.
"""

import argparse
import bisect
import html as html_mod
import math
import re
import subprocess
import sys
from itertools import product as iprod
from pathlib import Path
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .config import Config
from .igr_extractor import extract_igrs,_open_maybe_gzip


# Shared reverse-complement translation table.  Used by the motif-scan
# inner loop which builds two strands per sequence — recreating the
# translation table per row is measurable on large IGR sets.
_REVCOMP_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp_simple(s: str) -> str:
    """Fast reverse-complement using the module-level translation table."""
    return s.translate(_REVCOMP_TRANS)[::-1]


# =====================================================================
#  Motif scanning helpers (adapted from meme_scan.py)
# =====================================================================

_MOTIF_WIDTH = 6
# Spacer window between -35 and -10 motifs. Canonical σ⁷⁰ spacer is
# 17 ± 1 bp (Hawley & McClure 1983; Lisser & Margalit 1993); we
# accept 16-18 bp so the modal 17 bp case sits in the middle of the
# window with one nt of tolerance either side. Earlier versions of
# ProFinder used the looser 15-19 bp window, which captured more
# real promoters but also accepted edge cases where the geometry no
# longer supports σ⁷⁰ holoenzyme contacts at both motifs.
_SPACER_MIN = 16
_SPACER_MAX = 18
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _parse_meme_file(filepath):
    """Parse a MEME-format file. Returns {motif_name: [[pA, pC, pG, pT], ...]}.

    The motif name is the first token following ``MOTIF``. For
    ``all_unique_subgroups.meme`` these look like ``M001_m35``,
    ``M001_m10``, ``M002_m35`` etc — the prefix before the underscore
    identifies the subgroup, and the ``_m10`` / ``_m35`` suffix the
    element type.
    """
    motifs = {}
    current = None
    matrix = []
    in_matrix = False

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("MOTIF "):
                if current and matrix:
                    motifs[current] = matrix
                current = line.split()[1]
                matrix = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix:
                if not line:
                    in_matrix = False
                    continue
                vals = line.split()
                if len(vals) == 4:
                    try:
                        matrix.append([float(v) for v in vals])
                    except ValueError:
                        in_matrix = False
                else:
                    in_matrix = False

    if current and matrix:
        motifs[current] = matrix
    return motifs


_UNIFORM_BG = (0.25, 0.25, 0.25, 0.25)


def _compute_background_from_fasta(fasta_path):
    """Compute genome-wide A/C/G/T frequencies from a FASTA.

    Returns a 4-tuple ``(p_A, p_C, p_G, p_T)`` summed across all contigs,
    case-insensitive, ignoring ambiguous bases (N, IUPAC codes) and any
    non-letter characters. Falls back to uniform 0.25 if the FASTA is
    empty or unreadable. Used to calibrate PWM log-odds and the
    enumeration of k-mer probabilities for threshold / p-value lookup,
    per the comment in ``all_unique_subgroups.meme`` that "Per-organism
    A/C/G/T frequencies used in the scan are in per_organism_background.tsv".
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    try:
        with _open_maybe_gzip(fasta_path) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                s = str(rec.seq).upper()
                for b in "ACGT":
                    counts[b] += s.count(b)
    except Exception:
        return _UNIFORM_BG
    total = sum(counts.values())
    if total == 0:
        return _UNIFORM_BG
    return (counts["A"] / total, counts["C"] / total,
            counts["G"] / total, counts["T"] / total)


def _compute_low_complexity_mask(seq, window=12, min_entropy=1.4):
    """Return a list[bool] of length ``len(seq)`` marking motif-anchor
    positions whose surrounding context is low-complexity.

    For each position ``i``, the Shannon entropy of single-nucleotide
    frequencies is computed over a window of ``window`` bp centered on
    ``i`` (clipped at sequence boundaries). Positions with entropy
    below ``min_entropy`` bits are masked. Maximum entropy for a
    four-letter alphabet is 2.0 bits; homopolymers give 0.0, simple
    dinucleotide repeats give ~1.0, AT-only or GC-only tracts give
    ≤ 1.0, and well-mixed sequence typically exceeds 1.6 in a 12 bp
    window. The default cutoff of 1.4 catches homopolymers, simple
    repeats, and strongly biased tracts while leaving normal IGR
    composition untouched.

    Low-complexity masking before PWM scanning is standard practice
    (Wootton & Federhen, 1996; Morgulis et al., 2006); ProFinder
    previously applied no masking, leaving the -10 PWM (AT-rich
    consensus) prone to spurious calls inside poly-A and AT-repeat
    tracts common in bacterial IGRs.
    """
    n = len(seq)
    if n == 0:
        return []
    masked = [False] * n
    half = window // 2
    extra = window - half  # window may be odd
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + extra)
        win = seq[lo:hi]
        if len(win) < 4:
            masked[i] = True
            continue
        counts = {}
        for b in win:
            counts[b] = counts.get(b, 0) + 1
        total = len(win)
        ent = 0.0
        for c in counts.values():
            p = c / total
            ent -= p * math.log2(p)
        if ent < min_entropy:
            masked[i] = True
    return masked


def _freq_to_log_odds(freq_matrix, bg=_UNIFORM_BG, pseudocount=1e-6):
    """Convert a frequency matrix to log2-odds against ``bg``.

    ``bg`` is a 4-tuple (p_A, p_C, p_G, p_T). The pseudocount only
    matters when the input frequency is exact zero; the bundled
    ``all_unique_subgroups.meme`` has no zeros so it does not fire in
    practice. The same ``bg`` is reused by ``_enumerate_pwm_table`` to
    build the k-mer probability distribution against which thresholds
    and p-values are calibrated, so log-odds and significance share
    one background.
    """
    return [[math.log2(max(f, pseudocount) / max(bg[b], 1e-12))
             for b, f in enumerate(row)]
            for row in freq_matrix]


def _enumerate_pwm_table(log_odds_matrix, bg=_UNIFORM_BG):
    """Enumerate all 4^w k-mer scores under this PWM, return arrays
    suitable for background-aware threshold and p-value lookup.

    Returns ``(sorted_scores, cum_above)`` where ``sorted_scores`` is
    a list of all 4^w log-odds scores in ascending order and
    ``cum_above[i]`` is the cumulative probability — under ``bg`` — of
    a random k-mer scoring at least ``sorted_scores[i]``. By
    construction ``cum_above[0] == 1.0`` and ``cum_above[n] == 0.0``.
    With uniform background each k-mer contributes 1/4^w and this
    reduces to the prior count-based formulation; with non-uniform
    background, common k-mers (e.g. AT-rich hexamers in AT-rich
    genomes) carry more weight and the threshold for a given p-value
    is appropriately tightened.
    """
    w = len(log_odds_matrix)
    pairs = []
    for kmer in iprod(range(4), repeat=w):
        score = sum(log_odds_matrix[i][b] for i, b in enumerate(kmer))
        prob = 1.0
        for b in kmer:
            prob *= bg[b]
        pairs.append((score, prob))
    pairs.sort(key=lambda x: x[0])
    n = len(pairs)
    sorted_scores = [p[0] for p in pairs]
    cum_above = [0.0] * (n + 1)
    for i in range(n - 1, -1, -1):
        cum_above[i] = cum_above[i + 1] + pairs[i][1]
    return sorted_scores, cum_above


def _threshold_from_table(sorted_scores, cum_above, p_value):
    """Smallest score ``s`` such that P(random k-mer scores ≥ s | bg)
    is ≤ ``p_value``. Realised false-positive rate is therefore at
    most ``p_value``. If even the highest-scoring k-mer alone exceeds
    ``p_value`` (i.e. ``cum_above[n-1] > p_value``), returns a score
    one log-odds unit above the maximum (effectively unreachable).
    """
    n = len(sorted_scores)
    if n == 0:
        return float("inf")
    # cum_above[i] is monotonically decreasing in i. Find smallest i with
    # cum_above[i] <= p_value via binary search.
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if cum_above[mid] > p_value:
            lo = mid + 1
        else:
            hi = mid
    if lo >= n:
        return sorted_scores[-1] + 1.0
    return sorted_scores[lo]


def _pvalue_from_table(sorted_scores, cum_above, score):
    """Empirical p-value for ``score`` under the background:
    P(random k-mer scores ≥ score). For scores above every enumerated
    k-mer's score, returns the smallest single-k-mer probability in
    the table (the natural floor — that floor equals ``cum_above[n-1]``,
    the probability of the single highest-scoring k-mer)."""
    n = len(sorted_scores)
    if n == 0:
        return 1.0
    idx = bisect.bisect_left(sorted_scores, score)
    if idx >= n:
        return cum_above[n - 1]
    return cum_above[idx]


def _compute_score_threshold(log_odds_matrix, p_value, bg=_UNIFORM_BG):
    """Backwards-compatible wrapper: enumerate against ``bg``, then
    threshold. ``bg`` defaults to uniform 0.25 so existing callers
    still get the previous behaviour. Prefer enumerating once and
    reusing the table when running many lookups."""
    sorted_scores, cum_above = _enumerate_pwm_table(log_odds_matrix, bg)
    return _threshold_from_table(sorted_scores, cum_above, p_value)


def _score_kmer(seq, pos, log_odds_matrix):
    """Score a 6-mer at the given position. Returns None if non-ACGT bases present."""
    s = 0.0
    for j in range(_MOTIF_WIDTH):
        idx = _BASE_IDX.get(seq[pos + j])
        if idx is None:
            return None
        s += log_odds_matrix[j][idx]
    return s


class _MotifSet:
    """A collection of subgroup motifs for one element type (-10 or -35),
    each with its own log-odds matrix, background-aware score threshold,
    and cumulative-probability table for p-value lookup. Subgroup IDs
    (e.g. ``M001``) are used to test for "linked" hits where a -10 and
    a -35 come from the same subgroup."""

    def __init__(self, bg=_UNIFORM_BG):
        # list of (subgroup, log_odds_matrix, threshold, sorted_scores, cum_above)
        self.entries = []
        self.bg = bg

    def add(self, subgroup, freq_matrix, p_value):
        lom = _freq_to_log_odds(freq_matrix, self.bg)
        sorted_scores, cum_above = _enumerate_pwm_table(lom, self.bg)
        thresh = _threshold_from_table(sorted_scores, cum_above, p_value)
        self.entries.append((subgroup, lom, thresh, sorted_scores, cum_above))

    def best_hit(self, seq, pos):
        """Return (score, subgroup, pvalue) of the most-significant
        motif (lowest empirical p-value) that exceeds its threshold at
        this position, or None. Tie-break by higher raw score.

        Empirical p-values are computed under the same background
        ``self.bg`` that calibrated each PWM's threshold, so they are
        comparable across subgroups of differing information content
        and are the principled selection criterion. (Raw log-odds are
        not cross-subgroup comparable: a sharp PWM yields higher
        maxima than a flat PWM for the same biological motif.)
        """
        best_pv = None
        best = None
        for subgroup, lom, thresh, sorted_scores, cum_above in self.entries:
            s = _score_kmer(seq, pos, lom)
            if s is None or s < thresh:
                continue
            pv = _pvalue_from_table(sorted_scores, cum_above, s)
            if best_pv is None or pv < best_pv or (pv == best_pv and s > best[0]):
                best_pv = pv
                best = (s, subgroup, pv)
        return best

    def hit_for_subgroup(self, seq, pos, subgroup):
        """Return (score, pvalue) for a specific subgroup's motif at
        this position if it exceeds that subgroup's threshold, else
        None. Subgroups are unique per MotifSet so the loop terminates
        once a match is found."""
        for sub, lom, thresh, sorted_scores, cum_above in self.entries:
            if sub != subgroup:
                continue
            s = _score_kmer(seq, pos, lom)
            if s is not None and s >= thresh:
                pv = _pvalue_from_table(sorted_scores, cum_above, s)
                return (s, pv)
            return None
        return None


def _load_motifs_from_file(meme_path, p10, p35_strict, p35_relaxed,
                            bg=_UNIFORM_BG):
    """Load paired subgroup motifs from a single MEME file, calibrated
    against the supplied background ``bg`` (4-tuple of A/C/G/T
    frequencies).

    Expects motifs named ``<subgroup>_m10`` and ``<subgroup>_m35`` —
    e.g. ``M001_m10`` / ``M001_m35``. Returns
    ``(m10, m35_strict, m35_relaxed)`` MotifSets where each entry's
    source is the subgroup ID. The strict and relaxed -35 sets hold
    the same PWMs but with different score thresholds: strict is used
    for Path A and Path C, relaxed for Path B. All three sets are
    calibrated against the same ``bg``.
    """
    m10 = _MotifSet(bg=bg)
    m35_strict = _MotifSet(bg=bg)
    m35_relaxed = _MotifSet(bg=bg)

    meme_path = Path(meme_path)
    if not meme_path.is_file():
        print(f"  WARNING: MEME file not found: {meme_path}", file=sys.stderr)
        return m10, m35_strict, m35_relaxed

    motifs = _parse_meme_file(meme_path)
    if not motifs:
        print(f"  WARNING: no motifs parsed from {meme_path}", file=sys.stderr)
        return m10, m35_strict, m35_relaxed

    for name, matrix in motifs.items():
        # Expected naming: "<subgroup>_m10" or "<subgroup>_m35"
        if name.endswith("_m10"):
            subgroup = name[:-4]
            m10.add(subgroup, matrix, p10)
        elif name.endswith("_m35"):
            subgroup = name[:-4]
            m35_strict.add(subgroup, matrix, p35_strict)
            m35_relaxed.add(subgroup, matrix, p35_relaxed)

    print(f"  Loaded {len(m10.entries)} -10 and {len(m35_strict.entries)} -35 "
          f"subgroup motifs from {meme_path.name}")
    return m10, m35_strict, m35_relaxed


def _find_promoters_on_strand(seq, m10, m35_strict, m35_relaxed,
                                max_dist_to_cds_start=None):
    """Scan one strand for promoter candidates. Returns a list of result dicts.

    Before scoring, a low-complexity mask is computed over the strand
    (see ``_compute_low_complexity_mask``). Both -10 and -35 PWM
    lookups skip masked positions, so that AT-rich repeat tracts
    (which the -10 PWM is otherwise prone to call) do not contribute
    spurious hits.

    If ``max_dist_to_cds_start`` is provided, -10 hits whose right
    edge sits further than that many bases from the end of ``seq`` are
    dropped before the spacer search. This assumes ``seq`` is oriented
    5'→3' relative to the downstream gene so that the right end of
    ``seq`` is adjacent to the downstream CDS start — true for the
    ``+`` strand of CO_F/CO_R IGRs after step 2's orientation
    normalization, and for both strands of DP IGRs (each strand
    promotes one of the two flanking genes, and in each case the
    "downstream gene" sits immediately past the right end of the
    scanned strand). Distance from the right edge of the -10 to the
    CDS start = ``len(seq) - pos_10 - motif_width``.

    Each -10 hit is classified into one of four paths, in priority
    order. Paths A–C require a -35 hit with a 16-18 bp spacer
    (canonical σ⁷⁰ spacer is 17 ± 1); Path D requires only an
    extended -10:

        A  Linked -10 and strict -35 (same subgroup).
        B  Extended -10 ("TG" immediately upstream — at -2/-1 OR
           -3/-2 of the -10 hexamer) plus a linked relaxed -35 (same
           subgroup). Both anchor positions are accepted because the
           literature places extended -10 at either -15/-14 (canonical;
           "TG" at -3/-2 of a -10 anchored at -12 to -7) or
           -14/-13 (slipped; "TG" at -2/-1).
        C  Unlinked -10 and strict -35 (different subgroups).
        D  Extended -10 with no -35 in the 16-18 bp window.

    Anything else is no hit and is dropped.
    """
    w = _MOTIF_WIDTH
    max_pos = len(seq) - w

    # Per-position low-complexity mask. Both PWM scans consult this
    # mask before scoring, so masked anchor positions contribute
    # nothing to either the -10 hit list or the strict-by-pos -35
    # index. Skipping at the anchor is equivalent to masking the
    # entire 6-mer because the 6-mer extends from the anchor for
    # six bases — anchors inside a low-complexity tract describe
    # 6-mers that lie within that tract.
    mask = _compute_low_complexity_mask(seq)

    # Distance-to-CDS filter: when active, any -10 anchor with right
    # edge sitting more than ``max_dist_to_cds_start`` bases away from
    # the end of ``seq`` is rejected. min_pos_10 is the smallest
    # ``pos_10`` that satisfies the cutoff, so we can skip anchors
    # below it in one ``i < min_pos_10`` comparison. None disables.
    min_pos_10 = None
    if isinstance(max_dist_to_cds_start, (int, float)) and max_dist_to_cds_start >= 0:
        min_pos_10 = max(0, len(seq) - w - int(max_dist_to_cds_start))

    # Pre-scan all -10 hits
    hits_10 = []
    for i in range(max_pos + 1):
        if min_pos_10 is not None and i < min_pos_10:
            continue
        if i < len(mask) and mask[i]:
            continue
        hit = m10.best_hit(seq, i)
        if hit:
            # hit = (score, subgroup, pvalue)
            hits_10.append((i, hit[0], hit[1], hit[2]))

    if not hits_10:
        return []

    # Pre-scan strict -35 hits (best across subgroups) by position.
    # Relaxed -35 is looked up per-subgroup directly (see below), since
    # Path B requires the relaxed -35 to come from the same subgroup as
    # the -10 hit.
    strict_by_pos = {}
    for i in range(max_pos + 1):
        if i < len(mask) and mask[i]:
            continue
        hit_s = m35_strict.best_hit(seq, i)
        if hit_s:
            strict_by_pos[i] = hit_s  # (score, subgroup, pvalue)

    results = []
    for pos_10, score_10, subgroup_10, pvalue_10 in hits_10:
        # Accept "TG" at either of the two literature-reported anchor
        # positions for the extended -10. Both forms are documented in
        # the σ⁷⁰ literature; requiring exactly one anchor was
        # over-restrictive and excluded the canonical -15/-14 form
        # when the -10 hexamer is anchored at the textbook -12 to -7.
        ext10_at_m2 = pos_10 >= 2 and seq[pos_10 - 2:pos_10] == "TG"
        ext10_at_m3 = pos_10 >= 3 and seq[pos_10 - 3:pos_10 - 1] == "TG"
        has_ext10 = ext10_at_m2 or ext10_at_m3

        # Walk the 16-18 bp spacer window once, tracking:
        #   best_linked_strict   — same subgroup, strict threshold (Path A)
        #   best_linked_relaxed  — same subgroup, relaxed threshold (Path B)
        #   best_unlinked_strict — different subgroup, strict threshold (Path C)
        # Each candidate tuple is (score, subgroup, pvalue, pos, spacer).
        best_linked_strict = None
        best_linked_relaxed = None
        best_unlinked_strict = None
        for spacer in range(_SPACER_MIN, _SPACER_MAX + 1):
            p35 = pos_10 - w - spacer
            if p35 < 0:
                continue

            s_hit = strict_by_pos.get(p35)
            if s_hit is not None:
                s35, sub35, pv35 = s_hit
                if sub35 == subgroup_10:
                    if best_linked_strict is None or s35 > best_linked_strict[0]:
                        best_linked_strict = (s35, sub35, pv35, p35, spacer)
                else:
                    if best_unlinked_strict is None or s35 > best_unlinked_strict[0]:
                        best_unlinked_strict = (s35, sub35, pv35, p35, spacer)

            # Same-subgroup relaxed -35 lookup for Path B.
            r_hit = m35_relaxed.hit_for_subgroup(seq, p35, subgroup_10)
            if r_hit is not None:
                r_score, r_pv = r_hit
                if best_linked_relaxed is None or r_score > best_linked_relaxed[0]:
                    best_linked_relaxed = (r_score, subgroup_10, r_pv, p35, spacer)

        # Classify this -10 hit: A > B > C > D, else drop.
        if best_linked_strict is not None:
            path = "A"
            chosen_35 = best_linked_strict
        elif has_ext10 and best_linked_relaxed is not None:
            path = "B"
            chosen_35 = best_linked_relaxed
        elif best_unlinked_strict is not None:
            path = "C"
            chosen_35 = best_unlinked_strict
        elif has_ext10:
            path = "D"
            chosen_35 = None  # no -35 in the spacer window
        else:
            continue

        r = {
            "pos_10": pos_10,
            "seq_10": seq[pos_10:pos_10 + w],
            "score_10": round(score_10, 3),
            "pvalue_10": pvalue_10,
            "source_10": subgroup_10,
            "has_ext10": has_ext10,
            "path": path,
        }

        if chosen_35 is not None:
            s35, sub35, pv35, p35, spacer = chosen_35
            r["pos_35"] = p35
            r["seq_35"] = seq[p35:p35 + w]
            r["score_35"] = round(s35, 3)
            r["pvalue_35"] = pv35
            r["source_35"] = sub35
            r["spacer_len"] = spacer
            r["spacer_seq"] = seq[p35 + w:pos_10]
        else:
            r["pos_35"] = "."
            r["seq_35"] = "."
            r["score_35"] = "."
            r["pvalue_35"] = "."
            r["source_35"] = "."
            r["spacer_len"] = "."
            r["spacer_seq"] = "."

        results.append(r)

    return results


_MOTIF_COLUMNS = [
    "strand", "pos_10", "seq_10", "score_10", "pvalue_10", "source_10",
    "has_ext10", "pos_35", "seq_35", "score_35", "pvalue_35", "source_35",
    "spacer_len", "spacer_seq", "path",
]


def _neg_log10_p_combined(r):
    """Within-path sort key: -log10(p10) + -log10(p35) when both exist,
    else -log10(p10) alone (Path D). Higher = more significant."""
    pv10 = r.get("pvalue_10")
    pv35 = r.get("pvalue_35")
    if not isinstance(pv10, (int, float)) or pv10 <= 0:
        return float("-inf")
    val = -math.log10(pv10)
    if isinstance(pv35, (int, float)) and pv35 > 0:
        val += -math.log10(pv35)
    return val

_PATH_RANK = {"A": 0, "B": 1, "C": 2, "D": 3}


def _strands_for_orientation(orientation):
    """Which strand(s) of ``sequence_5p_to_3p`` can carry a real
    promoter for the IGR's downstream gene, given the PIGGY orientation
    label.

    The IGR has already been oriented 5'→3' relative to its associated
    downstream gene in step 2 (CO_R reverse-complemented; CO_F, DP,
    CONV kept on the + strand of the assembly). For tandem-cooriented
    IGRs (CO_F, CO_R), the downstream gene's promoter can only sit on
    the ``+`` strand of the oriented sequence — a ``-`` strand motif
    would describe a promoter pointing away from the marker, into a
    flanking gene that does not exist on that strand. For divergent
    IGRs (DP), each strand carries the promoter for one of the two
    diverging genes, so both strands are biologically meaningful.
    For convergent IGRs (CONV) there is no associated downstream
    promoter at all; we scan ``+`` only to mirror the CO_F default
    rather than double the work.

    Returns a list of strand labels to scan. Unknown / missing
    orientations fall back to both strands (preserves the original
    behavior for callers that don't supply orientation).
    """
    if orientation in ("CO_F", "CO_R", "CONV"):
        return ["+"]
    if orientation == "DP":
        return ["+", "-"]
    return ["+", "-"]


def _scan_sequences_for_motifs(df, m10, m35_strict, m35_relaxed,
                               seq_col="sequence_5p_to_3p",
                               id_col="igr_id",
                               orientation_col="orientation",
                               max_dist_to_cds_start=None):
    """Scan a DataFrame of sequences for -10/-35 promoter motifs.

    Strand selection is orientation-aware: see
    ``_strands_for_orientation``. For the manuscript's tandem-cooriented
    IGRs (CO_F / CO_R), only the ``+`` strand of the already-oriented
    ``sequence_5p_to_3p`` is scanned, because a ``-`` strand motif
    describes a promoter that cannot transcribe the IGR's downstream
    marker gene. If ``orientation_col`` is absent from ``df``, both
    strands are scanned (legacy behavior).

    ``max_dist_to_cds_start`` (optional) is passed through to
    ``_find_promoters_on_strand`` and constrains the -10 anchor to
    sit within that many bases of the downstream CDS start (= the
    right end of the scanned strand). ``None`` disables the filter.

    Returns two DataFrames:
        all_hits  — every motif hit found (multiple rows per sequence possible)
        best_hits — one row per sequence, keeping the best hit
                    (Path A > B > C > D, then highest combined
                    -log10(p-10) + -log10(p-35); -log10(p-10) alone
                    for Path D)
    """
    all_rows = []
    # igr_id -> (path_rank, neg_combined_neglogp, row_dict). Higher
    # -log10(p) is more significant, so we negate for min-comparison.
    best_per_seq = {}

    has_orientation = orientation_col in df.columns

    for _, row in df.iterrows():
        raw = row[seq_col]
        if not isinstance(raw, str) or not raw:
            continue
        raw = raw.upper()
        igr_id = row[id_col]

        orientation = row[orientation_col] if has_orientation else None
        strands_to_scan = _strands_for_orientation(orientation)
        strand_pairs = [(s, raw if s == "+" else _revcomp_simple(raw))
                        for s in strands_to_scan]

        for strand, seq in strand_pairs:
            for r in _find_promoters_on_strand(
                    seq, m10, m35_strict, m35_relaxed,
                    max_dist_to_cds_start=max_dist_to_cds_start):
                out = {id_col: igr_id, "strand": strand}
                for k in _MOTIF_COLUMNS:
                    if k != "strand":
                        out[k] = r.get(k, ".")
                all_rows.append(out)

                path_rank = _PATH_RANK.get(r["path"], 99)
                neglogp = _neg_log10_p_combined(r)
                prev = best_per_seq.get(igr_id)
                if prev is None or (path_rank, -neglogp) < (prev[0], prev[1]):
                    best_per_seq[igr_id] = (path_rank, -neglogp, out)

    all_hits = pd.DataFrame(all_rows) if all_rows else pd.DataFrame(columns=[id_col] + _MOTIF_COLUMNS)
    best_rows = [v[2] for v in best_per_seq.values()]
    best_hits = pd.DataFrame(best_rows) if best_rows else pd.DataFrame(columns=[id_col] + _MOTIF_COLUMNS)

    return all_hits, best_hits


# =====================================================================
#  Operon identification (carried over from multi-genome pipeline)
# =====================================================================

def _parse_gff_for_operons(gff_path: str):
    """Parse a Prokka GFF into a DataFrame of CDS features."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            rows.append({
                "seqid": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "attributes": cols[8],
            })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df.sort_values(["seqid", "start"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def _extract_gene_id(attributes: str) -> str:
    """Extract a gene identifier from GFF attribute string.

    Tries, in order: ID=, locus_tag=, gene=, Name= (stripping any
    trailing ' gene' or ' CDS' suffix from Geneious-style names).
    """
    for key in ("ID=", "locus_tag=", "gene=", "Name="):
        for attr in attributes.split(";"):
            if attr.startswith(key):
                val = attr[len(key):]
                # Geneious appends " gene" or " CDS" to Name values
                for suffix in (" gene", " CDS"):
                    if val.endswith(suffix):
                        val = val[:-len(suffix)]
                return val
    return ""


def _identify_operons(genes_df, max_internal_distance, min_flanking_distance):
    """Two-pass operon identification: proximity clustering then
    flanking-distance validation."""
    if genes_df.empty:
        return []

    # Pass 1: cluster consecutive genes within max_internal_distance
    raw_groups = []
    current_group = [0]
    for i in range(1, len(genes_df)):
        prev = genes_df.iloc[i - 1]
        curr = genes_df.iloc[i]
        same_contig = curr["seqid"] == prev["seqid"]
        close = (curr["start"] - prev["end"]) <= max_internal_distance
        if same_contig and close:
            current_group.append(i)
        else:
            raw_groups.append(current_group)
            current_group = [i]
    raw_groups.append(current_group)

    # Pass 2: validate flanking distances
    operons = []
    for group in raw_groups:
        first_i, last_i = group[0], group[-1]
        first_gene = genes_df.iloc[first_i]
        last_gene = genes_df.iloc[last_i]

        upstream_ok = True
        if first_i > 0:
            prev_gene = genes_df.iloc[first_i - 1]
            if prev_gene["seqid"] == first_gene["seqid"]:
                upstream_ok = (first_gene["start"] - prev_gene["end"]) >= min_flanking_distance

        downstream_ok = True
        if last_i < len(genes_df) - 1:
            next_gene = genes_df.iloc[last_i + 1]
            if next_gene["seqid"] == last_gene["seqid"]:
                downstream_ok = (next_gene["start"] - last_gene["end"]) >= min_flanking_distance

        if upstream_ok and downstream_ok:
            operons.append([genes_df.iloc[j] for j in group])

    return operons


# =====================================================================
#  HMM output parsing
# =====================================================================

_HMM_HEADER = [
    "target_name", "accession1", "query_name", "accession2",
    "full_sequence_evalue", "full_sequence_bitscore", "full_sequence_bias",
    "best_1_domain_evalue", "best_1_domain_score", "best_1_domain_bias",
    "exp", "reg", "glu", "ov", "env", "dom", "rep", "inc",
    "description_of_target",
]


def _parse_hmm_tblout(file_path: str) -> list:
    rows = []
    with open(file_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split(maxsplit=len(_HMM_HEADER) - 1)
            if len(parts) == len(_HMM_HEADER):
                rows.append(parts)
    return rows


# =====================================================================
#  Pipeline steps
# =====================================================================

def step01_run_prokka(cfg: Config, force: bool = False):
    """Annotate the input genome with Prokka."""
    if not force and cfg.gff_file.exists() and cfg.faa_file.exists():
        print("── Prokka output already exists, skipping ──")
        print("  Step 1 complete.\n")
        return

    print("── Running Prokka ──")
    cfg.prokka_dir.mkdir(parents=True, exist_ok=True)

    parts = []
    if cfg.conda_env_prokka:
        parts.append(
            f'eval "$(conda shell.bash hook 2>/dev/null)"; '
            f"conda activate {cfg.conda_env_prokka}; "
        )
    
    input_for_prokka = cfg.input_fasta

    if str(cfg.input_fasta).endswith(".gz"):
        input_for_prokka = cfg.output_dir / f"{cfg.input_fasta.stem}.decompressed.fna"
        with gzip.open(cfg.input_fasta, "rt") as fin, open(input_for_prokka, "w") as fout:
            fout.write(fin.read())

    parts.append(
        f"{cfg.prokka_bin}"
        f" --outdir {cfg.prokka_dir}"
        f" --prefix {cfg.prokka_prefix}"
        f" --kingdom {cfg.prokka_kingdom}"
        f" --cpus {cfg.threads}"
        f" --rfam"
        f" --force"
        f" {input_for_prokka}"
    )
    cmd = "set -euo pipefail; " + "".join(parts)
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[-500:]}")
        sys.exit(f"Prokka failed with return code {result.returncode}")

    print(f"  GFF -> {cfg.gff_file}")
    print(f"  FAA -> {cfg.faa_file}")
    print("  Step 1 complete.\n")


def step02_extract_igrs(cfg: Config, force: bool = False):
    """Extract intergenic regions from the Prokka GFF and genome FASTA."""
    if not force and cfg.igr_summary.exists() and cfg.igr_fasta.exists():
        print("── IGR files already exist, skipping ──")
        print("  Step 2 complete.\n")
        return

    print("── Extracting intergenic regions ──")
    cfg.igr_dir.mkdir(parents=True, exist_ok=True)

    # Use the Prokka .fna (nucleotide FASTA with contig IDs matching the GFF)
    # but only if it actually contains contig-level sequences. Per-CDS
    # nucleotide files (.ffn/.fnn) use gene locus tags as headers rather
    # than contig accessions, so they won't match the GFF seqid column.
    fasta_source = cfg.input_fasta
    if cfg.fna_file.exists():
        # Read the contig IDs from the GFF (first column of data lines)
        gff_contigs = set()
        with open(str(cfg.gff_file)) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                gff_contigs.add(line.split("\t", 1)[0])
        # Check if any GFF contig appears in the FNA headers
        fna_ids = {rec.id for rec in SeqIO.parse(str(cfg.fna_file), "fasta")}
        if gff_contigs & fna_ids:
            fasta_source = cfg.fna_file
        else:
            print(f"  NOTE: FNA file ({cfg.fna_file.name}) does not contain "
                  f"contig sequences matching the GFF — using genome FASTA instead.")

    print(f"  GFF:   {cfg.gff_file}")
    print(f"  FASTA: {fasta_source}")
    igr_df = extract_igrs(cfg.gff_file, fasta_source,
                          size_min=cfg.igr_size_min, size_max=cfg.igr_size_max)

    # Warn if most sequences are empty (contig ID mismatch)
    if not igr_df.empty:
        n_empty = igr_df["sequence"].isna().sum() + (igr_df["sequence"] == "").sum()
        if n_empty > 0:
            print(f"  WARNING: {n_empty}/{len(igr_df)} IGRs have empty sequences.")
            print("    This usually means contig IDs in the GFF don't match "
                  "the FASTA headers.")
            print("    If you supplied an FNA file in the batch table, make "
                  "sure it contains contig sequences (not per-CDS nucleotides).")

    # Add 5'→3' oriented sequence: reverse-complement CO_R, keep others as-is
    if not igr_df.empty:
        igr_df["sequence_5p_to_3p"] = igr_df.apply(
            lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
            else r["sequence"],
            axis=1,
        )

    igr_df.to_csv(cfg.igr_summary, sep="\t", index=False)
    print(f"  IGR summary ({len(igr_df)} regions) -> {cfg.igr_summary}")

    # Write IGR FASTA
    with open(cfg.igr_fasta, "w") as fh:
        for _, row in igr_df.iterrows():
            fh.write(f">{row['igr_id']}_{row['contig']}_{row['orientation']}"
                     f"_{row['left_gene']}_{row['right_gene']}\n")
            fh.write(f"{row['sequence']}\n")
    print(f"  IGR FASTA -> {cfg.igr_fasta}")
    print("  Step 2 complete.\n")


def step03_identify_operons(cfg: Config, force: bool = False):
    """Identify operons from the Prokka GFF."""
    if not force and cfg.operon_file.exists():
        print("── Operon file already exists, skipping ──")
        print("  Step 3 complete.\n")
        return

    print("── Identifying operons ──")
    genes_df = _parse_gff_for_operons(str(cfg.gff_file))
    if genes_df.empty:
        print("  WARNING: no CDS features found in GFF")
        pd.DataFrame().to_csv(cfg.operon_file, sep="\t", index=False)
        print("  Step 3 complete.\n")
        return

    operons = _identify_operons(genes_df, cfg.max_internal_distance,
                                cfg.min_flanking_distance)
    source = cfg.prokka_prefix

    rows = []
    for op_idx, operon in enumerate(operons, start=1):
        operon_id = f"operon{op_idx}"
        first_genes = set()
        if operon[0]["strand"] == "+":
            first_genes.add((0, "yes_1"))
        if operon[-1]["strand"] == "-":
            first_genes.add((len(operon) - 1, "yes_2"))

        for gi, gene in enumerate(operon):
            gene_id = _extract_gene_id(gene["attributes"])
            label = "No"
            for (idx, tag) in first_genes:
                if gi == idx:
                    label = tag
                    break
            rows.append({
                "Operon": operon_id,
                "SeqID": gene["seqid"],
                "Start": gene["start"],
                "End": gene["end"],
                "Attributes": gene["attributes"],
                "Gene": gene_id,
                "SourceFile": source,
                "FirstGene": label,
            })

    df = pd.DataFrame(rows)
    df.to_csv(cfg.operon_file, sep="\t", index=False)
    print(f"  Operons ({len(operons)} operons, {len(df)} gene rows) -> {cfg.operon_file}")
    print("  Step 3 complete.\n")


def step04_run_hmmsearch(cfg: Config, force: bool = False):
    """Run hmmsearch for every individual .hmm profile in the profiles directory."""
    if cfg.hmm_profiles_dir is None or not cfg.hmm_profiles_dir.is_dir():
        print("── HMM profiles directory not available ──")
        print(f"  hmm_profiles_dir: {cfg.hmm_profiles_dir}")
        print("  Provide --hmm-dir or reinstall to restore bundled profiles.")
        print("  Step 4 complete.\n")
        return

    hmm_files = sorted(cfg.hmm_profiles_dir.glob("*.hmm"))
    if not hmm_files:
        print(f"── No .hmm files found in {cfg.hmm_profiles_dir} ──")
        print("  Step 4 complete.\n")
        return

    print(f"  Using {len(hmm_files)} HMM profiles from: {cfg.hmm_profiles_dir}")

    tblout_dir = cfg.hmm_dir / "tblout"
    tblout_dir.mkdir(parents=True, exist_ok=True)

    # Check if all tblout files already exist
    if not force and all(
        (tblout_dir / (f.stem + ".tblout")).exists() for f in hmm_files
    ):
        print("── HMM output already exists for all profiles, skipping ──")
        print("  Step 4 complete.\n")
        return

    print(f"── Running hmmsearch ({len(hmm_files)} profiles) ──")
    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)

    failed = 0
    for i, hmm_path in enumerate(hmm_files, 1):
        out_file = tblout_dir / (hmm_path.stem + ".tblout")
        log_file = tblout_dir / (hmm_path.stem + ".log")

        if not force and out_file.exists():
            continue

        parts = []
        if cfg.conda_env_hmm:
            parts.append(
                f'eval "$(conda shell.bash hook 2>/dev/null)"; '
                f"conda activate {cfg.conda_env_hmm}; "
            )
        parts.append(
            f"{cfg.hmmsearch_bin}"
            f" --tblout {out_file}"
            f" -o {log_file}"
            f" --noali --cpu {cfg.threads}"
            f" {hmm_path} {cfg.faa_file}"
        )
        cmd = "set -euo pipefail; " + "".join(parts)
        result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  WARNING: hmmsearch failed for {hmm_path.name}: "
                  f"{result.stderr[-200:]}")
            failed += 1

        if i % 100 == 0 or i == len(hmm_files):
            print(f"  {i}/{len(hmm_files)} profiles searched...")

    if failed:
        print(f"  {failed} profile(s) failed during hmmsearch.")

    print("  Step 4 complete.\n")


def step05_filter_hmm(cfg: Config, force: bool = False):
    """Consolidate and filter HMM hits from individual profile searches."""
    if not force and cfg.hmm_filtered.exists():
        print("── Filtered HMM file already exists, skipping ──")
        print("  Step 5 complete.\n")
        return

    tblout_dir = cfg.hmm_dir / "tblout"
    tblout_files = sorted(tblout_dir.glob("*.tblout")) if tblout_dir.is_dir() else []

    if not tblout_files:
        print("── No HMM output to filter, skipping ──")
        print("  Step 5 complete.\n")
        return

    print(f"── Filtering HMM output ({len(tblout_files)} tblout files) ──")
    all_rows = []
    for f in tblout_files:
        all_rows.extend(_parse_hmm_tblout(str(f)))

    if not all_rows:
        print("  No HMM hits found.")
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_combined, sep="\t", index=False)
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_filtered, sep="\t", index=False)
        print("  Step 5 complete.\n")
        return

    df = pd.DataFrame(all_rows, columns=_HMM_HEADER)
    df.to_csv(cfg.hmm_combined, sep="\t", index=False)
    print(f"  Combined HMM hits: {len(df)}")

    df["full_sequence_bitscore"] = pd.to_numeric(df["full_sequence_bitscore"], errors="coerce")
    df = df[df["full_sequence_bitscore"] >= cfg.hmm_bitscore_min]

    # Keep ALL matches — a gene may match multiple profiles and all are retained.
    df.to_csv(cfg.hmm_filtered, sep="\t", index=False)
    n_genes = df["target_name"].nunique() if not df.empty else 0
    print(f"  Filtered HMM ({len(df)} hits across {n_genes} genes) -> {cfg.hmm_filtered}")
    print("  Step 5 complete.\n")


def step06_filter_operons_add_markers(cfg: Config, force: bool = False):
    """Filter operons and merge with HMM marker gene info."""
    if not force and cfg.operon_filtered_markers.exists():
        print("── Filtered operons with markers already exists, skipping ──")
        print("  Step 6 complete.\n")
        return

    print("── Filtering operons ──")
    df = pd.read_csv(cfg.operon_file, sep="\t")
    if df.empty:
        print("  No operons to filter.")
        df.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)
        print("  Step 6 complete.\n")
        return

    # Keep operons on a single contig
    single_contig = df.groupby("Operon").filter(
        lambda g: g["SeqID"].nunique() == 1
    )
    # Keep operons with at least one promoter-boundary gene
    has_first = single_contig.groupby("Operon").filter(
        lambda g: g["FirstGene"].isin(["yes_1", "yes_2"]).any()
    )
    # Remove operons with BOTH yes_1 and yes_2 (ambiguous)
    def _single_direction(g):
        labels = set(g["FirstGene"].unique())
        return not {"yes_1", "yes_2"}.issubset(labels)

    filtered = has_first.groupby("Operon").filter(_single_direction)
    filtered.to_csv(cfg.operon_filtered, sep="\t", index=False)
    print(f"  Filtered operons ({len(filtered)} rows) -> {cfg.operon_filtered}")

    # Merge with HMM markers if available
    if cfg.hmm_filtered.exists() and cfg.hmm_filtered.stat().st_size > 0:
        print("── Adding HMM marker info ──")
        try:
            hmm = pd.read_csv(cfg.hmm_filtered, sep="\t",
                               usecols=["target_name", "accession2"])
        except (pd.errors.EmptyDataError, ValueError):
            hmm = pd.DataFrame(columns=["target_name", "accession2"])

        merged = pd.merge(filtered, hmm, how="left",
                           left_on="Gene", right_on="target_name")
        valid_operons = merged.loc[
            merged["accession2"].notna() &
            (merged["accession2"] != "") &
            (merged["accession2"] != 0),
            "Operon"
        ].unique()
        merged = merged[merged["Operon"].isin(valid_operons)]
        merged.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)
        print(f"  Operons with markers ({len(merged)} rows) "
              f"-> {cfg.operon_filtered_markers}")

        if len(merged) == 0 and not hmm.empty and not filtered.empty:
            operon_genes = set(filtered["Gene"].dropna())
            hmm_genes = set(hmm["target_name"].dropna())
            overlap = operon_genes & hmm_genes
            print()
            print("  WARNING: 0 marker operons found despite HMM hits.")
            print(f"    Operon gene IDs (first 5): "
                  f"{sorted(operon_genes)[:5]}")
            print(f"    HMM target names (first 5): "
                  f"{sorted(hmm_genes)[:5]}")
            print(f"    Overlap: {len(overlap)} gene(s)")
            if not overlap:
                print("    The gene IDs in the GFF and the protein IDs in "
                      "the FAA appear to use different naming schemes.")
                print("    Check that both files come from the same "
                      "annotation run.")
    else:
        # No HMM data: use all filtered operons as-is
        print("  No HMM data available; using all filtered operons as markers")
        filtered["target_name"] = filtered["Gene"]
        filtered["accession2"] = "no_hmm"
        filtered.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)

    print("  Step 6 complete.\n")


def step07_match_igrs_to_markers(cfg: Config, force: bool = False):
    """Match IGRs to marker operon genes."""
    if not force and cfg.promoter_markers.exists() and cfg.promoter_markers_hmm.exists():
        print("── Promoter markers file already exists, skipping ──")
        print("  Step 7 complete.\n")
        return

    print("── Matching IGRs to marker operons ──")
    igr = pd.read_csv(cfg.igr_summary, sep="\t")
    try:
        markers = pd.read_csv(cfg.operon_filtered_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  No marker data available.")
        pd.DataFrame().to_csv(cfg.promoter_markers, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
        print("  Step 7 complete.\n")
        return

    marker_genes = set(markers["Gene"].dropna())

    # Match: an IGR is a promoter candidate if one of its flanking genes
    # is a marker gene, and the orientation is promoter-relevant.
    igr = igr.copy()
    igr["marker_match"] = "none"

    # CO_F: right_gene (downstream on + strand) is the marker
    co_f = (igr["orientation"] == "CO_F") & igr["right_gene"].isin(marker_genes)
    igr.loc[co_f, "marker_match"] = "CO_F"

    # CO_R: left_gene (downstream on - strand, since gene runs right to left) is the marker
    co_r = (igr["orientation"] == "CO_R") & igr["left_gene"].isin(marker_genes)
    igr.loc[co_r, "marker_match"] = "CO_R"

    # DP: either flanking gene could be a marker
    dp_left = (igr["orientation"] == "DP") & igr["left_gene"].isin(marker_genes)
    dp_right = (igr["orientation"] == "DP") & igr["right_gene"].isin(marker_genes)
    igr.loc[dp_left | dp_right, "marker_match"] = "DP"

    matched = igr[igr["marker_match"] != "none"].copy()
    matched.to_csv(cfg.promoter_markers, sep="\t", index=False)
    print(f"  Matched promoter markers ({len(matched)} rows) -> {cfg.promoter_markers}")

    # Build expanded variant with one row per IGR × HMM profile match.
    # Determine which gene column is the marker gene for each matched IGR.
    matched["marker_gene"] = matched.apply(
        lambda r: r["right_gene"] if r["marker_match"] == "CO_F" else r["left_gene"],
        axis=1,
    )

    # Load HMM filtered hits for the profile names
    hmm_profiles = pd.DataFrame()
    if cfg.hmm_filtered.exists():
        try:
            hmm_profiles = pd.read_csv(
                cfg.hmm_filtered, sep="\t",
                usecols=["target_name", "query_name"],
            )
        except (pd.errors.EmptyDataError, ValueError):
            pass

    if hmm_profiles.empty:
        matched["hmm_profile"] = "no_hmm"
        expanded = matched
    else:
        # Inner merge: only keep IGRs whose marker gene has a direct HMM
        # hit.  Step 6 keeps entire operons when any gene in the operon
        # matches, so some IGRs flank operon-member genes that have no HMM
        # hit themselves.  Those are excluded here.
        expanded = pd.merge(
            matched, hmm_profiles,
            how="inner",
            left_on="marker_gene", right_on="target_name",
        )
        expanded.rename(columns={"query_name": "hmm_profile"}, inplace=True)
        expanded.drop(columns=["target_name"], inplace=True, errors="ignore")

    expanded.drop(columns=["marker_gene"], inplace=True, errors="ignore")
    expanded.to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
    print(f"  Expanded markers with HMM profiles ({len(expanded)} rows) "
          f"-> {cfg.promoter_markers_hmm}")
    print("  Step 7 complete.\n")


def _reverse_complement(seq) -> str:
    if not isinstance(seq, str) or not seq:
        return ""
    return str(Seq(seq).reverse_complement())


def _write_fasta(df, output_path, short_header=False):
    """Write promoter FASTA from a DataFrame with sequence_5p_to_3p."""
    with open(output_path, "w") as fh:
        for _, row in df.iterrows():
            seq = row["sequence_5p_to_3p"]
            igr_id = row["igr_id"]
            orient = row["orientation"]
            if short_header:
                header = f">{igr_id}_{orient}"
            else:
                header = (f">{igr_id}_{row['contig']}_{orient}"
                          f"_{row['left_gene']}_{row['right_gene']}")
            fh.write(f"{header}\n{seq}\n")
    print(f"  FASTA ({len(df)} seqs, short={short_header}) -> {output_path}")


def _load_contigs(fasta_path):
    """Load contig sequences from a genome FASTA into a dict."""
    # return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}
    with _open_maybe_gzip(fasta_path) as fh:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
    


def _extract_cds_start(row, contigs, n_bp):
    """Return the first *n_bp* nucleotides of the CDS downstream of a
    promoter, oriented 5'→3' on the coding strand.

    For CO_F the downstream CDS is the right gene, sitting on the +
    strand immediately after the IGR.  Its first *n_bp* nt run left-to-
    right starting at ``igr_end + 1``.

    For CO_R the downstream CDS is the left gene, sitting on the −
    strand immediately before the IGR.  Its coding sequence starts at
    ``igr_start - 1`` and runs right-to-left.  We extract the last
    *n_bp* nt of the forward strand up to that position and reverse-
    complement them.

    Returns an empty string if the contig is missing or there aren't
    enough nucleotides available.
    """
    contig_seq = contigs.get(row["contig"], "")
    if not contig_seq:
        return ""

    orient = row["orientation"]
    if orient == "CO_F":
        # CDS starts at igr_end + 1 (1-based), on the + strand.
        cds_start_0 = row["end"]            # 0-based start (= igr_end in 1-based)
        fragment = contig_seq[cds_start_0 : cds_start_0 + n_bp]
    elif orient == "CO_R":
        # CDS ends at igr_start - 1 (1-based) on the + strand, but the
        # gene is on the − strand, so its coding sequence begins at the
        # high coordinate and runs leftward.
        cds_end_0 = row["start"] - 1        # 0-based position (= igr_start - 1 in 1-based)
        start_0 = max(0, cds_end_0 - n_bp)
        fragment = _reverse_complement(contig_seq[start_0 : cds_end_0])
    else:
        return ""

    return fragment


def _add_cds_column(df, contigs, n_bp):
    """Add a ``cds_start_seq`` column with the first *n_bp* nt of CDS,
    and a ``sequence_5p_to_3p_cds`` column with promoter + CDS joined.
    """
    df["cds_start_seq"] = df.apply(
        lambda r: _extract_cds_start(r, contigs, n_bp), axis=1)
    df["sequence_5p_to_3p_cds"] = df["sequence_5p_to_3p"] + df["cds_start_seq"]
    return df


# =====================================================================
#  Step 8 — Motif-based promoter filtering
#
#  Scans promoter sequences for -10 and -35 hexamer motifs using PWMs
#  from the bundled all_unique_subgroups.meme file, which contains
#  paired M###_m10 / M###_m35 subgroups.
#
#  Each -10 hit is classified into one of four paths. Paths A–C
#  require a -35 hit with a 16-18 bp spacer to the -10 hit (canonical
#  σ⁷⁰ spacer is 17 ± 1); Path D requires only an extended -10. -10
#  motifs whose right edge sits more than cfg.max_dist_to_cds_start
#  bp from the downstream CDS start (default 200) are dropped before
#  the spacer search.
#    A — linked -10 and strict -35 (same subgroup)
#    B — extended -10 (TG dinucleotide at -2/-1 or -3/-2) and a
#        linked relaxed -35 (same subgroup)
#    C — unlinked -10 and strict -35 (different subgroups)
#    D — extended -10 with no -35 in the 16-18 bp window
#  Anything else is regarded as no hit.
#
#  Each hit also carries an empirical p-value, computed per PWM as
#  the fraction of 4^w k-mers scoring at least as high. The best hit
#  per IGR and the HTML report rank within each path group by
#  -log10(p10) + -log10(p35), a cross-PWM-comparable joint
#  significance. Path D rows use -log10(p10) alone.
# =====================================================================

def step08_scan_motifs(cfg: Config, force: bool = False):
    """Scan promoter sequences for -10/-35 motifs and classify as A/B/C/D.

    Writes into ``motifs/``:
    * ``motif_hits_all.tsv``       — all hits for all promoter-orientation IGRs
    * ``motif_best_all.tsv``       — best hit per IGR (all promoters)
    * ``motif_hits_markers.tsv``   — all hits for marker promoters only
    * ``motif_best_markers.tsv``   — best hit per IGR (markers only)
    * ``promoter_markers_verified.tsv`` — marker promoters confirmed by motif scan

    Writes at the top level of ``output_dir``:
    * ``all_promoters_verified.fasta``    — all motif-confirmed promoters,
      CO_R sequences reverse-complemented so every record is 5'→3'
    * ``marker_promoters_verified.fasta`` — marker + motif-confirmed promoters,
      CO_R sequences reverse-complemented so every record is 5'→3'
    """
    if not force and cfg.motif_best_all.exists() and cfg.promoter_markers_verified.exists():
        print("── Motif scan results already exist, skipping ──")
        print("  Step 8 complete.\n")
        return

    # Locate MEME file
    meme_file = cfg.meme_file
    if meme_file is None or not Path(meme_file).is_file():
        bundled = Path(__file__).parent / "all_unique_subgroups.meme"
        if bundled.is_file():
            meme_file = bundled
        else:
            print("── No MEME file found, skipping motif scan ──")
            pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
            print("  Step 8 complete.\n")
            return

    # Compute per-genome A/C/G/T background. The bundled MEME file's
    # header explicitly notes that "Per-organism A/C/G/T frequencies
    # used in the scan are in per_organism_background.tsv" — that file
    # is not shipped, so we derive the background directly from the
    # input genome (Prokka .fna when available, otherwise the original
    # input). PWM log-odds and the 4^w k-mer enumeration used for
    # threshold / p-value lookup are calibrated against this same
    # background, so significance is internally consistent with the
    # composition of the genome being scanned.
    bg_source = cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta
    bg = _compute_background_from_fasta(bg_source)

    # Load paired subgroup motifs
    p10 = cfg.motif_p10
    p35 = cfg.motif_p35
    p35_relaxed = cfg.motif_p35_relaxed
    print(f"── Loading motifs from {Path(meme_file).name} "
          f"(p10={p10}, p35={p35}, p35_relaxed={p35_relaxed}) ──")
    print(f"  Background ({Path(bg_source).name}): "
          f"A={bg[0]:.3f} C={bg[1]:.3f} G={bg[2]:.3f} T={bg[3]:.3f} "
          f"(GC={bg[1] + bg[2]:.3f})")
    m10, m35_strict, m35_relaxed = _load_motifs_from_file(
        meme_file, p10, p35, p35_relaxed, bg=bg)

    if not m10.entries:
        print("  No -10 motifs loaded — cannot scan.")
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Report threshold ranges across subgroups.
    t10 = [e[2] for e in m10.entries]
    ts = [e[2] for e in m35_strict.entries]
    tr = [e[2] for e in m35_relaxed.entries]
    if t10:
        print(f"  -10 score thresholds: {min(t10):.3f} – {max(t10):.3f} "
              f"(p < {p10}, n={len(t10)} subgroups)")
    if ts:
        print(f"  -35 strict thresholds: {min(ts):.3f} – {max(ts):.3f} "
              f"(p < {p35}, n={len(ts)} subgroups)")
    if tr:
        print(f"  -35 relaxed thresholds: {min(tr):.3f} – {max(tr):.3f} "
              f"(p < {p35_relaxed}, n={len(tr)} subgroups)")

    # --- Scan ALL promoter-orientation IGRs ---
    try:
        all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
        all_igr = all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        all_igr = pd.DataFrame()

    if all_igr.empty:
        print("── No promoter-orientation IGRs to scan ──")
        for p in [cfg.motif_hits_all, cfg.motif_best_all,
                  cfg.motif_hits_markers, cfg.motif_best_markers]:
            pd.DataFrame().to_csv(p, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Orient sequences 5'→3'
    all_igr["sequence_5p_to_3p"] = all_igr.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    print(f"\n  Scanning {len(all_igr)} promoter-orientation IGRs "
          f"(orientation-aware: + strand only for CO_F/CO_R; "
          f"max_dist_to_cds_start={cfg.max_dist_to_cds_start} bp)...")
    all_hits, best_all = _scan_sequences_for_motifs(
        all_igr, m10, m35_strict, m35_relaxed,
        max_dist_to_cds_start=cfg.max_dist_to_cds_start)

    all_hits.to_csv(cfg.motif_hits_all, sep="\t", index=False)
    best_all.to_csv(cfg.motif_best_all, sep="\t", index=False)
    n_confirmed_all = len(best_all)
    print(f"  All promoters: {n_confirmed_all}/{len(all_igr)} confirmed by motif scan")
    print(f"    All hits -> {cfg.motif_hits_all}")
    print(f"    Best hits -> {cfg.motif_best_all}")

    # Write verified all-promoter FASTA
    confirmed_ids = set(best_all["igr_id"])
    confirmed_all_df = all_igr[all_igr["igr_id"].isin(confirmed_ids)].copy()
    _write_fasta(confirmed_all_df, cfg.all_promoters_verified_fasta, short_header=False)

    # --- Scan MARKER promoters ---
    marker_igr_ids = set()
    markers_df = pd.DataFrame()
    try:
        markers_df = pd.read_csv(cfg.promoter_markers, sep="\t")
        markers_df = markers_df[markers_df["orientation"].isin(["CO_F", "CO_R"])].copy()
        marker_igr_ids = set(markers_df["igr_id"])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pass

    if markers_df.empty:
        print("  No marker promoters to scan.")
        for p in [cfg.motif_hits_markers, cfg.motif_best_markers]:
            pd.DataFrame().to_csv(p, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Orient marker sequences 5'→3'
    if "sequence_5p_to_3p" not in markers_df.columns:
        markers_df["sequence_5p_to_3p"] = markers_df.apply(
            lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
            else r["sequence"],
            axis=1,
        )

    print(f"\n  Scanning {len(markers_df)} marker promoters "
          f"(orientation-aware: + strand only for CO_F/CO_R; "
          f"max_dist_to_cds_start={cfg.max_dist_to_cds_start} bp)...")
    marker_hits, best_markers = _scan_sequences_for_motifs(
        markers_df, m10, m35_strict, m35_relaxed,
        max_dist_to_cds_start=cfg.max_dist_to_cds_start)

    marker_hits.to_csv(cfg.motif_hits_markers, sep="\t", index=False)
    best_markers.to_csv(cfg.motif_best_markers, sep="\t", index=False)

    # Write verified markers: marker promoters that passed motif scan
    verified_ids = set(best_markers["igr_id"])
    verified_df = markers_df[markers_df["igr_id"].isin(verified_ids)].copy()

    # Merge best motif hit info onto verified markers
    motif_info = best_markers.set_index("igr_id")[
        ["strand", "pos_10", "seq_10", "score_10", "pvalue_10", "source_10",
         "has_ext10", "pos_35", "seq_35", "score_35", "pvalue_35", "source_35",
         "spacer_len", "spacer_seq", "path"]
    ].rename(columns=lambda c: f"motif_{c}")
    verified_df = verified_df.merge(motif_info, left_on="igr_id",
                                     right_index=True, how="left")

    verified_df.to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
    n_confirmed_markers = len(verified_df)
    print(f"  Marker promoters: {n_confirmed_markers}/{len(markers_df)} confirmed")
    print(f"    Marker hits -> {cfg.motif_hits_markers}")
    print(f"    Best marker hits -> {cfg.motif_best_markers}")
    print(f"    Verified markers -> {cfg.promoter_markers_verified}")

    # Write verified marker promoter FASTA
    _write_fasta(verified_df, cfg.marker_promoters_verified_fasta, short_header=False)

    print("  Step 8 complete.\n")


# =====================================================================
#  Diagnostic — sweep -10 / -35 PWM p-value thresholds
# =====================================================================
#
# Why a sweep: choosing motif_p10 and motif_p35 is consequential — the
# benchmark variants in the manuscript (AB / ABC / lowerps / samep /
# lower35) all differ only in these two (or three, including
# p35_relaxed) numbers, and produce yield differences spanning an
# order of magnitude. The defensible answer to "what should p10 / p35
# be?" depends on the genome and on what the user values (precision
# vs recall, path-A purity vs path-B/C recovery). This module records
# what the pipeline would have output at a wide range of (p10, p35)
# settings without re-running Prokka, hmmsearch, or operon ID.
#
# Three sweep modes, run by default:
#   matched   — set p10 = p35 = p for each swept value p
#   p10_only  — vary p10 across the swept values, hold p35 at the
#               cfg default (or whatever the run was launched with)
#   p35_only  — vary p35 across the swept values, hold p10 at the
#               cfg default
#
# In every mode the relaxed -35 threshold (used only for Path B) is
# pinned to `p35 * p35_relaxed_ratio` (default ratio 20, matching the
# pipeline's default 2.5e-3 strict / 5e-2 relaxed). Capped at 1.0 so
# extreme sweep points don't blow past valid probability.
#
# All sweep points share the per-genome A/C/G/T background, computed
# once and passed in, so the comparison across sweep points is
# internally consistent.
# =====================================================================

def _logspace_p_values(p_min, p_max, n):
    """Return `n` p-values log-spaced inclusive of both endpoints. For
    `n == 1`, returns just `p_max`."""
    if n <= 1:
        return [p_max]
    log_min = math.log10(p_min)
    log_max = math.log10(p_max)
    step = (log_max - log_min) / (n - 1)
    return [10.0 ** (log_min + i * step) for i in range(n)]


# Column layout for the per-sweep-point profinder_results.tsv. Kept in
# lockstep with step10_final_table's output so downstream tooling
# (e.g. Section 5 benchmarking) can ingest sweep tables identically.
_FINAL_TABLE_COLUMNS = [
    "igr_id", "contig", "start", "end", "length", "orientation",
    "associated_cds", "gene_name", "locus_tag", "product",
    "is_marker", "motif_confirmed", "motif_path", "motif_strand",
    "motif_pos_10", "motif_seq_10", "motif_score_10", "motif_pvalue_10",
    "motif_source_10", "motif_has_ext10",
    "motif_pos_35", "motif_seq_35", "motif_score_35", "motif_pvalue_35",
    "motif_source_35",
    "motif_spacer_len", "motif_spacer_seq",
    "sequence_5p_to_3p",
]

# The motif scan emits these per-IGR best-hit columns; they get the
# `motif_` prefix in the final table.
_FINAL_TABLE_MOTIF_COLS = [
    "strand", "pos_10", "seq_10", "score_10", "pvalue_10", "source_10",
    "has_ext10", "pos_35", "seq_35", "score_35", "pvalue_35", "source_35",
    "spacer_len", "spacer_seq", "path",
]


def _assemble_final_table(igr_df, best_hits_df, marker_ids, ann_map,
                           cds_bp=0, contigs=None):
    """Build a step-10-equivalent ``profinder_results.tsv``-shaped
    DataFrame from pre-loaded inputs and a per-IGR best-motif-hits
    DataFrame.

    ``igr_df`` is expected to contain ``igr_id``, ``contig``,
    ``start``, ``end``, ``length``, ``orientation``, ``left_gene``,
    ``right_gene``, and ``sequence`` columns (typical output of
    ``cfg.igr_summary``). ``associated_cds`` and ``sequence_5p_to_3p``
    are derived if absent. ``ann_map`` maps gene IDs to
    ``{product, gene, locus_tag}`` dicts (empty dict for no annotation).
    ``marker_ids`` is a set of IGR IDs flagged ``is_marker=True``.
    ``best_hits_df`` is the second return value of
    ``_scan_sequences_for_motifs`` (one row per motif-confirmed IGR),
    or any DataFrame with at least ``igr_id`` plus the columns listed
    in ``_FINAL_TABLE_MOTIF_COLS``. Caller is expected to filter
    ``best_hits_df`` upstream when only a path subset (e.g.
    ``path in {"A","B"}``) should be marked confirmed.

    The function mirrors ``step10_final_table``'s column layout so
    sweep outputs are drop-in compatible with downstream tooling that
    reads ``profinder_results.tsv``.
    """
    igr = igr_df.copy()
    if "associated_cds" not in igr.columns:
        igr["associated_cds"] = igr.apply(_get_associated_gene, axis=1)
    if "sequence_5p_to_3p" not in igr.columns:
        igr["sequence_5p_to_3p"] = igr.apply(
            lambda r: _reverse_complement(r["sequence"])
            if r["orientation"] == "CO_R" else r["sequence"],
            axis=1)

    igr["gene_name"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("gene", ""))
    igr["locus_tag"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("locus_tag", ""))
    igr["product"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("product", ""))
    igr["is_marker"] = (igr["igr_id"].isin(marker_ids)
                        if marker_ids else False)

    motif_map = {}
    if best_hits_df is not None and len(best_hits_df) > 0:
        for _, row in best_hits_df.iterrows():
            motif_map[row["igr_id"]] = {
                c: row.get(c, "") for c in _FINAL_TABLE_MOTIF_COLS}

    igr["motif_confirmed"] = igr["igr_id"].isin(set(motif_map.keys()))
    for c in _FINAL_TABLE_MOTIF_COLS:
        col = f"motif_{c}"
        igr[col] = igr["igr_id"].map(
            lambda g, _c=c: motif_map.get(g, {}).get(_c, ""))

    columns = list(_FINAL_TABLE_COLUMNS)
    if cds_bp > 0 and contigs is not None:
        igr = _add_cds_column(igr, contigs, cds_bp)
        columns.append("sequence_5p_to_3p_cds")

    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)
    return out


def _load_annotation_map(annotations_tsv):
    """Read ``cds_annotations.tsv`` into a {gene_id -> {gene, locus_tag,
    product}} dict. Returns empty dict if the file is missing or empty."""
    if not annotations_tsv.exists():
        return {}
    try:
        df = pd.read_csv(annotations_tsv, sep="\t")
    except (pd.errors.EmptyDataError, KeyError):
        return {}
    out = {}
    for _, row in df.iterrows():
        out[row["gene_id"]] = {
            "product": str(row.get("product", "")),
            "gene": str(row.get("gene", "")),
            "locus_tag": str(row.get("locus_tag", "")),
        }
    return out


def _sweep_table_filename(mode, p10, p35, p35_relaxed):
    """Return a deterministic filename for a sweep-point TSV that
    encodes mode + (p10, p35, p35_relaxed). Filename-safe and round-
    trips through ``str.split``: e.g.
    ``sweep_matched__p10_1.00e-03__p35_1.00e-03__p35r_2.00e-02.tsv``."""
    return (
        f"sweep_{mode}"
        f"__p10_{p10:.2e}"
        f"__p35_{p35:.2e}"
        f"__p35r_{p35_relaxed:.2e}"
        ".tsv"
    )


def _scan_counts(best_hits_df, marker_ids):
    """Summarize a best-hits DataFrame into per-path counts, both
    overall and restricted to marker-anchored IGRs."""
    out = {"n_total": 0, "n_path_A": 0, "n_path_B": 0,
           "n_path_C": 0, "n_path_D": 0,
           "n_marker": 0, "n_marker_path_A": 0, "n_marker_path_B": 0,
           "n_marker_path_C": 0, "n_marker_path_D": 0}
    if best_hits_df is None or len(best_hits_df) == 0:
        return out
    out["n_total"] = len(best_hits_df)
    path_counts = best_hits_df["path"].value_counts().to_dict()
    for k in ("A", "B", "C", "D"):
        out[f"n_path_{k}"] = int(path_counts.get(k, 0))
    if marker_ids:
        marker_mask = best_hits_df["igr_id"].isin(marker_ids)
        marker_df = best_hits_df[marker_mask]
        out["n_marker"] = int(len(marker_df))
        mpath = marker_df["path"].value_counts().to_dict()
        for k in ("A", "B", "C", "D"):
            out[f"n_marker_path_{k}"] = int(mpath.get(k, 0))
    return out


def _render_p_sweep_figure(sweep_df, out_path):
    """Optional 3-panel figure: one subplot per mode, x = swept p,
    y = stacked counts of marker-anchored hits by path (A / B / C / D).
    Returns False (and skips) if matplotlib is not available."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  matplotlib unavailable ({e}); skipping figure")
        return False

    modes = list(sweep_df["mode"].drop_duplicates())
    n_modes = len(modes)
    fig, axes = plt.subplots(1, max(1, n_modes),
                              figsize=(4.2 * max(1, n_modes), 3.5),
                              sharey=True, squeeze=False)
    axes = axes[0]
    path_colors = {"A": "#2e7d32", "B": "#e65100",
                    "C": "#c62828", "D": "#546e7a"}
    for ax, mode in zip(axes, modes):
        sub = sweep_df[sweep_df["mode"] == mode].sort_values("sweep_p")
        if len(sub) == 0:
            ax.set_title(f"{mode} (no data)")
            continue
        x = sub["sweep_p"].values
        bottom = [0] * len(sub)
        for path in ("A", "B", "C", "D"):
            y = sub[f"n_marker_path_{path}"].values
            ax.bar(range(len(sub)), y, bottom=bottom,
                   color=path_colors[path], label=f"path {path}",
                   width=0.85, edgecolor="white", linewidth=0.5)
            bottom = [b + v for b, v in zip(bottom, y)]
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels([f"{p:.1e}" for p in x],
                            rotation=45, ha="right", fontsize=7)
        ax.set_title(mode, fontsize=10)
        ax.set_xlabel("swept p-value")
        ax.set_xscale("linear")
    axes[0].set_ylabel("marker-anchored hits (count)")
    axes[-1].legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Figure -> {out_path}")
    return True


def diagnose_p_sweep(cfg: Config, sweep_min=1e-5, sweep_max=1e-1,
                      sweep_n=13, modes=("matched", "p10_only", "p35_only"),
                      p35_relaxed_ratio=20.0, force: bool = False):
    """Sweep -10 / -35 PWM p-value thresholds and record per-path
    promoter counts at each sweep point. Reuses cached outputs from
    steps 2 (igr_summary.tsv) and 7 (promoter_markers.tsv). Writes
    ``diagnostics/p_sweep.tsv`` and, when matplotlib is available,
    ``diagnostics/p_sweep.png``. Each sweep point reuses the same
    per-genome A/C/G/T background; only the score thresholds change.
    Outputs are deterministic for the same inputs and CLI flags."""
    import time

    out_tsv = cfg.diagnostics_dir / "p_sweep.tsv"
    out_png = cfg.diagnostics_dir / "p_sweep.png"
    if out_tsv.exists() and not force:
        print(f"── p-value sweep diagnostic already exists, skipping ──")
        print(f"  -> {out_tsv}")
        print("  Diagnostic complete.\n")
        return

    cfg.diagnostics_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = cfg.diagnostics_dir / "sweep_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Locate MEME file (same logic as step08)
    meme_file = cfg.meme_file
    if meme_file is None or not Path(meme_file).is_file():
        bundled = Path(__file__).parent / "all_unique_subgroups.meme"
        if bundled.is_file():
            meme_file = bundled
        else:
            print("── No MEME file found, skipping p-sweep diagnostic ──")
            return

    # Load inputs once
    try:
        all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
        all_igr = all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        all_igr = pd.DataFrame()
    if all_igr.empty:
        print("── No promoter-orientation IGRs to scan; "
              "run steps 1-2 first ──")
        return

    if "sequence_5p_to_3p" not in all_igr.columns:
        all_igr["sequence_5p_to_3p"] = all_igr.apply(
            lambda r: _reverse_complement(r["sequence"])
            if r["orientation"] == "CO_R" else r["sequence"], axis=1)
    if "associated_cds" not in all_igr.columns:
        all_igr["associated_cds"] = all_igr.apply(_get_associated_gene, axis=1)

    marker_ids = set()
    try:
        mdf = pd.read_csv(cfg.promoter_markers, sep="\t")
        mdf = mdf[mdf["orientation"].isin(["CO_F", "CO_R"])]
        marker_ids = set(mdf["igr_id"])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  WARNING: no marker IGR file (run step 7 first); "
              "marker counts will be 0.")

    # CDS annotations + contigs are loaded once and shared across every
    # sweep point. Both are p-value-independent.
    ann_map = _load_annotation_map(cfg.cds_annotations)
    contigs = (_load_contigs(cfg.input_fasta) if cfg.cds_bp > 0 else None)

    # Per-genome background (computed once, reused at every sweep point)
    bg_source = cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta
    bg = _compute_background_from_fasta(bg_source)

    sweep_ps = _logspace_p_values(sweep_min, sweep_max, sweep_n)
    print(f"── p-value sweep diagnostic ──")
    print(f"  IGRs to scan: {len(all_igr)}  (markers: {len(marker_ids)})")
    print(f"  Background ({Path(bg_source).name}): "
          f"A={bg[0]:.3f} C={bg[1]:.3f} G={bg[2]:.3f} T={bg[3]:.3f}")
    print(f"  Sweep: {sweep_n} log-spaced points in "
          f"[{sweep_min:.1e}, {sweep_max:.1e}]")
    print(f"  Modes: {', '.join(modes)}")
    print(f"  p35_relaxed_ratio = {p35_relaxed_ratio} "
          f"(p35_relaxed = ratio x p35, capped at 1.0)")
    print(f"  Spacer window: {_SPACER_MIN}-{_SPACER_MAX} bp; "
          f"max_dist_to_cds_start = {cfg.max_dist_to_cds_start} bp")
    print(f"  Per-sweep-point profinder_results.tsv files written to:")
    print(f"    {tables_dir}/")
    print(f"  Each table is filtered to motif_confirmed = True AND "
          f"motif_path in {{A, B}}.")

    records = []
    for mode in modes:
        print(f"\n  -- mode: {mode} --")
        for p in sweep_ps:
            if mode == "matched":
                p10, p35 = p, p
            elif mode == "p10_only":
                p10, p35 = p, cfg.motif_p35
            elif mode == "p35_only":
                p10, p35 = cfg.motif_p10, p
            else:
                print(f"    unknown mode {mode!r}, skipping")
                continue
            p35_relaxed = min(1.0, p35 * p35_relaxed_ratio)

            t0 = time.perf_counter()
            m10, m35s, m35r = _load_motifs_from_file(
                meme_file, p10, p35, p35_relaxed, bg=bg)
            _, best = _scan_sequences_for_motifs(
                all_igr, m10, m35s, m35r,
                max_dist_to_cds_start=cfg.max_dist_to_cds_start)
            elapsed = time.perf_counter() - t0

            counts = _scan_counts(best, marker_ids)

            # Build the path-A/B-only profinder_results.tsv for this
            # sweep point. Filtering best_hits BEFORE assembly ensures
            # motif_confirmed = False for any IGR whose best hit was
            # path C or D — those rows are then dropped from the table
            # to satisfy the user's "motif_confirmed = True AND
            # motif_path in {A, B}" requirement.
            if best is not None and len(best) > 0 and "path" in best.columns:
                best_ab = best[best["path"].isin({"A", "B"})].copy()
            else:
                best_ab = best.iloc[0:0] if best is not None else pd.DataFrame()
            final = _assemble_final_table(
                all_igr, best_ab, marker_ids, ann_map,
                cds_bp=cfg.cds_bp, contigs=contigs)
            final_filtered = final[final["motif_confirmed"]].copy()
            sweep_table_path = (tables_dir /
                                _sweep_table_filename(mode, p10, p35, p35_relaxed))
            final_filtered.to_csv(sweep_table_path, sep="\t", index=False)

            row = {"mode": mode, "sweep_p": p,
                    "p10": p10, "p35": p35, "p35_relaxed": p35_relaxed,
                    "bg_a": bg[0], "bg_c": bg[1], "bg_g": bg[2], "bg_t": bg[3],
                    "n_igrs_scanned": len(all_igr),
                    "elapsed_s": elapsed,
                    "sweep_table_path": str(sweep_table_path.relative_to(
                        cfg.output_dir)),
                    "n_in_sweep_table": int(len(final_filtered))}
            row.update(counts)
            records.append(row)
            print(f"    p={p:.2e}  p10={p10:.2e} p35={p35:.2e} "
                  f"p35r={p35_relaxed:.2e}  "
                  f"total={counts['n_total']} "
                  f"(A={counts['n_path_A']} B={counts['n_path_B']} "
                  f"C={counts['n_path_C']} D={counts['n_path_D']})  "
                  f"marker={counts['n_marker']}  "
                  f"AB_table={len(final_filtered)}  "
                  f"{elapsed:.1f}s")

    sweep_df = pd.DataFrame(records)
    cols = ["mode", "sweep_p", "p10", "p35", "p35_relaxed",
            "bg_a", "bg_c", "bg_g", "bg_t",
            "n_igrs_scanned",
            "n_total", "n_path_A", "n_path_B", "n_path_C", "n_path_D",
            "n_marker", "n_marker_path_A", "n_marker_path_B",
            "n_marker_path_C", "n_marker_path_D",
            "n_in_sweep_table", "sweep_table_path",
            "elapsed_s"]
    sweep_df = sweep_df[cols]
    sweep_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n  Sweep summary ({len(sweep_df)} rows) -> {out_tsv}")
    print(f"  Per-sweep-point tables ({len(sweep_df)} files) in "
          f"{tables_dir}/")

    try:
        _render_p_sweep_figure(sweep_df, out_png)
    except Exception as e:
        print(f"  Figure rendering failed: {e}")

    print("  Diagnostic complete.\n")


# =====================================================================
#  Step 9 — Annotate associated CDS from Prokka
# =====================================================================

def _parse_prokka_gff_annotations(gff_path, gene_ids=None):
    """Extract CDS annotations from a Prokka GFF.

    Parameters
    ----------
    gff_path : Path
        Prokka GFF3 file.
    gene_ids : set or None
        If provided, only return entries whose ID is in this set.
        If None, return all CDS entries.

    Returns
    -------
    dict
        Mapping gene_id -> {"product": str, "gene": str, "locus_tag": str}.
    """
    annotations = {}
    with open(str(gff_path)) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = cols[8]
            product = "hypothetical protein"
            gene_name = ""
            locus_tag = ""
            raw_id = ""
            for attr in attrs.split(";"):
                if attr.startswith("ID=") and not raw_id:
                    raw_id = attr[3:]
                elif attr.startswith("product="):
                    product = attr[8:]
                elif attr.startswith("gene="):
                    gene_name = attr[5:]
                elif attr.startswith("locus_tag="):
                    locus_tag = attr[10:]
            gene_id = _extract_gene_id(attrs)
            if gene_ids is None or gene_id in gene_ids:
                annotations[gene_id] = {
                    "product": product,
                    "gene": gene_name,
                    "locus_tag": locus_tag,
                }
    return annotations


def step09_annotate_cds(cfg: Config, force: bool = False):
    """Extract CDS annotations (product, gene name, locus tag) from Prokka
    GFF for ALL promoter-orientation IGRs."""
    annotation_tsv = cfg.cds_annotations
    if not force and annotation_tsv.exists():
        print("── CDS annotations already exist, skipping ──")
        print("  Step 9 complete.\n")
        return

    # Collect CDS gene IDs from ALL promoter-orientation IGRs.
    # For CO_F the associated CDS is right_gene; for CO_R it's left_gene.
    gene_ids = set()
    if cfg.igr_summary.exists():
        try:
            df = pd.read_csv(cfg.igr_summary, sep="\t")
            df = df[df["orientation"].isin(["CO_F", "CO_R"])]
            gene_ids.update(df.loc[df["orientation"] == "CO_F", "right_gene"].dropna())
            gene_ids.update(df.loc[df["orientation"] == "CO_R", "left_gene"].dropna())
        except (pd.errors.EmptyDataError, KeyError):
            pass

    if not gene_ids:
        print("── No promoter-associated CDS to annotate ──")
        pd.DataFrame(columns=["gene_id", "product", "gene", "locus_tag"]).to_csv(
            annotation_tsv, sep="\t", index=False)
        print("  Step 9 complete.\n")
        return

    print(f"── Annotating {len(gene_ids)} CDS from Prokka GFF ──")
    ann_map = _parse_prokka_gff_annotations(cfg.gff_file, gene_ids)

    rows = []
    for gid in sorted(gene_ids):
        info = ann_map.get(gid, {})
        rows.append({
            "gene_id": gid,
            "product": info.get("product", "hypothetical protein"),
            "gene": info.get("gene", ""),
            "locus_tag": info.get("locus_tag", ""),
        })
    result_df = pd.DataFrame(rows)
    result_df.to_csv(annotation_tsv, sep="\t", index=False)
    print(f"  Annotated {len(result_df)} CDS -> {annotation_tsv}")
    print("  Step 9 complete.\n")


# =====================================================================
#  Step 11 — HTML report
# =====================================================================

def _get_associated_gene(row):
    """Return the gene ID of the CDS immediately downstream of this promoter."""
    if row["orientation"] == "CO_F":
        return row.get("right_gene", "")
    elif row["orientation"] == "CO_R":
        return row.get("left_gene", "")
    return ""


def _build_motif_diagram_svg(seq_len, motif_hits, width=700, height=50):
    """Build an inline SVG showing motif positions along a promoter sequence.

    Returns an SVG string.
    """
    if seq_len == 0:
        return ""

    margin = 10
    track_y = 25
    track_h = 10
    bar_w = width - 2 * margin
    scale = bar_w / seq_len

    # Colour palette for motif elements
    db_colours = {
        "minus35": "#1565c0",
        "minus10": "#c62828",
        "ext10": "#ff8f00",
    }
    default_colour = "#9C27B0"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'style="background:#fafafa;border:1px solid #ddd;border-radius:4px;">',
        # Sequence track (grey bar)
        f'<rect x="{margin}" y="{track_y}" width="{bar_w}" height="{track_h}" '
        f'fill="#e0e0e0" rx="3"/>',
        # 5' and 3' labels
        f'<text x="{margin}" y="16" font-size="11" fill="#666">5\'</text>',
        f'<text x="{width - margin - 10}" y="16" font-size="11" fill="#666">3\'</text>',
        # Length label
        f'<text x="{width / 2}" y="16" font-size="10" fill="#999" text-anchor="middle">'
        f'{seq_len} bp</text>',
    ]

    for hit in motif_hits:
        start = hit.get("start", 0)
        stop = hit.get("stop", start + 10)
        db = hit.get("motif_database", "").lower()
        motif_name = hit.get("motif_id", "")
        colour = db_colours.get(db, default_colour)

        x = margin + (start - 1) * scale
        w = max((stop - start + 1) * scale, 3)  # at least 3px wide

        parts.append(
            f'<rect x="{x:.1f}" y="{track_y - 2}" width="{w:.1f}" height="{track_h + 4}" '
            f'fill="{colour}" opacity="0.7" rx="2">'
            f'<title>{html_mod.escape(motif_name)} ({db}) pos {start}-{stop}</title>'
            f'</rect>'
        )

    parts.append('</svg>')
    return "\n".join(parts)


def step11_generate_report(cfg: Config, force: bool = False):
    """Generate an HTML report combining promoters, CDS annotations,
    and motif hit positions with visual sequence diagrams."""
    if not force and cfg.report_html.exists():
        print("── HTML report already exists, skipping ──")
        print("  Step 11 complete.\n")
        return

    print("── Generating HTML report ──")

    # Load verified promoters (after motif scan).  Step 8 may have
    # written an empty placeholder file on a bailout path, so fall back
    # to promoter_markers.tsv if the verified file is empty or lacks the
    # expected columns.
    candidates = []
    if cfg.promoter_markers_verified.exists():
        candidates.append(cfg.promoter_markers_verified)
    if cfg.promoter_markers.exists():
        candidates.append(cfg.promoter_markers)

    promoter_df = None
    for source_file in candidates:
        try:
            df = pd.read_csv(source_file, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if "orientation" not in df.columns:
            continue
        df = df[df["orientation"].isin(["CO_F", "CO_R"])].copy()
        if not df.empty:
            promoter_df = df
            break

    if promoter_df is None or promoter_df.empty:
        print("  No verified promoter data available for report.")
        print("  Step 11 complete.\n")
        return

    # Count all promoter-orientation IGRs for reference
    n_all_promoters = 0
    if cfg.igr_summary.exists():
        try:
            all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
            n_all_promoters = len(all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])])
        except pd.errors.EmptyDataError:
            pass

    # Add associated gene column
    promoter_df["associated_gene"] = promoter_df.apply(_get_associated_gene, axis=1)

    # Orient sequences 5'→3'
    promoter_df["sequence_5p_to_3p"] = promoter_df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # Load CDS annotations (Prokka product names)
    annotation_map = {}
    annotation_tsv = cfg.cds_annotations
    if annotation_tsv.exists():
        try:
            ann_df = pd.read_csv(annotation_tsv, sep="\t")
            for _, row in ann_df.iterrows():
                annotation_map[row["gene_id"]] = str(row.get("product", ""))
        except (pd.errors.EmptyDataError, KeyError):
            pass

    # Build motif hit info from verified markers (motif_pos_10, motif_pos_35, etc.)
    # These columns were added by step 8 when writing promoter_markers_verified.tsv
    has_motif_cols = "motif_pos_10" in promoter_df.columns

    # Sort rows by motif path priority (A > B > C > D > other), then
    # within each group by the combined empirical significance of the
    # -10 and -35 hits, descending. We use -log10(p10) + -log10(p35)
    # where each p-value is the fraction of 4^w k-mers that score at
    # least as high as the observed hit under that PWM. The log-sum
    # is equivalent to -log10(p10 * p35), the joint significance under
    # independence, and is comparable across PWMs with differing
    # information content. Path D rows have no -35 and are ranked
    # within the D group by -log10(p10) alone.
    if has_motif_cols:
        def _path_rank(v):
            return {"A": 0, "B": 1, "C": 2, "D": 3}.get(str(v).strip(), 4)

        def _as_float(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        def _combined_neg_log10p(row):
            p10 = _as_float(row.get("motif_pvalue_10"))
            if p10 is None or p10 <= 0:
                return float("-inf")
            val = -math.log10(p10)
            p35 = _as_float(row.get("motif_pvalue_35"))
            if p35 is not None and p35 > 0:
                val += -math.log10(p35)
            return val

        promoter_df = promoter_df.assign(
            _path_rank=promoter_df["motif_path"].map(_path_rank),
            _score_key=promoter_df.apply(_combined_neg_log10p, axis=1),
        )
        promoter_df = promoter_df.sort_values(
            by=["_path_rank", "_score_key"], ascending=[True, False]
        ).drop(columns=["_path_rank", "_score_key"])

    # Build HTML
    n_promoters = len(promoter_df)
    n_annotated = sum(1 for _, r in promoter_df.iterrows()
                      if r["associated_gene"] in annotation_map
                      and annotation_map[r["associated_gene"]] != "hypothetical protein")
    n_with_motifs = sum(1 for _, r in promoter_df.iterrows()
                        if has_motif_cols and pd.notna(r.get("motif_pos_10"))
                        and str(r.get("motif_pos_10")) != ".")

    genome_name = cfg.input_fasta.stem

    # Paths for file links (relative to output_dir).
    all_verified_name = (cfg.all_promoters_verified_fasta.name
                         if cfg.all_promoters_verified_fasta.exists() else "")
    marker_verified_name = (cfg.marker_promoters_verified_fasta.name
                            if cfg.marker_promoters_verified_fasta.exists() else "")
    cds_ann_name = cfg.cds_annotations.name if cfg.cds_annotations.exists() else ""
    final_table_name = cfg.final_table.name if cfg.final_table.exists() else ""

    domain_label = "Archaea" if cfg.is_archaea else "Bacteria"

    html_parts = [_REPORT_HTML_HEAD.format(
        genome_name=html_mod.escape(genome_name),
        n_marker_promoters=n_promoters,
        n_all_promoters=n_all_promoters,
        n_annotated=n_annotated,
        n_with_motifs=n_with_motifs,
        all_verified_fasta=html_mod.escape(all_verified_name),
        marker_verified_fasta=html_mod.escape(marker_verified_name),
        cds_annotations_tsv=html_mod.escape(cds_ann_name),
        final_table=html_mod.escape(final_table_name),
        domain_label=domain_label,
    )]

    # Legend for motif diagram + path colour code
    html_parts.append("""
    <div class="legend">
        <strong>Motif path:</strong>
        <span class="legend-item"><span class="path path-a">A</span> linked &minus;10/&minus;35 (same subgroup), 15&ndash;19 bp</span>
        <span class="legend-item"><span class="path path-b">B</span> extended &minus;10 + linked relaxed &minus;35 (same subgroup), 15&ndash;19 bp</span>
        <span class="legend-item"><span class="path path-c">C</span> unlinked &minus;10/&minus;35 (different subgroups), 15&ndash;19 bp</span>
        <span class="legend-item"><span class="path path-d">D</span> extended &minus;10, no &minus;35 in window</span>
    </div>
    <div class="legend">
        <strong>Motif diagram:</strong>
        <span class="legend-item"><span class="legend-swatch" style="background:#1565c0"></span>&minus;35 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#c62828"></span>&minus;10 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#ff8f00"></span>Extended &minus;10 (TG)</span>
    </div>
    <p style="font-size:0.8rem; color:#777; margin-bottom:16px;">
        Rows are sorted by path priority (A &rarr; B &rarr; C), then by &minus;10 score.
    </p>
    """)

    # Table
    html_parts.append("""
    <table>
    <thead>
    <tr>
        <th>Promoter ID</th>
        <th>Contig</th>
        <th>Position</th>
        <th>Length</th>
        <th>Orientation</th>
        <th>Motif path</th>
        <th>&minus;10</th>
        <th>&minus;35</th>
        <th title="-log10(p-10) + -log10(p-35), i.e. -log10(p-10 × p-35). Higher = more significant. The within-path-group sort key.">Combined &minus;log<sub>10</sub>(p)</th>
        <th>Spacer</th>
        <th>Associated CDS</th>
        <th>Protein</th>
        <th>Sequence diagram</th>
        <th>Sequence (5'&rarr;3')</th>
    </tr>
    </thead>
    <tbody>
    """)

    for _, row in promoter_df.iterrows():
        igr_id = row["igr_id"]
        gene = row["associated_gene"]
        product = annotation_map.get(gene, "")

        # Build motif diagram from -10 and -35 positions
        motif_hits = []
        if has_motif_cols:
            pos_10 = row.get("motif_pos_10", ".")
            pos_35 = row.get("motif_pos_35", ".")
            if pos_10 != "." and pd.notna(pos_10):
                try:
                    p10 = int(float(pos_10))
                    motif_hits.append({
                        "motif_id": "-10",
                        "start": p10 + 1,  # 0-based to 1-based
                        "stop": p10 + _MOTIF_WIDTH,
                        "motif_database": "minus10",
                    })
                    # Extended -10 TG.  The column may be a Python bool
                    # (unread DataFrame), the string "True"/"False" (TSV
                    # round-trip), or missing altogether.
                    ext10_raw = row.get("motif_has_ext10", False)
                    if isinstance(ext10_raw, bool):
                        has_ext10_flag = ext10_raw
                    else:
                        has_ext10_flag = str(ext10_raw).strip().lower() == "true"
                    if has_ext10_flag:
                        motif_hits.append({
                            "motif_id": "ext-10 (TG)",
                            "start": p10 - 1,  # TG is 2 nt before -10
                            "stop": p10,
                            "motif_database": "ext10",
                        })
                except (ValueError, TypeError):
                    pass

            if pos_35 != "." and pd.notna(pos_35):
                try:
                    p35 = int(float(pos_35))
                    motif_hits.append({
                        "motif_id": "-35",
                        "start": p35 + 1,
                        "stop": p35 + _MOTIF_WIDTH,
                        "motif_database": "minus35",
                    })
                except (ValueError, TypeError):
                    pass

        svg = _build_motif_diagram_svg(row["length"], motif_hits)

        orient_class = "co-f" if row["orientation"] == "CO_F" else "co-r"
        seq = row.get("sequence_5p_to_3p", "")

        # Motif details
        motif_path = str(row.get("motif_path", "")).strip() if has_motif_cols else ""
        path_key = motif_path if motif_path in ("A", "B", "C", "D") else ""
        if path_key:
            path_cell = f'<span class="path path-{path_key.lower()}">{path_key}</span>'
            row_class = f"row-{path_key.lower()}"
        else:
            path_cell = '<span class="path path-none">&mdash;</span>'
            row_class = ""
        seq_10 = str(row.get("motif_seq_10", "")) if has_motif_cols else ""
        score_10 = row.get("motif_score_10", "") if has_motif_cols else ""
        source_10 = str(row.get("motif_source_10", "")) if has_motif_cols else ""
        seq_35 = str(row.get("motif_seq_35", "")) if has_motif_cols else ""
        score_35 = row.get("motif_score_35", "") if has_motif_cols else ""
        source_35 = str(row.get("motif_source_35", "")) if has_motif_cols else ""
        spacer_len = row.get("motif_spacer_len", "") if has_motif_cols else ""

        minus10_cell = f"{html_mod.escape(seq_10)}" if seq_10 and seq_10 != "." else '<span class="na">—</span>'
        if score_10 and str(score_10) != ".":
            minus10_cell += f' <small>({score_10}, {html_mod.escape(source_10)})</small>'

        minus35_cell = f"{html_mod.escape(seq_35)}" if seq_35 and seq_35 != "." else '<span class="na">—</span>'
        if score_35 and str(score_35) != ".":
            minus35_cell += f' <small>({score_35}, {html_mod.escape(source_35)})</small>'

        spacer_cell = str(spacer_len) if spacer_len and str(spacer_len) != "." else '<span class="na">—</span>'

        # Combined -log10(p-10) + -log10(p-35) — the within-group
        # sort key. Shows -log10(p-10) alone for Path D (no -35) so
        # rows stay rankable within the D group.
        def _to_float(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
        pvalue_10 = row.get("motif_pvalue_10", "") if has_motif_cols else ""
        pvalue_35 = row.get("motif_pvalue_35", "") if has_motif_cols else ""
        p10_f = _to_float(pvalue_10)
        p35_f = _to_float(pvalue_35)
        if p10_f is not None and p10_f > 0:
            combined_val = -math.log10(p10_f)
            if p35_f is not None and p35_f > 0:
                combined_val += -math.log10(p35_f)
            combined_cell = f"{combined_val:.2f}"
        else:
            combined_cell = '<span class="na">—</span>'

        product_cell = html_mod.escape(product) if product else '<span class="na">—</span>'

        html_parts.append(f"""
        <tr class="{row_class}">
            <td><code>{html_mod.escape(igr_id)}</code></td>
            <td>{html_mod.escape(str(row['contig']))}</td>
            <td>{row['start']:,}–{row['end']:,}</td>
            <td>{row['length']}</td>
            <td><span class="orient {orient_class}">{row['orientation']}</span></td>
            <td>{path_cell}</td>
            <td>{minus10_cell}</td>
            <td>{minus35_cell}</td>
            <td>{combined_cell}</td>
            <td>{spacer_cell}</td>
            <td><code>{html_mod.escape(str(gene))}</code></td>
            <td>{product_cell}</td>
            <td>{svg}</td>
            <td><code class="seq">{html_mod.escape(str(seq))}</code></td>
        </tr>
        """)

    html_parts.append("</tbody></table>")

    html_parts.append(_REPORT_HTML_FOOT)

    with open(cfg.report_html, "w") as fh:
        fh.write("\n".join(html_parts))
    print(f"  Report -> {cfg.report_html}")
    print("  Step 11 complete.\n")


# ── HTML template fragments ──────────────────────────────────────────

_REPORT_HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Promoter Report — {genome_name}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        color: #333; background: #fff; padding: 24px 32px;
        line-height: 1.5;
    }}
    h1 {{ font-size: 1.6rem; margin-bottom: 8px; }}
    h2 {{ font-size: 1.2rem; margin: 32px 0 12px; border-bottom: 2px solid #eee; padding-bottom: 6px; }}
    .summary {{
        display: flex; gap: 24px; margin: 16px 0 24px;
        flex-wrap: wrap;
    }}
    .summary .stat {{
        background: #f5f7fa; border-radius: 8px;
        padding: 12px 20px; min-width: 140px;
    }}
    .summary .stat .label {{ font-size: 0.8rem; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }}
    .summary .stat .value {{ font-size: 1.5rem; font-weight: 600; color: #1a73e8; }}
    .legend {{
        margin: 12px 0 20px; font-size: 0.85rem; color: #555;
    }}
    .legend-item {{ margin-left: 16px; }}
    .legend-swatch {{
        display: inline-block; width: 14px; height: 14px;
        border-radius: 3px; vertical-align: middle; margin-right: 4px;
    }}
    table {{
        border-collapse: collapse; width: 100%; font-size: 0.85rem;
        margin-bottom: 24px;
    }}
    th, td {{
        padding: 8px 10px; text-align: left;
        border-bottom: 1px solid #e8e8e8;
    }}
    th {{ background: #f5f7fa; font-weight: 600; position: sticky; top: 0; }}
    tr:hover {{ background: #f9fbfd; }}
    code {{ font-size: 0.82rem; background: #f0f0f0; padding: 1px 5px; border-radius: 3px; }}
    code.seq {{
        display: inline-block; max-width: 280px; overflow-x: auto;
        white-space: nowrap; font-size: 0.72rem; letter-spacing: 0.5px;
    }}
    .orient {{
        display: inline-block; padding: 2px 8px; border-radius: 4px;
        font-size: 0.78rem; font-weight: 600;
    }}
    .co-f {{ background: #e3f2fd; color: #1565c0; }}
    .co-r {{ background: #fce4ec; color: #c62828; }}
    .path {{
        display: inline-block; padding: 2px 10px; border-radius: 10px;
        font-size: 0.78rem; font-weight: 700; letter-spacing: 0.3px;
        min-width: 22px; text-align: center;
    }}
    .path-a {{ background: #e8f5e9; color: #2e7d32; }}
    .path-b {{ background: #fff8e1; color: #e65100; }}
    .path-c {{ background: #ffebee; color: #c62828; }}
    .path-d {{ background: #eceff1; color: #546e7a; }}
    .path-none {{ background: #f5f5f5; color: #999; }}
    tr.row-a {{ background: #f4faf4; }}
    tr.row-b {{ background: #fffcf2; }}
    tr.row-c {{ background: #fdf4f4; }}
    tr.row-d {{ background: #f6f8f9; }}
    tr.row-a:hover, tr.row-b:hover, tr.row-c:hover, tr.row-d:hover {{ background: #eef4ef; }}
    .na {{ color: #bbb; }}
    small {{ color: #888; }}
    footer {{
        margin-top: 40px; padding-top: 16px;
        border-top: 1px solid #eee; font-size: 0.78rem; color: #999;
    }}
</style>
</head>
<body>
<h1>Promoter identification report</h1>
<p style="color:#666; margin-bottom:4px;">Genome: <strong>{genome_name}</strong></p>

<div class="summary">
    <div class="stat"><div class="label">Verified marker promoters</div><div class="value">{n_marker_promoters}</div></div>
    <div class="stat"><div class="label">CDS annotated</div><div class="value">{n_annotated}</div></div>
    <div class="stat"><div class="label">With &minus;10/&minus;35 motifs</div><div class="value">{n_with_motifs}</div></div>
    <div class="stat"><div class="label">All promoter IGRs</div><div class="value">{n_all_promoters}</div></div>
</div>

<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    Domain: <strong>{domain_label}</strong> · Verification: <strong>&minus;10/&minus;35 motif scan</strong><br/>
    This report shows the <strong>{n_marker_promoters}</strong> marker-filtered promoters confirmed by &minus;10/&minus;35 motif scanning.
    A total of {n_all_promoters} promoter-orientation IGRs were identified in the genome.
</p>
<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    <strong>Downloads:</strong>
    <a href="{marker_verified_fasta}">marker_promoters_verified.fasta</a> &middot;
    <a href="{all_verified_fasta}">all_promoters_verified.fasta</a> &middot;
    <a href="{cds_annotations_tsv}">cds_annotations.tsv</a> &middot;
    <a href="{final_table}">profinder_results.tsv</a>
</p>

<h2>Verified marker promoters</h2>
"""

_REPORT_HTML_FOOT = """
<footer>
    Generated by ProFinder.
</footer>
</body>
</html>
"""

# =====================================================================
#  Step 10 — Final output table
# =====================================================================

def step10_final_table(cfg: Config, force: bool = False):
    """Build a comprehensive TSV table covering ALL promoter-orientation
    IGRs with every piece of metadata collected by the pipeline.

    Columns
    -------
    promoter_id           IGR identifier (igr_NNNNNN)
    contig                Source contig / scaffold
    start                 IGR start coordinate (1-based)
    end                   IGR end coordinate
    length                IGR length in bp
    orientation           CO_F or CO_R
    associated_cds        Prokka gene ID of downstream CDS
    gene_name             Gene name from Prokka GFF (gene= attribute)
    locus_tag             Locus tag from Prokka GFF (locus_tag= attribute)
    product               Protein product name from Prokka GFF
    is_marker             Whether this IGR flanks a marker-operon gene
    motif_confirmed       Whether a -10/-35 motif combination was found
    motif_path            Motif path classification (A, B, C, or D)
    motif_strand          Strand on which best motif hit was found
    motif_pos_10          Position of -10 element in scanned sequence
    motif_seq_10          Sequence of -10 element
    motif_score_10        log2-odds score of -10 hit
    motif_pvalue_10       Empirical p-value of -10 hit under its PWM
    motif_source_10       Subgroup ID of the best -10 hit (e.g. M001)
    motif_has_ext10       Whether extended -10 TG dinucleotide is present
    motif_pos_35          Position of -35 element
    motif_seq_35          Sequence of -35 element
    motif_score_35        log2-odds score of -35 hit
    motif_pvalue_35       Empirical p-value of -35 hit under its PWM
    motif_source_35       Subgroup ID of the best -35 hit (e.g. M001)
    motif_spacer_len      Spacer length between -35 and -10
    motif_spacer_seq      Spacer sequence
    sequence_5p_to_3p     Full promoter sequence oriented 5'→3'
    """
    if not force and cfg.final_table.exists():
        print("── Final output table already exists, skipping ──")
        print("  Step 10 complete.\n")
        return

    print("── Building final output table ──")

    # 1. Start from ALL promoter-orientation IGRs
    try:
        igr = pd.read_csv(cfg.igr_summary, sep="\t")
        igr = igr[igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        igr = pd.DataFrame()

    if igr.empty:
        print("  No promoter-orientation IGRs found.")
        pd.DataFrame().to_csv(cfg.final_table, sep="\t", index=False)
        print("  Step 10 complete.\n")
        return

    # 2. Associated CDS
    igr["associated_cds"] = igr.apply(_get_associated_gene, axis=1)

    # 3. Orient sequences 5'→3'
    igr["sequence_5p_to_3p"] = igr.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # 4. CDS annotations (gene name, locus tag, product)
    annotation_tsv = cfg.cds_annotations
    ann_map = {}   # gene_id -> {product, gene, locus_tag}
    if annotation_tsv.exists():
        try:
            ann_df = pd.read_csv(annotation_tsv, sep="\t")
            for _, row in ann_df.iterrows():
                ann_map[row["gene_id"]] = {
                    "product": str(row.get("product", "")),
                    "gene": str(row.get("gene", "")),
                    "locus_tag": str(row.get("locus_tag", "")),
                }
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["gene_name"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("gene", ""))
    igr["locus_tag"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("locus_tag", ""))
    igr["product"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("product", ""))

    # 5. Marker status
    marker_igr_ids = set()
    if cfg.promoter_markers.exists():
        try:
            m_df = pd.read_csv(cfg.promoter_markers, sep="\t")
            marker_igr_ids = set(m_df["igr_id"])
        except (pd.errors.EmptyDataError, KeyError):
            pass
    igr["is_marker"] = igr["igr_id"].isin(marker_igr_ids)

    # 6. Motif scan results (best hit per IGR from step 8)
    motif_cols = ["strand", "pos_10", "seq_10", "score_10", "pvalue_10", "source_10",
                  "has_ext10", "pos_35", "seq_35", "score_35", "pvalue_35", "source_35",
                  "spacer_len", "spacer_seq", "path"]
    motif_map = {}   # igr_id -> dict of motif columns
    if cfg.motif_best_all.exists():
        try:
            motif_df = pd.read_csv(cfg.motif_best_all, sep="\t")
            for _, row in motif_df.iterrows():
                motif_map[row["igr_id"]] = {c: row.get(c, "") for c in motif_cols}
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["motif_confirmed"] = igr["igr_id"].isin(set(motif_map.keys()))
    for c in motif_cols:
        col_name = f"motif_{c}"
        igr[col_name] = igr["igr_id"].map(
            lambda g, _c=c: motif_map.get(g, {}).get(_c, ""))

    # 7. CDS-extended sequences (optional)
    columns = [
        "igr_id", "contig", "start", "end", "length", "orientation",
        "associated_cds", "gene_name", "locus_tag", "product",
        "is_marker", "motif_confirmed", "motif_path", "motif_strand",
        "motif_pos_10", "motif_seq_10", "motif_score_10", "motif_pvalue_10",
        "motif_source_10", "motif_has_ext10",
        "motif_pos_35", "motif_seq_35", "motif_score_35", "motif_pvalue_35",
        "motif_source_35",
        "motif_spacer_len", "motif_spacer_seq",
        "sequence_5p_to_3p",
    ]

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        igr = _add_cds_column(igr, contigs, cfg.cds_bp)
        columns.append("sequence_5p_to_3p_cds")

    # 8. Select and order final columns
    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)

    out.to_csv(cfg.final_table, sep="\t", index=False)
    print(f"  Final table ({len(out)} rows) -> {cfg.final_table}")
    print("  Step 10 complete.\n")


STEPS = [
    (1,  "Run Prokka",                       step01_run_prokka),
    (2,  "Extract intergenic regions",       step02_extract_igrs),
    (3,  "Identify operons",                 step03_identify_operons),
    (4,  "Run hmmsearch",                    step04_run_hmmsearch),
    (5,  "Filter HMM output",                step05_filter_hmm),
    (6,  "Filter operons + add markers",     step06_filter_operons_add_markers),
    (7,  "Match IGRs to marker operons",     step07_match_igrs_to_markers),
    (8,  "Scan promoter motifs (-10/-35)",   step08_scan_motifs),
    (9,  "Annotate CDS (Prokka)",            step09_annotate_cds),
    (10, "Build final output table",         step10_final_table),
    (11, "Generate HTML report",             step11_generate_report),
]


def main():
    parser = argparse.ArgumentParser(
        description="ProFinder — bacterial and archaeal promoter identification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Single-sample and batch inputs are mutually exclusive.  Requirement
    # is enforced manually after parse_args so --list can run without
    # either (an argparse `required=` flag here would need to know about
    # --list before the parser has run).
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("-i", "--input", type=Path,
                             help="Input genome FASTA file (single-sample mode)")
    input_group.add_argument("--batch", type=Path,
                             help="TSV table for batch mode. Required columns: "
                                  "sample_id, fasta. Optional columns: gff, faa, "
                                  "fna. When gff and faa are provided, Prokka "
                                  "(step 1) is skipped for that sample.")
    parser.add_argument("-o", "--output", type=Path, default=Path("output"),
                        help="Output directory (default: output/). In batch mode "
                             "each sample gets a subdirectory.")
    parser.add_argument("--domain", choices=["bacteria", "archaea"],
                        default="bacteria",
                        help="Target domain (affects Prokka --kingdom). "
                             "Default: bacteria")
    parser.add_argument("--start", type=int, default=1,
                        help="First step to run (default: 1)")
    parser.add_argument("--end", type=int, default=max(s[0] for s in STEPS),
                        help="Last step to run (default: last)")
    parser.add_argument("--list", action="store_true",
                        help="List all steps and exit")
    parser.add_argument("--force", action="store_true",
                        help="Re-run steps even if output already exists")

    # Tool paths
    parser.add_argument("--prokka", default="prokka",
                        help="Path to prokka binary (default: prokka)")
    parser.add_argument("--hmmsearch", default="hmmsearch",
                        help="Path to hmmsearch binary (default: hmmsearch)")
    parser.add_argument("--hmm-dir", type=Path, default=None,
                        help="Directory containing individual .hmm profile files "
                             "(default: bundled profiles)")

    # Motif scanning
    parser.add_argument("--meme-file", type=Path, default=None,
                        help="Single MEME file with paired M###_m10 / "
                             "M###_m35 subgroup motifs "
                             "(default: bundled all_unique_subgroups.meme)")
    parser.add_argument("--p10", type=float, default=2.5e-3,
                        help="p-value threshold for -10 motif hits "
                             "(default: 0.0025)")
    parser.add_argument("--p35", type=float, default=2.5e-3,
                        help="Strict p-value threshold for -35 motif hits, "
                             "used for Path A and Path C (default: 0.0025)")
    parser.add_argument("--p35-relaxed", type=float, default=0.05,
                        help="Relaxed p-value threshold for -35 motif hits, "
                             "used only for Path B (extended -10) "
                             "(default: 0.05)")
    parser.add_argument("--max-dist-to-cds-start", type=int, default=200,
                        help="Maximum distance (bp) between the right edge "
                             "of the -10 motif and the downstream CDS start. "
                             "-10 hits further out are dropped before the "
                             "spacer search. Captures the bulk of literature "
                             "σ⁷⁰ leader lengths while filtering distal "
                             "motifs in long IGRs that are unlikely to "
                             "drive the downstream gene (default: 200)")

    # Conda environments
    parser.add_argument("--conda-prokka", default="",
                        help="Conda env for Prokka (blank = use current env)")
    parser.add_argument("--conda-hmm", default="",
                        help="Conda env for hmmsearch (blank = use current env)")

    # Parameters
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for external tools (default: 4)")
    parser.add_argument("--kingdom", default="Bacteria",
                        help="Prokka --kingdom (default: Bacteria)")
    parser.add_argument("--prefix", default="genome",
                        help="Prokka --prefix (default: genome)")
    parser.add_argument("--igr-min", type=int, default=81,
                        help="Minimum IGR length (default: 81)")
    parser.add_argument("--igr-max", type=int, default=1000,
                        help="Maximum IGR length (default: 1000)")
    parser.add_argument("--max-internal-dist", type=int, default=25,
                        help="Max distance between genes in an operon (default: 25)")
    parser.add_argument("--min-flanking-dist", type=int, default=75,
                        help="Min flanking distance for operon boundaries (default: 75)")
    parser.add_argument("--hmm-bitscore", type=float, default=25.0,
                        help="Minimum HMM bitscore (default: 25.0)")
    parser.add_argument("--cds-bp", type=int, default=0,
                        help="Number of CDS-start nucleotides to append to "
                             "promoter sequences in additional FASTA outputs "
                             "and the final table. 0 = disabled (default: 0)")

    # Diagnostic: -10 / -35 p-value sweep
    parser.add_argument("--diagnose-p-sweep", action="store_true",
                        help="Run the p-value sweep diagnostic after step 7 "
                             "instead of running step 8 onward. Sweeps -10 "
                             "and -35 PWM thresholds in three modes "
                             "(matched, p10_only, p35_only) and writes "
                             "diagnostics/p_sweep.tsv plus a figure.")
    parser.add_argument("--sweep-min", type=float, default=1e-5,
                        help="Smallest p-value in the sweep "
                             "(default: 1e-5)")
    parser.add_argument("--sweep-max", type=float, default=1e-1,
                        help="Largest p-value in the sweep "
                             "(default: 1e-1)")
    parser.add_argument("--sweep-n", type=int, default=13,
                        help="Number of log-spaced sweep points "
                             "(default: 13)")
    parser.add_argument("--sweep-modes", default="matched,p10_only,p35_only",
                        help="Comma-separated list of sweep modes "
                             "(default: matched,p10_only,p35_only). "
                             "Valid modes: matched | p10_only | p35_only.")
    parser.add_argument("--sweep-p35-relaxed-ratio", type=float, default=20.0,
                        help="In every sweep point, p35_relaxed is set to "
                             "this ratio x p35 (capped at 1.0). Default "
                             "ratio 20 matches the pipeline default of "
                             "2.5e-3 strict / 5e-2 relaxed.")

    args = parser.parse_args()

    if args.list:
        for num, name, _ in STEPS:
            print(f"  {num:2d}. {name}")
        sys.exit(0)

    # One of --input / --batch is required for any real run. We enforce
    # this here (rather than via argparse `required=`) so --list can run
    # without either.
    if args.input is None and args.batch is None:
        parser.error("one of the arguments -i/--input --batch is required")

    # When running in archaea mode, default Prokka kingdom to Archaea
    # unless the user explicitly provided --kingdom.
    kingdom = args.kingdom
    if args.domain == "archaea" and kingdom == "Bacteria":
        kingdom = "Archaea"

    # ── Batch mode ──────────────────────────────────────────────────
    if args.batch is not None:
        if not args.batch.exists():
            sys.exit(f"Batch file not found: {args.batch}")
        _run_batch(args, kingdom)
        return

    # ── Single-sample mode ──────────────────────────────────────────
    if not args.input.exists():
        sys.exit(f"Input file not found: {args.input}")

    cfg = _build_config(args, kingdom,
                        input_fasta=args.input.resolve(),
                        output_dir=args.output.resolve())
    _run_pipeline(cfg, args)

    print("\nPipeline complete.")


def _build_config(args, kingdom, *, input_fasta, output_dir, prokka_prefix=None):
    """Create a Config from CLI args, overriding paths as needed."""
    return Config(
        input_fasta=input_fasta,
        output_dir=output_dir,
        domain=args.domain,
        prokka_bin=args.prokka,
        hmmsearch_bin=args.hmmsearch,
        hmm_profiles_dir=args.hmm_dir,
        conda_env_prokka=args.conda_prokka,
        conda_env_hmm=args.conda_hmm,
        threads=args.threads,
        prokka_kingdom=kingdom,
        prokka_prefix=prokka_prefix or args.prefix,
        igr_size_min=args.igr_min,
        igr_size_max=args.igr_max,
        max_internal_distance=args.max_internal_dist,
        min_flanking_distance=args.min_flanking_dist,
        hmm_bitscore_min=args.hmm_bitscore,
        meme_file=args.meme_file,
        motif_p10=args.p10,
        motif_p35=args.p35,
        motif_p35_relaxed=args.p35_relaxed,
        max_dist_to_cds_start=args.max_dist_to_cds_start,
        cds_bp=args.cds_bp,
    )


def _run_pipeline(cfg, args):
    """Run pipeline steps on a single Config."""
    cfg.ensure_dirs()

    if cfg.cds_bp > 0 and cfg.cds_bp % 3 != 0:
        print()
        print("!" * 60)
        print("  WARNING: --cds-bp %d is not divisible by 3." % cfg.cds_bp)
        print()
        print("  The CDS extension will be %d nt, which is not a whole" % cfg.cds_bp)
        print("  number of codons. If you concatenate these promoter+CDS")
        print("  sequences upstream of a coding sequence, the downstream")
        print("  reading frame will be shifted. Use a multiple of 3")
        print("  (e.g. 90, 150, 300) to keep the CDS fragment in-frame.")
        print("!" * 60)
        print()

    if not args.force:
        print("Checkpoint mode: steps with existing output will be skipped.")
        print("Use --force to re-run all steps.\n")

    # Diagnostic mode: run the p-value-independent prerequisites
    # (Prokka, IGR extraction, operons, HMM, marker matching, CDS
    # annotation — steps 1-7 plus 9), skip the default-threshold motif
    # scan (step 8), and replace steps 10-11 with the p-value sweep.
    # Step 9 is included because the per-sweep-point profinder_results.tsv
    # files carry gene_name / locus_tag / product columns, which are
    # p-value-independent and worth computing once. Cached outputs from
    # any prior run (full or diagnostic) are reused.
    if getattr(args, "diagnose_p_sweep", False):
        DIAGNOSTIC_PREREQ_STEPS = {1, 2, 3, 4, 5, 6, 7, 9}
        for num, name, func in STEPS:
            if num not in DIAGNOSTIC_PREREQ_STEPS:
                continue
            if args.start <= num <= args.end:
                print(f"\n{'=' * 60}")
                print(f"  STEP {num}: {name}")
                print(f"{'=' * 60}\n")
                func(cfg, force=args.force)
        modes = tuple(m.strip() for m in args.sweep_modes.split(",") if m.strip())
        print(f"\n{'=' * 60}")
        print(f"  DIAGNOSTIC: -10 / -35 p-value sweep")
        print(f"{'=' * 60}\n")
        diagnose_p_sweep(
            cfg,
            sweep_min=args.sweep_min,
            sweep_max=args.sweep_max,
            sweep_n=args.sweep_n,
            modes=modes,
            p35_relaxed_ratio=args.sweep_p35_relaxed_ratio,
            force=args.force,
        )
        return

    for num, name, func in STEPS:
        if args.start <= num <= args.end:
            print(f"\n{'=' * 60}")
            print(f"  STEP {num}: {name}")
            print(f"{'=' * 60}\n")
            func(cfg, force=args.force)


def _link_prokka_files(cfg, row):
    """Symlink user-supplied Prokka files into the expected output layout.

    Returns the step number to start from: 2 if Prokka files were linked
    (skip step 1), or 1 if they were not provided.
    """
    gff_path = row.get("gff", "")
    faa_path = row.get("faa", "")

    if not gff_path or not faa_path:
        return 1  # no Prokka files supplied — run step 1 normally

    gff_src = Path(str(gff_path)).resolve()
    faa_src = Path(str(faa_path)).resolve()

    for label, src in [("gff", gff_src), ("faa", faa_src)]:
        if not src.exists():
            sys.exit(f"Batch table references missing {label} file: {src}")

    cfg.prokka_dir.mkdir(parents=True, exist_ok=True)

    # Derive the prefix from the GFF filename so downstream paths resolve.
    prefix = gff_src.stem  # e.g. "genome" from "genome.gff"
    cfg.prokka_prefix = prefix

    # Symlink GFF and FAA
    for src, suffix in [(gff_src, ".gff"), (faa_src, ".faa")]:
        dest = cfg.prokka_dir / f"{prefix}{suffix}"
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        dest.symlink_to(src)

    # Optional: FNA (nucleotide FASTA from Prokka)
    fna_path = row.get("fna", "")
    if fna_path:
        fna_src = Path(str(fna_path)).resolve()
        if fna_src.exists():
            dest = cfg.prokka_dir / f"{prefix}.fna"
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            dest.symlink_to(fna_src)

    return 2  # skip step 1


def _run_batch(args, kingdom):
    """Run the pipeline once per row in a batch TSV table."""
    batch_df = pd.read_csv(args.batch, sep="\t")

    required_cols = {"sample_id", "fasta"}
    missing = required_cols - set(batch_df.columns)
    if missing:
        sys.exit(f"Batch table is missing required columns: {', '.join(sorted(missing))}")

    n_samples = len(batch_df)
    print(f"Batch mode: {n_samples} sample(s) from {args.batch}\n")

    for idx, row in batch_df.iterrows():
        sample_id = str(row["sample_id"])
        fasta_path = Path(str(row["fasta"])).resolve()

        if not fasta_path.exists():
            sys.exit(f"FASTA not found for sample '{sample_id}': {fasta_path}")

        sample_output = args.output.resolve() / sample_id

        print("\n" + "#" * 60)
        print(f"  SAMPLE {idx + 1}/{n_samples}: {sample_id}")
        print("#" * 60 + "\n")

        cfg = _build_config(
            args, kingdom,
            input_fasta=fasta_path,
            output_dir=sample_output,
        )

        # If Prokka files are pre-supplied, symlink them and skip step 1.
        min_step = _link_prokka_files(cfg, row)
        effective_start = max(args.start, min_step)

        # Temporarily override start so _run_pipeline respects it.
        saved_start = args.start
        args.start = effective_start
        _run_pipeline(cfg, args)
        args.start = saved_start

        print(f"\n  Sample '{sample_id}' complete.\n")

    print(f"\nBatch complete. {n_samples} sample(s) processed.")


if __name__ == "__main__":
    main()
