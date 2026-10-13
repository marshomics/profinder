"""Core scoring primitives for the interval-level empirical-FDR redesign.

- parse the bundled MEME file into paired -10/-35 subgroup PWMs
- per-genome mononucleotide log-odds (kept as the *scorer*, per the design)
- numpy-vectorised position scoring of a whole motif stack over a sequence
"""

import math
import numpy as np

_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_MOTIF_WIDTH = 6


def parse_meme(path):
    """Parse a MEME file into {motif_name: (w x 4) freq ndarray}."""
    motifs = {}
    name = None
    rows = []
    in_matrix = False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("MOTIF "):
                if name and rows:
                    motifs[name] = np.asarray(rows, dtype=float)
                name = line.split()[1]
                rows = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix:
                parts = line.split()
                if len(parts) == 4:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        in_matrix = False
                else:
                    in_matrix = False
    if name and rows:
        motifs[name] = np.asarray(rows, dtype=float)
    return motifs


def paired_subgroups(motifs):
    """Return (subgroups, m10_freqs, m35_freqs) aligned by subgroup id."""
    subs = sorted({n[:-4] for n in motifs if n.endswith("_m10")}
                  & {n[:-4] for n in motifs if n.endswith("_m35")})
    m10 = np.stack([motifs[s + "_m10"] for s in subs])   # (S, 6, 4)
    m35 = np.stack([motifs[s + "_m35"] for s in subs])
    return subs, m10, m35


def freq_to_logodds(freqs, bg, pseudocount=1e-6):
    """(S,6,4) freq -> (S,6,4) log2-odds against bg (len-4)."""
    bg = np.asarray(bg, dtype=float)
    f = np.maximum(freqs, pseudocount)
    return (np.log2(f) - np.log2(bg[None, None, :])).astype(np.float32)


def encode(seq):
    """Sequence -> int8 array; non-ACGT -> -1."""
    a = np.full(len(seq), -1, dtype=np.int8)
    for b, i in _BASE_IDX.items():
        a[np.frombuffer(seq.encode(), dtype=np.uint8) == ord(b)] = i
    return a


def score_stack(code, pwm_stack):
    """Score every valid 6-mer start position against every PWM in the stack.

    code: int8 array (len L), -1 for non-ACGT.
    pwm_stack: (S, 6, 4) float32 log-odds.
    Returns (S, P) float32 scores, P = L-5; positions whose 6-mer contains a
    non-ACGT base are set to -inf.
    """
    L = code.shape[0]
    P = L - _MOTIF_WIDTH + 1
    if P <= 0:
        return np.empty((pwm_stack.shape[0], 0), dtype=np.float32)
    scores = np.zeros((pwm_stack.shape[0], P), dtype=np.float32)
    valid = np.ones(P, dtype=bool)
    for j in range(_MOTIF_WIDTH):
        col = code[j:j + P]                      # (P,) base index at offset j
        bad = col < 0
        safe = np.where(bad, 0, col)
        scores += pwm_stack[:, j, safe]          # (S, P) fancy index
        valid &= ~bad
    scores[:, ~valid] = -np.inf
    return scores


def score_scalar(seq, pos, lom):
    """Reference scalar scorer for one 6-mer (validation only)."""
    s = 0.0
    for j in range(_MOTIF_WIDTH):
        b = _BASE_IDX.get(seq[pos + j])
        if b is None:
            return None
        s += lom[j, b]
    return s
