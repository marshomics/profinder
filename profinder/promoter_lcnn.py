"""
PromoterLCNN — lightweight CNN predictor for bacterial promoter sequences.

Reimplemented from the PromoterLCNN model by occasumlux
(https://github.com/occasumlux/Promoters). The architecture and weights
are identical; this module removes the external class hierarchy so that
the pipeline only needs TensorFlow at runtime.

The model operates in two stages:
  1. Binary classifier (IsPromoter): promoter vs. non-promoter.
  2. Sigma-factor classifier (PromotersOnly): assigns one of six sigma
     subtypes (σ70, σ24, σ28, σ38, σ32, σ54) to each predicted promoter.

Input sequences must be exactly 81 nt.  The pipeline is responsible for
trimming or windowing before calling :func:`predict`.
"""

from enum import Enum
from pathlib import Path
from typing import List, Tuple

import numpy as np


# ── Sigma-factor labels ─────────────────────────────────────────────
class PromoterType(Enum):
    NON_PROMOTER = 1
    SIGMA_70 = 2
    SIGMA_24 = 3
    SIGMA_28 = 4
    SIGMA_38 = 5
    SIGMA_32 = 6
    SIGMA_54 = 7


# ── One-hot encoding ────────────────────────────────────────────────
_ENCODING = {"A": [1, 0, 0, 0], "T": [0, 1, 0, 0],
             "C": [0, 0, 1, 0], "G": [0, 0, 0, 1]}

_SEQ_LEN = 81


def _encode_sequences(seqs: List[str]) -> np.ndarray:
    """One-hot encode a list of 81-nt DNA strings -> (N, 81, 4) array."""
    encoded = []
    for i, seq in enumerate(seqs):
        upper = seq.upper()
        if len(upper) != _SEQ_LEN:
            raise ValueError(
                f"Sequence {i} has length {len(upper)} nt; "
                f"expected exactly {_SEQ_LEN} nt."
            )
        encoded.append([_ENCODING[nt] for nt in upper])
    return np.array(encoded, dtype=np.float32)


# ── Public API ──────────────────────────────────────────────────────
def load_models(is_promoter_dir: Path, promoters_only_dir: Path):
    """Load the two SavedModel directories and return a (model1, model2) tuple.

    Lazy-imports TensorFlow so that modules that don't call this function
    never pay the import cost.
    """
    import os
    # Suppress TensorFlow informational/warning logs (GPU probing, TensorRT,
    # oneDNN notices) that clutter pipeline output on CPU-only machines.
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    import tensorflow as tf  # noqa: delayed import
    tf.get_logger().setLevel("ERROR")
    # compile=False skips optimizer restoration — we only need forward passes,
    # and the original models were trained with tensorflow-addons' LazyAdam
    # which may not be installed at inference time.
    model1 = tf.keras.models.load_model(str(is_promoter_dir), compile=False)
    model2 = tf.keras.models.load_model(str(promoters_only_dir), compile=False)
    return model1, model2


def predict(
    models: Tuple,
    sequences: List[str],
) -> List[PromoterType]:
    """Run the two-stage PromoterLCNN cascade on *sequences*.

    Parameters
    ----------
    models : tuple
        (is_promoter_model, promoters_only_model) as returned by
        :func:`load_models`.
    sequences : list[str]
        DNA sequences, each exactly 81 nt.

    Returns
    -------
    list[PromoterType]
        One label per input sequence.
    """
    if not sequences:
        return []

    is_promoter_model, sigma_model = models
    n = len(sequences)

    # Encode
    encoded = _encode_sequences(sequences)  # (N, 81, 4)

    # Stage 1 — binary: class 0 = non-promoter, class ≥ 1 = promoter
    stage1_preds = is_promoter_model.predict(encoded, verbose=0)
    stage1_labels = stage1_preds.argmax(axis=1).ravel()

    results = [PromoterType.NON_PROMOTER] * n

    # Indices predicted as promoter (nonzero)
    promoter_idx = np.where(stage1_labels != 0)[0]
    if promoter_idx.size == 0:
        return results

    # Stage 2 — sigma subtype classification on promoter subset
    encoded_sub = encoded[promoter_idx]
    stage2_preds = sigma_model.predict(encoded_sub, verbose=0)
    stage2_labels = stage2_preds.argmax(axis=1).ravel()

    for j, idx in enumerate(promoter_idx):
        # +2 offset: PromoterType enum values start at 2 for SIGMA_70
        results[idx] = PromoterType(int(stage2_labels[j]) + 2)

    return results
