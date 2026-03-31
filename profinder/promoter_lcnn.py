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

_ZERO_VEC = [0, 0, 0, 0]  # ambiguous bases (N, R, Y, etc.)

_SEQ_LEN = 81


def _encode_sequences(seqs: List[str]) -> np.ndarray:
    """One-hot encode a list of 81-nt DNA strings -> (N, 81, 4) array.

    Ambiguous bases (N, R, Y, W, etc.) are encoded as the zero vector
    [0, 0, 0, 0] rather than raising an error.
    """
    encoded = []
    for i, seq in enumerate(seqs):
        upper = seq.upper()
        if len(upper) != _SEQ_LEN:
            raise ValueError(
                f"Sequence {i} has length {len(upper)} nt; "
                f"expected exactly {_SEQ_LEN} nt."
            )
        encoded.append([_ENCODING.get(nt, _ZERO_VEC) for nt in upper])
    return np.array(encoded, dtype=np.float32)


# ── Model loading ──────────────────────────────────────────────────

def _load_saved_model(path, tf):
    """Load a TensorFlow SavedModel directory, handling both Keras 2 and
    Keras 3 (TensorFlow ≥ 2.16).

    Keras 3 dropped support for ``tf.keras.models.load_model()`` on
    legacy SavedModel directories.  In that case we fall back to
    ``keras.layers.TFSMLayer``, which wraps the SavedModel as an
    inference-only callable layer.
    """
    try:
        return tf.keras.models.load_model(str(path), compile=False)
    except (ValueError, TypeError):
        # Keras 3+: use TFSMLayer for legacy SavedModel directories.
        return tf.keras.layers.TFSMLayer(
            str(path), call_endpoint="serving_default"
        )


# ── Public API ──────────────────────────────────────────────────────
def load_models(is_promoter_dir: Path, promoters_only_dir: Path):
    """Load the two SavedModel directories and return a (model1, model2) tuple.

    Lazy-imports TensorFlow so that modules that don't call this function
    never pay the import cost.  Works with both Keras 2
    (``tf.keras.models.load_model``) and Keras 3 / TF ≥ 2.16
    (``TFSMLayer`` fallback for legacy SavedModel directories).
    """
    import os
    # Suppress TensorFlow C++ runtime logs (CUDA probing, TensorRT,
    # oneDNN, absl warnings) BEFORE importing TF.  Level 3 = FATAL only.
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
    os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
    # Suppress absl logging that leaks through before TF's own logger
    # is configured (the "All log messages before absl::InitializeLog()"
    # warnings and cudart_stub notices).
    import logging
    logging.getLogger("tensorflow").setLevel(logging.ERROR)
    try:
        import absl.logging
        absl.logging.set_verbosity(absl.logging.ERROR)
        logging.getLogger("absl").setLevel(logging.ERROR)
    except ImportError:
        pass
    import tensorflow as tf  # noqa: delayed import
    tf.get_logger().setLevel("ERROR")

    model1 = _load_saved_model(is_promoter_dir, tf)
    model2 = _load_saved_model(promoters_only_dir, tf)
    return model1, model2


def _run_inference(model, encoded: np.ndarray) -> np.ndarray:
    """Run a forward pass, normalising the output to a plain ndarray.

    When loaded via ``tf.keras.models.load_model`` (Keras 2), calling
    ``model.predict(x)`` returns an ndarray directly.

    When loaded via ``TFSMLayer`` (Keras 3), calling the layer returns
    a dict mapping output-tensor names to tensors.  We extract the
    single output value and convert it to numpy.
    """
    # TFSMLayer is a layer, not a Model — it has no .predict() in the
    # Keras 2 sense.  Use __call__ which works for both.
    import tensorflow as tf

    result = model(tf.constant(encoded))

    if isinstance(result, dict):
        # TFSMLayer returns {"output_name": tensor}.  Grab the first
        # (and typically only) value.
        result = next(iter(result.values()))

    # Convert to numpy if it's still a tensor.
    if hasattr(result, "numpy"):
        result = result.numpy()

    return np.asarray(result)


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
    stage1_preds = _run_inference(is_promoter_model, encoded)
    stage1_labels = stage1_preds.argmax(axis=1).ravel()

    results = [PromoterType.NON_PROMOTER] * n

    # Indices predicted as promoter (nonzero)
    promoter_idx = np.where(stage1_labels != 0)[0]
    if promoter_idx.size == 0:
        return results

    # Stage 2 — sigma subtype classification on promoter subset
    encoded_sub = encoded[promoter_idx]
    stage2_preds = _run_inference(sigma_model, encoded_sub)
    stage2_labels = stage2_preds.argmax(axis=1).ravel()

    for j, idx in enumerate(promoter_idx):
        # +2 offset: PromoterType enum values start at 2 for SIGMA_70
        results[idx] = PromoterType(int(stage2_labels[j]) + 2)

    return results
