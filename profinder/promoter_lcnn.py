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

class _V1SessionModel:
    """Wraps a TF 1.x Session + SavedModel serving signature as a
    callable that behaves like a Keras model for inference.

    This is the most robust loader for legacy SavedModels that contain
    tensorflow-addons optimizer state.  The v1 loader treats the
    optimizer as opaque graph nodes and never tries to reconstruct
    Python optimizer objects, so it doesn't hit the ``add_slot``
    ``AttributeError`` that breaks both ``TFSMLayer`` and
    ``tf.saved_model.load()`` in TF 2.16+ / Keras 3.
    """

    def __init__(self, sess, input_tensor_name, output_tensor_name):
        self._sess = sess
        self._input_name = input_tensor_name
        self._output_name = output_tensor_name
        # Keep graph + session alive for the lifetime of this object.
        self._graph = sess.graph

    def __call__(self, x):
        if hasattr(x, "numpy"):
            x = x.numpy()
        return self._sess.run(self._output_name, {self._input_name: x})


def _load_via_v1_session(path, tf):
    """Load a SavedModel using the TF 1.x Session-based API.

    Returns a :class:`_V1SessionModel` whose ``__call__`` accepts a
    numpy array / tensor and returns a numpy array, matching the
    interface expected by :func:`_run_inference`.
    """
    graph = tf.Graph()
    sess = tf.compat.v1.Session(graph=graph)
    with graph.as_default():
        meta_graph = tf.compat.v1.saved_model.loader.load(
            sess, [tf.saved_model.SERVING], str(path),
        )
    sig = meta_graph.signature_def["serving_default"]
    input_name = list(sig.inputs.values())[0].name
    output_name = list(sig.outputs.values())[0].name
    return _V1SessionModel(sess, input_name, output_name)


def _load_saved_model(path, tf):
    """Load a TensorFlow SavedModel directory, handling Keras 2,
    Keras 3 (TF >= 2.16), and legacy SavedModels with optimizer state
    from tensorflow-addons.

    Four strategies are tried in order:

    1. ``tf.keras.models.load_model()`` — works on TF < 2.16.
    2. ``TFSMLayer`` — works on Keras 3 when the SavedModel has no
       problematic optimizer state.
    3. ``tf.saved_model.load()`` — works when the SavedModel has no
       optimizer state (extracts the serving signature).
    4. ``tf.compat.v1.saved_model.loader.load()`` — Session-based
       loading that completely avoids Python optimizer reconstruction.
       This is the catch-all for SavedModels trained with
       tensorflow-addons optimisers.
    """
    path_str = str(path)

    # --- Attempt 1: Keras 2-style load (works on TF < 2.16) -----------
    try:
        return tf.keras.models.load_model(path_str, compile=False)
    except (ValueError, TypeError, Exception):
        # Broad catch: Keras 3 raises various errors for legacy
        # SavedModel directories.
        pass

    # --- Attempt 2: TFSMLayer (Keras 3, clean SavedModels) ------------
    try:
        return tf.keras.layers.TFSMLayer(
            path_str, call_endpoint="serving_default"
        )
    except Exception:
        pass

    # --- Attempt 3: raw tf.saved_model.load ---------------------------
    try:
        obj = tf.saved_model.load(path_str)
        if hasattr(obj, "signatures") and "serving_default" in obj.signatures:
            return obj.signatures["serving_default"]
        if callable(obj):
            return obj
    except Exception:
        pass

    # --- Attempt 4: TF 1.x Session-based loader (most robust) --------
    return _load_via_v1_session(path, tf)


def _suppress_tf_logging():
    """Suppress TensorFlow's C++ and Python logging noise.

    Sets environment variables, configures Python loggers, and
    temporarily redirects stderr during ``import tensorflow`` to
    swallow the ``absl::InitializeLog`` and ``cudart_stub`` messages
    that are emitted by the C++ runtime before any Python-level
    logging is active.
    """
    import os
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
    os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"

    import logging
    logging.getLogger("tensorflow").setLevel(logging.ERROR)
    try:
        import absl.logging
        absl.logging.set_verbosity(absl.logging.ERROR)
        logging.getLogger("absl").setLevel(logging.ERROR)
    except ImportError:
        pass


def _import_tf_quietly():
    """Import TensorFlow while suppressing C++ stderr noise.

    The ``absl::InitializeLog()`` and ``cudart_stub.cc`` messages are
    written to stderr by the C++ runtime the instant the TF shared
    library is loaded, before any Python logger exists.  The only way
    to hide them is to redirect ``sys.stderr`` (file descriptor 2)
    during the import.
    """
    _suppress_tf_logging()

    import sys
    import os

    # Redirect fd 2 (C stderr) to /dev/null for the duration of import.
    _orig_fd = os.dup(2)
    _devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(_devnull, 2)
    try:
        import tensorflow as tf
    finally:
        os.dup2(_orig_fd, 2)
        os.close(_orig_fd)
        os.close(_devnull)

    tf.get_logger().setLevel("ERROR")
    return tf


# ── Public API ──────────────────────────────────────────────────────
def load_models(is_promoter_dir: Path, promoters_only_dir: Path):
    """Load the two SavedModel directories and return a (model1, model2) tuple.

    Lazy-imports TensorFlow so that modules that don't call this function
    never pay the import cost.  Works with Keras 2, Keras 3 / TF >= 2.16,
    and legacy SavedModels containing tensorflow-addons optimizer state.
    """
    tf = _import_tf_quietly()

    model1 = _load_saved_model(is_promoter_dir, tf)
    model2 = _load_saved_model(promoters_only_dir, tf)
    return model1, model2


def _run_inference(model, encoded: np.ndarray) -> np.ndarray:
    """Run a forward pass, normalising the output to a plain ndarray.

    Handles four loader flavours:

    * Keras ``Model`` — call directly, returns ndarray or tensor.
    * ``TFSMLayer`` — call directly, returns a dict of tensors.
    * ``ConcreteFunction`` (from ``signatures["serving_default"]``)
      — keyword arguments, returns a dict of tensors.
    * ``_V1SessionModel`` — takes a numpy array, returns a numpy array.
    """
    # _V1SessionModel already accepts numpy and returns numpy.
    if isinstance(model, _V1SessionModel):
        return np.asarray(model(encoded))

    import tensorflow as tf

    tensor_in = tf.constant(encoded)

    # ConcreteFunction expects keyword args.
    if hasattr(model, "structured_input_signature"):
        _, kwarg_spec = model.structured_input_signature
        if kwarg_spec:
            key = next(iter(kwarg_spec))
            result = model(**{key: tensor_in})
        else:
            result = model(tensor_in)
    else:
        result = model(tensor_in)

    if isinstance(result, dict):
        result = next(iter(result.values()))

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
