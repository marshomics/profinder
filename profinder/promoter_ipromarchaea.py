"""
iProm-Archaea — CNN predictor for archaeal promoter sequences.

Reimplemented from the iProm-Archaea model by Nazeer et al.
(https://github.com/PromoterTools/iPromArchaea). The architecture and
weights are identical; this module strips out the Flask web layer so
the pipeline can call the predictor directly.

The model is a single-stage binary classifier (promoter vs. non-promoter)
using 6-mer frequency encoding.  Input sequences should be ≥80 nt; the
pipeline is responsible for trimming or windowing before calling
:func:`predict`.

Unlike PromoterLCNN, there is no sigma-factor subtype classification —
archaea do not use bacterial sigma factors.
"""

from itertools import product as itertools_product
from pathlib import Path
from typing import List

import numpy as np


# ── Constants ────────────────────────────────────────────────────────

_K = 6
_MIN_SEQ_LEN = 80
_PREDICTION_THRESHOLD = 0.75

# All possible 6-mers in lexicographic ACGT order (4^6 = 4096 entries).
_POSSIBLE_KMERS = [''.join(p) for p in itertools_product('ACGT', repeat=_K)]


# ── 6-mer frequency encoding ────────────────────────────────────────

def _kmer_encode(sequences: List[str]) -> np.ndarray:
    """Encode sequences as normalised 6-mer frequency vectors.

    Parameters
    ----------
    sequences : list[str]
        DNA sequences (uppercase, ≥80 nt each).

    Returns
    -------
    np.ndarray
        Shape ``(N, 4096, 1)`` — ready for the Conv1D model.
    """
    vectors = []
    for seq in sequences:
        upper = seq.upper()
        kmers = [upper[i:i + _K] for i in range(len(upper) - _K + 1)]
        freq = {km: 0 for km in _POSSIBLE_KMERS}
        for km in kmers:
            if km in freq:
                freq[km] += 1
        total = len(kmers) if kmers else 1
        vec = [count / total for count in freq.values()]
        vectors.append(vec)
    arr = np.array(vectors, dtype=np.float32)
    # Add channel dimension for Conv1D: (N, 4096) -> (N, 4096, 1)
    return np.expand_dims(arr, axis=2)


# ── Model construction ───────────────────────────────────────────────

def _build_model():
    """Reconstruct the iProm-Archaea CNN architecture.

    Must match the original ``get_model()`` in app.py exactly so that
    the saved weights map onto the right layers.
    """
    import os
    import sys
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
    # Redirect C stderr to /dev/null during TF import to suppress
    # absl::InitializeLog() and cudart_stub.cc messages.
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
    from tensorflow.keras.layers import (Conv1D, MaxPooling1D, Dropout,
                                         Flatten, Dense, Input)
    from tensorflow.keras import regularizers, Model

    input_shape = (4 ** _K, 1)  # (4096, 1)
    inputs = Input(shape=input_shape)
    x = Conv1D(filters=16, kernel_size=2, activation='relu',
               kernel_regularizer=regularizers.l2(1e-5))(inputs)
    x = MaxPooling1D(pool_size=2, strides=2)(x)
    x = Conv1D(filters=16, kernel_size=2, activation='relu',
               kernel_regularizer=regularizers.l2(1e-5))(x)
    x = MaxPooling1D(pool_size=4, strides=2)(x)
    x = Flatten()(x)
    x = Dropout(0.22)(x)
    x = Dense(100, activation='relu',
              kernel_regularizer=regularizers.l2(1e-4))(x)
    x = Dropout(0.25)(x)
    out = Dense(1, activation='sigmoid')(x)

    model = Model(inputs=inputs, outputs=[out])
    model.compile(loss='binary_crossentropy', optimizer='adam',
                  metrics=['accuracy'])
    return model


# ── Public API ───────────────────────────────────────────────────────

def load_model(weights_path: Path):
    """Build the CNN and load pre-trained weights.

    Parameters
    ----------
    weights_path : Path
        Path to the ``.h5`` weights file (e.g. ``model_cnn.weights.h5``).

    Returns
    -------
    tensorflow.keras.Model
        Compiled model ready for inference.
    """
    model = _build_model()
    model.load_weights(str(weights_path))
    return model


def _split_windows(sequence: str, window: int = 100, min_len: int = 80) -> List[str]:
    """Split a sequence into non-overlapping windows, matching the
    original iProm-Archaea ``split_sequence()`` logic.

    Parameters
    ----------
    sequence : str
        Full-length DNA sequence.
    window : int
        Window size (default 100, matching the original tool).
    min_len : int
        Minimum acceptable window length (default 80).  Trailing
        fragments shorter than this are discarded.

    Returns
    -------
    list[str]
        Non-overlapping windows of *window* bp (the last window may be
        shorter, down to *min_len*).
    """
    windows = []
    for i in range(0, len(sequence), window):
        w = sequence[i:i + window]
        if len(w) >= min_len:
            windows.append(w)
    return windows


def _predict_windows(model, windows: List[str]) -> List[bool]:
    """Classify a batch of pre-split windows.

    Returns one bool per window: ``True`` = promoter.
    """
    if not windows:
        return []

    encoded = _kmer_encode(windows)
    raw_preds = model.predict(encoded, verbose=0)

    # Apply the same threshold logic as the original app.py:
    # values > 0.75 are rounded (to 1.0), then compared to 1.
    thresholded = np.where(raw_preds > _PREDICTION_THRESHOLD,
                           np.round(raw_preds), raw_preds)
    return [bool(p >= 1.0) for p in thresholded.ravel()]


def predict(model, sequences: List[str]) -> List[bool]:
    """Run iProm-Archaea on *sequences* using the original sliding-window
    approach.

    Each input sequence is split into non-overlapping 100 bp windows
    (matching the original tool's ``split_sequence``).  A sequence is
    classified as a **promoter** if *any* of its windows scores positive.
    This faithfully reproduces the scanning behaviour of the original
    iProm-Archaea web tool.

    Sequences shorter than 80 nt produce no valid windows and are
    classified as non-promoters.

    Parameters
    ----------
    model
        Keras model as returned by :func:`load_model`.
    sequences : list[str]
        DNA sequences of any length.

    Returns
    -------
    list[bool]
        ``True`` for predicted promoter, ``False`` for non-promoter.
    """
    if not sequences:
        return []

    # Build all windows across all sequences, tracking which sequence
    # each window belongs to.
    all_windows = []
    seq_indices = []       # parallel to all_windows: index into sequences
    for i, seq in enumerate(sequences):
        wins = _split_windows(seq.upper())
        for w in wins:
            all_windows.append(w)
            seq_indices.append(i)

    # Batch-predict all windows at once for efficiency.
    if all_windows:
        window_results = _predict_windows(model, all_windows)
    else:
        window_results = []

    # A sequence is a promoter if ANY of its windows scored positive.
    results = [False] * len(sequences)
    for idx, is_prom in zip(seq_indices, window_results):
        if is_prom:
            results[idx] = True

    return results
