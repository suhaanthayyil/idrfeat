"""Per-residue disorder annotation and IDR segment calling.

metapredict is the primary predictor and AIUPred the secondary one; AIUPred also supplies the
fold-on-binding propensity. Both are optional. When metapredict is not installed, segments are
called with ``heuristic_disorder_scores``, a transparent hydropathy-and-charge stand-in that
exists only so the pipeline runs offline. It is not a substitute for metapredict, and the
backend that defined each segment is recorded in the output. AIUPred models, when used, load
once into a process-level singleton.
"""

from __future__ import annotations

import numpy as np

from .constants import AROMATIC, DISORDER_PROMOTING, KD_NORM

_HEURISTIC_WINDOW = 11
_predictor = None


def segments_from_scores(
    scores: np.ndarray, threshold: float, min_len: int
) -> list[tuple[int, int]]:
    """Contiguous runs of scores >= threshold at least min_len long, 1-based inclusive."""
    above = np.asarray(scores) >= threshold
    segments = []
    start = None
    for i, hit in enumerate(above):
        if hit and start is None:
            start = i
        elif not hit and start is not None:
            if i - start >= min_len:
                segments.append((start + 1, i))
            start = None
    if start is not None and len(above) - start >= min_len:
        segments.append((start + 1, len(above)))
    return segments


def disordered(scores: np.ndarray, threshold: float) -> np.ndarray:
    """Boolean mask of residues with score >= threshold."""
    return np.asarray(scores) >= threshold


def metapredict_scores(seq: str) -> np.ndarray:
    """Per-residue metapredict disorder scores in [0, 1] for one sequence."""
    import metapredict as meta

    return np.asarray(meta.predict_disorder(seq, return_numpy=True), dtype=float)


def _get_predictor():
    global _predictor
    if _predictor is None:
        from aiupred import AIUPred

        _predictor = AIUPred(force_cpu=True)
    return _predictor


def aiupred_disorder_scores(seq: str) -> np.ndarray:
    """AIUPred per-residue disorder propensity in [0, 1]."""
    return np.asarray(_get_predictor().predict_disorder(seq), dtype=float)


def aiupred_binding_scores(seq: str) -> np.ndarray:
    """AIUPred per-residue binding (fold-on-binding) propensity in [0, 1]."""
    return np.asarray(_get_predictor().predict_binding(seq), dtype=float)


def _raw_propensity(residue: str) -> float:
    base = 1.0 - KD_NORM.get(residue, 0.5)
    if residue in DISORDER_PROMOTING:
        base += 0.2
    if residue in AROMATIC:
        base -= 0.2
    return min(1.0, max(0.0, base))


def heuristic_disorder_scores(seq: str) -> np.ndarray:
    """Deterministic fallback disorder score in [0, 1], window-smoothed.

    Hydrophilic, charged, and disorder-promoting residues score high; hydrophobic and aromatic
    residues score low. A stand-in for metapredict when it is not installed, not a replacement.
    """
    if not seq:
        return np.zeros(0, dtype=float)
    raw = np.array([_raw_propensity(c) for c in seq], dtype=float)
    half = _HEURISTIC_WINDOW // 2
    padded = np.pad(raw, half, mode="edge")
    kernel = np.ones(2 * half + 1) / (2 * half + 1)
    return np.convolve(padded, kernel, mode="valid")


def metapredict_available() -> bool:
    import importlib.util

    return importlib.util.find_spec("metapredict") is not None


def aiupred_available() -> bool:
    import importlib.util

    return importlib.util.find_spec("aiupred") is not None


def default_backend() -> str:
    """Primary disorder backend: metapredict when installed, else the heuristic."""
    return "metapredict" if metapredict_available() else "heuristic"


def primary_disorder(seq: str, backend: str | None = None) -> tuple[np.ndarray, str]:
    """Per-residue disorder scores from the chosen backend, with the backend name."""
    backend = backend or default_backend()
    if backend == "metapredict":
        return metapredict_scores(seq), backend
    if backend == "aiupred":
        return aiupred_disorder_scores(seq), backend
    if backend == "heuristic":
        return heuristic_disorder_scores(seq), backend
    raise ValueError(f"unknown disorder backend: {backend}")
