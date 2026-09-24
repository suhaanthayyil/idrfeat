"""Charge content and charge-patterning features for an IDR segment.

kappa follows Das and Pappu (2013): for each window size (5 and 6) the deviation of the local
charge asymmetry from the whole-sequence value is normalized by the same quantity for the
maximally segregated sequence of identical composition, and the two per-window ratios are
averaged. This matches the localCIDER reference implementation (Holehouse et al., 2017); the
equivalence is locked by a test. kappa is -1 when undefined, which happens when the sequence
carries fewer than two charge signs. scd is the sequence charge decoration of Sawle and Ghosh
(2015); it grows more negative as like charges cluster and opposite charges separate, and it
also matches localCIDER exactly.
"""

from __future__ import annotations

import math

from .constants import NEGATIVE, POSITIVE

_WINDOWS = (5, 6)


def charge_vector(seq: str) -> list[int]:
    """Per-residue charge as +1 (K, R), -1 (D, E), or 0."""
    return [1 if c in POSITIVE else -1 if c in NEGATIVE else 0 for c in seq]


def charge_features(seq: str) -> dict[str, float]:
    n = len(seq)
    npos = sum(1 for c in seq if c in POSITIVE)
    nneg = sum(1 for c in seq if c in NEGATIVE)
    return {
        "fcr": (npos + nneg) / n if n else 0.0,
        "ncpr": (npos - nneg) / n if n else 0.0,
        "net_charge": npos - nneg,
        "frac_positive": npos / n if n else 0.0,
        "frac_negative": nneg / n if n else 0.0,
        "kappa": kappa(seq),
        "scd": scd(seq),
    }


def kappa(seq: str) -> float:
    """Charge patterning from the localCIDER reference (Holehouse et al., 2017).

    Returns -1 when kappa is undefined (fewer than two charge signs, or too short), matching
    localCIDER. Returns NaN only when localCIDER is not installed, so the column records that
    the source was absent rather than reporting a wrong value.
    """
    charges = charge_vector(seq)
    if charges.count(1) == 0 or charges.count(-1) == 0 or len(charges) < min(_WINDOWS):
        return -1.0
    try:
        from localcider.sequenceParameters import SequenceParameters

        from .disorder import standardize_sequence
    except ImportError:
        return float("nan")
    return float(SequenceParameters(standardize_sequence(seq)).get_kappa())


def scd(seq: str) -> float:
    q = charge_vector(seq)
    n = len(q)
    if n < 2:
        return 0.0
    total = 0.0
    for m in range(1, n):
        qm = q[m]
        if qm == 0:
            continue
        for k in range(m):
            if q[k]:
                total += qm * q[k] * math.sqrt(m - k)
    return total / n
