"""Charge content and charge-patterning features for an IDR segment.

kappa follows Das and Pappu (2013): the deviation of the local charge asymmetry from the
whole-sequence value over sliding windows of size 5 and 6, normalized by the same quantity
for the maximally segregated sequence of identical composition. It is -1 when undefined,
matching localCIDER, which happens when the sequence carries fewer than two charge signs.
scd is the sequence charge decoration of Sawle and Ghosh (2015); it grows more negative as
like charges cluster and opposite charges separate.
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


def _sigma(npos: int, nneg: int, n: int) -> float:
    fcr = (npos + nneg) / n
    if fcr == 0.0:
        return 0.0
    ncpr = (npos - nneg) / n
    return ncpr * ncpr / fcr


def _delta(charges: list[int]) -> float:
    n = len(charges)
    sig_seq = _sigma(charges.count(1), charges.count(-1), n)
    accum = []
    for g in _WINDOWS:
        if g > n:
            continue
        total = 0.0
        windows = 0
        for i in range(n - g + 1):
            window = charges[i : i + g]
            total += (_sigma(window.count(1), window.count(-1), g) - sig_seq) ** 2
            windows += 1
        accum.append(total / windows)
    return sum(accum) / len(accum) if accum else 0.0


def kappa(seq: str) -> float:
    charges = charge_vector(seq)
    npos, nneg = charges.count(1), charges.count(-1)
    nzero = len(charges) - npos - nneg
    if npos == 0 or nneg == 0 or len(charges) < min(_WINDOWS):
        return -1.0
    segregated = [1] * npos + [0] * nzero + [-1] * nneg
    delta_max = _delta(segregated)
    if delta_max == 0.0:
        return -1.0
    return _delta(charges) / delta_max


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
