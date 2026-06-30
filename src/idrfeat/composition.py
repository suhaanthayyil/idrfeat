"""Amino acid composition and hydropathy features for an IDR segment.

All functions take a single amino acid string and return plain floats, so they work with no
external dependencies. Non-standard characters (X, U, gaps) are counted in the length but
contribute to no group, which matches how the downstream scales treat them.
"""

from __future__ import annotations

import math
from collections import Counter

from .constants import (
    AA,
    AROMATIC,
    CHARGED,
    DISORDER_PROMOTING,
    GLYCINE,
    KD_NORM,
    POLAR,
    PROLINE,
)


def aa_fractions(seq: str) -> dict[str, float]:
    """Fraction of each of the 20 standard amino acids, keyed by single letter."""
    n = len(seq)
    counts = Counter(seq)
    return {aa: (counts.get(aa, 0) / n if n else 0.0) for aa in AA}


def _group_fraction(seq: str, group: frozenset[str]) -> float:
    n = len(seq)
    if not n:
        return 0.0
    return sum(1 for c in seq if c in group) / n


def composition_features(seq: str) -> dict[str, float]:
    """Per-residue fractions plus grouped fractions used as IDR composition features."""
    feats = {f"aa_frac_{aa}": frac for aa, frac in aa_fractions(seq).items()}
    feats["frac_charged"] = _group_fraction(seq, CHARGED)
    feats["frac_polar"] = _group_fraction(seq, POLAR)
    feats["frac_aromatic"] = _group_fraction(seq, AROMATIC)
    feats["frac_proline"] = _group_fraction(seq, PROLINE)
    feats["frac_glycine"] = _group_fraction(seq, GLYCINE)
    feats["frac_disorder_promoting"] = _group_fraction(seq, DISORDER_PROMOTING)
    return feats


def mean_hydropathy(seq: str) -> float:
    """Mean normalized Kyte-Doolittle hydropathy in [0, 1] over scored residues."""
    vals = [KD_NORM[c] for c in seq if c in KD_NORM]
    return sum(vals) / len(vals) if vals else 0.0


def hydropathy_patterning(seq: str) -> float:
    """Blockiness of hydropathy along the segment.

    An SCD-style decoration computed on mean-centered normalized hydropathy: more negative
    when high- and low-hydropathy residues segregate into blocks, near zero when they are
    evenly interspersed or the sequence has uniform hydropathy.
    """
    vals = [KD_NORM[c] for c in seq if c in KD_NORM]
    n = len(vals)
    if n < 2:
        return 0.0
    mean = sum(vals) / n
    h = [v - mean for v in vals]
    total = sum(h[m] * h[k] * math.sqrt(m - k) for m in range(1, n) for k in range(m))
    return total / n
