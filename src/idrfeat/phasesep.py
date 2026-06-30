"""Sequence-based phase-separation features and curated-database membership.

These answer Carter's question on quantifying phase-separation propensity from sequence. The
charge-times-aromatic proxy is the literal product she guessed. ``llps_seq_score`` is a
transparent composite of the compositional signals most often linked to liquid-liquid phase
separation: aromatics, arginine, glycine, glutamine and asparagine, charge, and low complexity.
It is a heuristic, not a trained predictor; its value on the fitness labels is for phase 1 to
test, not to assume. Curated membership (PhaSePro, DrLLPS) is handled in ``annotations``.
"""

from __future__ import annotations

from .charge import charge_features
from .complexity import lowcomplexity_frac
from .composition import _group_fraction
from .constants import AROMATIC

_LLPS_RESIDUES_R = frozenset("R")
_LLPS_RESIDUES_G = frozenset("G")
_LLPS_RESIDUES_QN = frozenset("QN")


def charge_aromatic_proxy(seq: str) -> float:
    """Fraction of charged residues times fraction of aromatic residues."""
    return charge_features(seq)["fcr"] * _group_fraction(seq, AROMATIC)


def llps_seq_score(seq: str) -> float:
    """Mean of six normalized LLPS-promoting signals, in [0, 1].

    Components: aromatic fraction, arginine fraction, glycine fraction, glutamine plus
    asparagine fraction, fraction of charged residues, and low-complexity fraction.
    """
    if not seq:
        return 0.0
    components = [
        _group_fraction(seq, AROMATIC),
        _group_fraction(seq, _LLPS_RESIDUES_R),
        _group_fraction(seq, _LLPS_RESIDUES_G),
        _group_fraction(seq, _LLPS_RESIDUES_QN),
        charge_features(seq)["fcr"],
        lowcomplexity_frac(seq),
    ]
    return sum(components) / len(components)
