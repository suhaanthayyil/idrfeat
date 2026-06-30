"""Sequence-complexity features: k-mer content, low-complexity fraction, and runs.

The full k-mer count vectors (400 for k=2, 8000 for k=3) are available from ``kmer_spectrum``
for downstream use, but the per-segment table keeps only their richness and Shannon entropy so
the table stays tidy. Low complexity is flagged with a sliding-window entropy rule rather than
SEG so it needs no external tool and stays deterministic.
"""

from __future__ import annotations

import math
from collections import Counter

LOWCOMPLEXITY_WINDOW = 12
LOWCOMPLEXITY_ENTROPY_BITS = 2.0


def kmer_spectrum(seq: str, k: int) -> dict[str, int]:
    """Counts of every contiguous k-mer in ``seq``."""
    if k <= 0 or len(seq) < k:
        return {}
    return dict(Counter(seq[i : i + k] for i in range(len(seq) - k + 1)))


def _shannon_bits(counts: list[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c:
            p = c / total
            h -= p * math.log2(p)
    return h


def kmer_features(seq: str, ks: tuple[int, ...] = (1, 2, 3)) -> dict[str, float]:
    """Distinct count, Shannon entropy (bits), and normalized diversity per k."""
    feats: dict[str, float] = {}
    for k in ks:
        spectrum = kmer_spectrum(seq, k)
        windows = len(seq) - k + 1 if len(seq) >= k else 0
        distinct = len(spectrum)
        feats[f"kmer{k}_distinct"] = distinct
        feats[f"kmer{k}_entropy"] = _shannon_bits(list(spectrum.values()))
        feats[f"kmer{k}_diversity"] = distinct / windows if windows else 0.0
    return feats


def lowcomplexity_frac(
    seq: str,
    window: int = LOWCOMPLEXITY_WINDOW,
    entropy_bits: float = LOWCOMPLEXITY_ENTROPY_BITS,
) -> float:
    """Fraction of residues covered by a window whose residue entropy is below threshold."""
    n = len(seq)
    if n == 0:
        return 0.0
    w = min(window, n)
    flagged = [False] * n
    for i in range(n - w + 1):
        if _shannon_bits(list(Counter(seq[i : i + w]).values())) < entropy_bits:
            for j in range(i, i + w):
                flagged[j] = True
    return sum(flagged) / n


def longest_single_run(seq: str) -> int:
    """Length of the longest run of one identical residue."""
    best = run = 0
    prev = None
    for c in seq:
        run = run + 1 if c == prev else 1
        prev = c
        best = max(best, run)
    return best
