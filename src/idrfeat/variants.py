from __future__ import annotations

from .charge import charge_features
from .complexity import kmer_features, longest_single_run, lowcomplexity_frac
from .composition import composition_features, hydropathy_patterning, mean_hydropathy
from .phasesep import charge_aromatic_proxy, llps_seq_score


def segment_sequence_features(
    seg: str,
    ks: tuple[int, ...] = (1, 2, 3),
    lc_window: int = 12,
    lc_entropy: float = 2.0,
) -> dict[str, float]:
    f: dict[str, float] = {}
    f.update(composition_features(seg))
    f.update(charge_features(seg))
    f["hydropathy_mean"] = mean_hydropathy(seg)
    f["hydropathy_patterning"] = hydropathy_patterning(seg)
    f.update(kmer_features(seg, ks=ks))
    f["lowcomplexity_frac"] = lowcomplexity_frac(seg, window=lc_window, entropy_bits=lc_entropy)
    f["longest_single_run"] = float(longest_single_run(seg))
    f["charge_aromatic_proxy"] = charge_aromatic_proxy(seg)
    f["llps_seq_score"] = llps_seq_score(seg)
    return f


def apply_mutation(seq: str, pos: int, alt: str) -> str:
    if pos < 1 or pos > len(seq):
        raise ValueError("pos out of range")
    if len(alt) != 1:
        raise ValueError("alt must be a single residue")
    return seq[: pos - 1] + alt + seq[pos:]


def variant_delta(
    seq: str,
    pos: int,
    alt: str,
    start: int,
    end: int,
    ks: tuple[int, ...] = (1, 2, 3),
    lc_window: int = 12,
    lc_entropy: float = 2.0,
) -> dict[str, float]:
    if not (start <= pos <= end):
        raise ValueError("pos must fall inside the segment")
    kw = {"ks": ks, "lc_window": lc_window, "lc_entropy": lc_entropy}
    wt_seg = seq[start - 1 : end]
    mut_seg = apply_mutation(seq, pos, alt)[start - 1 : end]
    wt = segment_sequence_features(wt_seg, **kw)
    mu = segment_sequence_features(mut_seg, **kw)
    return {f"d_{k}": mu[k] - wt[k] for k in wt}
