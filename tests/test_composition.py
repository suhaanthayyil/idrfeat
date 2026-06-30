"""Tests for amino acid composition and hydropathy features."""

from __future__ import annotations

import math

from idrfeat.composition import (
    composition_features,
    hydropathy_patterning,
    mean_hydropathy,
)


def test_aa_fractions_sum_to_one() -> None:
    feats = composition_features("ACDEFGHIKLMNPQRSTVWY")
    total = sum(v for k, v in feats.items() if k.startswith("aa_frac_"))
    assert math.isclose(total, 1.0, abs_tol=1e-9)


def test_each_residue_fraction_is_correct() -> None:
    feats = composition_features("AAAACDEF")  # 4 A out of 8
    assert math.isclose(feats["aa_frac_A"], 0.5)
    assert math.isclose(feats["aa_frac_C"], 0.125)
    assert math.isclose(feats["aa_frac_Y"], 0.0)


def test_grouped_fractions() -> None:
    # D E K R are charged; F W Y aromatic; P proline; G glycine.
    feats = composition_features("DEKRFWYPG")  # length 9
    assert math.isclose(feats["frac_charged"], 4 / 9)
    assert math.isclose(feats["frac_aromatic"], 3 / 9)
    assert math.isclose(feats["frac_proline"], 1 / 9)
    assert math.isclose(feats["frac_glycine"], 1 / 9)


def test_disorder_promoting_fraction() -> None:
    # Disorder-promoting set is A R G Q S E K P. "ARGQSEKP" is all of them.
    feats = composition_features("ARGQSEKP")
    assert math.isclose(feats["frac_disorder_promoting"], 1.0)
    # I L V F W Y C M are not disorder-promoting.
    feats2 = composition_features("ILVFWYCM")
    assert math.isclose(feats2["frac_disorder_promoting"], 0.0)


def test_mean_hydropathy_normalized_range() -> None:
    # Isoleucine is the most hydrophobic Kyte-Doolittle residue (raw 4.5 -> 1.0).
    assert math.isclose(mean_hydropathy("I"), 1.0, abs_tol=1e-9)
    # Arginine is the most hydrophilic (raw -4.5 -> 0.0).
    assert math.isclose(mean_hydropathy("R"), 0.0, abs_tol=1e-9)
    mixed = mean_hydropathy("IR")
    assert 0.0 < mixed < 1.0


def test_hydropathy_patterning_blocky_vs_mixed() -> None:
    # Following the SCD convention, blocky hydropathy segregation is more negative than the
    # same residues evenly interspersed.
    clustered = hydropathy_patterning("IIIIISSSSS")
    alternating = hydropathy_patterning("ISISISISIS")
    assert clustered < alternating
    # A uniform-hydropathy sequence has no patterning.
    assert math.isclose(hydropathy_patterning("AAAAA"), 0.0, abs_tol=1e-9)
    assert hydropathy_patterning("") == 0.0


def test_empty_sequence_is_safe() -> None:
    feats = composition_features("")
    assert feats["aa_frac_A"] == 0.0
    assert feats["frac_charged"] == 0.0
    assert mean_hydropathy("") == 0.0
