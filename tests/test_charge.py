"""Tests for charge content and charge-patterning features."""

from __future__ import annotations

import math

from idrfeat.charge import charge_features, kappa, scd


def test_basic_charge_counts() -> None:
    # 2 positive (K, R), 2 negative (D, E), out of 8 residues.
    feats = charge_features("KRDEAAAA")
    assert math.isclose(feats["frac_positive"], 2 / 8)
    assert math.isclose(feats["frac_negative"], 2 / 8)
    assert math.isclose(feats["fcr"], 4 / 8)
    assert math.isclose(feats["ncpr"], 0.0)
    assert feats["net_charge"] == 0


def test_net_charge_sign() -> None:
    feats = charge_features("KKKR")  # all positive
    assert feats["net_charge"] == 4
    assert math.isclose(feats["ncpr"], 1.0)
    assert math.isclose(feats["fcr"], 1.0)


def test_kappa_segregated_high_mixed_low() -> None:
    segregated = kappa("KKKKKEEEEE")
    mixed = kappa("KEKEKEKEKE")
    assert segregated > 0.8
    assert mixed < 0.1
    assert segregated > mixed


def test_kappa_undefined_returns_minus_one() -> None:
    assert kappa("AAAAAAAAAA") == -1.0  # no charges
    assert kappa("KKKKKAAAAA") == -1.0  # only one charge sign


def test_scd_segregated_more_negative_than_alternating() -> None:
    segregated = scd("KKKKKEEEEE")
    alternating = scd("KEKEKEKEKE")
    assert segregated < alternating
    assert math.isclose(scd("AAAAAA"), 0.0)


def test_empty_sequence_is_safe() -> None:
    feats = charge_features("")
    assert feats["fcr"] == 0.0
    assert feats["net_charge"] == 0
    assert kappa("") == -1.0
    assert scd("") == 0.0
