"""Tests for the sequence-based phase-separation features."""

from __future__ import annotations

import math

from idrfeat.phasesep import charge_aromatic_proxy, llps_seq_score


def test_charge_aromatic_proxy_is_fcr_times_aromatic() -> None:
    # K D E are charged (FCR 0.5); F W Y aromatic (0.5). Product is 0.25.
    assert math.isclose(charge_aromatic_proxy("KDEFWY"), 0.25)
    assert math.isclose(charge_aromatic_proxy("AAAA"), 0.0)


def test_llps_score_higher_for_llps_like_sequence() -> None:
    llps_like = llps_seq_score("GYGRGRGRYGQNQNYGRGRG")
    ordered_like = llps_seq_score("ILVILVILVILVILVILVIL")
    assert llps_like > ordered_like
    assert 0.0 <= llps_like <= 1.0
    assert 0.0 <= ordered_like <= 1.0


def test_llps_score_empty_is_zero() -> None:
    assert llps_seq_score("") == 0.0
