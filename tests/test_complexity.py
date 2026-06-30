"""Tests for k-mer content, low-complexity fraction, and longest residue run."""

from __future__ import annotations

import math

from idrfeat.complexity import (
    kmer_features,
    kmer_spectrum,
    longest_single_run,
    lowcomplexity_frac,
)


def test_kmer_spectrum_counts() -> None:
    assert kmer_spectrum("AAAA", 1) == {"A": 4}
    assert kmer_spectrum("AAAA", 2) == {"AA": 3}
    assert kmer_spectrum("ABAB", 2) == {"AB": 2, "BA": 1}


def test_kmer_features_homopolymer_vs_diverse() -> None:
    homo = kmer_features("AAAAAAAA", ks=(1,))
    assert homo["kmer1_distinct"] == 1
    assert math.isclose(homo["kmer1_entropy"], 0.0)
    assert math.isclose(homo["kmer1_diversity"], 1 / 8)

    diverse = kmer_features("ACDE", ks=(1,))
    assert diverse["kmer1_distinct"] == 4
    assert math.isclose(diverse["kmer1_entropy"], 2.0)  # log2(4)
    assert math.isclose(diverse["kmer1_diversity"], 1.0)


def test_kmer_features_cover_requested_k() -> None:
    feats = kmer_features("ACDEFGACDEFG", ks=(1, 2, 3))
    for k in (1, 2, 3):
        assert f"kmer{k}_distinct" in feats
        assert f"kmer{k}_entropy" in feats
        assert f"kmer{k}_diversity" in feats


def test_lowcomplexity_fraction() -> None:
    assert math.isclose(lowcomplexity_frac("A" * 12), 1.0)
    assert math.isclose(lowcomplexity_frac("ACDEFGHIKLMN"), 0.0)


def test_longest_single_run() -> None:
    assert longest_single_run("AABBBBA") == 4
    assert longest_single_run("ABCDE") == 1
    assert longest_single_run("") == 0


def test_empty_sequence_is_safe() -> None:
    feats = kmer_features("", ks=(1, 2, 3))
    assert feats["kmer1_distinct"] == 0
    assert feats["kmer2_entropy"] == 0.0
    assert lowcomplexity_frac("") == 0.0
