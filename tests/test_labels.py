from __future__ import annotations

import math

import numpy as np
import pandas as pd

from idrfeat.labels import binarize, join_labels_to_idrs


def _features() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accession": ["P1", "P2", "P3"],
            "idr_index": [1, 1, 1],
            "idr_start": [1, 1, 1],
            "idr_end": [10, 10, 10],
        }
    )


def _labels() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accession": ["P1", "P2", "P3"],
            "start": [1, 50, 8],
            "end": [10, 60, 10],
            "log2fc": [-2.0, -3.0, -2.0],
        }
    )


def test_join_full_overlap_assigns_score() -> None:
    out = join_labels_to_idrs(_features(), _labels(), min_overlap=0.5)
    by_acc = dict(zip(out["accession"], out["fitness_log2fc"], strict=True))
    assert by_acc["P1"] == -2.0


def test_join_non_overlapping_is_nan() -> None:
    out = join_labels_to_idrs(_features(), _labels(), min_overlap=0.5)
    by_acc = dict(zip(out["accession"], out["fitness_log2fc"], strict=True))
    assert math.isnan(by_acc["P2"])


def test_join_below_min_overlap_is_nan() -> None:
    out = join_labels_to_idrs(_features(), _labels(), min_overlap=0.5)
    by_acc = dict(zip(out["accession"], out["fitness_log2fc"], strict=True))
    assert math.isnan(by_acc["P3"])


def test_binarize_cutoff_and_nan() -> None:
    out = binarize([-2.0, 0.0, np.nan], cutoff=-1.0)
    assert out[0] == 1.0
    assert out[1] == 0.0
    assert math.isnan(out[2])
