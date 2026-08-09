"""Loader and IDR join for the a functional screen fitness screen. Column names are configurable;
confirm the real schema, the region-to-IDR overlap rule, and the depletion cutoff with the data owner
before trusting the output."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_labels(
    path: str | Path,
    acc_col: str = "accession",
    start_col: str = "start",
    end_col: str = "end",
    score_col: str = "log2fc",
) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    out = df[[acc_col, start_col, end_col, score_col]].copy()
    out.columns = ["accession", "start", "end", "log2fc"]
    out["accession"] = out["accession"].astype(str)
    return out


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def join_labels_to_idrs(
    features_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    min_overlap: float = 0.5,
    agg: str = "mean",
) -> pd.DataFrame:
    reducer = {"mean": np.mean, "min": np.min, "max": np.max, "median": np.median}[agg]
    by_acc: dict[str, np.ndarray] = {}
    for acc, sub in labels_df.groupby("accession"):
        by_acc[str(acc)] = sub[["start", "end", "log2fc"]].to_numpy(dtype=float)
    scores = []
    for row in features_df.itertuples(index=False):
        seg_len = row.idr_end - row.idr_start + 1
        regions = by_acc.get(str(row.accession))
        if regions is None or seg_len <= 0:
            scores.append(np.nan)
            continue
        vals = [
            val
            for r_start, r_end, val in regions
            if _overlap(row.idr_start, row.idr_end, int(r_start), int(r_end)) / seg_len
            >= min_overlap
        ]
        scores.append(float(reducer(vals)) if vals else np.nan)
    out = features_df.copy()
    out["fitness_log2fc"] = scores
    return out


def binarize(log2fc, cutoff: float = -1.0) -> np.ndarray:
    s = pd.Series(log2fc, dtype=float)
    label = (s <= cutoff).astype(float)
    label[s.isna()] = np.nan
    return label.to_numpy()
