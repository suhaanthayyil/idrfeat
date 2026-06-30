"""End-to-end test of the `idrfeat features` command on the bundled example."""

from __future__ import annotations

import pandas as pd

from idrfeat.cli import main
from idrfeat.features import TABLE_COLUMNS
from idrfeat.io import example_fasta_path


def test_features_command_writes_table(tmp_path) -> None:
    out = tmp_path / "features.parquet"
    rc = main(
        [
            "features",
            "--fasta",
            str(example_fasta_path()),
            "--out",
            str(out),
            "--backend",
            "heuristic",
        ]
    )
    assert rc == 0
    assert out.exists()
    assert out.with_suffix(".csv").exists()
    df = pd.read_parquet(out)
    assert list(df.columns) == TABLE_COLUMNS
    assert len(df) >= 1
    assert (df["disorder_backend"] == "heuristic").all()
