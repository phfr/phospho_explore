#!/usr/bin/env python3
"""
Alternative rank lists: per condition, split **up** (Log2FC > 0) vs **down** (Log2FC < 0).

Site score = Log2FC × (−log10 p), same as repo-root `generate_rank_lists.py`.
- **Up** list: per gene, **max** score among sites with Log2FC > 0; sort descending.
- **Down** list: per gene, **min** score among sites with Log2FC < 0 (most negative);
  **RankScore = −min** so strongest down-regulation sorts to the top (descending .rnk).

Writes up to 16 files under `test_updown_split/output/` (8 conditions × up/down).
Empty direction lists are skipped (no file) with a console warning.

Does not modify the original `generate_rank_lists.py` or root `output/*.rnk`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Same condition mapping as ../generate_rank_lists.py
CONDITION_COLUMNS: list[tuple[str, str, str]] = [
    ("αCD3_5min", "Log2FC_αCD3_5min", "-Log_p-value_αCD3_5min"),
    ("αCD3_αCD226_5min", "Log2FC_αCD3_αCD226_5min", "-Log_p-value_αCD3_αCD226_5min"),
    ("αCD3_αICOS_5min", "Log2FC_αCD3_αICOS_5min", "-Log_p-value_αCD3_αICOS_5min"),
    ("αCD3_αCD2_5min", "Log2FC_αCD3_αCD2_5min)", "-Log_p-value_αCD3_αCD2_5min"),
    ("αCD3_10min", "Log2FC_αCD3_10min", "-Log_p-value_αCD3_10min"),
    (
        "αCD3_αCD226_10min",
        "Log2FC_αCD3_αCD226_10min",
        "-Log_p-value_αCD3+αCD226_10min",
    ),
    ("αCD3_αICOS_10min", "Log2FC_αCD3_αICOS_10min", "-Log_p-value_αCD3_αICOS_10min"),
    (
        "αCD3_αCD2_10min",
        "Log2FC_αCD3_αCD2_10min",
        "-Log_p-value_αCD3+αCD2_10min",
    ),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-rank",
        action="store_true",
        help="Write gene symbols only (same order as ranked two-column mode).",
    )
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    repo = here.parent
    tsv_path = repo / "data.tsv"
    out_dir = here / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not tsv_path.is_file():
        raise SystemExit(f"Missing {tsv_path} (data.tsv stays at repo root).")

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if "geneid" not in df.columns:
        raise KeyError("Expected column 'geneid' in data.tsv")
    df["geneid"] = df["geneid"].astype(str).str.strip()

    summary: list[tuple[str, int]] = []

    for experiment, log2_col, p_col in CONDITION_COLUMNS:
        missing = [c for c in (log2_col, p_col) if c not in df.columns]
        if missing:
            raise KeyError(f"Missing columns for {experiment}: {missing}")

        log2 = pd.to_numeric(df[log2_col], errors="coerce")
        neglogp = pd.to_numeric(df[p_col], errors="coerce")
        score = log2 * neglogp

        tmp = pd.DataFrame({"geneid": df["geneid"], "log2": log2, "score": score}).dropna(
            subset=["score"]
        )

        # --- up: Log2FC > 0, max site score per gene ---
        up_df = tmp[tmp["log2"] > 0]
        up_path = out_dir / f"{experiment}_up.rnk"
        if up_df.empty:
            print(f"  WARN {up_path.name}: no sites with Log2FC > 0 — skipped")
        else:
            up_series = up_df.groupby("geneid", sort=False)["score"].max()
            up_out = up_series.rename("RankScore").reset_index()
            up_out = up_out.replace([np.inf, -np.inf], np.nan).dropna(subset=["RankScore"])
            up_out = up_out.sort_values("RankScore", ascending=False)
            write_df = up_out[["geneid"]] if args.skip_rank else up_out
            write_df.to_csv(up_path, sep="\t", header=False, index=False)
            summary.append((up_path.name, len(up_out)))

        # --- down: Log2FC < 0, min (most negative) site score; rank by -min descending ---
        down_df = tmp[tmp["log2"] < 0]
        down_path = out_dir / f"{experiment}_down.rnk"
        if down_df.empty:
            print(f"  WARN {down_path.name}: no sites with Log2FC < 0 — skipped")
        else:
            mn = down_df.groupby("geneid", sort=False)["score"].min()
            down_out = (-mn).rename("RankScore").reset_index()
            down_out = down_out.replace([np.inf, -np.inf], np.nan).dropna(subset=["RankScore"])
            down_out = down_out.sort_values("RankScore", ascending=False)
            write_df_d = down_out[["geneid"]] if args.skip_rank else down_out
            write_df_d.to_csv(down_path, sep="\t", header=False, index=False)
            summary.append((down_path.name, len(down_out)))

    print(f"Wrote {len(summary)} lists under {out_dir}:")
    for name, n in summary:
        print(f"  {name}: {n} genes")


if __name__ == "__main__":
    main()
