#!/usr/bin/env python3
"""Build GSEA-style .rnk gene lists from phosphorylation TSV (per-condition rank scores)."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

# (output experiment label, Log2FC column, -log10 p-value column)
CONDITION_COLUMNS: list[tuple[str, str, str]] = [
    ("αCD3_5min", "Log2FC_αCD3_5min", "-Log_p-value_αCD3_5min"),
    ("αCD3_αCD226_5min", "Log2FC_αCD3_αCD226_5min", "-Log_p-value_αCD3_αCD226_5min"),
    ("αCD3_αICOS_5min", "Log2FC_αCD3_αICOS_5min", "-Log_p-value_αCD3_αICOS_5min"),
    ("αCD3_αCD2_5min", "Log2FC_αCD3_αCD2_5min)", "-Log_p-value_αCD3_αCD2_5min"),
    ("αCD3_10min", "Log2FC_αCD3_10min", "-Log_p-value_αCD3_10min"),
    ("αCD3_αCD226_10min", "Log2FC_αCD3_αCD226_10min", "-Log_p-value_αCD3+αCD226_10min"),
    ("αCD3_αICOS_10min", "Log2FC_αCD3_αICOS_10min", "-Log_p-value_αCD3_αICOS_10min"),
    ("αCD3_αCD2_10min", "Log2FC_αCD3_αCD2_10min", "-Log_p-value_αCD3+αCD2_10min"),
]


def aggregate_max_abs_by_gene(df: pd.DataFrame, score_col: str) -> pd.Series:
    """For each geneid, keep the site-level score with the largest absolute value."""

    def _max_abs(series: pd.Series) -> float:
        s = series.dropna()
        if s.empty:
            return math.nan
        return float(s.loc[s.abs().idxmax()])

    return df.groupby("geneid", sort=False)[score_col].agg(_max_abs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-rank",
        action="store_true",
        help="Write only gene symbols (omit rank score column); order is unchanged.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    tsv_path = root / "data.tsv"
    out_dir = root / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if "geneid" not in df.columns:
        raise KeyError("Expected column 'geneid' in data.tsv")

    df["geneid"] = df["geneid"].astype(str).str.strip()

    summary: list[tuple[str, int]] = []

    for experiment, log2_col, p_col in CONDITION_COLUMNS:
        missing = [c for c in (log2_col, p_col) if c not in df.columns]
        if missing:
            raise KeyError(f"Missing columns for {experiment}: {missing}")

        rank_col = f"RankScore_{experiment}"
        log2 = pd.to_numeric(df[log2_col], errors="coerce")
        neglogp = pd.to_numeric(df[p_col], errors="coerce")
        df[rank_col] = log2 * neglogp

        per_gene = aggregate_max_abs_by_gene(df[["geneid", rank_col]], rank_col)
        out = per_gene.rename("RankScore").reset_index()
        out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["RankScore"])
        out = out.sort_values("RankScore", ascending=False)

        rnk_path = out_dir / f"{experiment}.rnk"
        write_df = out[["geneid"]] if args.skip_rank else out
        write_df.to_csv(rnk_path, sep="\t", header=False, index=False)
        summary.append((experiment, len(out)))

    print("Unique genes per ranked list (after NaN/inf removal):")
    for experiment, n in summary:
        print(f"  {experiment}: {n}")


if __name__ == "__main__":
    main()
