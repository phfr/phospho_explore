#!/usr/bin/env python3
"""
Find genes with opposing phosphorylation sites (positive vs negative site scores) per condition.

Uses site-level Log2FC and -log10 p from data.tsv; score = Log2FC * (-log10 p) (same spirit as rank lists).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from _conditions import CONDITIONS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--data-tsv", type=Path, default=root / "data.tsv")
    parser.add_argument(
        "--out-tsv",
        type=Path,
        default=root / "output" / "analyze" / "multiphos_switches.tsv",
    )
    parser.add_argument(
        "--min-sites",
        type=int,
        default=2,
        help="Minimum phosphorylation sites per gene to consider.",
    )
    parser.add_argument(
        "--pos-min",
        type=float,
        default=2.0,
        help="A site counts as strongly positive if signed score >= this.",
    )
    parser.add_argument(
        "--neg-max",
        type=float,
        default=-2.0,
        help="A site counts as strongly negative if signed score <= this.",
    )
    parser.add_argument(
        "--min-neglogp",
        type=float,
        default=0.0,
        help="Optional: require -log10 p at each pole to be at least this (0=off).",
    )
    args = parser.parse_args()

    path = args.data_tsv.resolve()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "geneid" not in df.columns or "id" not in df.columns:
        raise SystemExit("data.tsv must contain geneid and id columns.")

    df["geneid"] = df["geneid"].astype(str).str.strip()

    rows: list[dict] = []
    for spec in tqdm(list(CONDITIONS), desc="Conditions", unit="cond"):
        log2 = pd.to_numeric(df[spec.log2fc_col], errors="coerce")
        nlp = pd.to_numeric(df[spec.neglogp_col], errors="coerce")
        if args.min_neglogp > 0:
            mask_cov = nlp >= args.min_neglogp
        else:
            mask_cov = pd.Series(True, index=df.index)
        score = log2 * nlp
        sub = pd.DataFrame(
            {
                "geneid": df["geneid"],
                "site": df["id"],
                "score": score,
                "log2fc": log2,
                "neglogp": nlp,
            }
        )
        sub = sub[mask_cov]
        groups = list(sub.groupby("geneid", sort=False))
        for gid, gdf in tqdm(groups, desc=spec.condition_id, leave=False, unit="gene"):
            gdf = gdf.dropna(subset=["score"])
            if len(gdf) < args.min_sites:
                continue
            best_pos_idx = gdf["score"].idxmax()
            best_neg_idx = gdf["score"].idxmin()
            sp, sn = float(gdf.loc[best_pos_idx, "score"]), float(gdf.loc[best_neg_idx, "score"])
            if sp < args.pos_min or sn > args.neg_max:
                continue
            rows.append(
                {
                    "condition": spec.condition_id,
                    "geneid": gid,
                    "n_sites": len(gdf),
                    "best_pos_site": gdf.loc[best_pos_idx, "site"],
                    "best_pos_score": sp,
                    "best_pos_log2fc": float(gdf.loc[best_pos_idx, "log2fc"]),
                    "best_neg_site": gdf.loc[best_neg_idx, "site"],
                    "best_neg_score": sn,
                    "best_neg_log2fc": float(gdf.loc[best_neg_idx, "log2fc"]),
                    "score_span": sp - sn,
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["condition", "score_span"], ascending=[True, False])
    outp = args.out_tsv.resolve()
    outp.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outp, sep="\t", index=False)
    print(f"Wrote {outp} ({len(out)} gene × condition switches)")


if __name__ == "__main__":
    main()
