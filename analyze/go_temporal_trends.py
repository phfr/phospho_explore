#!/usr/bin/env python3
"""
Classify GO terms into Early / Late / Sustained / Absent per co-stim arm (5 min vs 10 min).

Reads go_bp_prerank_long.tsv from run_go_prerank_batch.py.
Default significance is nominal p<=0.05 (use --threshold-mode fdr for FDR-based calls).

Writes the full table plus a non-Absent-only TSV for quick review. Optional scatter plots
(5 min vs 10 min -log10 p) are the usual way to *see* Early vs Late vs Sustained clusters.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _conditions import TEMPORAL_ARMS

PATTERN_ORDER = ("Early", "Late", "Sustained", "Absent")
PATTERN_COLORS = {
    "Early": "#0072B2",
    "Late": "#D55E00",
    "Sustained": "#009E73",
    "Absent": "#BDBDBD",
}


def neglog10_values(s: pd.Series) -> np.ndarray:
    """Map p-values / FDR to -log10; non-finite or missing → 0 (no evidence)."""
    a = pd.to_numeric(s, errors="coerce").to_numpy(dtype=float)
    out = np.zeros_like(a, dtype=float)
    m = np.isfinite(a) & (a > 0)
    out[m] = -np.log10(np.clip(a[m], 1e-300, 1.0))
    return out


def _scatter_xy_max(sub: pd.DataFrame, c5: str, c10: str, thr_line: float) -> float:
    x = neglog10_values(sub[c5])
    y = neglog10_values(sub[c10])
    return max(float(np.nanmax(x)), float(np.nanmax(y)), thr_line * 1.05, 0.5)


def draw_temporal_scatter_ax(
    ax: plt.Axes,
    sub: pd.DataFrame,
    val_col: str,
    c5: str,
    c10: str,
    thr_line: float,
    *,
    xy_max: float | None = None,
    title: str | None = None,
    show_legend: bool = True,
) -> None:
    x = neglog10_values(sub[c5])
    y = neglog10_values(sub[c10])
    for pat in PATTERN_ORDER:
        m = (sub["pattern"].to_numpy() == pat)
        if not m.any():
            continue
        ax.scatter(
            x[m],
            y[m],
            c=PATTERN_COLORS[pat],
            s=14 if pat == "Absent" else 22,
            alpha=0.35 if pat == "Absent" else 0.75,
            edgecolors="none",
            label=pat,
            rasterized=True,
        )
    ax.axvline(thr_line, color="0.35", ls="--", lw=0.9)
    ax.axhline(thr_line, color="0.35", ls="--", lw=0.9)
    ax.set_xlabel(f"-log10({val_col}) 5 min", fontsize=8)
    ax.set_ylabel(f"-log10({val_col}) 10 min", fontsize=8)
    if title:
        ax.set_title(title, fontsize=9)
    lo = 0.0
    hi = xy_max if xy_max is not None else _scatter_xy_max(sub, c5, c10, thr_line)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    if show_legend:
        ax.legend(loc="upper left", fontsize=7, framealpha=0.9)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--long-tsv",
        type=Path,
        default=root / "output" / "analyze" / "go_bp_prerank_long.tsv",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        default=root / "output" / "analyze" / "go_temporal_classification.tsv",
    )
    parser.add_argument(
        "--threshold-mode",
        choices=("fdr", "pval"),
        default="pval",
        help="Significance column: fdr (FDR q-value) or pval (nominal NOM p-val from prerank).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="Cutoff for significance (default: nominal p<=0.05). For --threshold-mode fdr, try 0.25–1.0; "
        "strict FDR (e.g. 0.25) often leaves almost all terms 'Absent' with prerank.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Write per-arm bar charts of pattern counts under output/analyze/figures/",
    )
    parser.add_argument(
        "--plot-scatter",
        action="store_true",
        help="Per-arm scatter: -log10(sig) at 5 min (x) vs 10 min (y); shows Early/Late/Sustained quadrants.",
    )
    parser.add_argument(
        "--plot-scatter-facet",
        action="store_true",
        help="Single 2×2 figure with one temporal scatter per co-stim arm (shared axis limits).",
    )
    parser.add_argument(
        "--no-nonabsent-tsv",
        action="store_true",
        help="Do not write <out-stem>_nonabsent.tsv (Early/Late/Sustained rows only).",
    )
    args = parser.parse_args()

    long_path = args.long_tsv.resolve()
    if not long_path.is_file():
        raise SystemExit(f"Missing {long_path}")

    df = pd.read_csv(long_path, sep="\t")
    label_col = "term_label" if "term_label" in df.columns else "Term"
    val_col = "fdr" if args.threshold_mode == "fdr" and "fdr" in df.columns else "pval"
    if val_col not in df.columns:
        raise SystemExit(f"Need column {val_col} in long TSV")

    go_map = (
        df[[label_col, "go_id"]].dropna(subset=["go_id"]).drop_duplicates(subset=label_col).set_index(label_col)["go_id"]
        if "go_id" in df.columns
        else None
    )

    def is_sig(v: float) -> bool:
        return bool(pd.notna(v) and v <= args.threshold)

    rows: list[dict] = []
    for arm in TEMPORAL_ARMS:
        d5 = df[df["condition"] == arm.t5_condition_id].set_index(label_col)[val_col]
        d10 = df[df["condition"] == arm.t10_condition_id].set_index(label_col)[val_col]
        terms = d5.index.union(d10.index)
        for term in terms:
            v5 = d5.get(term, float("nan"))
            v10 = d10.get(term, float("nan"))
            s5 = is_sig(float(v5))
            s10 = is_sig(float(v10))
            if s5 and s10:
                pat = "Sustained"
            elif s5 and not s10:
                pat = "Early"
            elif (not s5) and s10:
                pat = "Late"
            else:
                pat = "Absent"
            gid = go_map.get(term) if go_map is not None else None
            rows.append(
                {
                    "arm": arm.arm_id,
                    "go_id": gid,
                    "term_label": term,
                    "pattern": pat,
                    val_col + "_5min": v5,
                    val_col + "_10min": v10,
                }
            )

    out = pd.DataFrame(rows)
    out_path = args.out_tsv.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path} ({len(out)} rows)")

    if not args.no_nonabsent_tsv:
        na_path = out_path.with_name(out_path.stem + "_nonabsent.tsv")
        c5, c10 = val_col + "_5min", val_col + "_10min"
        sub_na = out[out["pattern"] != "Absent"].copy()
        sub_na["_min_sig"] = sub_na[[c5, c10]].min(axis=1)
        sub_na = sub_na.sort_values(["arm", "pattern", "_min_sig"], ascending=[True, True, True]).drop(
            columns=["_min_sig"]
        )
        sub_na.to_csv(na_path, sep="\t", index=False)
        print(f"Wrote {na_path} ({len(sub_na)} non-Absent rows) — browse Early/Late/Sustained here")
        try:
            rel_na = na_path.relative_to(root)
        except ValueError:
            rel_na = na_path
        print(
            f"  All 8 phospho conditions: python analyze/cross_condition_go_heatmap.py "
            f'--terms-tsv "{rel_na}" --out-png output/analyze/figures/go_nonabsent_crosscondition_heatmap.png '
            f"--dotplot"
        )

    fig_dir = root / "output" / "analyze" / "figures"
    if args.plot or args.plot_scatter or args.plot_scatter_facet:
        fig_dir.mkdir(parents=True, exist_ok=True)

    if args.plot:
        for arm in TEMPORAL_ARMS:
            sub = out[out["arm"] == arm.arm_id]
            counts = sub["pattern"].value_counts()
            fig, ax = plt.subplots(figsize=(5, 3))
            ser = counts.reindex(PATTERN_ORDER).fillna(0).astype(int)
            ser.plot(kind="bar", ax=ax, color="steelblue")
            ax.set_title(f"GO temporal patterns: {arm.arm_id}")
            ax.set_ylabel("GO terms")
            fig.tight_layout()
            fp = fig_dir / f"go_temporal_{arm.arm_id}.png"
            fig.savefig(fp, bbox_inches="tight")
            plt.close(fig)
            print(f"  plot {fp}")

    if args.plot_scatter or args.plot_scatter_facet:
        c5 = val_col + "_5min"
        c10 = val_col + "_10min"
        thr_line = -np.log10(max(args.threshold, 1e-300))

        if args.plot_scatter_facet:
            xy_max = max(
                _scatter_xy_max(out[out["arm"] == arm.arm_id], c5, c10, thr_line) for arm in TEMPORAL_ARMS
            )
            fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
            for ax, arm in zip(np.ravel(axes), TEMPORAL_ARMS):
                sub = out[out["arm"] == arm.arm_id]
                draw_temporal_scatter_ax(
                    ax,
                    sub,
                    val_col,
                    c5,
                    c10,
                    thr_line,
                    xy_max=xy_max,
                    title=arm.arm_id,
                    show_legend=False,
                )
            leg_patches = [mpatches.Patch(color=PATTERN_COLORS[p], label=p) for p in PATTERN_ORDER]
            fig.legend(
                handles=leg_patches,
                loc="upper right",
                bbox_to_anchor=(0.99, 0.99),
                fontsize=8,
                framealpha=0.95,
            )
            fig.suptitle(f"GO temporal (5 vs 10 min), threshold={args.threshold} ({val_col})", fontsize=11, y=1.02)
            fig.tight_layout()
            fp = fig_dir / "go_temporal_scatter_facet.png"
            fig.savefig(fp, bbox_inches="tight", dpi=160)
            plt.close(fig)
            print(f"  scatter facet {fp}")

        if args.plot_scatter:
            for arm in TEMPORAL_ARMS:
                sub = out[out["arm"] == arm.arm_id]
                fig, ax = plt.subplots(figsize=(5.2, 5))
                draw_temporal_scatter_ax(
                    ax,
                    sub,
                    val_col,
                    c5,
                    c10,
                    thr_line,
                    title=f"GO temporal: {arm.arm_id} (vs threshold {args.threshold})",
                    show_legend=True,
                )
                fig.tight_layout()
                fp = fig_dir / f"go_temporal_scatter_{arm.arm_id}.png"
                fig.savefig(fp, bbox_inches="tight", dpi=160)
                plt.close(fig)
                print(f"  scatter {fp}")


if __name__ == "__main__":
    main()
