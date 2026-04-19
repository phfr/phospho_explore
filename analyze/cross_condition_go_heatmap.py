#!/usr/bin/env python3
"""
Build a GO term × condition heatmap from go_bp_prerank_long.tsv (-log10 p-value in cells).

Use --terms-tsv to restrict rows (e.g. go_temporal_classification_nonabsent.tsv) so every
column is still all experimental conditions, but rows are only temporally interesting GO sets.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from _conditions import CONDITIONS


def render_crosscondition_dotplot(
    wide_plot: pd.DataFrame,
    nes_plot: pd.DataFrame | None,
    out_path: Path,
    cap_neglog10: float,
    *,
    title: str,
    figwidth: float,
    figheight: float,
) -> None:
    """Bubble grid: x = condition, y = GO term, size = -log10(p), color = NES when available."""
    n_term, n_cond = wide_plot.shape
    vals = np.clip(np.nan_to_num(wide_plot.to_numpy(dtype=float), nan=0.0), 0.0, cap_neglog10)
    yy, xx = np.meshgrid(np.arange(n_term), np.arange(n_cond), indexing="ij")
    flat_x = xx.ravel()
    flat_y = (n_term - 1 - yy).ravel()
    s_norm = np.clip(vals.ravel() / cap_neglog10, 0.0, 1.0)
    s_pts = 12.0 + s_norm * 200.0

    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    if nes_plot is not None and nes_plot.shape == wide_plot.shape:
        nes_m = nes_plot.to_numpy(dtype=float)
        flat_c = nes_m.ravel()
        finite = np.isfinite(flat_c)
        abs_n = float(np.nanpercentile(np.abs(nes_m[np.isfinite(nes_m)]), 95)) if np.any(np.isfinite(nes_m)) else 1.0
        abs_n = max(abs_n, 0.5)
        if (~finite).any():
            ax.scatter(
                flat_x[~finite],
                flat_y[~finite],
                s=np.clip(s_pts[~finite], 6.0, 35.0),
                c="#999999",
                alpha=0.5,
                edgecolors="none",
                rasterized=True,
            )
        if finite.any():
            sc = ax.scatter(
                flat_x[finite],
                flat_y[finite],
                c=flat_c[finite],
                s=s_pts[finite],
                cmap="RdBu_r",
                vmin=-abs_n,
                vmax=abs_n,
                edgecolors="0.25",
                linewidths=0.15,
                alpha=0.92,
                rasterized=True,
            )
            fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.02, label="NES")
    else:
        sc = ax.scatter(
            flat_x,
            flat_y,
            c=vals.ravel(),
            s=s_pts,
            cmap="YlOrRd",
            vmin=0.0,
            vmax=cap_neglog10,
            edgecolors="0.25",
            linewidths=0.15,
            alpha=0.9,
            rasterized=True,
        )
        fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.02, label="-log10(nominal p)")

    ax.set_xticks(np.arange(n_cond))
    ax.set_xticklabels(list(wide_plot.columns), rotation=55, ha="right", fontsize=8)
    ax.set_yticks(np.arange(n_term))
    yfs = max(2.5, min(7.0, 520.0 / max(n_term, 1)))
    ax.set_yticklabels(list(wide_plot.index), fontsize=yfs)
    ax.set_xlim(-0.5, n_cond - 0.5)
    ax.set_ylim(-0.5, n_term - 0.5)
    ax.grid(True, axis="both", color="0.88", lw=0.35, zorder=0)
    ax.set_title(title + " (marker size ∝ -log10 p)")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", dpi=140)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--long-tsv",
        type=Path,
        default=root / "output" / "analyze" / "go_bp_prerank_long.tsv",
        help="Long table from run_go_prerank_batch.py",
    )
    parser.add_argument(
        "--out-png",
        type=Path,
        default=root / "output" / "analyze" / "figures" / "go_crosscondition_heatmap.png",
    )
    parser.add_argument(
        "--out-matrix",
        type=Path,
        default=root / "output" / "analyze" / "figures" / "go_crosscondition_matrix.tsv",
    )
    parser.add_argument(
        "--fdr-any",
        type=float,
        default=0.25,
        help="Keep GO rows significant in at least one condition (FDR ≤ this). Use --pval-any for p-value filter.",
    )
    parser.add_argument(
        "--pval-any",
        type=float,
        default=0.0,
        help="If >0, keep rows where min nominal p-value across conditions ≤ this (overrides FDR filter when set).",
    )
    parser.add_argument(
        "--cap-neglog10",
        type=float,
        default=20.0,
        help="Cap -log10(p) for color scale (default 20).",
    )
    parser.add_argument(
        "--clustermap",
        action="store_true",
        help="Use seaborn.clustermap instead of heatmap (reorders rows/cols).",
    )
    parser.add_argument(
        "--terms-tsv",
        type=Path,
        default=None,
        help="TSV with a term_label column; keep only those GO rows (union of labels). "
        "Skips --fdr-any / --pval-any row filtering. Typical: go_temporal_classification_nonabsent.tsv",
    )
    parser.add_argument("--figwidth", type=float, default=14.0)
    parser.add_argument("--figheight", type=float, default=24.0, help="Max figure height (rows are shrunk if needed).")
    parser.add_argument(
        "--dotplot",
        action="store_true",
        help="Also write a bubble/dot grid (size = -log10 p, color = NES) next to the heatmap style.",
    )
    parser.add_argument(
        "--out-dot-png",
        type=Path,
        default=None,
        help="Output path for dot plot (default: same stem as --out-png with _dotplot.png).",
    )
    parser.add_argument(
        "--dotplot-figwidth",
        type=float,
        default=None,
        help="Figure width for dot plot only (defaults to --figwidth).",
    )
    parser.add_argument(
        "--dotplot-figheight",
        type=float,
        default=None,
        help="Figure height for dot plot only (defaults scaled like heatmap).",
    )
    args = parser.parse_args()

    long_path: Path = args.long_tsv.resolve()
    if not long_path.is_file():
        raise SystemExit(f"Missing long TSV: {long_path}. Run analyze/run_go_prerank_batch.py first.")

    df = pd.read_csv(long_path, sep="\t")
    if "pval" not in df.columns:
        raise SystemExit("Expected column 'pval' in long TSV.")
    col_order = [c.condition_id for c in CONDITIONS]
    df = df[df["condition"].isin(col_order)]

    pv = np.asarray(pd.to_numeric(df["pval"], errors="coerce"), dtype=float)
    pv = np.clip(pv, 1e-300, np.inf)
    df["neglog10p"] = -np.log10(pv)
    df["neglog10p"] = df["neglog10p"].clip(upper=args.cap_neglog10)

    label_col = "term_label" if "term_label" in df.columns else "Term"
    id_col = "go_id" if "go_id" in df.columns else None

    wide = df.pivot_table(
        index=label_col,
        columns="condition",
        values="neglog10p",
        aggfunc="max",
    )
    wide = wide.reindex(columns=[c for c in col_order if c in wide.columns])

    if args.terms_tsv is not None:
        tpath = args.terms_tsv.resolve()
        if not tpath.is_file():
            raise SystemExit(f"Missing --terms-tsv: {tpath}")
        lab = pd.read_csv(tpath, sep="\t")
        if "term_label" not in lab.columns:
            raise SystemExit("--terms-tsv must contain a term_label column.")
        keep = lab["term_label"].dropna().unique()
        hit = wide.index.intersection(pd.Index(keep))
        wide_plot = wide.loc[hit].sort_index()
        if wide_plot.empty:
            raise SystemExit("No overlap between --terms-tsv term_label values and long TSV terms.")
    elif "fdr" in df.columns and args.pval_any <= 0:
        fdr_wide = df.pivot_table(index=label_col, columns="condition", values="fdr", aggfunc="min")
        fdr_wide = fdr_wide.reindex(columns=[c for c in col_order if c in fdr_wide.columns])
        mask_sig = (fdr_wide <= args.fdr_any).any(axis=1)
        wide_plot = wide.loc[mask_sig]
    elif args.pval_any > 0:
        pwide = df.pivot_table(index=label_col, columns="condition", values="pval", aggfunc="min")
        pwide = pwide.reindex(columns=[c for c in col_order if c in pwide.columns])
        mask_sig = (pwide <= args.pval_any).any(axis=1)
        wide_plot = wide.loc[mask_sig]
    else:
        wide_plot = wide

    if wide_plot.empty:
        raise SystemExit("No rows after significance filter; loosen --fdr-any or --pval-any.")

    nes_plot: pd.DataFrame | None = None
    if "NES" in df.columns:
        nes_wide = df.pivot_table(index=label_col, columns="condition", values="NES", aggfunc="first")
        nes_wide = nes_wide.reindex(columns=[c for c in col_order if c in nes_wide.columns])
        nes_plot = nes_wide.reindex(index=wide_plot.index, columns=wide_plot.columns)

    out_png: Path = args.out_png.resolve()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_matrix: Path = args.out_matrix.resolve()
    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    wide_plot.to_csv(out_matrix, sep="\t")

    plt.rcParams.update({"figure.dpi": 120})
    h = min(args.figheight, 0.22 * len(wide_plot) + 3)
    w = args.figwidth

    if args.clustermap:
        g = sns.clustermap(
            wide_plot.fillna(0.0),
            cmap="YlOrRd",
            figsize=(w, h),
            row_cluster=True,
            col_cluster=True,
            cbar_kws={"label": "-log10(nominal p)"},
        )
        g.savefig(out_png, bbox_inches="tight")
        plt.close(g.fig)
    else:
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(
            wide_plot.fillna(0.0),
            cmap="YlOrRd",
            ax=ax,
            cbar_kws={"label": "-log10(nominal p)"},
        )
        title = "GO BP: -log10(nominal p), selected terms" if args.terms_tsv else "GO Biological Process: -log10(nominal p) from prerank"
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(out_png, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote {out_png} ({wide_plot.shape[0]} GO terms × {wide_plot.shape[1]} conditions)")
    print(f"Matrix: {out_matrix}")

    if args.dotplot:
        dot_path = (args.out_dot_png or out_png.with_name(out_png.stem + "_dotplot.png")).resolve()
        dw = args.dotplot_figwidth or args.figwidth
        dh = args.dotplot_figheight or min(args.figheight, 0.11 * len(wide_plot) + 3.5)
        dtitle = "GO BP: selected terms" if args.terms_tsv else "GO BP: prerank across conditions"
        render_crosscondition_dotplot(
            wide_plot,
            nes_plot,
            dot_path,
            args.cap_neglog10,
            title=dtitle,
            figwidth=dw,
            figheight=dh,
        )
        print(f"Wrote dot plot {dot_path}")


if __name__ == "__main__":
    main()
