#!/usr/bin/env python3
"""
Phospho-level views for genes in leading-edge lists: gene × condition heatmap (rank score)
and per-condition volcano panels (sites restricted to those genes).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from _conditions import CONDITIONS


def find_repo_root(start: Path) -> Path:
    r = start.resolve()
    for anc in [r, *r.parents]:
        if (anc / "generate_rank_lists.py").is_file():
            return anc
    raise SystemExit("Cannot find project root (expected generate_rank_lists.py in a parent directory).")


def aggregate_max_abs_by_gene(df: pd.DataFrame, score_col: str) -> pd.Series:
    def _max_abs(series: pd.Series) -> float:
        s = series.dropna()
        if s.empty:
            return math.nan
        return float(s.loc[s.abs().idxmax()])

    return df.groupby("geneid", sort=False)[score_col].agg(_max_abs)


def load_manifest_genes(manifest_path: Path, root_dir: Path) -> tuple[set[str], list[str]]:
    man = pd.read_csv(manifest_path, sep="\t")
    genes: set[str] = set()
    for rel in man["file"].astype(str):
        p = Path(rel) if Path(rel).is_absolute() else (root_dir / rel)
        if not p.is_file():
            p2 = manifest_path.parent / Path(rel).name
            p = p2 if p2.is_file() else p
        if not p.is_file():
            continue
        for line in p.read_text(encoding="utf-8").splitlines():
            g = line.strip()
            if g:
                genes.add(g)
    return genes, list(man["file"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--manifest", type=Path, default=root / "output" / "analyze" / "leading_edge" / "manifest.tsv")
    parser.add_argument("--data-tsv", type=Path, default=root / "data.tsv")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge_followup",
    )
    parser.add_argument("--max-genes-heatmap", type=int, default=200, help="Cap rows by variance across conditions.")
    args = parser.parse_args()

    manifest = args.manifest.resolve()
    if not manifest.is_file():
        raise SystemExit(f"Missing manifest: {manifest}")
    root_dir = find_repo_root(manifest)

    genes, _ = load_manifest_genes(manifest, root_dir)
    if not genes:
        raise SystemExit("No genes found from manifest file paths.")

    data_path = args.data_tsv.resolve()
    if not data_path.is_file():
        raise SystemExit(f"Missing data.tsv: {data_path}")

    df = pd.read_csv(data_path, sep="\t", low_memory=False)
    if "geneid" not in df.columns:
        raise SystemExit("data.tsv needs geneid")
    df["geneid"] = df["geneid"].astype(str).str.strip()

    rank_cols: dict[str, pd.Series] = {}
    for c in CONDITIONS:
        if c.log2fc_col not in df.columns or c.neglogp_col not in df.columns:
            raise SystemExit(f"Missing columns for {c.condition_id}")
        log2 = pd.to_numeric(df[c.log2fc_col], errors="coerce")
        neg = pd.to_numeric(df[c.neglogp_col], errors="coerce")
        rs = log2 * neg
        rank_cols[c.condition_id] = aggregate_max_abs_by_gene(df.assign(_rs=rs), "_rs")

    mat = pd.DataFrame(rank_cols).reindex(index=sorted(genes))
    mat = mat.replace([np.inf, -np.inf], np.nan).dropna(how="all").fillna(0.0)

    if len(mat) > args.max_genes_heatmap:
        v = mat.var(axis=1).fillna(0)
        mat = mat.loc[v.nlargest(args.max_genes_heatmap).index]

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out_dir / "leading_edge_gene_rankscore_matrix.tsv", sep="\t")

    plt.rcParams.update({"figure.dpi": 120})
    h = min(28.0, 0.18 * len(mat) + 3.0)
    g = sns.clustermap(
        mat,
        cmap="RdBu_r",
        center=0.0,
        figsize=(10, h),
        row_cluster=True,
        col_cluster=True,
        cbar_kws={"label": "RankScore (max |LFC×−log10p| per gene)"},
    )
    g.fig.suptitle("Leading-edge genes: phospho rank score across conditions", y=1.02, fontsize=11)
    g.savefig(out_dir / "leading_edge_phospho_rankscore_clustermap.png", bbox_inches="tight")
    plt.close(g.fig)
    print(f"Wrote {out_dir / 'leading_edge_phospho_rankscore_clustermap.png'}")

    # Volcano grid: sites for genes in full union (not capped)
    genes_all, _ = load_manifest_genes(manifest, root_dir)
    sub = df[df["geneid"].isin(genes_all)].copy()
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), sharex=False, sharey=False)
    axes = np.ravel(axes)
    for ax, c in zip(axes, CONDITIONS):
        log2 = pd.to_numeric(sub[c.log2fc_col], errors="coerce")
        neg = pd.to_numeric(sub[c.neglogp_col], errors="coerce")
        m = log2.notna() & neg.notna()
        ax.scatter(log2[m], neg[m], s=8, alpha=0.35, c="0.35", edgecolors="none", rasterized=True)
        ax.axhline(1.3, color="0.7", ls="--", lw=0.6)
        ax.axvline(0, color="0.7", ls="-", lw=0.5)
        ax.set_title(c.condition_id, fontsize=8)
        ax.set_xlabel("Log2FC", fontsize=7)
        ax.set_ylabel("−log10 p", fontsize=7)
    fig.suptitle("Sites in leading-edge gene union (all lists in manifest)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "leading_edge_phospho_volcano_facets.png", bbox_inches="tight", dpi=140)
    plt.close(fig)
    print(f"Wrote {out_dir / 'leading_edge_phospho_volcano_facets.png'}")


if __name__ == "__main__":
    main()
