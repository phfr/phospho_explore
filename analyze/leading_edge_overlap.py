#!/usr/bin/env python3
"""
Overlap structure across leading-edge gene lists: membership matrix, pairwise Jaccard,
and top shared intersections (no Cytoscape/STRING workflow).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_sets(leading_dir: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for p in sorted(leading_dir.glob("*.txt")):
        genes: set[str] = set()
        for line in p.read_text(encoding="utf-8").splitlines():
            g = line.strip()
            if g:
                genes.add(g)
        if genes:
            sets[p.name] = genes
    return sets


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--leading-edge-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge_followup",
    )
    parser.add_argument("--max-sets-jaccard", type=int, default=60, help="If more lists, keep largest by gene count.")
    args = parser.parse_args()

    led = args.leading_edge_dir.resolve()
    if not led.is_dir():
        raise SystemExit(f"Missing directory: {led}")

    sets_map = load_sets(led)
    if len(sets_map) < 2:
        raise SystemExit("Need at least two non-empty .txt lists in leading-edge dir.")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_genes = sorted(set().union(*sets_map.values()))
    names = list(sets_map.keys())
    if len(names) > args.max_sets_jaccard:
        sizes = {n: len(sets_map[n]) for n in names}
        names = sorted(names, key=lambda x: sizes[x], reverse=True)[: args.max_sets_jaccard]
        sets_map = {k: sets_map[k] for k in names}

    inc = pd.DataFrame(0, index=all_genes, columns=list(sets_map.keys()), dtype=np.int8)
    for col, gset in sets_map.items():
        inc.loc[inc.index.intersection(pd.Index(list(gset))), col] = 1
    inc.to_csv(out_dir / "leading_edge_membership_matrix.tsv", sep="\t")
    print(f"Wrote membership matrix {inc.shape}")

    names = list(sets_map.keys())
    n = len(names)
    jac = np.ones((n, n), dtype=float)
    inter = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = sets_map[names[i]], sets_map[names[j]]
            u = len(a | b)
            k = len(a & b)
            inter[i, j] = inter[j, i] = k
            jac[i, j] = jac[j, i] = (k / u) if u else 0.0

    jdf = pd.DataFrame(jac, index=names, columns=names)
    jdf.to_csv(out_dir / "leading_edge_pairwise_jaccard.tsv", sep="\t")

    rows: list[dict] = []
    for i in range(n):
        for j in range(i + 1, n):
            rows.append(
                {
                    "set_a": names[i],
                    "set_b": names[j],
                    "intersection": int(inter[i, j]),
                    "jaccard": jac[i, j],
                }
            )
    pair_df = pd.DataFrame(rows).sort_values("intersection", ascending=False)
    pair_df.to_csv(out_dir / "leading_edge_pairwise_overlap_long.tsv", sep="\t", index=False)

    plt.rcParams.update({"figure.dpi": 120})
    w = min(22.0, 0.22 * n + 4)
    fig, ax = plt.subplots(figsize=(w, w))
    sns.heatmap(jdf, cmap="YlOrRd", vmin=0, vmax=1, ax=ax, square=True, cbar_kws={"label": "Jaccard"})
    ax.set_title("Leading-edge gene lists: pairwise Jaccard")
    fig.tight_layout()
    fig.savefig(out_dir / "leading_edge_jaccard_heatmap.png", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_dir / 'leading_edge_jaccard_heatmap.png'}")

    top = pair_df.head(25)
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(top))
    ax.barh(y, top["intersection"].to_numpy(), color="steelblue")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ∩ {b}"[:100] for a, b in zip(top["set_a"], top["set_b"])], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Intersection size")
    ax.set_title("Top 25 pairwise intersections of leading-edge lists")
    fig.tight_layout()
    fig.savefig(out_dir / "leading_edge_top_intersections.png", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_dir / 'leading_edge_top_intersections.png'}")


if __name__ == "__main__":
    main()
