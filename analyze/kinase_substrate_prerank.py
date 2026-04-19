#!/usr/bin/env python3
"""
Gene-level kinase–substrate enrichment: gseapy.prerank with a kinase GMT (KSEA-like summary).

Classic KSEA is often site-weighted; this uses substrate gene sets vs your ranked gene list.
See analyze/README.md for building the GMT from PhosphoSitePlus.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import gseapy as gp
import pandas as pd
from tqdm import tqdm

from _conditions import CONDITIONS
from go_utils import prepare_rnk_for_gseapy

DEFAULT_KINASE_GMT = Path("data/gene_sets/kinase_substrates.gmt")
SAMPLE_KINASE_GMT = Path("data/gene_sets/kinase_substrates_sample.gmt")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--rnk-dir", type=Path, default=root / "output")
    parser.add_argument("--gmt-path", type=Path, default=None, help="Kinase GMT (default: data/gene_sets/kinase_substrates.gmt if present).")
    parser.add_argument(
        "--out-long",
        type=Path,
        default=root / "output" / "analyze" / "ksea_prerank_long.tsv",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=root / "output" / "analyze" / "kinase_prerank",
    )
    parser.add_argument("--permutation-num", type=int, default=1000)
    parser.add_argument("--min-size", type=int, default=5)
    parser.add_argument("--max-size", type=int, default=400)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--heatmap-top",
        type=int,
        default=25,
        help="Top kinases per condition by -log10(p) for summary heatmap (0=skip plot).",
    )
    parser.add_argument(
        "--heatmap-cap-neglog10",
        type=float,
        default=12.0,
        help="Cap -log10(p) in the PNG heatmap only (nominal p=0 -> 300 would hide weaker signals).",
    )
    args = parser.parse_args()

    if args.gmt_path:
        gmt = args.gmt_path.resolve()
    else:
        g1 = (root / DEFAULT_KINASE_GMT).resolve()
        g2 = (root / SAMPLE_KINASE_GMT).resolve()
        gmt = g1 if g1.is_file() else g2
    if not gmt.is_file():
        raise SystemExit(
            f"Kinase GMT not found: {gmt}\n"
            "Create one with: python analyze/build_kinase_gmt.py --input your_psp_table.tsv ...\n"
            f"Or keep the bundled sample at {SAMPLE_KINASE_GMT}."
        )

    rnk_dir = args.rnk_dir.resolve()
    out_root = args.out_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    expected = {c.condition_id for c in CONDITIONS}

    all_rnk = sorted(rnk_dir.glob("*.rnk"))
    to_run = [p for p in all_rnk if p.stem in expected]
    for p in all_rnk:
        if p.stem not in expected:
            print(f"skip {p.name}")

    parts: list[pd.DataFrame] = []
    for rnk_path in tqdm(to_run, desc="Kinase prerank", unit="rnk"):
        stem = rnk_path.stem
        cond_dir = out_root / stem
        if cond_dir.exists():
            shutil.rmtree(cond_dir)
        cond_dir.mkdir(parents=True, exist_ok=True)
        prep, mode = prepare_rnk_for_gseapy(rnk_path, cond_dir / "_input_for_gseapy.rnk")
        if mode == "rank_fallback":
            tqdm.write(
                f"WARN {rnk_path.name}: single-column .rnk — using positional ranks; "
                "use two-column gene<TAB>RankScore.",
                file=sys.stderr,
            )
        pre = gp.prerank(
            rnk=str(prep),
            gene_sets=str(gmt),
            outdir=str(cond_dir),
            permutation_num=args.permutation_num,
            min_size=args.min_size,
            max_size=args.max_size,
            threads=args.threads,
            no_plot=True,
            verbose=False,
        )
        df = pre.res2d.copy()
        df.insert(0, "condition", stem)
        df.rename(
            columns={k: v for k, v in {"NOM p-val": "pval", "FDR q-val": "fdr"}.items() if k in df.columns},
            inplace=True,
        )
        if "pval" not in df.columns and "P-value" in df.columns:
            df.rename(columns={"P-value": "pval"}, inplace=True)
        df["kinase"] = df["Term"]
        parts.append(df)
        tqdm.write(f"OK {stem}: {len(df)} kinases ({mode})")

    if not parts:
        raise SystemExit("No kinase prerank results.")

    long = pd.concat(parts, ignore_index=True)
    out_long = args.out_long.resolve()
    out_long.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(out_long, sep="\t", index=False)
    print(f"Wrote {out_long}")

    if args.heatmap_top > 0:
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns

        pv = np.asarray(pd.to_numeric(long["pval"], errors="coerce"), dtype=float)
        pv = np.clip(pv, 1e-300, np.inf)
        long["neglog10p"] = -np.log10(pv)
        tops = []
        for cond in [c.condition_id for c in CONDITIONS]:
            sub = long[long["condition"] == cond].nlargest(args.heatmap_top, "neglog10p")
            tops.append(sub)
        topdf = pd.concat(tops)
        kinases = sorted(topdf["kinase"].unique())
        mat = topdf.pivot_table(index="kinase", columns="condition", values="neglog10p", aggfunc="max")
        mat = mat.reindex(index=[k for k in kinases if k in mat.index])
        col_order = [c.condition_id for c in CONDITIONS]
        mat = mat[[c for c in col_order if c in mat.columns]]
        cap = float(args.heatmap_cap_neglog10)
        mat_show = mat.fillna(0).clip(upper=cap)
        fig, ax = plt.subplots(figsize=(14, max(4, 0.25 * len(mat))))
        sns.heatmap(
            mat_show,
            cmap="viridis",
            ax=ax,
            vmin=0.0,
            vmax=cap,
            cbar_kws={"label": f"-log10(nominal p), capped at {cap}"},
        )
        ax.set_title("Kinase substrate sets (prerank): top kinases per condition (color capped)")
        fig_dir = root / "output" / "analyze" / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        fp = fig_dir / "kinase_crosscondition_heatmap.png"
        fig.tight_layout()
        fig.savefig(fp, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {fp}")


if __name__ == "__main__":
    main()
