#!/usr/bin/env python3
"""
For each leading-edge gene list, test overlap with kinase substrate sets (GMT)
using Fisher's exact test vs the gene universe from data.tsv.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.stats import fisher_exact
except ImportError as e:
    raise SystemExit("scipy is required (pip install scipy)") from e


def find_repo_root(start: Path) -> Path:
    r = start.resolve()
    for anc in [r, *r.parents]:
        if (anc / "generate_rank_lists.py").is_file():
            return anc
    raise SystemExit("Cannot find project root.")


def parse_gmt(path: Path) -> dict[str, set[str]]:
    kin: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = {p.strip().upper() for p in parts[2:] if p.strip()}
            if name and genes:
                kin[name] = genes
    return kin


def collect_leading_sets(leading_dir: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for p in sorted(leading_dir.glob("*.txt")):
        g = {line.strip().upper() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()}
        if g:
            out[p.name] = g
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--leading-edge-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge",
    )
    parser.add_argument("--data-tsv", type=Path, default=root / "data.tsv")
    parser.add_argument(
        "--gmt-path",
        type=Path,
        default=None,
        help="Kinase GMT (default: kinase_substrates.gmt or sample GMT).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge_followup",
    )
    parser.add_argument("--fdr-alpha", type=float, default=0.1, help="BH-FDR within each leading-edge list.")
    args = parser.parse_args()

    led = args.leading_edge_dir.resolve()
    if not led.is_dir():
        raise SystemExit(f"Missing {led}")

    data_path = args.data_tsv.resolve()
    if not data_path.is_file():
        raise SystemExit(f"Missing {data_path}")

    gmt_path = args.gmt_path
    if gmt_path is None:
        g1 = (root / "data/gene_sets/kinase_substrates.gmt").resolve()
        g2 = (root / "data/gene_sets/kinase_substrates_sample.gmt").resolve()
        gmt_path = g1 if g1.is_file() else g2
    else:
        gmt_path = gmt_path.resolve()
    if not gmt_path.is_file():
        raise SystemExit(f"Kinase GMT not found: {gmt_path}")

    universe = (
        pd.read_csv(data_path, sep="\t", usecols=["geneid"], low_memory=False)["geneid"]
        .astype(str)
        .str.strip()
        .str.upper()
        .dropna()
        .unique()
    )
    U = set(universe)
    n_u = len(U)

    kinases = parse_gmt(gmt_path)
    for k, s in kinases.items():
        kinases[k] = {g for g in s if g in U}

    le_sets = collect_leading_sets(led)
    if not le_sets:
        raise SystemExit("No non-empty leading-edge .txt files.")

    rows: list[dict] = []
    for le_name, Lraw in le_sets.items():
        L = {g for g in Lraw if g in U}
        n_l = len(L)
        if n_l < 1:
            continue
        for kin_name, S in kinases.items():
            S_in = {g for g in S if g in U}
            n_s = len(S_in)
            if n_s < 3:
                continue
            a = len(L & S_in)
            b = n_l - a
            c = n_s - a
            d = n_u - n_l - n_s + a
            if d < 0:
                d = 0
            _, p_two = fisher_exact([[a, b], [c, d]])
            rows.append(
                {
                    "leading_edge_list": le_name,
                    "kinase": kin_name,
                    "overlap": a,
                    "lead_n": n_l,
                    "substrate_n": n_s,
                    "universe_n": n_u,
                    "fisher_p_two_sided": p_two,
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        raise SystemExit("No overlap rows (check gene symbol overlap with GMT / data.tsv).")

    out["fdr_bh"] = out.groupby("leading_edge_list", group_keys=False)["fisher_p_two_sided"].transform(
        lambda s: pd.Series(_benjamini_hochberg_fdr(s.astype(float).to_numpy()), index=s.index)
    )
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "leading_edge_kinase_fisher.tsv"
    out.sort_values(["leading_edge_list", "fisher_p_two_sided"]).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path}")

    sig = out[out["fdr_bh"] <= args.fdr_alpha].sort_values(["leading_edge_list", "fdr_bh"])
    sig.to_csv(out_dir / "leading_edge_kinase_fisher_fdr.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'leading_edge_kinase_fisher_fdr.tsv'} ({len(sig)} rows at FDR≤{args.fdr_alpha})")


def _benjamini_hochberg_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    sorted_p = p[order]
    q = sorted_p * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


if __name__ == "__main__":
    main()
