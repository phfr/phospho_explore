#!/usr/bin/env python3
"""
Second-layer enrichment (Enrichr) on the union of genes from leading-edge lists.

Requires network access. Writes result tables under output/analyze/leading_edge_followup/enrichr/.
"""

from __future__ import annotations

import argparse
import shutil
import tempfile
from pathlib import Path

import pandas as pd

try:
    import gseapy as gp
except ImportError as e:  # pragma: no cover
    raise SystemExit("gseapy is required") from e


def find_repo_root(start: Path) -> Path:
    r = start.resolve()
    for anc in [r, *r.parents]:
        if (anc / "generate_rank_lists.py").is_file():
            return anc
    raise SystemExit("Cannot find project root.")


def collect_genes_from_manifest(manifest_path: Path, root_dir: Path) -> list[str]:
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
            g = line.strip().upper()
            if g:
                genes.add(g)
    return sorted(genes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--manifest", type=Path, default=root / "output" / "analyze" / "leading_edge" / "manifest.tsv")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge_followup" / "enrichr",
    )
    parser.add_argument(
        "--libraries",
        type=str,
        default="GO_Molecular_Function_2021,Reactome_2022,KEGG_2021_Human",
        help="Comma-separated Enrichr library names.",
    )
    args = parser.parse_args()

    manifest = args.manifest.resolve()
    if not manifest.is_file():
        raise SystemExit(f"Missing manifest: {manifest}")
    root_dir = find_repo_root(manifest)
    genes = collect_genes_from_manifest(manifest, root_dir)
    if len(genes) < 3:
        raise SystemExit("Too few genes for Enrichr (need ≥3).")

    libs = [x.strip() for x in args.libraries.split(",") if x.strip()]
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    tmp = Path(tempfile.mkdtemp(prefix="enrichr_"))
    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=libs,
            organism="human",
            outdir=str(tmp),
            no_plot=True,
        )
    except Exception as e:
        shutil.rmtree(tmp, ignore_errors=True)
        raise SystemExit(f"Enrichr failed (network or API): {e}") from e

    shutil.rmtree(tmp, ignore_errors=True)

    if enr is None or enr.results is None or enr.results.empty:
        print("Enrichr returned no results.")
        return

    res = enr.results.copy()
    res.to_csv(out_dir / "enrichr_combined_union.tsv", sep="\t", index=False)
    print(f"Wrote {out_dir / 'enrichr_combined_union.tsv'} ({len(res)} rows)")


if __name__ == "__main__":
    main()
