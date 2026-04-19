#!/usr/bin/env python3
"""
Run gseapy.prerank (GSEA preranked) for each .rnk against GO Biological Process gene sets.

Writes:
  - output/analyze/go_prerank/<condition>/  (per-run gseapy outputs)
  - output/analyze/go_bp_prerank_long.tsv   (concatenated long table)

Gene sets: place MSigDB `data/gene_sets/c5.go.bp.*.Hs.symbols.gmt` (C5 GO BP, human symbols), or pass
--gmt-path, or use --enrichr-library GO_Biological_Process_2021 (network).
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
from go_utils import parse_go_id_from_term, prepare_rnk_for_gseapy, term_label_without_go

def discover_msigdb_go_bp_gmt(project_root: Path) -> Path | None:
    """Prefer any data/gene_sets/c5.go.bp.*.Hs.symbols.gmt (newest name wins), else legacy v2023.1 path."""
    gene_dir = project_root / "data" / "gene_sets"
    if not gene_dir.is_dir():
        return None
    found = sorted(gene_dir.glob("c5.go.bp.*.Hs.symbols.gmt"))
    if found:
        return found[-1].resolve()
    legacy = (gene_dir / "c5.go.bp.v2023.1.Hs.symbols.gmt").resolve()
    return legacy if legacy.is_file() else None


def truncate_gmt(src: Path, max_sets: int, dest: Path) -> Path:
    lines: list[str] = []
    with src.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= max_sets:
                break
            lines.append(line.rstrip("\n"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--rnk-dir",
        type=Path,
        default=root / "output",
        help="Directory containing .rnk files (default: project output/)",
    )
    parser.add_argument(
        "--gmt-path",
        type=Path,
        default=None,
        help="Path to MSigDB GO BP GMT (e.g. data/gene_sets/c5.go.bp.v2026.1.Hs.symbols.gmt).",
    )
    parser.add_argument(
        "--enrichr-library",
        type=str,
        default="",
        help="If set (and --gmt-path omitted), use this Enrichr library name, e.g. GO_Biological_Process_2021.",
    )
    parser.add_argument(
        "--out-long",
        type=Path,
        default=root / "output" / "analyze" / "go_bp_prerank_long.tsv",
        help="Output long-format TSV path.",
    )
    parser.add_argument(
        "--out-prerank-root",
        type=Path,
        default=root / "output" / "analyze" / "go_prerank",
        help="Directory for per-condition gseapy output folders.",
    )
    parser.add_argument("--permutation-num", type=int, default=1000)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--max-sets",
        type=int,
        default=0,
        help="If >0, only the first N gene sets from the GMT are used (smoke tests).",
    )
    parser.add_argument(
        "--conditions",
        type=str,
        default="",
        help="Comma-separated condition ids to run (default: all 8 matching .rnk stems).",
    )
    args = parser.parse_args()

    rnk_dir: Path = args.rnk_dir.resolve()
    out_long: Path = args.out_long.resolve()
    out_root: Path = args.out_prerank_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    gmt_path: Path | None = args.gmt_path.resolve() if args.gmt_path else None
    enrichr = (args.enrichr_library or "").strip()

    if gmt_path is None and not enrichr:
        candidate = discover_msigdb_go_bp_gmt(root)
        if candidate is not None and candidate.is_file():
            gmt_path = candidate
            print(f"Using local GO BP GMT: {gmt_path.relative_to(root)}")
        else:
            enrichr = "GO_Biological_Process_2021"
            print(
                "No data/gene_sets/c5.go.bp.*.Hs.symbols.gmt found; using Enrichr "
                "GO_Biological_Process_2021 (requires network). "
                "Add MSigDB C5 GO BP symbols GMT under data/gene_sets/ or pass --gmt-path."
            )

    if gmt_path is not None and enrichr:
        raise SystemExit("Use only one of --gmt-path or --enrichr-library.")

    trunc_gmt_path: Path | None = None
    gene_sets: str
    if gmt_path is not None:
        if not gmt_path.is_file():
            raise SystemExit(f"GMT not found: {gmt_path}")
        if args.max_sets > 0:
            trunc_gmt_path = out_root / "_truncated_gmt_input.gmt"
            truncate_gmt(gmt_path, args.max_sets, trunc_gmt_path)
            gene_sets = str(trunc_gmt_path)
        else:
            gene_sets = str(gmt_path)
    else:
        gene_sets = enrichr

    expected = {c.condition_id for c in CONDITIONS}
    only = {s.strip() for s in args.conditions.split(",") if s.strip()}
    if only:
        expected = expected & only
    all_rnk = sorted(rnk_dir.glob("*.rnk"))
    for p in all_rnk:
        if p.stem not in expected:
            print(f"skip (unknown condition id): {p.name}")
    rnk_files = [p for p in all_rnk if p.stem in expected]
    if not rnk_files:
        raise SystemExit(f"No matching .rnk files in {rnk_dir} for selected conditions.")

    long_parts: list[pd.DataFrame] = []

    for rnk_path in tqdm(rnk_files, desc="GO prerank", unit="rnk"):
        stem = rnk_path.stem
        cond_dir = out_root / stem
        if cond_dir.exists():
            shutil.rmtree(cond_dir)
        cond_dir.mkdir(parents=True, exist_ok=True)

        prep, mode = prepare_rnk_for_gseapy(rnk_path, cond_dir / "_input_for_gseapy.rnk")
        if mode == "rank_fallback":
            tqdm.write(
                f"WARN {rnk_path.name}: single-column .rnk — using positional rank scores "
                f"(not column-2 RankScore). Prefer two-column gene<TAB>score.",
                file=sys.stderr,
            )
        pre = gp.prerank(
            rnk=str(prep),
            gene_sets=gene_sets,
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
        ren = {
            "NOM p-val": "pval",
            "FDR q-val": "fdr",
            "FWER p-val": "fwer",
            "Tag %": "tag_pct",
            "Gene %": "gene_pct",
        }
        df.rename(columns={k: v for k, v in ren.items() if k in df.columns}, inplace=True)
        df["go_id"] = df["Term"].map(parse_go_id_from_term)
        df["term_label"] = df["Term"].map(term_label_without_go)
        long_parts.append(df)
        tqdm.write(f"OK {stem}: {len(df)} gene sets in res2d ({mode})")

    if not long_parts:
        raise SystemExit("No matching .rnk files for known condition ids.")

    long = pd.concat(long_parts, ignore_index=True)
    out_long.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(out_long, sep="\t", index=False)
    print(f"Wrote {out_long} ({len(long)} rows)")
    if trunc_gmt_path is not None and trunc_gmt_path.exists():
        trunc_gmt_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
