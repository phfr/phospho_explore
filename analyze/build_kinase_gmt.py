#!/usr/bin/env python3
"""
Build a GMT file of kinase -> substrate gene sets from a tabular kinase–substrate table.

Supports plain TSV/CSV and gzip (e.g. PhosphoSitePlus Kinase_Substrate_Dataset.gz from
resources/phosphosite/). For the official PhosphoSite Kinase_Substrate file, use
--psp-kinase-substrate to skip the header preamble rows.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Kinase–substrate TSV/CSV or .gz.")
    parser.add_argument(
        "--out-gmt",
        type=Path,
        default=Path("data/gene_sets/kinase_substrates.gmt"),
        help="Output GMT path.",
    )
    parser.add_argument(
        "--kinase-col",
        type=str,
        default="KINASE",
        help="Column name for kinase gene symbol (PhosphoSite: KINASE).",
    )
    parser.add_argument(
        "--substrate-col",
        type=str,
        default="SUB_GENE",
        help="Column name for substrate gene symbol (PhosphoSite: SUB_GENE).",
    )
    parser.add_argument(
        "--sep",
        type=str,
        default="\t",
        help="Field separator (default tab). Use ',' for CSV.",
    )
    parser.add_argument(
        "--skip-rows",
        type=int,
        default=-1,
        help="Rows to skip before header (default: -1 = auto for Kinase_Substrate_Dataset name, else 0).",
    )
    parser.add_argument(
        "--psp-kinase-substrate",
        action="store_true",
        help="PhosphoSite Kinase_Substrate_Dataset format: skip first 3 preamble lines.",
    )
    parser.add_argument(
        "--human-only",
        action="store_true",
        help="Keep rows where SUB_ORGANISM and KIN_ORGANISM are human (requires those columns).",
    )
    args = parser.parse_args()

    inp = args.input.resolve()
    sep = "\t" if args.sep.lower() in ("tab", "\\t", "\t") else args.sep
    skip = args.skip_rows
    if skip < 0:
        if args.psp_kinase_substrate or "Kinase_Substrate" in inp.name:
            skip = 3
        else:
            skip = 0

    comp = "infer" if str(inp).lower().endswith(".gz") else None
    df = pd.read_csv(
        inp,
        sep=sep,
        comment="#",
        low_memory=False,
        skiprows=skip,
        compression=comp,
        encoding="latin-1",
    )
    if args.human_only:
        for col in ("SUB_ORGANISM", "KIN_ORGANISM"):
            if col not in df.columns:
                raise SystemExit(f"--human-only requires column {col!r}")
        df = df[(df["SUB_ORGANISM"].astype(str).str.lower() == "human") & (df["KIN_ORGANISM"].astype(str).str.lower() == "human")]

    if args.kinase_col not in df.columns or args.substrate_col not in df.columns:
        raise SystemExit(
            f"Columns not found. Have: {list(df.columns)[:40]}...\n"
            f"Expected --kinase-col {args.kinase_col!r} and --substrate-col {args.substrate_col!r}"
        )

    buckets: dict[str, set[str]] = defaultdict(set)
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Kinase→substrate", unit="row"):
        k = str(row[args.kinase_col]).strip()
        s = str(row[args.substrate_col]).strip()
        if not k or k.lower() in ("nan", "none") or not s or s.lower() in ("nan", "none"):
            continue
        buckets[k].add(s)

    out = args.out_gmt.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for kin, genes in sorted(buckets.items()):
        if len(genes) < 3:
            continue
        glist = "\t".join(sorted(genes))
        lines.append(f"{kin}\tPhosphoSite_Kinase_Substrate\t{glist}")

    out.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    print(f"Wrote {out} ({len(lines)} kinases, min 3 substrates each)")


if __name__ == "__main__":
    main()
