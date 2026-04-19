#!/usr/bin/env python3
"""
Export leading-edge genes (from gseapy prerank) for selected GO terms — STRING-friendly lists.

Reads output/analyze/go_bp_prerank_long.tsv (must include Lead_genes column).
"""

from __future__ import annotations

import argparse
import urllib.parse
from pathlib import Path

import pandas as pd

from _conditions import CONDITIONS


def string_network_url(genes: list[str], species: int = 9606) -> str:
    """Build STRING 'multiple proteins' search URL (length may truncate for huge lists)."""
    body = "\r".join(sorted(set(genes)))
    q = urllib.parse.urlencode({"identifiers": body, "species": str(species)})
    return f"https://string-db.org/cgi/network?{q}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--long-tsv",
        type=Path,
        default=root / "output" / "analyze" / "go_bp_prerank_long.tsv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "output" / "analyze" / "leading_edge",
    )
    parser.add_argument(
        "--go-ids",
        type=str,
        default="",
        help="Comma-separated GO IDs, e.g. GO:0045087,GO:0007411",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=0,
        help="If >0 and --go-ids empty, export top N terms per condition by ascending pval.",
    )
    parser.add_argument(
        "--max-url-genes",
        type=int,
        default=80,
        help="Max genes embedded in STRING URL (STRING GET limits); larger lists still get .txt files.",
    )
    parser.add_argument("--species", type=int, default=9606, help="NCBI taxonomy id for STRING (9606 = human).")
    args = parser.parse_args()

    long_path = args.long_tsv.resolve()
    if not long_path.is_file():
        raise SystemExit(f"Missing {long_path}")

    df = pd.read_csv(long_path, sep="\t")
    if "Lead_genes" not in df.columns:
        raise SystemExit("Long TSV must contain 'Lead_genes' from prerank output.")

    go_ids = {g.strip() for g in args.go_ids.split(",") if g.strip()}
    col_order = [c.condition_id for c in CONDITIONS]

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict] = []

    if go_ids:
        if "go_id" not in df.columns:
            raise SystemExit("go_id column required when using --go-ids")
        sel = df[df["go_id"].isin(go_ids)]
    elif args.top_n > 0:
        parts = []
        for cond in col_order:
            sub = df[df["condition"] == cond].copy()
            sub = sub.sort_values("pval", ascending=True).head(args.top_n)
            parts.append(sub)
        sel = pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0]
    else:
        raise SystemExit("Provide --go-ids or --top-n > 0.")

    for _, row in sel.iterrows():
        cond = row["condition"]
        gid = row.get("go_id", "")
        if pd.isna(gid) or str(gid).strip().lower() in ("", "nan", "none"):
            gid = ""
        label = row.get("term_label", row.get("Term", ""))
        raw = str(row["Lead_genes"] or "")
        genes = [g.strip() for g in raw.replace(",", ";").split(";") if g.strip()]
        safe = str(gid or label or "term").replace(":", "_").replace("/", "_")[:120]
        fname = f"{cond}__{safe}.txt"
        fpath = out_dir / fname
        fpath.write_text("\n".join(genes) + ("\n" if genes else ""), encoding="utf-8")
        url_genes = genes[: args.max_url_genes]
        url = string_network_url(url_genes, species=args.species) if url_genes else ""
        try:
            rel_file = str(fpath.relative_to(root))
        except ValueError:
            rel_file = str(fpath)
        manifest_rows.append(
            {
                "condition": cond,
                "go_id": gid,
                "term_label": label,
                "n_lead_genes": len(genes),
                "file": rel_file,
                "string_url_truncated": url[:500] + ("..." if len(url) > 500 else ""),
            }
        )

    man = pd.DataFrame(manifest_rows)
    man_path = out_dir / "manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)
    print(f"Wrote {len(manifest_rows)} gene lists under {out_dir}")
    print(f"Manifest: {man_path}")


if __name__ == "__main__":
    main()
