"""Helpers for reading ranked lists and GO term strings."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

_GO_ID_RE = re.compile(r"\(GO:(\d+)\)\s*$")


def parse_go_id_from_term(term: str) -> str | None:
    """Extract 'GO:nnnnnnn' from strings like '... (GO:0032270)'."""
    m = _GO_ID_RE.search(term.strip())
    if not m:
        return None
    return f"GO:{m.group(1)}"


def term_label_without_go(term: str) -> str:
    """Strip trailing ' (GO:xxxx)' if present."""
    t = term.strip()
    m = _GO_ID_RE.search(t)
    if m:
        return t[: m.start()].rstrip()
    return t


def prepare_rnk_for_gseapy(rnk_path: Path, dest_path: Path) -> tuple[Path, str]:
    """
    Ensure a two-column tab file: gene, numeric score (gseapy requirement).

    **Two-column .rnk (gene + score):** uses column 0 as gene symbol and column 1 as the
    numeric ranking metric (e.g. RankScore). **Descending positional ranks are never used**
    in this mode.

    **Single-column .rnk (gene only):** falls back to descending positional ranks
    (n..1) so gseapy still has a metric — avoid this if you have RankScore in column 2.
    """
    df = pd.read_csv(rnk_path, sep="\t", header=None, comment="#", dtype=str)
    df = df.apply(lambda c: c.str.strip() if c.dtype == object else c)
    df = df.dropna(how="all")

    ncols = df.shape[1]
    if ncols >= 2:
        genes = df.iloc[:, 0].astype(str)
        scores = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        out = pd.DataFrame({"gene": genes, "score": scores})
        out = out.dropna(subset=["score"])
        if out.empty:
            raise ValueError(
                f"{rnk_path}: second column has no numeric scores after parsing. "
                "Expected gene<TAB>RankScore (or other numeric metric)."
            )
        mode = "gene_plus_score"
    else:
        genes = df.iloc[:, 0].astype(str)
        n = len(genes)
        out = pd.DataFrame({"gene": genes, "score": range(n, 0, -1)})
        mode = "rank_fallback"

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest_path, sep="\t", header=False, index=False)
    return dest_path, mode
