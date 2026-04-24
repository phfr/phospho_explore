"""Load and normalize `data/ptm.txt` for joins with `data.tsv` phospho rows.

`ptm.txt` is large (~400k+ rows). Use `dtype=str` for text columns and optional
categoricals for `ptm_type` / `source` after load to reduce memory.

Join keys:
  - Site-level: match `site_key(substrate_UniProtAC, site)` to
    `site_key(proteinid, hl3d)` from phospho TSV (UniProt + residue, uppercased).
  - Gene-level: `substrate_genename` vs `geneid` (symbols, uppercased).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PTM_COLUMNS: list[str] = [
    "ptm_type",
    "source",
    "substrate_UniProtAC",
    "substrate_genename",
    "organism",
    "site",
    "enzyme_UniProtAC",
    "enzyme_genename",
    "note",
    "pmid",
]


def site_key(uniprot: str, site: str) -> str:
    """Canonical key for merging phospho `proteinid` + `hl3d` with PTM substrate + site."""
    u = str(uniprot).strip().upper()
    s = str(site).strip().upper()
    return f"{u}|{s}" if u and s else ""


def load_ptm(path: str | Path, *, categoricals: bool = True) -> pd.DataFrame:
    """Read tab-separated `ptm.txt` with fixed column names (see `data/ptm.info.txt`)."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(p)
    df = pd.read_csv(
        p,
        sep="\t",
        header=None,
        names=PTM_COLUMNS,
        dtype=str,
        na_filter=False,
        low_memory=False,
    )
    if categoricals:
        for col in ("ptm_type", "source"):
            if col in df.columns:
                df[col] = df[col].astype("category")
    return df


def normalize_gene_symbol(s: str) -> str:
    return str(s).strip().upper()


def filter_human(df: pd.DataFrame) -> pd.DataFrame:
    """Keep Homo sapiens rows (handles 'Homo sapiens' and 'Homo sapiens (Human)')."""
    o = df["organism"].astype(str).str.lower()
    mask = o.str.contains("sapiens", na=False)
    return df.loc[mask].copy()


def filter_ptm_type(df: pd.DataFrame, ptm_type: str) -> pd.DataFrame:
    t = ptm_type.strip().upper()
    pt = df["ptm_type"].astype(str).str.upper()
    return df.loc[pt == t].copy()


def add_site_key_ptm(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_site_key"] = [
        site_key(u, s) for u, s in zip(out["substrate_UniProtAC"], out["site"], strict=False)
    ]
    return out


def add_site_key_phospho(df: pd.DataFrame, uniprot_col: str = "proteinid", site_col: str = "hl3d") -> pd.DataFrame:
    out = df.copy()
    out["_site_key"] = [site_key(u, s) for u, s in zip(out[uniprot_col], out[site_col], strict=False)]
    return out


def n_pmids(pmids: str) -> int:
    """Count distinct PMIDs in a comma-separated string."""
    if not pmids or not str(pmids).strip():
        return 0
    parts = [x.strip() for x in str(pmids).split(",") if x.strip() and x.strip().isdigit()]
    return len(set(parts))


def has_enzyme(row: pd.Series) -> bool:
    e1 = str(row.get("enzyme_UniProtAC", "") or "").strip()
    e2 = str(row.get("enzyme_genename", "") or "").strip()
    return bool(e1 or e2)
