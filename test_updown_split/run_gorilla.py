#!/usr/bin/env python3
"""
Submit ranked .rnk lists to GOrilla (Technion) and save Biological Process GO tree PNGs.

Same behavior as repo-root `run_gorilla.py`; defaults use **this folder**:
  --input-dir  test_updown_split/output/
  --output-dir test_updown_split/gorilla_bp/

Depends on: pip install requests tqdm
Service: https://cbl-gorilla.cs.technion.ac.il/
"""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import requests
from tqdm import tqdm

GORILLA_BASE = "https://cbl-gorilla.cs.technion.ac.il"
SERVLET_URL = f"{GORILLA_BASE}/servlet/GOrilla"

SPECIES_CHOICES = {
    "human": "HOMO_SAPIENS",
    "mouse": "MUS_MUSCULUS",
    "rat": "RATTUS_NORVEGICUS",
    "fly": "DROSOPHILA_MELANOGASTER",
    "worm": "CAENORHABDITIS_ELEGANS",
    "yeast": "SACCHAROMYCES_CEREVISIAE",
    "zebrafish": "DANIO_RERIO",
    "arabidopsis": "ARABIDOPSIS_THALIANA",
}


def read_ranked_list(path: Path, genes_only: bool) -> str:
    """Return newline-separated list for GOrilla (rank order preserved)."""
    rows: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        rows.append(line.split("\t", 1)[0] if genes_only else line)
    return "\n".join(rows)


def extract_run_id(results_url: str) -> str | None:
    m = re.search(r"/GOrilla/([^/]+)/", results_url)
    return m.group(1) if m else None


def wait_for_results(
    session: requests.Session, start_url: str, poll_s: float, timeout_s: float
) -> requests.Response:
    deadline = time.monotonic() + timeout_s
    last: requests.Response | None = None
    while time.monotonic() < deadline:
        last = session.get(start_url, timeout=120, allow_redirects=True)
        if "Calculating Enrichment" not in last.text:
            return last
        time.sleep(poll_s)
    raise TimeoutError(f"GOrilla still calculating after {timeout_s}s: {start_url}")


def submit_gorilla(
    session: requests.Session,
    *,
    ranked_text: str,
    species_value: str,
    pvalue_thresh: str,
    fast_mode: bool,
) -> tuple[requests.Response, str]:
    """POST job; return (results_response, final_results_url)."""
    fields: dict[str, tuple[None, str]] = {
        "application": (None, "gorilla"),
        "species": (None, species_value),
        "run_mode": (None, "mhg"),
        "target_set": (None, ranked_text),
        "db": (None, "proc"),  # Biological process
        "pvalue_thresh": (None, pvalue_thresh),
        "run_gogo_button": (None, "Search Enriched GO terms"),
    }
    if fast_mode:
        fields["fast_mode"] = (None, "on")

    post = session.post(SERVLET_URL, files=fields, timeout=300, allow_redirects=False)
    post.raise_for_status()
    if post.status_code != 302:
        raise RuntimeError(f"Unexpected POST status {post.status_code} (expected 302)")
    loc = post.headers.get("Location")
    if not loc:
        raise RuntimeError("POST missing Location header")
    if loc.startswith("http"):
        start_url = loc
    else:
        start_url = f"{GORILLA_BASE}/servlet/{loc}"

    results = wait_for_results(session, start_url, poll_s=0.5, timeout_s=600)
    return results, results.url


def download_go_png(session: requests.Session, results_url: str) -> bytes:
    run_id = extract_run_id(results_url)
    if not run_id:
        raise RuntimeError(f"Could not parse GOrilla run id from URL: {results_url}")
    png_url = f"{GORILLA_BASE}/GOrilla/{run_id}/GO.png"
    r = session.get(png_url, timeout=300)
    if r.status_code == 404:
        raise FileNotFoundError(f"No GO tree PNG at {png_url} (often means no enrichment).")
    r.raise_for_status()
    if "image" not in (r.headers.get("content-type") or "").lower():
        raise RuntimeError(f"Expected PNG, got {r.headers.get('content-type')} from {png_url}")
    return r.content


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=here / "output",
        help="Directory containing .rnk files (default: test_updown_split/output/)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=here / "gorilla_bp",
        help="Where to write PNGs (default: test_updown_split/gorilla_bp/)",
    )
    parser.add_argument(
        "--species",
        choices=sorted(SPECIES_CHOICES.keys()),
        default="human",
        help="Organism (default: human / Homo sapiens)",
    )
    parser.add_argument(
        "--pvalue",
        default="0.001",
        choices=[
            "0.001",
            "0.0001",
            "0.00001",
            "0.000001",
            "0.0000001",
            "0.00000001",
            "0.000000001",
            "0.0000000001",
            "0.00000000001",
        ],
        help="GOrilla p-value cutoff (default: 10^-3)",
    )
    parser.add_argument(
        "--genes-only",
        action="store_true",
        help="Submit only the first column of each .rnk line (gene symbols); default sends full lines (gene + rank score).",
    )
    parser.add_argument(
        "--no-fast-mode",
        action="store_true",
        help="Omit fast mode (slower GOrilla run). Default: fast mode ON (equivalent to the site's checkbox).",
    )
    args = parser.parse_args()

    input_dir: Path = args.input_dir.resolve()
    out_dir: Path = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rnk_files = sorted(input_dir.glob("*.rnk"))
    if not rnk_files:
        raise SystemExit(f"No .rnk files under {input_dir}")

    species_value = SPECIES_CHOICES[args.species]
    fast_mode = not args.no_fast_mode

    session = requests.Session()
    session.headers.update(
        {"User-Agent": "phospho_mhg_anno/test_updown_split/run_gorilla.py (requests)"}
    )

    print(f"Input: {input_dir} ({len(rnk_files)} .rnk files)")
    print(f"Output PNGs: {out_dir}")
    print(
        f"Species: {args.species} ({species_value}), ontology: Biological process (proc), fast_mode={fast_mode}"
    )

    for path in tqdm(rnk_files, desc="GOrilla", unit="rnk"):
        ranked = read_ranked_list(path, genes_only=args.genes_only)
        n_lines = ranked.count("\n") + (1 if ranked else 0)
        stem = path.stem
        dest = out_dir / f"{stem}.png"
        try:
            results, final_url = submit_gorilla(
                session,
                ranked_text=ranked,
                species_value=species_value,
                pvalue_thresh=args.pvalue,
                fast_mode=fast_mode,
            )
            title_m = re.search(r"<title>\s*([^<]+?)\s*</title>", results.text, re.I | re.S)
            title = title_m.group(1).strip() if title_m else "?"

            if "No GO Enrichment" in results.text or title.lower().startswith("no go"):
                print(f"  {stem}: no enrichment — skipped PNG ({n_lines} lines submitted)")
                continue

            png_bytes = download_go_png(session, final_url)
            dest.write_bytes(png_bytes)
            print(f"  {stem}: OK — {dest.name} ({len(png_bytes)} bytes, {n_lines} lines)")
        except (TimeoutError, requests.RequestException, FileNotFoundError, RuntimeError, OSError) as exc:
            print(f"  {stem}: FAILED — {exc}")


if __name__ == "__main__":
    main()
