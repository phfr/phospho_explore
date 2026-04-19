#!/usr/bin/env python3
"""
Run leading-edge follow-up analyses (phospho plots, overlap, Enrichr, kinase Fisher).

Excludes Cytoscape/STRING-only workflows; Enrichr is best-effort if the network fails.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    py = sys.executable

    def run(cmd: list[str], title: str, *, required: bool = True) -> None:
        print(f"\n=== {title} ===", flush=True)
        r = subprocess.run(cmd, cwd=root)
        if r.returncode != 0 and required:
            raise SystemExit(r.returncode)
        if r.returncode != 0:
            print(f"(non-fatal exit {r.returncode})", flush=True)

    run([py, "analyze/leading_edge_phospho_plots.py"], "Phospho heatmap + volcano facets", required=True)
    run([py, "analyze/leading_edge_overlap.py"], "List overlap (Jaccard + top intersections)", required=True)
    run([py, "analyze/leading_edge_enrichr_second.py"], "Enrichr (second layer)", required=False)
    run([py, "analyze/leading_edge_kinase_overlap.py"], "Kinase substrate overlap (Fisher)", required=True)
    print("\n=== leading_edge_followup_all done ===", flush=True)


if __name__ == "__main__":
    main()
