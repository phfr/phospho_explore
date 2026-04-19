#!/usr/bin/env python3
"""
Run the full analysis pipeline from the repo root.

Order: rank lists → GO prerank → GO heatmap (all terms) → temporal GO + figures →
non-absent heatmap + dotplot (your sizing) → leading edge (top 5 per condition) →
leading-edge follow-up (phospho / overlap / Enrichr / kinase Fisher) →
kinase prerank (if a kinase GMT exists) → multiphos switches (if data.tsv exists).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def run(cmd: list[str], title: str) -> None:
    print(f"\n=== {title} ===", flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    py = sys.executable

    run([py, "generate_rank_lists.py"], "Rank lists (.rnk from data.tsv)")
    run([py, "analyze/run_go_prerank_batch.py"], "GO prerank (gseapy, all conditions)")
    run(
        [
            py,
            "analyze/cross_condition_go_heatmap.py",
            "--out-matrix",
            "output/analyze/go_crosscondition_matrix.tsv",
        ],
        "GO cross-condition heatmap (default row filter)",
    )
    run(
        [
            py,
            "analyze/go_temporal_trends.py",
            "--plot",
            "--plot-scatter",
            "--plot-scatter-facet",
        ],
        "GO temporal classification + bar / scatter / facet figures",
    )
    run(
        [
            py,
            "analyze/cross_condition_go_heatmap.py",
            "--terms-tsv",
            "output/analyze/go_temporal_classification_nonabsent.tsv",
            "--out-png",
            "output/analyze/figures/go_nonabsent_crosscondition_heatmap.png",
            "--out-matrix",
            "output/analyze/go_nonabsent_crosscondition_matrix.tsv",
            "--dotplot",
            "--dotplot-figheight",
            "28",
            "--figwidth",
            "6",
        ],
        "GO non-absent cross-condition heatmap + dotplot",
    )
    run([py, "analyze/leading_edge_export.py", "--top-n", "5"], "Leading edge export (top 5 per condition)")
    run([py, "analyze/leading_edge_followup_all.py"], "Leading-edge follow-up (phospho, overlap, Enrichr, kinase)")

    gmt_full = ROOT / "data/gene_sets/kinase_substrates.gmt"
    gmt_sample = ROOT / "data/gene_sets/kinase_substrates_sample.gmt"
    if gmt_full.is_file() or gmt_sample.is_file():
        run([py, "analyze/kinase_substrate_prerank.py"], "Kinase–substrate prerank")
    else:
        print("\n=== Skipping kinase prerank (no kinase GMT under data/gene_sets/) ===", flush=True)

    if (ROOT / "data.tsv").is_file():
        run([py, "analyze/multiphos_switches.py"], "Multi-phos switches")
    else:
        print("\n=== Skipping multiphos_switches (no data.tsv) ===", flush=True)

    print("\n=== Done ===", flush=True)


if __name__ == "__main__":
    main()
