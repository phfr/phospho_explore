# Analyze scripts

Pipeline order:

1. **Rank lists** (repo root): `python generate_rank_lists.py` — writes `output/*.rnk` (gene + `RankScore`; gseapy needs two columns).
2. **GO prerank (all downstream GO analyses)**: `python analyze/run_go_prerank_batch.py`  
   - Preferred: put MSigDB **`data/gene_sets/c5.go.bp.*.Hs.symbols.gmt`** (any version) — auto-detected — or pass `--gmt-path`.  
   - If no such file is found, defaults to Enrichr `GO_Biological_Process_2021` (downloads; needs network).  
   - Outputs: `output/analyze/go_bp_prerank_long.tsv`, `output/analyze/go_prerank/<condition>/`.
3. **Cross-condition heatmap**: `python analyze/cross_condition_go_heatmap.py`
4. **Temporal GO patterns** (5 min vs 10 min per arm): `python analyze/go_temporal_trends.py [--plot]`
5. **Leading edge → STRING lists**: `python analyze/leading_edge_export.py --go-ids GO:0007411` or `--top-n 5`
6. **Leading-edge follow-up** (after step 5): `python analyze/leading_edge_followup_all.py`  
   - Phospho **clustermap** + **volcano facets** for union genes (`leading_edge_followup/`).  
   - **Overlap**: membership TSV, pairwise Jaccard, top intersection bar chart.  
   - **Enrichr** (union genes, needs network): `enrichr/enrichr_combined_union.tsv` — tune `--libraries` if rate-limited.  
   - **Kinase overlap**: Fisher vs kinase GMT substrates → `leading_edge_kinase_fisher.tsv`.
7. **PhosphoSite data (optional):** `python scripts/download_phosphosite_resources.py` → `resources/phosphosite/*.gz` ([mirror](https://rescued.omnipathdb.org/phosphosite/)).
8. **Kinase / substrate prerank** (gene-level KSEA-like):  
   `python analyze/build_kinase_gmt.py --input resources/phosphosite/Kinase_Substrate_Dataset.gz --psp-kinase-substrate --human-only`  
   then `python analyze/kinase_substrate_prerank.py` (or use `data/gene_sets/kinase_substrates_sample.gmt` for a smoke test).  
   Exploratory notebook: `notebooks/kinase_comprehensive_analysis.ipynb`.
9. **Multi-phos switches** (site-level): `python analyze/multiphos_switches.py`

Install dependencies:

```bash
pip install -r requirements-analyze.txt
```

Notes:

- **Ranked lists:** With the usual **two-column** `.rnk` (gene + `RankScore`), gseapy uses **that score column** only — never positional ranks. If you pass a **gene-only** one-column list, `go_utils.prepare_rnk_for_gseapy` falls back to descending ranks and prints a warning.
- **STRING** URLs in the leading-edge manifest may be truncated for very long gene lists; full lists are always in the `.txt` files.
