# External resources

## PhosphoSitePlus archives (Omnipath mirror)

Download curated tables from [rescued.omnipathdb.org/phosphosite/](https://rescued.omnipathdb.org/phosphosite/) into `resources/phosphosite/`:

```bash
pip install tqdm
python scripts/download_phosphosite_resources.py
```

Typical files:

- `Kinase_Substrate_Dataset.gz` — kinase → substrate relationships (use with `analyze/build_kinase_gmt.py`).
- `Phosphorylation_site_dataset.gz` — site-level phosphorylation records.

Cite PhosphoSitePlus / CST per the license text inside each file.

`resources/phosphosite/.gitignore` ignores `*.gz` so large archives are not committed by mistake; remove that line if you want them in git.
