#!/usr/bin/env python3
"""
Download PhosphoSitePlus-style archives from the Omnipath rescued mirror into resources/.

Source: https://rescued.omnipathdb.org/phosphosite/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.request import Request, urlopen

try:
    from tqdm import tqdm
except ImportError:
    print("Install tqdm: pip install tqdm", file=sys.stderr)
    raise

BASE = "https://rescued.omnipathdb.org/phosphosite/"

DEFAULT_FILES = (
    "Kinase_Substrate_Dataset.gz",
    "Phosphorylation_site_dataset.gz",
)


def download(url: str, dest: Path, chunk: int = 1 << 20) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "phospho_mhg_anno/download_phosphosite_resources.py"})
    with urlopen(req, timeout=120) as resp:
        total = resp.headers.get("Content-Length")
        total_n = int(total) if total and total.isdigit() else None
        with open(dest, "wb") as out, tqdm(
            total=total_n,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=dest.name,
        ) as bar:
            while True:
                block = resp.read(chunk)
                if not block:
                    break
                out.write(block)
                bar.update(len(block))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "resources" / "phosphosite",
        help="Directory for downloaded .gz files.",
    )
    parser.add_argument(
        "--files",
        nargs="*",
        default=list(DEFAULT_FILES),
        help=f"Filenames under phosphosite/ (default: {DEFAULT_FILES}).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the file already exists with non-zero size.",
    )
    args = parser.parse_args()

    out_dir: Path = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for name in args.files:
        url = BASE + name
        dest = out_dir / name
        if dest.is_file() and dest.stat().st_size > 0 and not args.force:
            print(f"skip (exists): {dest}")
            continue
        print(f"fetch {url}")
        download(url, dest)

    readme = out_dir / "README.txt"
    readme.write_text(
        "Downloaded from https://rescued.omnipathdb.org/phosphosite/\n"
        "PhosphoSitePlus data are subject to the PhosphoSitePlus license; cite accordingly.\n"
        "Use with: python analyze/build_kinase_gmt.py --input resources/phosphosite/Kinase_Substrate_Dataset.gz ...\n",
        encoding="utf-8",
    )
    print(f"Done. Files in {out_dir}")


if __name__ == "__main__":
    main()
