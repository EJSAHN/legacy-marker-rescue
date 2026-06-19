#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import write_tsv


def main():
    ap = argparse.ArgumentParser(description="Create a source-data manifest for archive outputs.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    rows = []
    for path in sorted(outdir.glob("*.tsv")):
        rows.append({"file": str(path), "role": "derived workflow output"})
    write_tsv(outdir / "source_data_manifest.tsv", rows)

if __name__ == "__main__":
    main()
