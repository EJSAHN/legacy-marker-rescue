#!/usr/bin/env python
from __future__ import annotations
import argparse, random
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import write_tsv
from legacy_marker_rescue.random_primer import random_primer_set


def main():
    ap = argparse.ArgumentParser(description="Generate GC-matched random RAPD primer sets for exploratory analyses.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    n = int(cfg.get("random_primer_null", {}).get("n_sets", 100))
    seed = int(cfg.get("random_primer_null", {}).get("seed", 1))
    primers = cfg["rapd"]["primers"]
    rng = random.Random(seed)
    rows = []
    for i in range(1, n+1):
        ps = random_primer_set(primers, rng)
        row = {"set_id": f"random_set_{i:04d}"}
        row.update(ps)
        rows.append(row)
    write_tsv(outdir / "random_primer_sets.tsv", rows)

if __name__ == "__main__":
    main()
