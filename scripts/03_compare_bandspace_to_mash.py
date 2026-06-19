#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import read_tsv, write_tsv
from legacy_marker_rescue.distance import pairwise_band_distances, pearson, spearman


def key(a, b):
    return tuple(sorted([a, b]))


def main():
    ap = argparse.ArgumentParser(description="Compare RAPD band-space distances with Mash genome distances.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    matrix = read_tsv(outdir / "rapd_band_matrix.tsv")
    rapd_pairs = pairwise_band_distances(matrix)
    mash_rows = read_tsv(cfg["genomes"]["mash_dist"])
    mash = {}
    for r in mash_rows:
        a = Path(r.get("query", r.get("assembly_a", ""))).name
        b = Path(r.get("reference", r.get("assembly_b", ""))).name
        d = r.get("distance", r.get("mash_distance", ""))
        if a and b and d != "":
            mash[key(a, b)] = float(d)
    joined = []
    for r in rapd_pairs:
        k = key(r["assembly_a"], r["assembly_b"])
        if k in mash:
            row = dict(r)
            row["mash_distance"] = mash[k]
            joined.append(row)
    x = [float(r["rapd_distance"]) for r in joined]
    y = [float(r["mash_distance"]) for r in joined]
    summary = [{"n_pairs": len(joined), "pearson_r": pearson(x, y), "spearman_r": spearman(x, y)}]
    write_tsv(outdir / "rapd_vs_mash_pairwise.tsv", joined)
    write_tsv(outdir / "rapd_vs_mash_summary.tsv", summary)

if __name__ == "__main__":
    main()
