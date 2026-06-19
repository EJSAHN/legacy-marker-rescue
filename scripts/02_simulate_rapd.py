#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import read_tsv, write_tsv, ensure_dir
from legacy_marker_rescue.rapd import simulate_rapd_for_fasta, build_band_matrix


def main():
    ap = argparse.ArgumentParser(description="Simulate in-silico RAPD from genome FASTA files.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = ensure_dir(cfg["project"]["output_dir"])
    inv = read_tsv(cfg["genomes"]["fasta_inventory"])
    primers = cfg["rapd"]["primers"]
    mm = int(cfg["rapd"]["max_mismatch_per_site"])
    max_bp = int(cfg["rapd"]["max_amplicon_bp"])
    tol = float(cfg["rapd"]["band_relative_tolerance"])
    rows = []
    for r in inv:
        rows.extend(simulate_rapd_for_fasta(r["fasta_path"], primers, mm, max_bp))
    matrix, bands = build_band_matrix(rows, tol)
    write_tsv(Path(outdir) / "rapd_amplicons.tsv", rows)
    write_tsv(Path(outdir) / "rapd_band_matrix.tsv", matrix)
    write_tsv(Path(outdir) / "rapd_band_summary.tsv", bands)

if __name__ == "__main__":
    main()
