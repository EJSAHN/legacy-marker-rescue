#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.io import write_tsv


def main():
    ap = argparse.ArgumentParser(description="Build a clean genome FASTA inventory.")
    ap.add_argument("--genome-dir", required=True, help="Directory containing genome FASTA files or NCBI Datasets output.")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()
    root = Path(args.genome_dir)
    rows = []
    for path in root.rglob("*.fna"):
        name = path.name.lower()
        if "cds_from_genomic" in name:
            continue
        if "genomic" not in name:
            continue
        rows.append({"assembly": path.stem, "fasta_path": str(path), "bytes": path.stat().st_size})
    write_tsv(args.out, rows, ["assembly", "fasta_path", "bytes"])

if __name__ == "__main__":
    main()
