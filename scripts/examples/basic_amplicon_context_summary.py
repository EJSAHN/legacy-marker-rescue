#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import read_tsv, write_tsv
from legacy_marker_rescue.sequence import gc_fraction


def main():
    ap = argparse.ArgumentParser(description="Summarize basic sequence context for predicted RAPD amplicons.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    amplicons = read_tsv(outdir / "rapd_amplicons.tsv")
    # This example summarizes predicted amplicon sizes. Manuscript-level context summaries are validated from the source-data workbook.
    sizes = [int(r["amplicon_size_bp"]) for r in amplicons]
    summary = [{
        "n_amplicons": len(sizes),
        "mean_amplicon_size_bp": round(sum(sizes) / len(sizes), 3) if sizes else 0,
        "min_amplicon_size_bp": min(sizes) if sizes else 0,
        "max_amplicon_size_bp": max(sizes) if sizes else 0,
    }]
    write_tsv(outdir / "rapd_context_summary.tsv", summary)

if __name__ == "__main__":
    main()
