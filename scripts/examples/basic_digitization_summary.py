#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import read_tsv, write_tsv


def main():
    ap = argparse.ArgumentParser(description="Build a basic digitization summary for digitized legacy segments.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    segs = read_tsv(outdir / "legacy_detected_segments.tsv")
    sizes = [float(s["kb"]) for s in segs]
    summary = [{
        "n_segments": len(segs),
        "min_kb": min(sizes) if sizes else "",
        "max_kb": max(sizes) if sizes else "",
    }]
    write_tsv(outdir / "legacy_digitization_qc_summary.tsv", summary)

if __name__ == "__main__":
    main()
