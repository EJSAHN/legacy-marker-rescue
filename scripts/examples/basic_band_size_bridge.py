#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import read_tsv, write_tsv
from legacy_marker_rescue.bridge import match_legacy_to_modern


def main():
    ap = argparse.ArgumentParser(description="Bridge legacy band-size bins to modern in-silico RAPD band-size bins.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(cfg["project"]["output_dir"])
    legacy = [{"legacy_band": r["segment_id"], "kb": r["kb"]} for r in read_tsv(outdir / "legacy_detected_segments.tsv")]
    modern = read_tsv(outdir / "rapd_band_summary.tsv")
    modern2 = []
    for r in modern:
        modern2.append({"band": r["band"], "bp": r["mean_size_bp"]})
    rows = match_legacy_to_modern(legacy, modern2, float(cfg["bridge"]["tolerance_kb"]))
    write_tsv(outdir / "legacy_modern_bridge.tsv", rows)

if __name__ == "__main__":
    main()
