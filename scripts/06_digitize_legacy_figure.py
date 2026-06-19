#!/usr/bin/env python
from __future__ import annotations
import argparse
from pathlib import Path
from legacy_marker_rescue.config import load_config
from legacy_marker_rescue.io import write_tsv
from legacy_marker_rescue.digitize import render_pdf_page, crop_image, detect_vertical_segments, x_to_kb


def main():
    ap = argparse.ArgumentParser(description="Digitize a configured legacy banding figure into segment coordinates.")
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    lf = cfg["legacy_figure"]
    outdir = Path(cfg["project"]["output_dir"])
    image = render_pdf_page(lf["pdf"], int(lf["page_index"]), dpi=300)
    crop = crop_image(image, lf["crop_box_pixels"])
    segs = detect_vertical_segments(crop)
    axis = lf["axis_kb"]
    for i, s in enumerate(segs, 1):
        s["segment_id"] = f"segment_{i:05d}"
        s["kb"] = round(x_to_kb(float(s["x_mid"]), axis["x0_pixel"], axis["x25_pixel"], axis.get("min_kb", 0.0), axis.get("max_kb", 2.5)), 4)
    write_tsv(outdir / "legacy_detected_segments.tsv", segs)

if __name__ == "__main__":
    main()
