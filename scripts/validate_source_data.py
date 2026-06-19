#!/usr/bin/env python
"""Validate manuscript source-data summaries.

The script recomputes headline numerical summaries from the multisheet
source-data workbook rather than copying values from summary cells. It is meant
as a compact public check that the archived display-item data support the
reported manuscript numbers.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import mean, pstdev

try:
    from openpyxl import load_workbook
except ImportError as exc:  # pragma: no cover
    raise SystemExit("openpyxl is required. Install with: pip install openpyxl") from exc


def as_float(x):
    if x is None or x == "":
        return None
    return float(x)


def read_sheet(path: Path, sheet: str) -> list[dict]:
    wb = load_workbook(path, read_only=True, data_only=True)
    if sheet not in wb.sheetnames:
        raise KeyError(f"Sheet not found: {sheet}")
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    if not rows:
        return []
    headers = [str(h).strip() if h is not None else f"Column {i+1}" for i, h in enumerate(rows[0])]
    out = []
    for row in rows[1:]:
        if row is None or not any(v is not None for v in row):
            continue
        out.append({headers[i]: row[i] if i < len(row) else None for i in range(len(headers))})
    return out


def pearson(x, y):
    pairs = [(float(a), float(b)) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 2:
        return float("nan")
    xs, ys = zip(*pairs)
    mx, my = mean(xs), mean(ys)
    num = sum((a - mx) * (b - my) for a, b in pairs)
    denx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    deny = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (denx * deny) if denx and deny else float("nan")


def ranks(values):
    indexed = sorted((float(v), i) for i, v in enumerate(values))
    out = [0.0] * len(indexed)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][0] == indexed[i][0]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            out[indexed[k][1]] = rank
        i = j
    return out


def spearman(x, y):
    pairs = [(float(a), float(b)) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 2:
        return float("nan")
    xs, ys = zip(*pairs)
    return pearson(ranks(xs), ranks(ys))


def near(name, observed, expected, tol=1e-3):
    if abs(float(observed) - float(expected)) > tol:
        raise AssertionError(f"{name}: observed {observed!r}, expected {expected!r}")
    return (name, observed)


def count_where(rows, column, value):
    return sum(1 for r in rows if str(r.get(column, "")).strip() == value)


def main():
    ap = argparse.ArgumentParser(description="Recompute manuscript-level source-data summaries from Supplementary Data S1.")
    ap.add_argument("--workbook", default="source_data/Supplementary_Data_S1.xlsx", help="Path to the source-data workbook.")
    ap.add_argument("--strict", action="store_true", help="Fail if any check differs beyond the tolerance.")
    args = ap.parse_args()
    wb_path = Path(args.workbook)
    if not wb_path.exists():
        raise FileNotFoundError(wb_path)

    checks = []

    pairs = read_sheet(wb_path, "WGS distance pairs")
    rapd = [as_float(r["RAPD band distance"]) for r in pairs]
    mash = [as_float(r["Mash distance"]) for r in pairs]
    checks.append(near("WGS pairwise comparisons", len(pairs), 253, 0))
    checks.append(near("RAPD-vs-Mash Pearson r", pearson(rapd, mash), 0.755981579, 5e-6))
    checks.append(near("RAPD-vs-Mash Spearman r", spearman(rapd, mash), 0.63628683, 5e-6))

    loo = read_sheet(wb_path, "WGS leave-one-out")
    min_sp = min(as_float(r["Spearman r"]) for r in loo)
    checks.append(near("Minimum leave-one-assembly-out Spearman r", min_sp, 0.521533, 5e-6))

    ctx = read_sheet(wb_path, "Genome context summary")
    def metric_row(metric):
        for r in ctx:
            if r.get("Grouping variable") == "overall" and r.get("Group value") == "all" and r.get("Metric") == metric:
                return r
        raise KeyError(metric)
    gene = metric_row("Fraction overlapping annotated genes")
    contig10 = metric_row("Fraction within 10 kb of contig end")
    gc = metric_row("Mean GC fraction")
    checks.append(near("Gene overlap observed", gene["Observed value"], 0.7690783, 5e-6))
    checks.append(near("Gene overlap matched-random mean", gene["Random mean"], 0.65776016, 5e-6))
    checks.append(near("Contig-end 10 kb observed", contig10["Observed value"], 0.36555024, 5e-6))
    checks.append(near("Mean GC observed", gc["Observed value"], 0.55227209, 5e-6))

    rsets = read_sheet(wb_path, "Random primer sets")
    random_sets = [r for r in rsets if str(r.get("Set type")) == "GC-matched random"]
    historical = [r for r in rsets if str(r.get("Set type")) == "Historical"]
    checks.append(near("Accepted random primer sets", len(random_sets), 25, 0))
    checks.append(near("Historical total amplicons", historical[0]["Total amplicons"], 14630, 0))
    checks.append(near("Historical non-duplicate Spearman r", historical[0]["Spearman r (non-duplicate pairs)"], 0.2954562386, 5e-6))
    checks.append(near("Historical fixed-core band fraction", historical[0]["Fixed-core band fraction"], 0.8723404255, 5e-6))
    checks.append(near("Historical polymorphic band fraction", historical[0]["Polymorphic band fraction"], 0.04255319149, 5e-6))

    gsum = read_sheet(wb_path, "Guthrie summary agreement")
    global_rows = [r for r in gsum if r.get("Denominator") == "Global denominator"]
    diffs = [as_float(r["Digitized shared-locus (%)"]) - as_float(r["Published shared-locus (%)"]) for r in global_rows]
    rmse = math.sqrt(mean(d * d for d in diffs))
    loci = max(int(r["Composite loci"]) for r in global_rows)
    checks.append(near("Guthrie composite loci", loci, 31, 0))
    checks.append(near("Guthrie published-summary RMSE", rmse, 2.492207, 5e-4))

    btol = read_sheet(wb_path, "Guthrie bridge tolerance")
    b005 = min(btol, key=lambda r: abs(as_float(r["Band-size tolerance (kb)"]) - 0.05))
    checks.append(near("Guthrie legacy bands", b005["Legacy bands"], 31, 0))
    checks.append(near("Guthrie bands with any modern bridge", b005["Legacy bands with modern match"], 27, 0))
    checks.append(near("Guthrie any-bridge fraction", b005["Legacy match fraction"], 0.870968, 5e-6))

    bclass = read_sheet(wb_path, "Guthrie bridge classes")
    checks.append(near("Guthrie polymorphic bridge bins", count_where(bclass, "Bridge class", "Polymorphic bridge"), 6, 0))
    checks.append(near("Guthrie core-only bridge bins", count_where(bclass, "Bridge class", "Core-only bridge"), 21, 0))
    checks.append(near("Guthrie no-modern-match bins", count_where(bclass, "Bridge class", "No modern match"), 4, 0))

    den_segs = read_sheet(wb_path, "Denoyes detected segments")
    opc02 = [r for r in den_segs if r.get("Panel") == "OPC02"]
    opa13 = [r for r in den_segs if r.get("Panel") == "OPA13"]
    checks.append(near("Denoyes detected segments", len(den_segs), 204, 0))
    checks.append(near("Denoyes OPC02 detected segments", len(opc02), 87, 0))
    checks.append(near("Denoyes OPA13 detected segments", len(opa13), 117, 0))

    lane = read_sheet(wb_path, "Denoyes lane audit")
    checks.append(near("Denoyes retained nonzero sample lanes", count_where(lane, "Lane status", "Retained sample lane"), 33, 0))
    checks.append(near("Denoyes manual-check sample lanes", count_where(lane, "Lane status", "Manual-check sample lane"), 5, 0))

    print("Source-data validation passed")
    for name, value in checks:
        if isinstance(value, float):
            print(f"- {name}: {value:.6g}")
        else:
            print(f"- {name}: {value}")


if __name__ == "__main__":
    main()
