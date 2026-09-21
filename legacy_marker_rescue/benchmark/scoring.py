"""One-to-one object-localization evaluation against published foreground masks."""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import linear_sum_assignment

def overlap_edges(predictions: list[dict], labels: np.ndarray, radius: float) -> list[dict[int, float]]:
    """Candidate-center to foreground-component edges in native coordinates.

    At radius zero, a center must fall in a foreground pixel. Positive radii use
    distance to foreground pixel centers; the containing pixel remains eligible.
    No learned or reference-dependent expansion is used.
    """
    if radius < 0 or not math.isfinite(radius):
        raise ValueError('Invalid tolerance')
    h, w = labels.shape
    rows = []
    for p in predictions:
        x = float(p['x_px'])
        y = float(p['y_px'])
        if not math.isfinite(x + y) or not (0 <= x < w and 0 <= y < h):
            raise ValueError('Prediction outside native image')
        hit = {}
        iy, ix = (int(math.floor(y)), int(math.floor(x)))
        component = int(labels[iy, ix])
        if component:
            hit[component] = 0.0
        if radius:
            l = max(0, int(math.floor(x - radius - 0.5)))
            r = min(w, int(math.ceil(x + radius + 0.5)))
            t = max(0, int(math.floor(y - radius - 0.5)))
            b = min(h, int(math.ceil(y + radius + 0.5)))
            yy, xx = np.mgrid[t:b, l:r]
            distance = np.sqrt((xx + 0.5 - x) ** 2 + (yy + 0.5 - y) ** 2)
            local = labels[t:b, l:r]
            for cid in np.unique(local[(distance <= radius) & (local > 0)]):
                d = float(distance[local == cid].min())
                cid = int(cid)
                hit[cid] = min(hit.get(cid, float('inf')), d)
        rows.append(hit)
    return rows

def assign_components(edges: list[dict[int, float]], n_reference: int) -> list[tuple[int, int, float]]:
    n = len(edges)
    m = n_reference
    if n > 5000 or m > 5000 or n * (m + n) > 25000000:
        raise ValueError('MATCHING_SAFETY_LIMIT')
    if not n or not m:
        return []
    costs = np.full((n, m + n), 2.0, dtype=float)
    costs[:, m:] = 1.0
    max_d = max((v for e in edges for v in e.values()), default=0.0)
    for i, es in enumerate(edges):
        for cid, d in es.items():
            if not 1 <= cid <= m or d < 0 or (not math.isfinite(d)):
                raise ValueError('Invalid matching edge')
            costs[i, cid - 1] = d / (max_d + 1) / (2 * (min(n, m) + 1))
    ii, jj = linear_sum_assignment(costs)
    return [(int(i), int(j) + 1, float(edges[int(i)][int(j) + 1])) for i, j in zip(ii, jj) if j < m and costs[i, j] < 1.0]

def evaluate(predictions: list[dict], labels: np.ndarray, refs: list[dict], radius: float) -> dict:
    matched = assign_components(overlap_edges(predictions, labels, radius), len(refs))
    n = len(predictions)
    m = len(refs)
    k = len(matched)
    return {'tolerance_native_px': radius, 'predictions': n, 'reference_components': m, 'matched': k, 'unmatched_predictions': n - k, 'unmatched_reference_components': m - k, 'precision_vs_components': k / n if n else None, 'recall_vs_components': k / m if m else None, 'f1_vs_components': 2 * k / (n + m) if n + m else None, 'matches': [{'prediction_index': i, 'component_id': c, 'distance_to_foreground_px': d} for i, c, d in matched], 'unmatched_prediction_indices': sorted(set(range(n)) - {i for i, _, _ in matched}), 'unmatched_component_ids': sorted(set(range(1, m + 1)) - {c for _, c, _ in matched})}

def summaries(rows: list[dict]) -> list[dict]:
    """Keep collections/methods separate; no band-level independence assumption."""
    out = []
    groups = {}
    for r in rows:
        groups.setdefault((r['collection'], r['method'], r['tolerance_native_px']), []).append(r)
    for (collection, method, tol), g in sorted(groups.items()):
        p = sum((x['predictions'] for x in g))
        n = sum((x['reference_components'] for x in g))
        k = sum((x['matched'] for x in g))
        vals = [x['f1_vs_components'] for x in g if x['f1_vs_components'] is not None]
        out.append({'collection': collection, 'method': method, 'tolerance_native_px': tol, 'evaluated_images': len(g), 'unique_raw_image_hashes': len({x['raw_image_sha256'] for x in g}), 'predictions': p, 'reference_components': n, 'matched': k, 'pooled_precision_vs_components': k / p if p else None, 'pooled_recall_vs_components': k / n if n else None, 'pooled_f1_vs_components': 2 * k / (p + n) if p + n else None, 'median_image_f1': float(np.median(vals)) if vals else None, 'q25_image_f1': float(np.percentile(vals, 25)) if vals else None, 'q75_image_f1': float(np.percentile(vals, 75)) if vals else None, 'minimum_image_f1': min(vals) if vals else None, 'maximum_image_f1': max(vals) if vals else None, 'confidence_intervals_calculated': False})
    return out
