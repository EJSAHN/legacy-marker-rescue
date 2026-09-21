"""One-to-one agreement for independently completed band annotations."""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from .common import read_tsv, write_tsv, write_json
BAND_COLUMNS = ['image_id', 'image_sha256', 'lane_id', 'band_id', 'x_px', 'y_px', 'reviewer', 'status']
LANE_COLUMNS = ['image_id', 'image_sha256', 'lane_id', 'lane_type', 'review_status', 'reviewer', 'notes']

def match_bands(reference, predicted, tolerance):
    """Maximize valid one-to-one matches, then minimize absolute vertical error."""
    if tolerance <= 0:
        raise ValueError('Tolerance must be positive and fixed before evaluation')
    r = np.asarray(reference, dtype=float)
    p = np.asarray(predicted, dtype=float)
    if not np.isfinite(r).all() or not np.isfinite(p).all():
        raise ValueError('Nonfinite band coordinates')
    if not len(r) or not len(p):
        return []
    cost = np.abs(r[:, None] - p[None, :])
    valid = cost <= tolerance
    penalized = np.where(valid, cost, (min(len(r), len(p)) + 1) * (tolerance + 1))
    ri, pi = linear_sum_assignment(penalized)
    return [(int(a), int(b), float(cost[a, b])) for a, b in zip(ri, pi) if valid[a, b]]

def evaluate(reference_lanes, reference_bands, predictions, tolerance, out):
    lanes = read_tsv(reference_lanes)
    refs = read_tsv(reference_bands)
    preds = read_tsv(predictions)
    if not lanes:
        raise ValueError('Reference lanes are empty; independent annotation is required')
    keys = {}
    for lane in lanes:
        key = (lane['image_id'], lane['lane_id'])
        if key in keys:
            raise ValueError('Duplicate reference lane')
        if not lane.get('reviewer') or len(lane.get('image_sha256', '')) != 64:
            raise ValueError('Reference provenance is incomplete')
        keys[key] = lane
    for group in (refs, preds):
        seen = set()
        for row in group:
            key = (row['image_id'], row['lane_id'])
            uid = (row['image_id'], row['band_id'])
            if uid in seen:
                raise ValueError('Duplicate band identifier')
            seen.add(uid)
            if key in keys and row['image_sha256'] != keys[key]['image_sha256']:
                raise ValueError('Image coordinate frames differ')
    rows = []
    pairs = []
    tp = fp = fn = 0
    for key, lane in keys.items():
        if lane['lane_type'] != 'sample' or lane['review_status'] != 'complete':
            continue
        rr = [r for r in refs if (r['image_id'], r['lane_id']) == key]
        if any((r['status'] != 'confirmed' for r in rr)):
            raise ValueError('Uncertain bands in a complete lane; mark lane unreadable or adjudicate')
        pp = [r for r in preds if (r['image_id'], r['lane_id']) == key]
        matches = match_bands([float(r['y_px']) for r in rr], [float(p['y_px']) for p in pp], tolerance)
        a = len(matches)
        b = len(pp) - a
        c = len(rr) - a
        tp += a
        fp += b
        fn += c
        rows.append({'image_id': key[0], 'lane_id': key[1], 'TP': a, 'FP': b, 'FN': c, 'precision': a / (a + b) if a + b else np.nan, 'recall': a / (a + c) if a + c else np.nan, 'f1': 2 * a / (2 * a + b + c) if 2 * a + b + c else np.nan})
        for i, j, error in matches:
            pairs.append({'image_id': key[0], 'lane_id': key[1], 'reference_band': rr[i]['band_id'], 'predicted_band': pp[j]['band_id'], 'absolute_vertical_error_px': error})
    if not rows:
        raise ValueError('No independently reviewed complete sample lanes')
    report = {'scored_sample_lanes': len(rows), 'unscored_reference_lanes': len(lanes) - len(rows), 'TP': tp, 'FP': fp, 'FN': fn, 'precision': tp / (tp + fp) if tp + fp else np.nan, 'recall': tp / (tp + fn) if tp + fn else np.nan, 'f1': 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else np.nan, 'tolerance_px': tolerance, 'scope': 'image annotation agreement on fully reviewed lanes; no genotype-level accuracy claim'}
    write_tsv(out / 'per_lane_accuracy.tsv', rows)
    write_tsv(out / 'band_matches.tsv', pairs)
    write_json(out / 'accuracy_summary.json', report)
    return report
