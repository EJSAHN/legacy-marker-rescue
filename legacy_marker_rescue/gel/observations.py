"""Structural annotation checks and one-to-one ordered point matching.

Coordinates are native embedded-image pixels, with the origin at the top left.
Reference annotations must be entered by a human; this module does not generate
reference lanes, bands, labels, or biological genotype calls.
"""
from __future__ import annotations
import math
from typing import Any
PANELS = ('A_OPC02', 'B_OPA13')
ROLES = {'sample', 'marker', 'gap', 'unresolved'}
LABEL_SOURCES = {'printed_here', 'shared_header', 'unresolved'}

def finite(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} must be numeric')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f'{name} is not finite')
    return value

def inside(x: float, y: float, box: list[float]) -> bool:
    return box[0] <= x < box[2] and box[1] <= y < box[3]

def overlap(a: list[float], b: list[float]) -> bool:
    return min(a[2], b[2]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[1], b[1])

def validate_box(box: Any, bounds: list[float], name: str) -> list[float]:
    if not isinstance(box, list) or len(box) != 4:
        raise ValueError(f'{name} must have four coordinates')
    box = [finite(x, name) for x in box]
    if box[2] <= box[0] or box[3] <= box[1]:
        raise ValueError(f'{name} has zero or reversed extent')
    if box[0] < bounds[0] or box[1] < bounds[1] or box[2] > bounds[2] or (box[3] > bounds[3]):
        raise ValueError(f'{name} lies outside its panel')
    return box

def validate_review(review: dict, protocol: dict, protocol_sha256: str, require_final: bool=True) -> dict:
    """Reject malformed or incomplete input rather than manufacture accuracy."""
    if not isinstance(review, dict) or review.get('schema_version') != '1.0':
        raise ValueError('Unknown review schema')
    if review.get('dataset_id') != protocol['dataset_id']:
        raise ValueError('Wrong dataset')
    if review.get('protocol_sha256') != protocol_sha256:
        raise ValueError('Review protocol hash does not match this package')
    if review.get('image_sha256') != protocol['image_sha256']:
        raise ValueError('Image hash does not match this package')
    if review.get('width_px') != protocol['width_px'] or review.get('height_px') != protocol['height_px']:
        raise ValueError('Wrong image dimensions or coordinate frame')
    if not isinstance(review.get('session_id'), str) or not review['session_id']:
        raise ValueError('Review session identifier is missing')
    if require_final and review.get('status') != 'complete':
        raise ValueError('Review is a draft; no accuracy score was computed')
    person = review.get('reviewer', {})
    if require_final:
        if not str(person.get('name', '')).strip():
            raise ValueError('Reviewer name is missing')
        if person.get('prior_algorithm_exposure') not in {'yes', 'no', 'uncertain'}:
            raise ValueError('Declare prior exposure to algorithm outputs')
        if person.get('annotation_method') != 'manual_image_reading':
            raise ValueError('This evaluation requires manual image annotations')
        for key in ('direct_image_reading', 'no_prediction_overlay_used', 'coverage_checked'):
            if review.get('attestations', {}).get(key) is not True:
                raise ValueError(f'Missing completion declaration: {key}')
        if not review.get('completed_at'):
            raise ValueError('Completion time is missing')
    panels = review.get('panels')
    if not isinstance(panels, dict) or set(panels) != set(PANELS):
        raise ValueError('Both named panels must be present')
    ids: set[str] = set()
    counts = {'lanes': 0, 'sample_lanes': 0, 'readable_sample_lanes': 0, 'definite_bands': 0, 'uncertain_bands': 0, 'exclusion_regions': 0}
    for pid in PANELS:
        panel = panels[pid]
        if require_final and panel.get('coverage_complete') is not True:
            raise ValueError(f'{pid}: panel coverage is not complete')
        lanes = panel.get('lanes', [])
        bands = panel.get('bands', [])
        regions = panel.get('exclusions', [])
        if not all((isinstance(x, list) for x in (lanes, bands, regions))):
            raise ValueError(f'{pid}: invalid annotation lists')
        if require_final and (not lanes):
            raise ValueError(f'{pid}: no lanes were annotated')
        lm = {}
        for lane in lanes:
            lid = lane.get('id')
            if not isinstance(lid, str) or not lid or lid in ids:
                raise ValueError('Missing or duplicate annotation identifier')
            ids.add(lid)
            box = validate_box(lane.get('box'), protocol['panel_bounds'][pid], 'Lane rectangle')
            if lane.get('role') not in ROLES:
                raise ValueError('Unknown lane role')
            if lane.get('assessment') not in {'readable', 'unreadable'}:
                raise ValueError('Unknown lane assessment')
            if lane.get('label_source') not in LABEL_SOURCES:
                raise ValueError('Unknown lane label source')
            if lane.get('label_source') != 'unresolved' and (not str(lane.get('label', '')).strip()):
                raise ValueError('Enter the printed label or set label source to unresolved')
            if require_final and lane.get('reviewed') is not True:
                raise ValueError(f'{pid} {lid}: lane has not been marked reviewed')
            if require_final and (lane['assessment'] == 'unreadable' or lane['role'] == 'unresolved') and (not str(lane.get('note', '')).strip()):
                raise ValueError('Explain unreadable or unresolved lanes in the lane note')
            for other in lm.values():
                if overlap(box, other['box']):
                    raise ValueError(f'{pid}: lane rectangles overlap; redraw their boundaries')
            lm[lid] = lane
            counts['lanes'] += 1
            if lane['role'] == 'sample':
                counts['sample_lanes'] += 1
                if lane['assessment'] == 'readable':
                    counts['readable_sample_lanes'] += 1
        for band in bands:
            bid = band.get('id')
            if not isinstance(bid, str) or not bid or bid in ids:
                raise ValueError('Missing or duplicate annotation identifier')
            ids.add(bid)
            lane = lm.get(band.get('lane_id'))
            if lane is None:
                raise ValueError('Band refers to a missing lane')
            x, y = (finite(band.get('x'), 'Band x'), finite(band.get('y'), 'Band y'))
            if not inside(x, y, lane['box']):
                raise ValueError('Band center is outside its assigned lane')
            if band.get('certainty') not in {'definite', 'uncertain'}:
                raise ValueError('Band certainty is not specified')
            if lane['assessment'] == 'unreadable' and band['certainty'] == 'definite':
                raise ValueError('Do not place definite bands in an unreadable lane')
            if lane['role'] in {'gap', 'unresolved'} and band['certainty'] == 'definite':
                raise ValueError('Resolve the lane role before entering definite bands')
            counts['definite_bands' if band['certainty'] == 'definite' else 'uncertain_bands'] += 1
        for region in regions:
            rid = region.get('id')
            if not isinstance(rid, str) or not rid or rid in ids:
                raise ValueError('Missing or duplicate annotation identifier')
            ids.add(rid)
            validate_box(region.get('box'), protocol['panel_bounds'][pid], 'Exclusion rectangle')
            if not str(region.get('reason', '')).strip():
                raise ValueError('Exclusion regions require a reason')
            counts['exclusion_regions'] += 1
        for lane in lanes:
            lb = [b for b in bands if b['lane_id'] == lane['id']]
            definite = [b for b in lb if b['certainty'] == 'definite']
            uncertain = [b for b in lb if b['certainty'] == 'uncertain']
            if require_final and lane['role'] == 'sample' and (lane['assessment'] == 'readable') and (not lb) and (lane.get('no_visible_bands_confirmed') is not True):
                raise ValueError('An empty readable sample lane requires an explicit no-visible-bands confirmation')
            if lane.get('no_visible_bands_confirmed') is True and lb:
                raise ValueError('An empty-lane confirmation contradicts existing band marks')
            positions = [(b['x'], b['y']) for b in lb]
            if len(positions) != len(set(positions)):
                raise ValueError('Duplicate reference band centers in one lane')
            for b in definite:
                if any((inside(b['x'], b['y'], r['box']) for r in regions)):
                    raise ValueError('A definite band lies in an excluded region')
                if any((abs(b['y'] - u['y']) <= protocol['primary_y_tolerance_px'] for u in uncertain)):
                    raise ValueError('Definite and uncertain band neighborhoods overlap')
    if require_final and (not counts['readable_sample_lanes']):
        raise ValueError('No readable sample lanes are available')
    return counts

def match_1d(reference_y: list[float], prediction_y: list[float], tolerance: float) -> list[tuple[int, int, float]]:
    """Maximum-cardinality, minimum-L1 matching of ordered 1-D points.

The dynamic program uses the noncrossing optimal alignment of sorted points.
Ties are resolved deterministically. Returned indices refer to original lists.
"""
    tolerance = finite(tolerance, 'Tolerance')
    if tolerance < 0:
        raise ValueError('Tolerance must be nonnegative')
    ref = sorted(((finite(v, 'Reference position'), i) for i, v in enumerate(reference_y)))
    pred = sorted(((finite(v, 'Prediction position'), i) for i, v in enumerate(prediction_y)))
    n, m = (len(ref), len(pred))
    counts = [[0] * (m + 1) for _ in range(n + 1)]
    costs = [[0.0] * (m + 1) for _ in range(n + 1)]
    actions = [[''] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            choices = [(counts[i - 1][j], costs[i - 1][j], 'r'), (counts[i][j - 1], costs[i][j - 1], 'p')]
            delta = abs(ref[i - 1][0] - pred[j - 1][0])
            if delta <= tolerance:
                choices.append((counts[i - 1][j - 1] + 1, costs[i - 1][j - 1] + delta, 'm'))
            best = min(choices, key=lambda z: (-z[0], z[1], {'m': 0, 'r': 1, 'p': 2}[z[2]]))
            counts[i][j], costs[i][j], actions[i][j] = best
    out = []
    i, j = (n, m)
    while i and j:
        action = actions[i][j]
        if action == 'm':
            out.append((ref[i - 1][1], pred[j - 1][1], abs(ref[i - 1][0] - pred[j - 1][0])))
            i -= 1
            j -= 1
        elif action == 'r':
            i -= 1
        else:
            j -= 1
    return list(reversed(out))
