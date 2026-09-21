"""Describe candidate agreement with approximate, non-exhaustive reader marks.

Unmatched predictions are review candidates, not demonstrated false positives.
The reader's marks are never changed, completed, or used to fit the detector.
"""
from __future__ import annotations
from collections import Counter
import math
from .observations import match_1d, inside

def canonical_legacy(rows: list[dict]) -> list[dict]:
    return [{'candidate_id': r['prediction_id'], 'panel_id': r['panel_id'], 'x_px': float(r['x_px']), 'y_px': float(r['y_px']), 'label': r['archived_label'], 'evidence_tier': 'legacy_unclassified'} for r in rows]

def compare_snapshot(review: dict, candidates: list[dict], tolerance: float, mask_mode: str='recorded_masks', uncertain_radius: float=3.0) -> dict:
    if mask_mode not in {'recorded_masks', 'no_masks_diagnostic'}:
        raise ValueError('Unknown mask mode')
    if tolerance < 0 or not math.isfinite(tolerance) or uncertain_radius < 0 or (not math.isfinite(uncertain_radius)):
        raise ValueError('Invalid localization tolerance')
    ids = set()
    for row in candidates:
        cid = row['candidate_id']
        if not cid or cid in ids:
            raise ValueError('Candidate identifiers must be unique')
        ids.add(cid)
        if row['panel_id'] not in review['panels'] or not all((math.isfinite(float(row[k])) for k in ('x_px', 'y_px'))):
            raise ValueError('Bad candidate coordinate frame')
    audits, matches, unmatched_marks, lane_rows = ([], [], [], [])
    for pid, panel in review['panels'].items():
        lanes = panel['lanes']
        eligible = {l['id']: [] for l in lanes}
        definite_by_lane = {l['id']: [b for b in panel['bands'] if b['lane_id'] == l['id'] and b['certainty'] == 'definite'] for l in lanes}
        for row in [r for r in candidates if r['panel_id'] == pid]:
            x, y = (float(row['x_px']), float(row['y_px']))
            audit = {'candidate_id': row['candidate_id'], 'panel_id': pid, 'x_px': x, 'y_px': y, 'detector_label': row.get('label', ''), 'evidence_tier': row.get('evidence_tier', ''), 'manual_lane_id': '', 'manual_label': '', 'matched_mark_id': '', 'absolute_y_error_px': None, 'inside_recorded_exclusion': False, 'near_uncertain_mark': False, 'disposition': ''}
            containing = [l for l in lanes if inside(x, y, l['box'])]
            if not inside(x, y, [0, 0, review['width_px'], review['height_px']]):
                audit['disposition'] = 'outside_native_image'
            elif not containing:
                audit['disposition'] = 'outside_supplied_lane_rectangles'
            elif len(containing) > 1:
                raise ValueError('Ambiguous overlapping lane assignment')
            else:
                lane = containing[0]
                lid = lane['id']
                audit.update(manual_lane_id=lid, manual_label=lane['label'])
                audit['inside_recorded_exclusion'] = any((inside(x, y, r['box']) for r in panel['exclusions']))
                audit['near_uncertain_mark'] = any((b['certainty'] == 'uncertain' and b['lane_id'] == lid and (abs(y - b['y']) <= uncertain_radius) for b in panel['bands']))
                if lane['role'] != 'sample':
                    audit['disposition'] = 'non_sample_' + lane['role']
                elif lane['assessment'] != 'readable':
                    audit['disposition'] = 'whole_lane_unreadable'
                elif mask_mode == 'recorded_masks' and audit['inside_recorded_exclusion']:
                    audit['disposition'] = 'recorded_unreadable_region'
                elif mask_mode == 'recorded_masks' and audit['near_uncertain_mark']:
                    audit['disposition'] = 'uncertain_mark_neighborhood'
                else:
                    audit['disposition'] = 'candidate_without_matching_mark'
                    eligible[lid].append(audit)
            audits.append(audit)
        for lane in sorted(lanes, key=lambda l: l['box'][0]):
            if lane['role'] != 'sample' or lane['assessment'] != 'readable':
                continue
            refs = definite_by_lane[lane['id']]
            preds = eligible[lane['id']]
            paired = match_1d([r['y'] for r in refs], [p['y_px'] for p in preds], tolerance)
            used = set()
            for ri, pi, err in paired:
                used.add(ri)
                p = preds[pi]
                b = refs[ri]
                p.update(disposition='matched_definite_mark', matched_mark_id=b['id'], absolute_y_error_px=err)
                matches.append({'panel_id': pid, 'lane_id': lane['id'], 'label': lane['label'], 'mark_id': b['id'], 'candidate_id': p['candidate_id'], 'mark_y_px': b['y'], 'candidate_y_px': p['y_px'], 'absolute_y_error_px': err})
            for i, b in enumerate(refs):
                if i not in used:
                    unmatched_marks.append({'panel_id': pid, 'lane_id': lane['id'], 'label': lane['label'], 'mark_id': b['id'], 'x_px': b['x'], 'y_px': b['y'], 'disposition': 'marked_band_without_candidate_at_this_tolerance'})
            n, m, t = (len(refs), len(preds), len(paired))
            lane_rows.append({'panel_id': pid, 'lane_id': lane['id'], 'label': lane['label'], 'n_definite_marks': n, 'n_eligible_candidates': m, 'n_matched': t, 'n_unmatched_marks': n - t, 'n_unmatched_candidates': m - t, 'marked_band_recovery_fraction': t / n if n else None, 'matched_fraction_of_candidates': t / m if m else None})
    summaries = []
    for pid in [*review['panels'], 'pooled']:
        rows = [r for r in lane_rows if pid == 'pooled' or r['panel_id'] == pid]
        aa = [r for r in audits if pid == 'pooled' or r['panel_id'] == pid]
        n = sum((r['n_definite_marks'] for r in rows))
        m = sum((r['n_eligible_candidates'] for r in rows))
        t = sum((r['n_matched'] for r in rows))
        summaries.append({'panel_id': pid, 'tolerance_px': tolerance, 'mask_mode': mask_mode, 'n_sample_lanes': len(rows), 'n_definite_marks': n, 'n_eligible_candidates': m, 'n_matched': t, 'n_unmatched_marks': n - t, 'n_unmatched_candidates': m - t, 'marked_band_recovery_fraction': t / n if n else None, 'matched_fraction_of_candidates': t / m if m else None, 'n_all_candidates': len(aa), 'n_not_evaluated_candidates': len(aa) - m, 'reference_coverage_assumed_complete': False, 'unmatched_candidate_interpretation': 'unreviewed_not_false_positive'})
    pooled = summaries[-1]
    if len(audits) != len(candidates) or len(matches) != pooled['n_matched'] or len(unmatched_marks) != pooled['n_unmatched_marks']:
        raise AssertionError('Comparison count conservation failed')
    return {'summary': summaries, 'lanes': lane_rows, 'matches': matches, 'unmatched_marks': unmatched_marks, 'candidate_audit': audits, 'dispositions': dict(Counter((r['disposition'] for r in audits)))}

def covered_area(review: dict, uncertain_radius: float=3.0) -> list[dict]:
    """Native pixel-center area within supplied sample ROIs, not whole-figure coverage."""
    import numpy as np
    rows = []
    for pid, panel in review['panels'].items():
        for lane in panel['lanes']:
            if lane['role'] != 'sample' or lane['assessment'] != 'readable':
                continue
            from .geometry import pixel_range
            l, r = pixel_range(lane['box'][0], lane['box'][2])
            t, b = pixel_range(lane['box'][1], lane['box'][3])
            yy, xx = np.mgrid[t:b, l:r]
            xx = xx + 0.5
            yy = yy + 0.5
            ex = np.zeros(xx.shape, dtype=bool)
            un = ex.copy()
            for reg in panel['exclusions']:
                x0, y0, x1, y1 = reg['box']
                ex |= (xx >= x0) & (xx < x1) & (yy >= y0) & (yy < y1)
            for mark in panel['bands']:
                if mark['lane_id'] == lane['id'] and mark['certainty'] == 'uncertain':
                    un |= np.abs(yy - mark['y']) <= uncertain_radius
            total = ex.size
            rows.append({'panel_id': pid, 'lane_id': lane['id'], 'label': lane['label'], 'n_sample_roi_pixels': total, 'n_exclusion_pixels': int(ex.sum()), 'n_uncertain_neighborhood_pixels': int(un.sum()), 'n_union_mask_pixels': int((ex | un).sum()), 'masked_roi_fraction': float((ex | un).mean()) if total else None, 'exhaustive_sample_coverage_claimed': False})
    return rows
