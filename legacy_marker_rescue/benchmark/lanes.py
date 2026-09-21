"""Image-only lane proposals and a separately labelled reference-geometry control.

All proposals span the FULL image height. Source masks never enter the primary
image-only proposal function. It is a fixed new adapter, not part of the older
manually lane-guided method and not a learned lane detector.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import grey_opening, gaussian_filter1d, binary_closing, label
from scipy.signal import find_peaks
from legacy_marker_rescue.gel.profile import detect_array

def runs(active: np.ndarray) -> list[tuple[int, int]]:
    a = np.r_[False, np.asarray(active, dtype=bool), False].astype(np.int8)
    return list(zip(np.flatnonzero(np.diff(a) == 1).tolist(), np.flatnonzero(np.diff(a) == -1).tolist()))

def image_envelopes(gray: np.ndarray, config: dict) -> tuple[list[tuple[int, int]], dict]:
    if gray.ndim != 2 or not np.isfinite(gray).all():
        raise ValueError('Invalid image')
    bg = grey_opening(gray, size=(config['vertical_opening_px'], 1), mode='reflect')
    score = gaussian_filter1d(np.mean(np.maximum(gray - bg, 0), axis=0), config['x_smoothing_sigma_px'], mode='reflect')
    threshold = max(config['minimum_projection_signal'], config['projection_fraction_of_p99'] * float(np.percentile(score, 99)))
    active = score >= threshold
    n = config['close_gap_px']
    active = np.pad(active, (n + 1, n + 1))
    active = binary_closing(active, structure=np.ones(n + 1, dtype=bool))[n + 1:-n - 1]
    boxes = []
    for left, right in runs(active):
        if right - left < config['minimum_envelope_width_px']:
            continue
        segment = score[left:right]
        seeds, _ = find_peaks(segment, distance=config['split_peak_distance_px'], prominence=max(config['minimum_projection_signal'], config['split_peak_prominence_fraction'] * float(segment.max())))
        if len(seeds) < 2:
            boxes.append((left, right))
            continue
        seeds = list(seeds + left)
        cuts = [left]
        for a, b in zip(seeds, seeds[1:]):
            cuts.append(int(a + np.argmin(score[a:b + 1])))
        cuts.append(right)
        for a, b in zip(cuts, cuts[1:]):
            if b - a >= config['minimum_envelope_width_px']:
                boxes.append((a, b))
    if len(boxes) > config['maximum_lanes']:
        raise ValueError('LANE_PROPOSAL_SAFETY_LIMIT')
    return (boxes, {'lane_source': 'image_only_vertical_contrast_projection', 'n_proposals': len(boxes), 'projection_threshold': float(threshold), 'projection': score.tolist(), 'full_height': True, 'reference_mask_used': False})

def reference_x_envelopes(mask: np.ndarray) -> list[tuple[int, int]]:
    """Oracle-like diagnostic, not independent automatic lane detection.

    Uses only foreground x support; never crops in y or passes band centers to
    the detector. Lanes whose foreground x projections touch are merged. Blank
    lanes cannot be supplied by this control and all reference objects are still
    evaluated, even if an ROI is too small for the detector.
    """
    return runs(mask.any(axis=0))

def make_geometry(gray: np.ndarray, envelopes: list[tuple[int, int]], identifier: str) -> dict:
    h, w = gray.shape
    return {'schema_version': '1.0', 'coordinate_system': 'native_pixel_edges', 'image_sha256': identifier, 'width_px': w, 'height_px': h, 'panels': {'gel': {'lanes': [{'id': f'window_{i + 1:04d}', 'box': [int(x0), 0, int(x1), h], 'label': f'window_{i + 1:04d}', 'label_source': 'unresolved', 'role': 'sample', 'assessment': 'readable'} for i, (x0, x1) in enumerate(envelopes)]}}}

def extract_candidates(gray: np.ndarray, envelopes: list[tuple[int, int]], detector_config: dict, identifier: str, sx=1.0, sy=1.0) -> tuple[list[dict], list[dict]]:
    results = detect_array(gray, make_geometry(gray, envelopes, identifier), detector_config, 1.0)
    candidates = []
    summary = []
    for r in results:
        for c in r['candidates']:
            row = dict(c)
            row['role'] = 'unclassified_gel_window'
            row['flags'] = ['clipped_working_copy_pixels' if f == 'saturated_source_pixels' else f for f in row['flags']]
            row['x_px'] /= sx
            row['y_px'] /= sy
            for key in ('support_x0_px', 'support_x1_px'):
                row[key] = float(row[key]) / sx
            for key in ('half_prominence_y0_px', 'half_prominence_y1_px', 'peak_width_px'):
                row[key] /= sy
            candidates.append(row)
        summary.append({k: v for k, v in r.items() if k not in {'profile', 'maxima', 'candidates'}})
    return (candidates, summary)
