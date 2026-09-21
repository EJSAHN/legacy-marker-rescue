"""Extract bright/dark band candidates from pixels inside supplied lanes.

This module has no reference-annotation input and never reads reader marks,
uncertain neighborhoods, or exclusion regions. Detection does not force a band
count, snap to a click, or borrow bands from another lane.
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d, grey_opening
from scipy.signal import find_peaks
from .geometry import validate_geometry, pixel_range
from .io_utils import digest

def load_gray(path: Path) -> np.ndarray:
    from PIL import Image
    with Image.open(path) as image:
        if image.mode not in ('L', 'RGB', 'RGBA'):
            raise ValueError('This detector accepts 8-bit grayscale/RGB images only')
        if image.mode == 'RGBA' and np.asarray(image)[..., 3].min() < 255:
            raise ValueError('Transparent source pixels are not supported')
        return np.array(image.convert('L'), dtype=np.float64)

def check_config(c: dict) -> None:
    if c.get('schema_version') != '1.0' or c.get('polarity') not in {'bright', 'dark'}:
        raise ValueError('Invalid detector schema/polarity')
    fraction_keys = ('center_width_fraction', 'candidate_dynamic_fraction', 'strong_dynamic_fraction', 'minimum_candidate_support_fraction', 'minimum_strong_support_fraction', 'saturation_flag_fraction')
    for k in fraction_keys:
        v = c[k]
        if isinstance(v, bool) or not isinstance(v, (float, int)) or (not math.isfinite(v)) or (not 0 < v <= 1):
            raise ValueError(f'Invalid fraction: {k}')
    integer_keys = ('minimum_lane_width_px', 'minimum_lane_height_px', 'background_opening_size_px', 'minimum_peak_distance_px', 'edge_guard_px')
    for k in integer_keys:
        if isinstance(c[k], bool) or not isinstance(c[k], int) or c[k] < 1:
            raise ValueError(f'Invalid positive integer: {k}')
    if c['background_opening_size_px'] % 2 != 1:
        raise ValueError('Background window must be odd')
    for k in ('smoothing_sigma_px', 'candidate_absolute_prominence', 'candidate_noise_multiplier', 'strong_absolute_prominence', 'strong_noise_multiplier', 'noise_floor', 'column_support_noise_multiplier', 'column_support_floor', 'minimum_peak_width_px', 'maximum_strong_peak_width_px', 'saturation_threshold', 'primary_threshold_multiplier'):
        v = c[k]
        if isinstance(v, bool) or not isinstance(v, (float, int)) or (not math.isfinite(v)) or (v <= 0):
            raise ValueError(f'Invalid positive number: {k}')
    if not c['sensitivity_threshold_multipliers'] or any((not isinstance(x, (int, float)) or isinstance(x, bool) or (not math.isfinite(x)) or (x <= 0) for x in c['sensitivity_threshold_multipliers'])):
        raise ValueError('Invalid sensitivity multipliers')
    if c['minimum_candidate_support_fraction'] > c['minimum_strong_support_fraction']:
        raise ValueError('Strong support threshold below candidate threshold')

def robust_noise(signal: np.ndarray, floor: float) -> float:
    d = np.diff(signal)
    if not len(d):
        return float(floor)
    return max(float(floor), float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0)))

def analyze_lane(gray: np.ndarray, lane: dict, panel: str, c: dict, multiplier: float=1.0) -> dict:
    """Median column-wise top-hat profile, followed by deterministic peak screening."""
    x0, y0, x1, y1 = [float(v) for v in lane['box']]
    left, right = pixel_range(x0, x1)
    top, bottom = pixel_range(y0, y1)
    result = {'panel_id': panel, 'lane_id': lane['id'], 'label': lane['label'], 'role': lane['role'], 'box': [x0, y0, x1, y1], 'candidates': [], 'maxima': [], 'profile': []}
    if lane['role'] not in {'sample', 'marker'}:
        result['status'] = 'ROLE_NOT_SCORABLE'
        return result
    if lane['assessment'] != 'readable':
        result['status'] = 'READER_MARKED_WHOLE_LANE_UNREADABLE'
        return result
    if right - left < c['minimum_lane_width_px'] or bottom - top < c['minimum_lane_height_px']:
        result['status'] = 'ROI_TOO_SMALL'
        return result
    width = x1 - x0
    center = (x0 + x1) / 2
    il, ir = pixel_range(center - width * c['center_width_fraction'] / 2, center + width * c['center_width_fraction'] / 2)
    il, ir = (max(left, il), min(right, ir))
    if ir - il < 3:
        result['status'] = 'CENTER_STRIP_TOO_SMALL'
        return result
    source = gray[top:bottom, il:ir]
    oriented = source if c['polarity'] == 'bright' else 255.0 - source
    background = grey_opening(oriented, size=(c['background_opening_size_px'], 1), mode='reflect')
    residual_columns = oriented - background
    raw = np.median(oriented, axis=1)
    background_profile = np.median(background, axis=1)
    residual = np.median(residual_columns, axis=1)
    smooth = gaussian_filter1d(residual, c['smoothing_sigma_px'], mode='reflect')
    column_smooth = gaussian_filter1d(residual_columns, c['smoothing_sigma_px'], axis=0, mode='reflect')
    noise = robust_noise(smooth, c['noise_floor'])
    dynamic = max(0.0, float(np.percentile(raw, 95) - np.percentile(raw, 5)))
    threshold = float(multiplier * max(c['candidate_absolute_prominence'], c['candidate_noise_multiplier'] * noise, c['candidate_dynamic_fraction'] * dynamic))
    strong_threshold = float(multiplier * max(c['strong_absolute_prominence'], c['strong_noise_multiplier'] * noise, c['strong_dynamic_fraction'] * dynamic))
    column_cutoff = max(c['column_support_floor'], c['column_support_noise_multiplier'] * noise)
    peaks, prop = find_peaks(smooth, prominence=(None, None), width=(None, None), plateau_size=(None, None))
    rows = []
    for i, peak in enumerate(peaks):
        peak = int(peak)
        prominence = float(prop['prominences'][i])
        peak_width = float(prop['widths'][i])
        support = float(np.mean(column_smooth[peak] >= column_cutoff))
        saturation = float(np.mean(oriented[peak] >= c['saturation_threshold']))
        flags = []
        if peak_width > c['maximum_strong_peak_width_px']:
            flags.append('broad_profile_feature')
        if saturation >= c['saturation_flag_fraction']:
            flags.append('saturated_source_pixels')
        if int(prop['plateau_sizes'][i]) > 2:
            flags.append('flat_profile_peak')
        edge = min(peak, len(smooth) - 1 - peak) < c['edge_guard_px']
        if edge:
            flags.append('roi_edge')
        y = top + peak + 0.5
        row = {'candidate_id': f"{panel}:{lane['id']}:y{top + peak}", 'panel_id': panel, 'lane_id': lane['id'], 'label': lane['label'], 'role': lane['role'], 'x_px': float((il + ir) / 2), 'y_px': float(y), 'support_x0_px': int(il), 'support_x1_px': int(ir), 'half_prominence_y0_px': float(top + prop['left_ips'][i] + 0.5), 'half_prominence_y1_px': float(top + prop['right_ips'][i] + 0.5), 'peak_width_px': peak_width, 'prominence': prominence, 'profile_height': float(smooth[peak]), 'noise_mad_estimate': noise, 'candidate_threshold': threshold, 'strong_threshold': strong_threshold, 'horizontal_support_fraction': support, 'saturated_fraction': saturation, 'flags': flags, 'evidence_tier': 'rejected', 'reason': '', 'accepted': False}
        if edge:
            row['reason'] = 'roi_edge'
        elif prominence < threshold:
            row['reason'] = 'below_prominence_threshold'
        elif peak_width < c['minimum_peak_width_px']:
            row['reason'] = 'too_narrow'
        elif support < c['minimum_candidate_support_fraction']:
            row['reason'] = 'limited_horizontal_support'
        else:
            row['accepted'] = True
            row['evidence_tier'] = 'strong_profile' if prominence >= strong_threshold and support >= c['minimum_strong_support_fraction'] and (not flags) else 'review_candidate'
            row['reason'] = 'profile_criteria_met'
        rows.append(row)
    kept = []
    for row in sorted([r for r in rows if r['accepted']], key=lambda r: (-r['prominence'], r['y_px'])):
        if any((abs(row['y_px'] - other['y_px']) < c['minimum_peak_distance_px'] for other in kept)):
            row.update(accepted=False, evidence_tier='rejected', reason='close_to_stronger_profile_peak')
        else:
            kept.append(row)
    result['maxima'] = sorted(rows, key=lambda r: r['y_px'])
    result['candidates'] = sorted(kept, key=lambda r: r['y_px'])
    for j in range(len(smooth)):
        result['profile'].append({'panel_id': panel, 'lane_id': lane['id'], 'label': lane['label'], 'y_px': float(top + j + 0.5), 'raw_median_intensity': float(raw[j]), 'background_median_intensity': float(background_profile[j]), 'median_residual': float(residual[j]), 'smoothed_residual': float(smooth[j]), 'candidate_prominence_threshold': threshold, 'strong_prominence_threshold': strong_threshold})
    result.update(status='PROFILE_PROCESSED', native_roi_indices=[left, top, right, bottom], center_strip_indices=[il, top, ir, bottom], noise_mad_estimate=noise, dynamic_range_p95_p5=dynamic, candidate_threshold=threshold, strong_threshold=strong_threshold, n_sampled_rows=bottom - top, n_sampled_columns=ir - il, n_candidates=len(kept), n_strong=sum((r['evidence_tier'] == 'strong_profile' for r in kept)), n_review=sum((r['evidence_tier'] == 'review_candidate' for r in kept)), n_rejected_maxima=sum((not r['accepted'] for r in rows)), edge_guard_px=c['edge_guard_px'])
    return result

def detect_array(gray: np.ndarray, geometry: dict, config: dict, multiplier: float=1.0) -> list[dict]:
    validate_geometry(geometry)
    check_config(config)
    if gray.shape != (geometry['height_px'], geometry['width_px']) or gray.ndim != 2:
        raise ValueError('Pixel array and geometry dimensions differ')
    if not np.isfinite(gray).all() or gray.min() < 0 or gray.max() > 255:
        raise ValueError('Only finite 0-255 image values are supported')
    if not math.isfinite(multiplier) or multiplier <= 0:
        raise ValueError('Invalid threshold multiplier')
    return [analyze_lane(gray, lane, pid, config, multiplier) for pid, panel in geometry['panels'].items() for lane in sorted(panel['lanes'], key=lambda l: (l['box'][0], l['id']))]

def detect_file(image_path: Path, geometry: dict, config: dict, multiplier: float=1.0) -> list[dict]:
    if digest(image_path) != geometry['image_sha256']:
        raise ValueError('The source image hash does not match the native geometry')
    return detect_array(load_gray(image_path), geometry, config, multiplier)
