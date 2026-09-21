"""Lane-guided image analysis and partial-reader correspondence."""
from __future__ import annotations
from pathlib import Path
import copy
from .io_utils import read_json, write_json, write_table, read_table, digest
from .geometry import strip_geometry
from .observations import validate_review
from .detect import run as detect
from .comparison import compare_snapshot, covered_area, canonical_legacy


def run(data: Path, config_dir: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    image = data / 'source_image.png'
    observations = data / 'observations.json'
    review = read_json(observations)
    protocol = data / 'review_protocol.json'
    validate_review(review, read_json(protocol), digest(protocol), require_final=False)
    if digest(image) != review['image_sha256']: raise ValueError('Image and reader coordinate frames differ')
    scope = read_json(data / 'reader_scope.json')
    if scope['annotation_sha256'] != digest(observations): raise ValueError('Reader scope identifies a different observation file')
    geometry = strip_geometry(review)
    geometry_path = out / 'lane_geometry.json'
    write_json(geometry_path, geometry)
    altered = copy.deepcopy(review)
    for panel in altered['panels'].values():
        panel['bands'] = []; panel['exclusions'] = []; panel['coverage_complete'] = True
    if strip_geometry(altered) != geometry: raise ValueError('Band or exclusion data entered detector geometry')
    config = read_json(config_dir / 'detection.json')
    comparison_config = read_json(config_dir / 'comparison.json')
    variants = {'primary': config['primary_threshold_multiplier']}
    variants.update({f'sensitivity_{x:g}'.replace('.', 'p'): x for x in config['sensitivity_threshold_multipliers']})
    candidates = {}
    for name, multiplier in variants.items():
        dest = out / 'detector_runs' / name
        print(f'DETECTION: {name}', flush=True)
        detect(image, geometry_path, config_dir / 'detection.json', dest, multiplier)
        candidates[name] = read_json(dest / 'candidates.json')
    previous = canonical_legacy(read_table(data / 'previous_predictions.tsv'))
    all_candidates = {'legacy': previous, **candidates}
    output = out / 'comparison'; output.mkdir()
    summaries, queue = [], []
    for name, values in all_candidates.items():
        for mode in comparison_config['exclusion_modes']:
            for tolerance in [comparison_config['original_primary_y_tolerance_px'], *comparison_config['sensitivity_y_tolerances_px']]:
                result = compare_snapshot(review, values, tolerance, mode, comparison_config['uncertain_neighborhood_px'])
                summaries.extend({'method': name, **row} for row in result['summary'])
                stem = f'{name}_{tolerance:g}px_{mode}'
                for field in ('lanes', 'matches', 'unmatched_marks', 'candidate_audit'):
                    write_table(output / (stem + '_' + field + '.tsv'), result[field])
                if name == 'primary' and tolerance == 3 and mode == 'recorded_masks':
                    for row in result['candidate_audit']:
                        if row['manual_lane_id'] and row['disposition'] not in {'matched_definite_mark', 'non_sample_marker'}:
                            queue.append({'panel_id': row['panel_id'], 'lane_id': row['manual_lane_id'], 'label': row['manual_label'],
                                          'item_id': row['candidate_id'], 'x_px': row['x_px'], 'y_px': row['y_px'],
                                          'reason': row['disposition'], 'action': 'review_source_image'})
                    for row in result['unmatched_marks']:
                        queue.append({'panel_id': row['panel_id'], 'lane_id': row['lane_id'], 'label': row['label'],
                                      'item_id': row['mark_id'], 'x_px': row['x_px'], 'y_px': row['y_px'],
                                      'reason': row['disposition'], 'action': 'review_source_image'})
    write_table(output / 'agreement_summary.tsv', summaries)
    write_table(output / 'discrepancy_queue.tsv', queue)
    write_table(output / 'mask_area_within_sample_rois.tsv', covered_area(review, comparison_config['uncertain_neighborhood_px']))
    counts = {}
    for panel_id in geometry['panels']:
        rows = [r for r in candidates['primary'] if r['panel_id'] == panel_id]
        counts[panel_id] = {'sample': sum(r['role'] == 'sample' for r in rows), 'marker': sum(r['role'] == 'marker' for r in rows)}
    status = {'status': 'COMPLETE', 'candidate_counts': counts, 'image_sha256': digest(image),
              'observations_sha256': digest(observations), 'manual_lane_geometry_used': True,
              'manual_band_positions_used_to_detect': False, 'complete_reference_assumed': False,
              'interpretation': 'Correspondence with partial reader observations; unmatched candidates are not established false positives.'}
    write_json(out / 'STATUS.json', status)
    return status
