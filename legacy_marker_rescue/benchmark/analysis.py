"""External image evaluation with explicit reference-geometry controls."""
from __future__ import annotations
import numpy as np
from .imaging import rescale_height
from .lanes import image_envelopes, extract_candidates

def auto_detect(gray: np.ndarray, method: str, benchmark: dict, config: dict, identifier: str) -> dict:
    """Primary prediction boundary: this function has NO mask/reference argument."""
    if method == 'image_lanes_native':
        working = gray
        sx = sy = 1.0
    elif method == 'image_lanes_height256':
        working, sx, sy = rescale_height(gray, benchmark['preprocessing']['sensitivity_height_px'])
    else:
        raise ValueError('Not an image-only method')
    envelopes, props = image_envelopes(working, benchmark['image_lane_proposal'])
    candidates, lanes = extract_candidates(working, envelopes, config, identifier, sx, sy)
    return {'method': method, 'reference_mask_used': False, 'working_shape': list(working.shape), 'scale_x': sx, 'scale_y': sy, 'native_envelopes': [[l / sx, 0, r / sx, gray.shape[0]] for l, r in envelopes], 'lane_proposal': props, 'candidates': candidates, 'lane_summaries': lanes}

import hashlib
import json
import traceback
import warnings
import zipfile
from pathlib import Path
from ..gel.io_utils import write_json, write_table, digest
from .acquire import retrieve, pair_members
from .imaging import decode_image, decode_mask, component_reference
from .lanes import reference_x_envelopes
from .scoring import evaluate, summaries


def process_pair(z: zipfile.ZipFile, pair: dict, collection: str,
                 benchmark: dict, detector_config: dict, out: Path) -> tuple[list[dict], dict]:
    folder = out / 'images' / pair['id']
    folder.mkdir(parents=True, exist_ok=False)
    row = {'collection': collection, **pair, 'status': 'PROCESSING', 'methods': {}}
    scores = []
    try:
        image_bytes = z.read(pair['image_member'])
        image_sha = hashlib.sha256(image_bytes).hexdigest()
        row['raw_image_sha256'] = image_sha
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            gray, properties = decode_image(image_bytes, benchmark['preprocessing'])
        row['decoding_warnings'] = [{'category': w.category.__name__, 'message': str(w.message)} for w in captured]
        row['image_properties'] = properties
        for method in ('image_lanes_native', 'image_lanes_height256'):
            try:
                predicted = auto_detect(gray, method, benchmark, detector_config, image_sha)
                file = folder / (method + '_predictions.json')
                write_json(file, predicted)
                row['methods'][method] = {'status': 'PREDICTED', 'file': file.relative_to(out).as_posix(),
                                          'sha256_before_reading_reference': digest(file)}
            except Exception as error:
                row['methods'][method] = {'status': 'PREDICTION_FAILED', 'error': type(error).__name__ + ': ' + str(error)}
        mask_bytes = z.read(pair['mask_member'])
        row['raw_mask_sha256'] = hashlib.sha256(mask_bytes).hexdigest()
        mask, mask_info = decode_mask(mask_bytes, gray.shape)
        labels, references = component_reference(mask)
        row['mask_properties'] = mask_info
        row['reference_components'] = len(references)
        row['foreground_pixels'] = int(mask.sum())
        write_table(folder / 'reference_components.tsv', references,
                    ['component_id', 'x0', 'x1', 'y0', 'y1', 'x_center', 'y_center', 'area_px', 'touches_image_edge'])
        diagnostic = 'reference_x_envelopes_native'
        try:
            envelopes = reference_x_envelopes(mask)
            candidates, lane_summaries = extract_candidates(gray, envelopes, detector_config, image_sha)
            predicted = {'method': diagnostic, 'reference_mask_used': True,
                         'reference_information_used': 'foreground_x_envelopes_only',
                         'working_shape': list(gray.shape), 'scale_x': 1., 'scale_y': 1.,
                         'native_envelopes': [[l, 0, r, gray.shape[0]] for l, r in envelopes],
                         'candidates': candidates, 'lane_summaries': lane_summaries}
            file = folder / (diagnostic + '_predictions.json')
            write_json(file, predicted)
            row['methods'][diagnostic] = {'status': 'PREDICTED', 'file': file.relative_to(out).as_posix()}
        except Exception as error:
            row['methods'][diagnostic] = {'status': 'PREDICTION_FAILED', 'error': type(error).__name__ + ': ' + str(error)}
        tolerances = [benchmark['evaluation']['primary_tolerance_native_px'], *benchmark['evaluation']['sensitivity_tolerances_native_px']]
        for method, entry in row['methods'].items():
            if entry['status'] != 'PREDICTED': continue
            file = out / entry['file']
            if 'sha256_before_reading_reference' in entry and digest(file) != entry['sha256_before_reading_reference']:
                raise RuntimeError('Predictions changed after reference access')
            predicted = json.loads(file.read_text(encoding='utf-8'))
            metric_rows = []
            try:
                for tolerance in tolerances:
                    result = evaluate(predicted['candidates'], labels, references, tolerance)
                    write_json(folder / f'{method}_matches_{tolerance}px.json', result)
                    metrics = {k: v for k, v in result.items() if k not in {'matches', 'unmatched_prediction_indices', 'unmatched_component_ids'}}
                    metric_rows.append({'collection': collection, 'image_id': pair['id'], 'raw_image_sha256': image_sha,
                                        'original_name': pair['original_name'], 'method': method,
                                        'reference_geometry_used': method == diagnostic, **metrics})
                scores.extend(metric_rows)
                entry['status'] = 'EVALUATED'
            except Exception as error:
                entry['status'] = 'EVALUATION_FAILED'; entry['error'] = type(error).__name__ + ': ' + str(error)
        row['provenance'] = {'image_only_predictions_committed_before_mask_read': True, 'pixels_for_prediction': 'normalized_source_image', 'reference_changed': False, 'detector_tuned': False, 'manual_Denoyes_marks_used': False}
        row['status'] = 'COMPLETE' if all(v['status'] == 'EVALUATED' for v in row['methods'].values()) else 'PARTIAL_METHOD_FAILURE'
    except Exception as error:
        row.update(status='FAILED', error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
    write_json(folder / 'IMAGE_STATUS.json', row)
    return scores, row


def run(benchmark: dict, detector_config: dict, cache: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    sources = out / 'sources'; sources.mkdir()
    available, collections, planned = [], [], []
    for collection in benchmark['collections']:
        entry = {'collection': collection['id'], 'expected_image_pairs': collection['expected_image_pairs'], 'role': collection['role']}
        try:
            path, receipt = retrieve(collection, benchmark['record_id'], cache, cache)
            write_json(sources / (collection['id'] + '_receipt.json'), receipt)
            pairs, notes = pair_members(path)
            entry.update(status='PAIRED', pair_count=len(pairs), pairing_issues=notes,
                         source_count_matches=len(pairs) == collection['expected_image_pairs'], archive_sha256=receipt['sha256'])
            available.append((collection, path, pairs))
            planned.extend({'collection': collection['id'], **p} for p in pairs)
        except Exception as error:
            entry.update(status='ACQUISITION_OR_PAIRING_FAILED', error=type(error).__name__ + ': ' + str(error))
        collections.append(entry)
    write_json(sources / 'selected_images_before_scoring.json', planned)
    write_table(sources / 'selected_images_before_scoring.tsv', planned)
    write_json(sources / 'collection_status.json', collections)
    statuses, metric_rows = [], []
    index = 0
    for collection, path, pairs in available:
        with zipfile.ZipFile(path) as archive:
            for pair in pairs:
                index += 1
                print(f'EXTERNAL {index}/{len(planned)}: {collection["id"]} / {pair["original_name"]}', flush=True)
                scores, status = process_pair(archive, pair, collection['id'], benchmark, detector_config, out)
                statuses.append(status); metric_rows.extend(scores)
                write_json(out / 'image_status.json', statuses)
    write_json(out / 'per_image_metrics.json', metric_rows)
    write_table(out / 'per_image_metrics.tsv', metric_rows)
    aggregate = summaries(metric_rows)
    write_json(out / 'summary.json', aggregate); write_table(out / 'summary.tsv', aggregate)
    complete = (all(x['status'] == 'PAIRED' and x.get('source_count_matches') and not x.get('pairing_issues') for x in collections)
                and len(statuses) == sum(x['expected_image_pairs'] for x in benchmark['collections'])
                and all(x['status'] == 'COMPLETE' for x in statuses))
    status = {'status': 'COMPLETE' if complete else 'INCOMPLETE', 'collections': collections,
              'images': len(statuses), 'reference_used_for_primary_prediction': False,
              'interpretation': 'Center-to-foreground correspondence. Reference-geometry and scaled-image results are separate diagnostics.'}
    write_json(out / 'STATUS.json', status)
    if not complete: raise RuntimeError('External dataset processing incomplete; inspect the per-image status records')
    return status
