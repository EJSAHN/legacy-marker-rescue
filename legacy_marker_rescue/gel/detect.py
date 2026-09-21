"""Pixel-based candidate extraction in supplied lane regions."""
from pathlib import Path
from .io_utils import read_json, write_json, write_table, digest
from .profile import detect_file

def run(image: Path, geometry: Path, config: Path, output: Path, multiplier: float=1.0):
    output.mkdir(parents=True, exist_ok=False)
    results = detect_file(image, read_json(geometry), read_json(config), multiplier)
    candidates = [r for lane in results for r in lane['candidates']]
    write_json(output / 'candidates.json', candidates)
    write_table(output / 'band_candidates.tsv', candidates, list(candidates[0]) if candidates else ['candidate_id', 'panel_id', 'lane_id', 'label', 'role', 'x_px', 'y_px', 'evidence_tier'])
    maxima = [r for lane in results for r in lane['maxima']]
    write_table(output / 'all_profile_maxima.tsv', maxima)
    write_table(output / 'lane_profiles.tsv', [r for lane in results for r in lane['profile']])
    summaries = [{k: v for k, v in lane.items() if k not in {'profile', 'candidates', 'maxima'}} for lane in results]
    write_json(output / 'lane_summaries.json', summaries)
    write_table(output / 'lane_summaries.tsv', summaries)
    write_json(output / 'DETECTION_STATUS.json', {'status': 'CANDIDATE_EXTRACTION_COMPLETE', 'source_image_sha256': digest(image), 'geometry_sha256': digest(geometry), 'config_sha256': digest(config), 'threshold_multiplier': multiplier, 'n_candidates': len(candidates), 'manual_band_coordinates_used': False, 'manual_exclusions_used_to_detect': False, 'manual_lane_geometry_used': True, 'reference_coverage_assumed_complete': False})
    for pid in read_json(geometry)['panels']:
        rows = [r for r in candidates if r['panel_id'] == pid and r['role'] == 'sample']
        print(f"{pid}: {len(rows)} sample candidates; {sum((r['evidence_tier'] == 'strong_profile' for r in rows))} strong-profile candidates", flush=True)
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', type=Path, required=True)
    parser.add_argument('--geometry', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threshold-multiplier', type=float, default=1.0)
    args = parser.parse_args()
    run(args.image, args.geometry, args.config, args.output_dir, args.threshold_multiplier)
