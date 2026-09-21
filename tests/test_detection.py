"""Software controls only; synthetic arrays are not study observations."""
from __future__ import annotations
import copy
import itertools
import json
import math
from pathlib import Path
import random
import tempfile
import unittest
import numpy as np
from legacy_marker_rescue.gel.geometry import strip_geometry, validate_geometry, pixel_range
from legacy_marker_rescue.gel.io_utils import read_json, write_json, write_table, read_table, digest
from legacy_marker_rescue.gel.profile import detect_array, detect_file, load_gray, check_config, robust_noise
from legacy_marker_rescue.gel.observations import match_1d
from legacy_marker_rescue.gel.comparison import compare_snapshot, covered_area
ROOT = Path(__file__).resolve().parents[1]
CONFIG = read_json(ROOT / 'config/detection.json')

def lane(lid='L1', box=None):
    return {'id': lid, 'box': box or [10.0, 10.0, 40.0, 180.0], 'label': 'sample', 'label_source': 'printed_here', 'role': 'sample', 'assessment': 'readable'}

def geometry(lanes=None):
    return {'schema_version': '1.0', 'coordinate_system': 'native_pixel_edges', 'image_sha256': 'test-fixture', 'width_px': 100, 'height_px': 200, 'panels': {'A': {'lanes': lanes or [lane()]}}}

def picture(ys=(60, 110), amplitudes=None, sigma=2.0):
    image = np.ones((200, 100), dtype=float) * 30
    amplitudes = amplitudes or [100.0] * len(ys)
    for y, a in zip(ys, amplitudes):
        image[:, 10:40] += a * np.exp(-0.5 * ((np.arange(200) - y) / sigma) ** 2)[:, None]
    return np.clip(image, 0, 255)

def review_fixture():
    return {'width_px': 100, 'height_px': 200, 'panels': {'A': {'lanes': [lane()], 'bands': [{'id': 'b1', 'lane_id': 'L1', 'x': 25.0, 'y': 60.5, 'certainty': 'definite'}], 'exclusions': [], 'coverage_complete': False}}, 'status': 'draft', 'reviewer': {'name': 'test', 'prior_algorithm_exposure': 'yes'}}

class ProfileTests(unittest.TestCase):

    def test_two_isolated_bands(self):
        r = detect_array(picture(), geometry(), CONFIG)[0]
        self.assertEqual(len(r['candidates']), 2)
        self.assertEqual([x['y_px'] for x in r['candidates']], [60.5, 110.5])

    def test_blank_has_no_band(self):
        self.assertEqual(detect_array(np.zeros((200, 100)), geometry(), CONFIG)[0]['candidates'], [])

    def test_uniform_bright_not_a_band(self):
        self.assertEqual(detect_array(np.ones((200, 100)) * 255, geometry(), CONFIG)[0]['candidates'], [])

    def test_line_outside_lane_not_seen(self):
        im = np.zeros((200, 100))
        im[60:63, 60:90] = 255
        self.assertEqual(detect_array(im, geometry(), CONFIG)[0]['candidates'], [])

    def test_narrow_vertical_streak_not_band(self):
        im = np.zeros((200, 100))
        im[:, 24:26] = 250
        self.assertEqual(detect_array(im, geometry(), CONFIG)[0]['candidates'], [])

    def test_edge_column_artifact_not_used(self):
        im = np.zeros((200, 100))
        im[59:62, 10:12] = 255
        self.assertEqual(detect_array(im, geometry(), CONFIG)[0]['candidates'], [])

    def test_reversed_polarity_equivalent(self):
        config = copy.deepcopy(CONFIG)
        config['polarity'] = 'dark'
        a = detect_array(picture(), geometry(), CONFIG)[0]['candidates']
        b = detect_array(255 - picture(), geometry(), config)[0]['candidates']
        self.assertEqual([r['y_px'] for r in a], [r['y_px'] for r in b])

    def test_positive_gain_fixed_noise_floor(self):
        r = detect_array(picture(), geometry(), CONFIG)[0]['candidates']
        self.assertTrue(all((x['prominence'] > 0 and x['horizontal_support_fraction'] >= 0.65 for x in r)))

    def test_saturation_is_flagged(self):
        im = picture((60,), [300], sigma=3)
        r = detect_array(im, geometry(), CONFIG)[0]['candidates']
        self.assertTrue(any(('saturated_source_pixels' in x['flags'] for x in r)))
        self.assertTrue(all((x['evidence_tier'] == 'review_candidate' for x in r)))

    def test_reproducible_noise(self):
        im = picture() + np.random.default_rng(123).normal(0, 0.8, (200, 100))
        self.assertEqual(detect_array(im, geometry(), CONFIG), detect_array(im, geometry(), CONFIG))

    def test_input_pixels_unchanged(self):
        im = picture()
        before = im.copy()
        detect_array(im, geometry(), CONFIG)
        np.testing.assert_array_equal(im, before)

    def test_label_does_not_change_positions(self):
        g = geometry()
        a = detect_array(picture(), g, CONFIG)[0]['candidates']
        g['panels']['A']['lanes'][0]['label'] = 'unrelated label'
        b = detect_array(picture(), g, CONFIG)[0]['candidates']
        self.assertEqual([(r['x_px'], r['y_px'], r['prominence']) for r in a], [(r['x_px'], r['y_px'], r['prominence']) for r in b])

    def test_translation_preserves_native_coordinates(self):
        a = detect_array(picture(), geometry(), CONFIG)[0]['candidates']
        im = np.roll(picture(), 5, axis=0)
        g = geometry()
        g['panels']['A']['lanes'][0]['box'] = [10.0, 15.0, 40.0, 185.0]
        b = detect_array(im, g, CONFIG)[0]['candidates']
        self.assertEqual([r['y_px'] + 5 for r in a], [r['y_px'] for r in b])

    def test_distance_suppression(self):
        c = copy.deepcopy(CONFIG)
        c['minimum_peak_distance_px'] = 100
        self.assertEqual(len(detect_array(picture(), geometry(), c)[0]['candidates']), 1)

    def test_minimum_lane_size_is_reported(self):
        r = detect_array(picture(), geometry([lane(box=[10.0, 10.0, 13.0, 180.0])]), CONFIG)[0]
        self.assertEqual(r['status'], 'ROI_TOO_SMALL')

    def test_marker_kept_separate(self):
        l = lane()
        l['role'] = 'marker'
        self.assertTrue(all((r['role'] == 'marker' for r in detect_array(picture(), geometry([l]), CONFIG)[0]['candidates'])))

    def test_unreadable_lane_not_reported_zero_absence(self):
        l = lane()
        l['assessment'] = 'unreadable'
        r = detect_array(picture(), geometry([l]), CONFIG)[0]
        self.assertEqual(r['status'], 'READER_MARKED_WHOLE_LANE_UNREADABLE')

    def test_gap_not_analyzed(self):
        l = lane()
        l['role'] = 'gap'
        self.assertEqual(detect_array(picture(), geometry([l]), CONFIG)[0]['status'], 'ROLE_NOT_SCORABLE')

    def test_nan_image_rejected(self):
        im = picture()
        im[0, 0] = np.nan
        with self.assertRaises(ValueError):
            detect_array(im, geometry(), CONFIG)

    def test_image_shape_mismatch_rejected(self):
        with self.assertRaises(ValueError):
            detect_array(np.zeros((20, 20)), geometry(), CONFIG)

    def test_every_local_maximum_is_accounted_for(self):
        r = detect_array(picture(), geometry(), CONFIG)[0]
        self.assertEqual(len(r['maxima']), r['n_candidates'] + r['n_rejected_maxima'])

    def test_sensitive_threshold_not_selected_from_annotations(self):
        self.assertEqual(CONFIG['primary_threshold_multiplier'], 1.0)
        self.assertEqual(CONFIG['sensitivity_threshold_multipliers'], [0.75, 1.25])

    def test_noise_floor(self):
        self.assertEqual(robust_noise(np.ones(20), 0.5), 0.5)

    def test_even_window_rejected(self):
        c = copy.deepcopy(CONFIG)
        c['background_opening_size_px'] = 30
        with self.assertRaises(ValueError):
            check_config(c)

    def test_bad_fraction_rejected(self):
        c = copy.deepcopy(CONFIG)
        c['center_width_fraction'] = 2
        with self.assertRaises(ValueError):
            check_config(c)

    def test_negative_multiplier_rejected(self):
        with self.assertRaises(ValueError):
            detect_array(picture(), geometry(), CONFIG, -1)

class GeometryTests(unittest.TestCase):

    def test_half_open_pixel_center_conversion(self):
        self.assertEqual(pixel_range(10.5, 20.5), (10, 20))
        self.assertEqual(pixel_range(10.6, 20.4), (11, 20))

    def test_geometry_disallows_reference_lists(self):
        g = geometry()
        g['panels']['A']['bands'] = []
        with self.assertRaises(ValueError):
            validate_geometry(g)

    def test_overlap_rejected(self):
        with self.assertRaises(ValueError):
            validate_geometry(geometry([lane(), lane('L2', [20.0, 10.0, 45.0, 180.0])]))

    def test_out_of_bounds_rejected(self):
        with self.assertRaises(ValueError):
            validate_geometry(geometry([lane(box=[-1.0, 10.0, 40.0, 180.0])]))

    def test_duplicate_id_rejected(self):
        with self.assertRaises(ValueError):
            validate_geometry(geometry([lane(), lane(box=[50.0, 10.0, 80.0, 180.0])]))

    def test_marks_and_masks_do_not_enter_geometry(self):
        r = review_fixture()
        r['image_sha256'] = 'test-fixture'
        original = copy.deepcopy(r)
        g1 = strip_geometry(r)
        r['panels']['A']['bands'][0]['y'] = 70.5
        r['panels']['A']['exclusions'] = [{'box': [12, 20, 30, 50], 'reason': 'rough'}]
        r['panels']['A']['lanes'][0]['reviewed'] = True
        g2 = strip_geometry(r)
        self.assertEqual(g1, g2)
        self.assertEqual(detect_array(picture(), g1, CONFIG), detect_array(picture(), g2, CONFIG))
        self.assertEqual(original['panels']['A']['bands'][0]['y'], 60.5)

class AgreementTests(unittest.TestCase):

    def candidate(self, y=60.5, x=25, cid='p1'):
        return {'candidate_id': cid, 'panel_id': 'A', 'x_px': x, 'y_px': y, 'label': 'ignored'}

    def test_exact_match(self):
        r = compare_snapshot(review_fixture(), [self.candidate()], 3)
        self.assertEqual(r['summary'][-1]['n_matched'], 1)

    def test_unmarked_candidate_not_called_false_positive(self):
        r = compare_snapshot(review_fixture(), [self.candidate(100)], 3)
        self.assertEqual(r['candidate_audit'][0]['disposition'], 'candidate_without_matching_mark')
        self.assertNotIn('"false_positive":', json.dumps(r))
        self.assertNotIn('f1', json.dumps(r))

    def test_input_review_is_not_changed(self):
        r = review_fixture()
        b = copy.deepcopy(r)
        compare_snapshot(r, [self.candidate()], 3)
        self.assertEqual(r, b)

    def test_one_prediction_matches_at_most_one_mark(self):
        r = review_fixture()
        r['panels']['A']['bands'].append({'id': 'b2', 'lane_id': 'L1', 'x': 25.0, 'y': 62.0, 'certainty': 'definite'})
        c = compare_snapshot(r, [self.candidate()], 3)
        self.assertEqual(c['summary'][-1]['n_matched'], 1)

    def test_two_predictions_match_one_mark_once(self):
        c = compare_snapshot(review_fixture(), [self.candidate(), self.candidate(61, cid='p2')], 3)
        self.assertEqual(c['summary'][-1]['n_matched'], 1)

    def test_tolerance_boundary(self):
        self.assertEqual(compare_snapshot(review_fixture(), [self.candidate(63.5)], 3)['summary'][-1]['n_matched'], 1)
        self.assertEqual(compare_snapshot(review_fixture(), [self.candidate(63.51)], 3)['summary'][-1]['n_matched'], 0)

    def test_exclusions_are_reported_not_deleted(self):
        r = review_fixture()
        r['panels']['A']['exclusions'] = [{'box': [10, 90, 40, 120], 'reason': 'smear'}]
        c = compare_snapshot(r, [self.candidate(100)], 3)
        self.assertEqual(c['candidate_audit'][0]['disposition'], 'recorded_unreadable_region')
        self.assertEqual(c['summary'][-1]['n_all_candidates'], 1)
        self.assertEqual(c['summary'][-1]['n_eligible_candidates'], 0)

    def test_no_mask_diagnostic_keeps_excluded_candidate(self):
        r = review_fixture()
        r['panels']['A']['exclusions'] = [{'box': [10, 90, 40, 120], 'reason': 'smear'}]
        c = compare_snapshot(r, [self.candidate(100)], 3, 'no_masks_diagnostic')
        self.assertEqual(c['summary'][-1]['n_eligible_candidates'], 1)

    def test_uncertain_mark_is_not_confirmed(self):
        r = review_fixture()
        r['panels']['A']['bands'][0]['certainty'] = 'uncertain'
        c = compare_snapshot(r, [self.candidate()], 3)
        self.assertEqual(c['summary'][-1]['n_definite_marks'], 0)
        self.assertIsNone(c['summary'][-1]['marked_band_recovery_fraction'])

    def test_outside_lane_kept_in_audit(self):
        c = compare_snapshot(review_fixture(), [self.candidate(x=70)], 3)
        self.assertEqual(c['candidate_audit'][0]['disposition'], 'outside_supplied_lane_rectangles')

    def test_label_not_used_to_force_matching(self):
        c = self.candidate()
        c['label'] = 'not_the_sample'
        self.assertEqual(compare_snapshot(review_fixture(), [c], 3)['summary'][-1]['n_matched'], 1)

    def test_counter_conservation(self):
        p = [self.candidate(), self.candidate(100, cid='p2'), self.candidate(120, x=80, cid='p3')]
        s = compare_snapshot(review_fixture(), p, 3)['summary'][-1]
        self.assertEqual(s['n_eligible_candidates'], s['n_matched'] + s['n_unmatched_candidates'])
        self.assertEqual(s['n_definite_marks'], s['n_matched'] + s['n_unmatched_marks'])
        self.assertEqual(s['n_all_candidates'], s['n_eligible_candidates'] + s['n_not_evaluated_candidates'])

    def test_duplicate_prediction_id_rejected(self):
        with self.assertRaises(ValueError):
            compare_snapshot(review_fixture(), [self.candidate(), self.candidate()], 3)

    def test_mask_area_union_not_double_counted(self):
        r = review_fixture()
        reg = {'box': [10, 20, 40, 30], 'reason': 'smear'}
        r['panels']['A']['exclusions'] = [reg, copy.deepcopy(reg)]
        a = covered_area(r)[0]
        self.assertEqual(a['n_exclusion_pixels'], 300)
        self.assertEqual(a['n_union_mask_pixels'], 300)

    def test_malformed_mask_mode_rejected(self):
        with self.assertRaises(ValueError):
            compare_snapshot(review_fixture(), [], 3, 'wrong')

class MatchingTests(unittest.TestCase):

    def test_empty(self):
        self.assertEqual(match_1d([], [], 3), [])

    def test_invalid_tolerance(self):
        with self.assertRaises(ValueError):
            match_1d([1], [1], -1)

    def test_non_greedy_assignment(self):
        self.assertEqual(len(match_1d([0, 4], [3, 7], 3)), 2)

    def test_no_match_not_forced(self):
        self.assertEqual(match_1d([0], [10], 3), [])

    def test_unsorted_indices_preserved(self):
        m = match_1d([10, 1], [1, 10], 0)
        self.assertEqual(sorted(((i, j) for i, j, e in m)), [(0, 1), (1, 0)])

    def test_against_bruteforce_small_cases(self):
        rng = random.Random(20260921)
        for _ in range(40):
            a = [rng.randrange(10) for _ in range(3)]
            b = [rng.randrange(10) for _ in range(3)]
            tol = 2
            options = [(0, 0.0)]
            for k in range(1, 4):
                for ii in itertools.combinations(range(3), k):
                    for jj in itertools.permutations(range(3), k):
                        costs = [abs(a[i] - b[j]) for i, j in zip(ii, jj)]
                        if max(costs) <= tol:
                            options.append((k, sum(costs)))
            best = min(options, key=lambda z: (-z[0], z[1]))
            m = match_1d(a, b, tol)
            self.assertEqual((len(m), sum((e for _, _, e in m))), best)

class IOTests(unittest.TestCase):

    def test_json_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.json'
            write_json(p, {'label': 'sample', 'value': 1})
            self.assertEqual(read_json(p), {'label': 'sample', 'value': 1})

    def test_nonfinite_json_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                write_json(Path(d) / 'a.json', {'a': float('nan')})

    def test_table_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.tsv'
            write_table(p, [{'a': 'x', 'b': 2}])
            self.assertEqual(read_table(p), [{'a': 'x', 'b': '2'}])

    def test_image_hash_mismatch(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.png'
            Image.fromarray(picture().astype('uint8')).save(p)
            with self.assertRaises(ValueError):
                detect_file(p, geometry(), CONFIG)

    def test_image_file_coordinate_roundtrip(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.png'
            Image.fromarray(picture().astype('uint8')).save(p)
            g = geometry()
            g['image_sha256'] = digest(p)
            self.assertEqual([r['y_px'] for r in detect_file(p, g, CONFIG)[0]['candidates']], [60.5, 110.5])

    def test_transparency_rejected(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.png'
            Image.new('RGBA', (10, 10), (0, 0, 0, 0)).save(p)
            with self.assertRaises(ValueError):
                load_gray(p)
if __name__ == '__main__':
    unittest.main(verbosity=2)
