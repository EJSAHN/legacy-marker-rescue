"""Synthetic software tests only; no external scientific scores are fabricated."""
from __future__ import annotations
import copy
import hashlib
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import numpy as np
from PIL import Image
from legacy_marker_rescue.benchmark.acquire import validate_member, verify_archive, pair_members, retrieve
from legacy_marker_rescue.benchmark.imaging import decode_image, decode_mask, component_reference, rescale_height
from legacy_marker_rescue.benchmark.lanes import image_envelopes, reference_x_envelopes, make_geometry, extract_candidates, runs
from legacy_marker_rescue.benchmark.scoring import overlap_edges, assign_components, evaluate, summaries
from legacy_marker_rescue.benchmark.analysis import auto_detect, process_pair
ROOT = Path(__file__).resolve().parents[1]
B = json.loads((ROOT / 'config/benchmark.json').read_text())
C = json.loads((ROOT / 'config/detection.json').read_text())

def encode(arr, fmt='TIFF'):
    f = io.BytesIO()
    Image.fromarray(arr).save(f, format=fmt)
    return f.getvalue()

def fixture():
    im = np.full((200, 160), 25.0, dtype=float)
    m = np.zeros_like(im, dtype=np.uint8)
    for x0, x1 in [(20, 48), (85, 113)]:
        for y in [45, 85, 140]:
            im[:, x0:x1] += 180 * np.exp(-0.5 * ((np.arange(200) - y) / 2) ** 2)[:, None]
            m[y - 3:y + 4, x0:x1] = 1
    return (im.astype('uint8'), m)

def samplezip(path):
    a, m = fixture()
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr('gels/images/one.tif', encode(a))
        z.writestr('gels/masks/one.tif', encode(m))
        z.writestr('gels/test_images/two.png', encode(a, 'PNG'))
        z.writestr('gels/test_masks/two.tif', encode(m))

class SourceArchiveTests(unittest.TestCase):

    def test_pairing_all_partitions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            samplezip(p)
            pairs, issues = pair_members(p)
            self.assertEqual(len(pairs), 2)
            self.assertFalse(issues)
            self.assertEqual({x['source_partition'] for x in pairs}, {'train', 'test'})

    def test_png_pairs_to_tif(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            samplezip(p)
            pairs, _ = pair_members(p)
            self.assertTrue(any((x['image_member'].endswith('.png') and x['mask_member'].endswith('.tif') for x in pairs)))

    def test_missing_pair_kept_as_issue(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            with zipfile.ZipFile(p, 'w') as z:
                z.writestr('g/images/a.tif', b'x')
            pairs, issues = pair_members(p)
            self.assertFalse(pairs)
            self.assertEqual(issues[0]['status'], 'IMAGE_WITHOUT_MASK')

    def test_mask_without_image_is_issue(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            with zipfile.ZipFile(p, 'w') as z:
                z.writestr('g/masks/a.tif', b'x')
            pairs, issues = pair_members(p)
            self.assertEqual(issues[0]['status'], 'MASK_WITHOUT_IMAGE')

    def test_unknown_directory_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            with zipfile.ZipFile(p, 'w') as z:
                z.writestr('g/unknown/a.tif', b'x')
            pairs, issues = pair_members(p)
            self.assertEqual(issues[0]['status'], 'UNRECOGNIZED_IMAGE_OR_MASK_DIRECTORY')

    def test_ambiguous_pair_fails(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'data.zip'
            with zipfile.ZipFile(p, 'w') as z:
                z.writestr('g/images/a.tif', b'x')
                z.writestr('g/images/a.png', b'x')
            with self.assertRaises(ValueError):
                pair_members(p)

    def test_unsafe_names(self):
        for name in ['../a', '/a', 'drive:/a', 'a\\b', 'x/../../a', 'a:b', 'a\x00b']:
            with self.subTest(name=name), self.assertRaises(ValueError):
                validate_member(name)

    def test_safe_unicode_name(self):
        validate_member('g/images/gel alpha β.tif')

    def test_archive_checksum(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.zip'
            samplezip(p)
            md5 = hashlib.md5(p.read_bytes()).hexdigest()
            self.assertEqual(verify_archive(p, md5)['md5'], md5)
            with self.assertRaises(ValueError):
                verify_archive(p, '0' * 32)

    def test_verified_cache_no_network(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'x.zip'
            samplezip(p)
            c = {'file': 'x.zip', 'md5': hashlib.md5(p.read_bytes()).hexdigest()}
            with patch('legacy_marker_rescue.benchmark.acquire.urlopen', side_effect=AssertionError('network')):
                _, r = retrieve(c, 'test', Path(d))
                self.assertEqual(r['acquisition'], 'verified_cache')

    def test_downloads_duplicate_name(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            dl = p / 'Downloads'
            dl.mkdir()
            q = dl / 'x (2).zip'
            samplezip(q)
            c = {'file': 'x.zip', 'md5': hashlib.md5(q.read_bytes()).hexdigest()}
            with patch('legacy_marker_rescue.benchmark.acquire.urlopen', side_effect=AssertionError('network')):
                _, r = retrieve(c, 'test', p / 'cache', dl)
                self.assertEqual(r['acquisition'], 'verified_downloads_copy')

    def test_network_response_verified(self):

        class Response(io.BytesIO):
            headers = {'Content-Type': 'application/zip'}
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            samplezip(p / 'a.zip')
            data = (p / 'a.zip').read_bytes()
            c = {'file': 'x.zip', 'md5': hashlib.md5(data).hexdigest()}
            with patch('legacy_marker_rescue.benchmark.acquire.urlopen', return_value=Response(data)):
                dest, r = retrieve(c, 'test', p / 'cache')
                self.assertEqual(dest.read_bytes(), data)
                self.assertEqual(r['acquisition'], 'https_download')

    def test_network_failure_is_not_success(self):
        with tempfile.TemporaryDirectory() as d, patch('legacy_marker_rescue.benchmark.acquire.urlopen', side_effect=OSError('offline')), patch('legacy_marker_rescue.benchmark.acquire.time.sleep'):
            with self.assertRaises(RuntimeError):
                retrieve({'file': 'x.zip', 'md5': '0' * 32}, 'test', Path(d))

class DecodeTests(unittest.TestCase):

    def test_8bit(self):
        a, _ = fixture()
        g, info = decode_image(encode(a), B['preprocessing'])
        self.assertEqual(g.shape, a.shape)
        self.assertFalse(info['inverted_from_image_median'])

    def test_16bit(self):
        a, _ = fixture()
        g, info = decode_image(encode(a.astype('uint16') * 128), B['preprocessing'])
        self.assertEqual(g.max(), 255.0)
        self.assertIn(info['original_mode'], ['I;16', 'I;16L', 'I'])

    def test_dark_auto_polarity(self):
        a, _ = fixture()
        g, info = decode_image(encode(255 - a), B['preprocessing'])
        h, _ = decode_image(encode(a), B['preprocessing'])
        self.assertTrue(info['inverted_from_image_median'])
        np.testing.assert_allclose(g, h, atol=1e-12)

    def test_flat_not_invented_signal(self):
        g, info = decode_image(encode(np.full((30, 40), 200, dtype='uint8')), B['preprocessing'])
        self.assertTrue(info['constant_image'])
        self.assertEqual(g.sum(), 0.0)

    def test_percentile_fallback_not_division_zero(self):
        a = np.zeros((100, 100), dtype='uint8')
        a[1, 1] = 255
        g, info = decode_image(encode(a), B['preprocessing'])
        self.assertTrue(info['normalization_fallback_minmax'])
        self.assertEqual(g[1, 1], 255)

    def test_nonfinite_image_rejected(self):
        a = np.zeros((10, 10), dtype='float32')
        a[1, 1] = np.nan
        with self.assertRaises(ValueError):
            decode_image(encode(a), B['preprocessing'])

    def test_transparent_image_rejected(self):
        with self.assertRaises(ValueError):
            decode_image(encode(np.zeros((10, 10, 4), dtype='uint8'), 'PNG'), B['preprocessing'])

    def test_multiframe_not_first_frame(self):
        f = io.BytesIO()
        a = Image.fromarray(np.zeros((10, 10), dtype='uint8'))
        a.save(f, format='TIFF', save_all=True, append_images=[a])
        with self.assertRaises(ValueError):
            decode_image(f.getvalue(), B['preprocessing'])

    def test_mask01(self):
        _, m = fixture()
        a, info = decode_mask(encode(m), m.shape)
        np.testing.assert_equal(a, m != 0)
        self.assertEqual(info['foreground_value'], 1)

    def test_mask255(self):
        _, m = fixture()
        a, info = decode_mask(encode(m * 255), m.shape)
        np.testing.assert_equal(a, m != 0)
        self.assertEqual(info['foreground_value'], 255)

    def test_palette_mask_uses_indices(self):
        a = np.array([[0, 1], [1, 0]], dtype='uint8')
        im = Image.fromarray(a).convert('P')
        im.putpalette([255, 255, 255, 100, 20, 20] + [0] * 762)
        f = io.BytesIO()
        im.save(f, format='PNG')
        m, _ = decode_mask(f.getvalue(), (2, 2))
        np.testing.assert_array_equal(m, a)

    def test_wrong_mask_dimensions_fail(self):
        with self.assertRaises(ValueError):
            decode_mask(encode(np.zeros((10, 10), dtype='uint8')), (11, 10))

    def test_unknown_mask_encoding_fail(self):
        with self.assertRaises(ValueError):
            decode_mask(encode(np.array([[0, 2]], dtype='uint8')), (1, 2))

    def test_rgb_mask_fail(self):
        with self.assertRaises(ValueError):
            decode_mask(encode(np.zeros((10, 10, 3), dtype='uint8')), (10, 10))

    def test_mask_empty_allowed(self):
        m, _ = decode_mask(encode(np.zeros((10, 10), dtype='uint8')), (10, 10))
        labels, refs = component_reference(m)
        self.assertEqual(refs, [])

    def test_single_pixel_not_pruned(self):
        m = np.zeros((10, 10), bool)
        m[5, 5] = True
        _, r = component_reference(m)
        self.assertEqual(len(r), 1)
        self.assertEqual(r[0]['area_px'], 1)

    def test_touching_pixels_one_component(self):
        m = np.zeros((10, 10), bool)
        m[5, 5] = m[6, 6] = True
        _, r = component_reference(m)
        self.assertEqual(len(r), 1)

    def test_edge_component_counted(self):
        m = np.zeros((10, 10), bool)
        m[0, 0] = True
        _, r = component_reference(m)
        self.assertTrue(r[0]['touches_image_edge'])

    def test_rescale_coordinate_factors(self):
        a = np.zeros((100, 250))
        g, sx, sy = rescale_height(a, 256)
        self.assertEqual(g.shape, (256, 640))
        self.assertEqual((sx, sy), (2.56, 2.56))

class LaneAdapterTests(unittest.TestCase):

    def test_image_only_has_no_mask_argument(self):
        import inspect
        self.assertNotIn('mask', inspect.signature(auto_detect).parameters)

    def test_detect_function_does_not_read_reference(self):
        a, _ = fixture()
        gray, _ = decode_image(encode(a), B['preprocessing'])
        with patch('legacy_marker_rescue.benchmark.analysis.decode_mask', side_effect=AssertionError('mask read')):
            r = auto_detect(gray, 'image_lanes_native', B, C, 'fixture')
            self.assertFalse(r['reference_mask_used'])

    def test_no_lanes_flat(self):
        self.assertFalse(image_envelopes(np.zeros((100, 200)), B['image_lane_proposal'])[0])

    def test_two_synthetic_lanes(self):
        a, _ = fixture()
        boxes, _ = image_envelopes(a.astype(float), B['image_lane_proposal'])
        self.assertEqual(len(boxes), 2)

    def test_proposals_nonoverlap_fullheight(self):
        a, _ = fixture()
        boxes, _ = image_envelopes(a.astype(float), B['image_lane_proposal'])
        g = make_geometry(a, boxes, 'test')
        lanes = g['panels']['gel']['lanes']
        for r in lanes:
            self.assertEqual(r['box'][1], 0)
            self.assertEqual(r['box'][3], 200)
        for l, r in zip(lanes, lanes[1:]):
            self.assertLessEqual(l['box'][2], r['box'][0])

    def test_native_reference_position_not_input(self):
        a, _ = fixture()
        first = auto_detect(a.astype(float), 'image_lanes_native', B, C, 'test')
        second = auto_detect(a.astype(float), 'image_lanes_native', B, C, 'test')
        self.assertEqual(first, second)

    def test_reference_x_geometry_no_y_windows(self):
        _, m = fixture()
        r = reference_x_envelopes(m)
        g = make_geometry(m, r, 't')
        self.assertEqual(len(r), 2)
        self.assertTrue(all((x['box'][1] == 0 and x['box'][3] == m.shape[0] for x in g['panels']['gel']['lanes'])))

    def test_reference_y_changes_keep_envelopes(self):
        _, m = fixture()
        self.assertEqual(reference_x_envelopes(m), reference_x_envelopes(np.roll(m, 10, axis=0)))

    def test_oracle_merges_overlapping_x_not_guessed_lanes(self):
        m = np.zeros((100, 100), bool)
        m[10:12, 20:50] = 1
        m[70:74, 45:75] = 1
        self.assertEqual(reference_x_envelopes(m), [(20, 75)])

    def test_small_reference_roi_status_not_discarded(self):
        a = np.zeros((50, 50))
        p, rows = extract_candidates(a, [(3, 4)], C, 'test')
        self.assertFalse(p)
        self.assertEqual(rows[0]['status'], 'ROI_TOO_SMALL')

    def test_bright_blank_uniform_no_bands(self):
        a = np.full((200, 100), 255.0)
        p, rows = extract_candidates(a, [(10, 40)], C, 'test')
        self.assertFalse(p)

    def test_native_candidate_positions(self):
        a, _ = fixture()
        r = auto_detect(a.astype(float), 'image_lanes_native', B, C, 'test')
        self.assertEqual(len(r['candidates']), 6)

    def test_scale_output_back_to_original_bounds(self):
        a, _ = fixture()
        r = auto_detect(a.astype(float), 'image_lanes_height256', B, C, 'test')
        self.assertTrue(all((0 <= p['y_px'] < 200 and 0 <= p['x_px'] < 160 for p in r['candidates'])))

class MatchingTests(unittest.TestCase):

    def setUp(self):
        self.mask = np.zeros((50, 70), bool)
        self.mask[10:15, 10:25] = 1
        self.mask[30:35, 40:55] = 1
        self.lab, self.refs = component_reference(self.mask)

    def test_two_matches(self):
        r = evaluate([{'x_px': 18.0, 'y_px': 12.0}, {'x_px': 48.0, 'y_px': 32.0}], self.lab, self.refs, 3)
        self.assertEqual(r['matched'], 2)

    def test_single_prediction_one_reference(self):
        self.assertEqual(len(assign_components([{1: 1.0, 2: 1.0}], 2)), 1)

    def test_two_predictions_one_reference(self):
        self.assertEqual(len(assign_components([{1: 0.0}, {1: 0.0}], 1)), 1)

    def test_max_cardinality_not_greedy(self):
        a = assign_components([{1: 0.0, 2: 2.0}, {1: 0.0}], 2)
        self.assertEqual(len(a), 2)

    def test_no_match_not_forced(self):
        self.assertEqual(assign_components([{}, {}], 2), [])

    def test_false_extra_kept_in_denominator(self):
        p = [{'x_px': 18.0, 'y_px': 12.0}, {'x_px': 60.0, 'y_px': 40.0}]
        r = evaluate(p, self.lab, self.refs, 0)
        self.assertEqual(r['predictions'], 2)
        self.assertEqual(r['unmatched_predictions'], 1)

    def test_missing_lane_refs_still_denominator(self):
        r = evaluate([{'x_px': 18.0, 'y_px': 12.0}], self.lab, self.refs, 0)
        self.assertEqual(r['unmatched_reference_components'], 1)

    def test_tolerance_is_native_euclidean(self):
        p = [{'x_px': 18.5, 'y_px': 17.5}]
        self.assertEqual(evaluate(p, self.lab, self.refs, 3)['matched'], 1)
        self.assertEqual(evaluate(p, self.lab, self.refs, 2)['matched'], 0)

    def test_zero_tolerance_inside_area_not_centroid(self):
        self.assertEqual(evaluate([{'x_px': 10.1, 'y_px': 10.1}], self.lab, self.refs, 0)['matched'], 1)

    def test_empty_predictions_not_perfect(self):
        r = evaluate([], self.lab, self.refs, 3)
        self.assertEqual(r['recall_vs_components'], 0.0)
        self.assertEqual(r['f1_vs_components'], 0.0)
        self.assertIsNone(r['precision_vs_components'])

    def test_empty_truth_not_perfect(self):
        r = evaluate([], np.zeros((10, 10), int), [], 3)
        self.assertIsNone(r['f1_vs_components'])

    def test_invalid_prediction_fails(self):
        with self.assertRaises(ValueError):
            evaluate([{'x_px': -1, 'y_px': 4}], self.lab, self.refs, 3)

    def test_no_mask_mutation(self):
        before = self.lab.copy()
        evaluate([{'x_px': 11, 'y_px': 12}], self.lab, self.refs, 3)
        np.testing.assert_equal(self.lab, before)

    def test_f1_formula(self):
        r = evaluate([{'x_px': 18.0, 'y_px': 12.0}], self.lab, self.refs, 0)
        self.assertAlmostEqual(r['f1_vs_components'], 2 / 3)

    def test_exhaustive_small_graphs(self):
        for bits in itertools.product([0, 1], repeat=6):
            edges = [{j + 1: float((i + j) % 3) for j in range(3) if bits[i * 3 + j]} for i in range(2)]
            best = 0
            for a, b in itertools.product(range(4), repeat=2):
                if a and a not in edges[0] or (b and b not in edges[1]) or (a and a == b):
                    continue
                best = max(best, bool(a) + bool(b))
            self.assertEqual(len(assign_components(edges, 3)), best)

class IntegrationTests(unittest.TestCase):

    def test_synthetic_pair_whole_pipeline(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            samplezip(p / 'input.zip')
            with zipfile.ZipFile(p / 'input.zip') as z:
                pairs, _ = pair_members(p / 'input.zip')
                metrics, row = process_pair(z, pairs[0], 'test', B, C, p / 'result')
            self.assertEqual(row['status'], 'COMPLETE')
            self.assertEqual(len(metrics), 9)
            self.assertTrue(row['provenance']['image_only_predictions_committed_before_mask_read'])
            self.assertTrue(all((r['reference_components'] == 6 for r in metrics)))

    def test_second_mask_cannot_change_primary_predictions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            a, m = fixture()
            files = []
            for i, mask in enumerate([m, np.roll(m, 19, axis=0)]):
                zpath = p / f'{i}.zip'
                with zipfile.ZipFile(zpath, 'w') as z:
                    z.writestr('x/images/a.tif', encode(a))
                    z.writestr('x/masks/a.tif', encode(mask))
                with zipfile.ZipFile(zpath) as z:
                    pair = pair_members(zpath)[0][0]
                    process_pair(z, pair, 'test', B, C, p / f'out{i}')
                files.append(p / f'out{i}' / 'images' / pair['id'] / 'image_lanes_native_predictions.json')
            self.assertEqual(files[0].read_bytes(), files[1].read_bytes())

    def test_bad_mask_records_failure_and_primary_predictions(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            a, m = fixture()
            with zipfile.ZipFile(p / 'a.zip', 'w') as z:
                z.writestr('x/images/a.tif', encode(a))
                z.writestr('x/masks/a.tif', encode(m * 2))
            with zipfile.ZipFile(p / 'a.zip') as z:
                pair = pair_members(p / 'a.zip')[0][0]
                metrics, row = process_pair(z, pair, 'test', B, C, p / 'out')
            self.assertEqual(row['status'], 'FAILED')
            self.assertFalse(metrics)
            self.assertTrue((p / 'out' / 'images' / pair['id'] / 'image_lanes_native_predictions.json').exists())

    def test_summary_separates_methods(self):
        rows = []
        for method in ['image_lanes_native', 'reference_x_envelopes_native']:
            rows.append({'collection': 'x', 'method': method, 'tolerance_native_px': 3, 'predictions': 2, 'reference_components': 3, 'matched': 1, 'f1_vs_components': 0.4, 'raw_image_sha256': 'x'})
        r = summaries(rows)
        self.assertEqual(len(r), 2)
        self.assertFalse(r[0]['confidence_intervals_calculated'])
if __name__ == '__main__':
    unittest.main(verbosity=2)
