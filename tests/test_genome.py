from __future__ import annotations
import itertools, json, tempfile, unittest, sys
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from legacy_marker_rescue.genome.stats import jaccard, correlation, permutation_family, make_permutations
from legacy_marker_rescue.genome.sequence import SiteScanner, revcomp, pair_positions, cluster_sizes, pooled_matrix, training_frozen_matrix, fixed_log_matrix, select_sizes, scan_assembly
from legacy_marker_rescue.genome.common import metadata_groups, accession, write_tsv, read_tsv
from legacy_marker_rescue.genome.annotations import match_bands, evaluate
import tests.sequence_reference as old

class IdentityTests(unittest.TestCase):

    def test_filename_accession(self):
        self.assertEqual(accession('example/GCA_000149035.1_M1_001_genomic.fna'), 'GCA_000149035.1')

    def test_ambiguous_identifier_fails(self):
        with self.assertRaises(ValueError):
            accession('GCA_000149035.1+GCA_001951205.1')

    def test_strain_across_distinct_biosamples_is_grouped(self):
        inv = [{'assembly': 'GCA_2.1', 'biosample': 'SAM2', 'strain': 'M1.001', 'isolate': ''}, {'assembly': 'GCA_1.1', 'biosample': 'SAM1', 'strain': '', 'isolate': 'M1.001'}, {'assembly': 'GCA_3.1', 'biosample': 'SAM3', 'strain': '', 'isolate': ''}]
        rows, reps = metadata_groups(inv)
        self.assertEqual(len(reps), 2)
        self.assertEqual(rows[0]['analysis_group'], rows[1]['analysis_group'])
        self.assertIn(1, reps)

    def test_unknown_labels_do_not_merge(self):
        inv = [{'assembly': f'GCA_{i}.1', 'strain': 'unknown', 'isolate': '', 'biosample': ''} for i in range(4)]
        self.assertEqual(len(metadata_groups(inv)[1]), 4)

class MatrixTests(unittest.TestCase):

    def test_jaccard_ignores_shared_absence(self):
        b = np.array([[1, 0, 1, 0], [1, 1, 0, 0], [0, 0, 0, 0]])
        d = jaccard(b)
        self.assertAlmostEqual(d[0, 1], 2 / 3)
        self.assertEqual(d[0, 2], 1)
        self.assertTrue(np.allclose(d, d.T, equal_nan=True))

    def test_two_empty_profiles_are_undefined_not_identical(self):
        d = jaccard(np.zeros((2, 3), dtype=int))
        self.assertTrue(np.isnan(d[0, 1]))

    def test_nonbinary_fails(self):
        with self.assertRaises(ValueError):
            jaccard([[1, 2], [0, 1]])

    def test_ties_match_scipy(self):
        x = np.array([1, 1, 2, 4, 5, 5])
        y = np.array([5, 3, 3, 2, 1, 1])
        self.assertAlmostEqual(correlation(x, y), spearmanr(x, y).statistic, places=13)
        self.assertAlmostEqual(correlation(x, y, 'pearson'), pearsonr(x, y).statistic, places=13)

    def test_constant_undefined(self):
        self.assertTrue(np.isnan(correlation([1, 1, 1], [1, 2, 3])))

    def test_label_permutation_matches_direct_reference(self):
        rng = np.random.default_rng(17)
        a = rng.random((6, 6))
        a = (a + a.T) / 2
        np.fill_diagonal(a, 0)
        b = rng.random((6, 6))
        b = (b + b.T) / 2
        np.fill_diagonal(b, 0)
        other = b * 0.5 + a * 0.5
        perms = make_permutations(6, 25, 31)
        result = permutation_family(np.array([a, other]), b, perms)
        ix = np.triu_indices(6, 1)
        for i, p in enumerate(perms):
            direct = b[np.ix_(p, p)]
            self.assertAlmostEqual(result['spearman']['permuted'][i, 0], spearmanr(a[ix], direct[ix]).statistic, places=13)
            self.assertAlmostEqual(result['pearson']['permuted'][i, 0], pearsonr(a[ix], direct[ix]).statistic, places=13)

    def test_adjusted_p_not_smaller(self):
        rng = np.random.default_rng(2)
        z = rng.random((4, 7, 7))
        z = (z + z.transpose(0, 2, 1)) / 2
        for a in z:
            np.fill_diagonal(a, 0)
        r = permutation_family(z[:3], z[3], make_permutations(7, 99, 88))
        for x in r.values():
            self.assertTrue(np.all(x['p_search_adjusted'] >= x['p_two_sided'] - 1e-12))

    def test_permutation_identity_and_seed(self):
        p = make_permutations(6, 20, 10)
        self.assertTrue(np.array_equal(p, make_permutations(6, 20, 10)))
        self.assertTrue(all((sorted(x) == list(range(6)) for x in p)))

    def test_block_permutation_preserves_blocks(self):
        blocks = ['a', 'a', 'b', 'b', 'c']
        p = make_permutations(5, 50, 72, blocks)
        for perm in p:
            self.assertEqual(perm[4], 4)
            self.assertTrue(all((blocks[i] == blocks[j] for i, j in enumerate(perm))))

class SequenceTests(unittest.TestCase):
    primers = ['TCGGGAGGGT', 'ACGCTCAAAC', 'AGAATCGGGG']

    def test_vectorized_matches_archived_regex_and_scalar(self):
        rng = np.random.default_rng(281)
        seq = ''.join(rng.choice(list('ACGTN'), 5000))
        for primer in self.primers:
            seq += primer + 'N' + revcomp(primer) + 'ACGT'
        scanner = SiteScanner(self.primers, 1, 37)
        self.assertEqual(scanner.scan(seq), scanner.reference_scan(seq))
        for i, p in enumerate(self.primers):
            for side, q in enumerate((p, revcomp(p))):
                vs = old.generate_variants(q, 1)
                pat = old.compile_variant_regex(vs)
                expected = [(pos, mm) for pos, mm, s in old.find_sites(seq, pat, vs)]
                self.assertEqual(scanner.scan(seq)[i][side], expected)

    def test_one_mismatch_detected_two_not(self):
        p = 'TCGGGAGGGT'
        scanner = SiteScanner([p], 1, 13)
        one = 'A' + p[1:]
        two = 'AA' + p[2:]
        self.assertEqual(scanner.scan(one)[0][0], [(0, 1)])
        self.assertEqual(scanner.scan(two)[0][0], [])

    def test_chunks_do_not_duplicate_hits(self):
        p = self.primers[0]
        seq = ('N' + p) * 20
        for chunk in (1, 9, 10, 11, 20, 100):
            s = SiteScanner([p], 1, chunk)
            self.assertEqual(s.scan(seq), s.reference_scan(seq))

    def test_ambiguous_sequence_not_treated_as_A(self):
        self.assertEqual(SiteScanner(['AAAAAAAAAA']).scan('NNNNNNNNNN')[0][0], [])

    def test_inward_and_outward_have_distinct_geometry(self):
        p = self.primers[0]
        scanner = SiteScanner([p], 0)
        inward = p + 'N' * 100 + revcomp(p)
        outward = revcomp(p) + 'N' * 100 + p
        f, r = scanner.scan(inward)[0]
        self.assertEqual(len(list(pair_positions(f, r))), 1)
        f, r = scanner.scan(outward)[0]
        self.assertEqual(len(list(pair_positions(f, r))), 0)
        self.assertEqual(len(list(pair_positions(r, f))), 1)

    def test_reverse_complement_preserves_inward_product(self):
        p = self.primers[0]
        scanner = SiteScanner([p], 0)
        seq = 'N' * 3 + p + 'A' * 115 + revcomp(p) + 'N' * 7
        f, r = scanner.scan(seq)[0]
        a = list(pair_positions(f, r, min_bp=100))
        f, r = scanner.scan(revcomp(seq))[0]
        b = list(pair_positions(f, r, min_bp=100))
        self.assertEqual([x[2] for x in a], [x[2] for x in b])

    def test_pair_enumeration_matches_archived_both_orientation(self):
        f = [(0, 0, 'x'), (30, 1, 'x'), (160, 0, 'x')]
        r = [(15, 1, 'y'), (120, 0, 'y'), (170, 1, 'y')]
        expected = old.pair_sites(f, r, 10, 100, 180)
        observed = []
        for label, left, right in [('F_to_RC', f, r), ('RC_to_F', r, f)]:
            for start, end, size, lmm, rmm in pair_positions([(p, m) for p, m, s in left], [(p, m) for p, m, s in right], max_bp=180, min_bp=100):
                observed.append((start, end, size, label, lmm, rmm))
        self.assertEqual(sorted((x[:6] for x in expected)), sorted(observed))

    def test_minimum_length_applied(self):
        pairs = list(pair_positions([(0, 0)], [(1, 0), (89, 0), (90, 0), (100, 0)], min_bp=100, max_bp=110))
        self.assertEqual([x[2] for x in pairs], [100, 110])

    def test_cache_checked_by_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fa = p / 'genomic.fna'
            fa.write_text('>c\n' + self.primers[0] + 'N' * 100 + revcomp(self.primers[0]) + '\n')
            cfg = {'max_amplicon_bp': 3000, 'max_mismatches_per_site': 1, 'max_sites_per_primer_contig': 20000, 'max_amplicons_per_set_assembly': 8000, 'reference_scan_bases_per_assembly': 20000}
            scanner = SiteScanner([self.primers[0]], 1, 300)
            a, info, cached = scan_assembly(fa, scanner, ['test'], cfg, p / 'cache', lambda x: None)
            b, other, cached2 = scan_assembly(fa, scanner, ['test'], cfg, p / 'cache', lambda x: None)
            self.assertFalse(cached)
            self.assertTrue(cached2)
            self.assertTrue(np.array_equal(a, b))
            self.assertEqual(info['total_bases'], 120)

class BinningTests(unittest.TestCase):

    def test_greedy_clustering_matches_archived(self):
        sizes = [100, 100, 101, 110, 118, 200, 209, 210, 250, 251, 990]
        for tol in (0.03, 0.05, 0.1):
            mapping, centers, bins = cluster_sizes(sizes, tol)
            ref = old.cluster_band_sizes(sizes, tol)
            for s, i in mapping.items():
                self.assertEqual(ref[s], f'B{i + 1:03d}_{int(centers[i])}bp')

    def test_heldout_sizes_cannot_change_training_features(self):
        records = [{'primer1': [100, 200]}, {'primer1': [101, 220]}, {'primer1': [999]}]
        x, first = training_frozen_matrix(records, [0, 1], 0.05)
        records[2]['primer1'] = [301, 350, 1000, 1001, 9999]
        y, second = training_frozen_matrix(records, [0, 1], 0.05)
        self.assertTrue(np.array_equal(x[:2], y[:2]))
        self.assertEqual(x.shape[1], y.shape[1])

    def test_unassigned_test_bands_recorded(self):
        x, audit = training_frozen_matrix([{'p': [100]}, {'p': [200]}, {'p': [700]}], [0, 1], 0.05)
        self.assertEqual(audit[-1]['unassigned_products'], 1)

    def test_fixed_bins_training_distances_invariant_to_test_changes(self):
        a = [{'p': [100, 202]}, {'p': [101, 250]}, {'p': [1000]}]
        d = jaccard(fixed_log_matrix(a, 0.05))
        a[2] = {'p': [700, 1800, 2200]}
        e = jaccard(fixed_log_matrix(a, 0.05))
        self.assertAlmostEqual(d[0, 1], e[0, 1])

    def test_orientation_filter_is_explicit(self):
        arr = np.array([[0, 0, 0, 150, 150, 0, 0, 0], [0, 0, 200, 350, 150, 1, 0, 0], [0, 0, 400, 430, 30, 0, 0, 0]])
        both = select_sizes([arr], {'p': 0}, 'both', 100, 3000)
        inward = select_sizes([arr], {'p': 0}, 'inward', 100, 3000)
        self.assertEqual(both[0]['p'], [150, 150])
        self.assertEqual(inward[0]['p'], [150])

class AnnotationTests(unittest.TestCase):

    def test_one_prediction_cannot_match_two_bands(self):
        self.assertEqual(len(match_bands([10, 11], [10.5], 2)), 1)

    def test_assignment_maximizes_cardinality(self):
        matches = match_bands([10, 13], [11, 8], 3)
        self.assertEqual(len(matches), 2)

    def test_no_overlap_is_zero_not_forced(self):
        self.assertEqual(match_bands([10], [30], 2), [])

    def test_empty_reference_not_invented(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            for name in ('lanes', 'ref', 'pred'):
                write_tsv(p / (name + '.tsv'), [], ['image_id', 'lane_id'])
            with self.assertRaises(ValueError):
                evaluate(p / 'lanes.tsv', p / 'ref.tsv', p / 'pred.tsv', 5, p / 'out')

    def test_wrong_coordinate_frame_fails(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            lane = {'image_id': 'a', 'lane_id': '1', 'image_sha256': 'a' * 64, 'reviewer': 'r', 'lane_type': 'sample', 'review_status': 'complete'}
            ref = {'image_id': 'a', 'lane_id': '1', 'image_sha256': 'b' * 64, 'band_id': 'b1', 'status': 'confirmed', 'y_px': 100}
            write_tsv(p / 'lane.tsv', [lane])
            write_tsv(p / 'ref.tsv', [ref])
            write_tsv(p / 'pred.tsv', [], list(ref))
            with self.assertRaises(ValueError):
                evaluate(p / 'lane.tsv', p / 'ref.tsv', p / 'pred.tsv', 5, p / 'out')
if __name__ == '__main__':
    unittest.main()
