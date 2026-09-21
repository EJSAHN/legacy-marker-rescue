from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest
import numpy as np
from legacy_marker_rescue.context.annotation import GeneIndex, parse_gff, choose_gff
from legacy_marker_rescue.context.metrics import sequence_metrics, interval_metrics, sample_start, iteration_rng, aggregate_comparisons, summarize_intervals
from legacy_marker_rescue.context.io_utils import read_json, write_json, read_tsv, write_tsv, fasta_records, sha256, checked_relative
from legacy_marker_rescue.context.binning import cluster_sizes
from legacy_marker_rescue.context.analysis import sequence_site_audit, cache_is_valid

class GeneQueryTests(unittest.TestCase):

    def test_half_open_boundaries(self):
        idx = GeneIndex([(10, 20, 'g')])
        self.assertEqual(idx.query(0, 10, 0)['overlaps_annotated_gene'], 0)
        self.assertEqual(idx.query(20, 25, 0)['overlaps_annotated_gene'], 0)
        self.assertEqual(idx.query(19, 20, 0)['overlaps_annotated_gene'], 1)

    def test_nested_long_gene_not_lost(self):
        intervals = [(0, 10000, 'long')] + [(i, i + 1, str(i)) for i in range(100, 5000, 2)]
        self.assertIn('long', GeneIndex(intervals).overlapping_ids(9000, 9001))

    def test_exact_query_matches_bruteforce(self):
        rng = random.Random(17)
        intervals = [(s, s + rng.randrange(1, 400), 'g' + str(i)) for i, s in enumerate((rng.randrange(2000) for _ in range(500)))]
        idx = GeneIndex(intervals)
        for _ in range(150):
            s = rng.randrange(2400)
            e = s + rng.randrange(1, 100)
            flank = rng.randrange(100)
            direct = {i for gs, ge, i in intervals if gs < e and ge > s}
            near = {i for gs, ge, i in intervals if gs < e + flank and ge > max(0, s - flank)}
            distance = min((0 if gs < e and ge > s else s - ge if ge <= s else gs - e for gs, ge, _ in intervals))
            q = idx.query(s, e, flank)
            self.assertEqual(idx.overlapping_ids(s, e), direct)
            self.assertEqual(q['nearest_gene_distance_bp'], distance)
            self.assertEqual(q['gene_count_within_flank'], len(near))

    def test_multipart_feature_counted_once(self):
        q = GeneIndex([(1, 10, 'same'), (20, 30, 'same')]).query(5, 25, 0)
        self.assertEqual(q['gene_count_within_flank'], 1)

    def test_empty_gene_set_is_distinct_from_no_annotation(self):
        self.assertEqual(GeneIndex([]).query(0, 10, 0)['overlaps_annotated_gene'], 0)
        self.assertIsNone(GeneIndex([]).query(0, 10, 0)['nearest_gene_distance_bp'])

    def test_invalid_coordinates_rejected(self):
        with self.assertRaises(ValueError):
            GeneIndex([(9, 2, 'bad')])
        with self.assertRaises(ValueError):
            GeneIndex([]).query(10, 5, 0)

class AnnotationTests(unittest.TestCase):

    def parse(self, text, lengths=None):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'genomic.gff'
            p.write_text(text)
            return parse_gff(p, lengths or {'chr1': 100, 'chr2': 100}, 'GCA_000000001.1')

    def test_one_based_conversion(self):
        index, _, features, status = self.parse('chr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
        self.assertEqual(features[0]['start_0based'], 0)
        self.assertEqual(features[0]['end_exclusive'], 10)
        self.assertEqual(status['annotation_status'], 'USABLE')

    def test_missing_contig_not_invented_gene_absence(self):
        index, coverage, _, _ = self.parse('chr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
        self.assertNotIn('chr2', index)
        self.assertEqual(coverage['chr2'], 'annotation_coverage_unknown')

    def test_full_declared_contig_with_no_genes(self):
        index, coverage, _, _ = self.parse('##sequence-region chr2 1 100\nchr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
        self.assertIn('chr2', index)
        self.assertEqual(coverage['chr2'], 'declared_contig_no_gene_features')

    def test_transcripts_not_called_genes(self):
        index, _, _, status = self.parse('chr1\tX\tmRNA\t1\t10\t.\t+\t.\tID=t1\n')
        self.assertEqual(status['annotation_status'], 'NO_GENE_FEATURES')
        self.assertFalse(index)

    def test_mismatched_contig_rejected(self):
        index, _, _, status = self.parse('wrong\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
        self.assertFalse(index)
        self.assertEqual(status['annotation_status'], 'INVALID_GFF')

    def test_out_of_bounds_not_clipped(self):
        _, _, _, status = self.parse('chr1\tX\tgene\t1\t101\t.\t+\t.\tID=g1\n')
        self.assertEqual(status['annotation_status'], 'INVALID_GFF')

    def test_build_accession_mismatch(self):
        _, _, _, status = self.parse('#!genome-build-accession NCBI_Assembly:GCA_000000002.1\nchr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
        self.assertEqual(status['annotation_status'], 'INVALID_GFF')

    def test_comments_and_embedded_fasta(self):
        _, _, features, status = self.parse('##gff-version 3\nchr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n##FASTA\n>chr1\nACGT\n')
        self.assertEqual(len(features), 1)
        self.assertEqual(status['annotation_status'], 'USABLE')

    def test_crlf_gff_is_readable(self):
        _, _, _, status = self.parse('##gff-version 3\r\nchr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\r\n')
        self.assertEqual(status['annotation_status'], 'USABLE')

    def test_conflicting_gff_requires_mapping(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            f = p / 'genomic.fna'
            f.write_text('>x\nACGT\n')
            (p / 'a.gff').write_text('a')
            (p / 'b.gff').write_text('b')
            self.assertEqual(choose_gff(p, f, 'GCA_000000001.1')[1], 'AMBIGUOUS_GFF_FILES')
            self.assertEqual(choose_gff(p, f, 'GCA_000000001.1', {'GCA_000000001.1': 'a.gff'})[1], 'EXPLICIT_MAPPING')

    def test_missing_gff_reported(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            f = p / 'x.fna'
            f.write_text('>x\nACGT\n')
            self.assertEqual(choose_gff(p, f, 'GCA_000000001.1')[1], 'MISSING_GFF')

    def test_gzip_gff(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'a.gff.gz'
            with gzip.open(p, 'wt') as f:
                f.write('chr1\tX\tgene\t1\t10\t.\t+\t.\tID=g1\n')
            self.assertEqual(parse_gff(p, {'chr1': 100}, 'GCA_000000001.1')[3]['annotation_status'], 'USABLE')

class SequenceAndSamplingTests(unittest.TestCase):

    def test_gc_denominators(self):
        result = sequence_metrics('GCCANN')
        self.assertEqual(result['gc_fraction_all_bases'], 0.5)
        self.assertEqual(result['gc_fraction_acgt'], 0.75)
        self.assertEqual(result['n_fraction'], 2 / 6)

    def test_all_ambiguous_not_gc_zero_acgt(self):
        r = sequence_metrics('NNRY')
        self.assertIsNone(r['gc_fraction_acgt'])
        self.assertIsNone(r['shannon_entropy_acgt'])
        self.assertEqual(r['ambiguous_fraction'], 1)

    def test_ambiguous_base_breaks_homopolymer(self):
        self.assertEqual(sequence_metrics('AAANAANN')['longest_acgt_run_bp'], 3)

    def test_interval_bounds(self):
        with self.assertRaises(ValueError):
            interval_metrics('ACGT', 0, 5, None, {'contig_end_windows_bp': []})

    def test_equal_length_sampling_has_one_start(self):
        self.assertEqual(sample_start(100, 100, random.Random(1)), 0)

    def test_invalid_sampling_rejected(self):
        with self.assertRaises(ValueError):
            sample_start(99, 100, random.Random(1))

    def test_last_start_included(self):

        class Last:

            def randrange(self, n):
                return n - 1
        self.assertEqual(sample_start(110, 100, Last()), 10)

    def test_sampling_seed_reproducible(self):
        a, b = (iteration_rng(123, 'GCA_000000001.1', 1), iteration_rng(123, 'GCA_000000001.1', 1))
        self.assertEqual([a.randrange(1000) for _ in range(100)], [b.randrange(1000) for _ in range(100)])

    def test_seed_changes_by_assembly_and_iteration(self):
        vals = [iteration_rng(123, acc, it).randrange(2 ** 40) for acc, it in [('a', 1), ('a', 2), ('b', 1)]]
        self.assertEqual(len(set(vals)), 3)

    def test_primer_sites_verified(self):
        product = {'primer_sequence': 'ACGCTCAAAC', 'contig': 'c', 'start_0based': 0, 'end_exclusive': 30, 'left_mismatches': 0, 'right_mismatches': 0}
        sequence_site_audit({'c': 'ACGCTCAAAC' + 'A' * 10 + 'GTTTGAGCGT'}, [product])
        with self.assertRaises(ValueError):
            sequence_site_audit({'c': 'TCGCTCAAAC' + 'A' * 10 + 'GTTTGAGCGT'}, [product])

    def test_gene_NA_vs_annotated_zero(self):
        c = {'gene_flank_bp': 5, 'contig_end_windows_bp': [5]}
        self.assertIsNone(interval_metrics('A' * 30, 5, 15, None, c)['overlaps_annotated_gene'])
        self.assertEqual(interval_metrics('A' * 30, 5, 15, GeneIndex([]), c)['overlaps_annotated_gene'], 0)

class AggregationTests(unittest.TestCase):

    def rows(self):
        rows = []
        for acc, n, observed, random_val in [('a', 1, 1.0, 0.0), ('b', 9, 0.0, 1.0)]:
            for it in range(3):
                rows.extend(summarize_intervals([{'metric': observed if it == 0 else random_val}] * n, acc, it, ['metric']))
        return rows

    def test_interval_and_assembly_weights_differ(self):
        result = aggregate_comparisons(self.rows(), {'a', 'b'}, 2, ['metric'])
        self.assertAlmostEqual(result[0]['observed'], 0.1)
        self.assertAlmostEqual(result[2]['observed'], 0.5)
        self.assertIsNone(result[2]['p_value'])

    def test_representative_scope_excludes_nonrepresentative(self):
        result = aggregate_comparisons(self.rows(), {'a'}, 2, ['metric'])
        self.assertEqual(result[1]['observed'], 1.0)
        self.assertEqual(result[1]['assemblies_selected'], 1)

    def test_missing_metric_not_filled_zero(self):
        rows = []
        for it in range(3):
            rows.extend(summarize_intervals([{'metric': None}], 'a', it, ['metric']))
        result = aggregate_comparisons(rows, {'a'}, 2, ['metric'])
        self.assertIsNone(result[0]['observed'])
        self.assertEqual(result[0]['nonmissing_observed_intervals'], 0)

    def test_denominator_mismatch_detected_for_gene_metric(self):
        rows = self.rows()
        rows[1]['nonmissing_intervals'] = 0
        with self.assertRaises(ValueError):
            aggregate_comparisons(rows, {'a', 'b'}, 2, ['metric'])

    def test_variable_acgt_denominators_are_reported(self):
        rows = []
        for it in range(3):
            rows.extend(summarize_intervals([{'gc_fraction_acgt': 0.5}, {'gc_fraction_acgt': 0.5 if it == 0 else None}], 'a', it, ['gc_fraction_acgt']))
        result = aggregate_comparisons(rows, {'a'}, 2, ['gc_fraction_acgt'])
        self.assertEqual(result[0]['nonmissing_observed_intervals'], 2)
        self.assertEqual(result[0]['random_nonmissing_interval_min'], 1)
        self.assertEqual(result[0]['observed'], 0.5)

    def test_nonfinite_metric_fails(self):
        with self.assertRaises(ValueError):
            summarize_intervals([{'metric': float('nan')}], 'a', 0, ['metric'])

class IOAndIntegrityTests(unittest.TestCase):

    def test_no_traversal(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                checked_relative(Path(d), '../outside')

    def test_table_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'table.tsv'
            write_tsv(p, [{'a': 1, 'b': None}])
            self.assertEqual(read_tsv(p), [{'a': '1', 'b': 'NA'}])
            self.assertNotIn(b'\r', p.read_bytes())

    def test_duplicate_fasta_id_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'x.fna'
            p.write_text('>a\nACGT\n>a\nACGT\n')
            with self.assertRaises(ValueError):
                list(fasta_records(p))

    def test_ragged_tsv_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'x.tsv'
            p.write_text('a\tb\n1\n')
            with self.assertRaises(ValueError):
                read_tsv(p)

    def test_cache_corruption_invalidates(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p / 'a.txt').write_text('original')
            write_json(p / 'CACHE.json', {'signature': 'test', 'files': {'a.txt': sha256(p / 'a.txt')}})
            self.assertTrue(cache_is_valid(p, 'test'))
            (p / 'a.txt').write_text('changed')
            self.assertFalse(cache_is_valid(p, 'test'))

    def test_new_cache_signature_invalidates(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            write_json(p / 'CACHE.json', {'signature': 'one', 'files': {}})
            self.assertFalse(cache_is_valid(p, 'two'))

    def test_bin_rule_stable_for_duplicate_sizes(self):
        self.assertEqual(cluster_sizes([100, 101, 101, 101, 200], 0.05), cluster_sizes([100, 101, 200], 0.05))
if __name__ == '__main__':
    unittest.main()
