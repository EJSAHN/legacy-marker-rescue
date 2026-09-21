"""End-to-end tests on artificial inputs, distinct from the study observations."""
from __future__ import annotations
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from legacy_marker_rescue.cli import main, ROOT
from legacy_marker_rescue.genome.common import write_tsv, sha256
from legacy_marker_rescue.genome.sequence import SiteScanner, scan_assembly, revcomp, select_sizes, pooled_matrix
from legacy_marker_rescue.genome.stats import jaccard
from legacy_marker_rescue.io import resolve_path
from tests.test_external import samplezip


def fixture(root):
    """Small artificial genomes and gel pairs; no research results are generated here."""
    cfg = root / 'config'; shutil.copytree(ROOT / 'config', cfg)
    config = json.loads((cfg / 'genome.json').read_text())
    config.update(permutations=9, holdout_splits=2, holdout_permutations=9,
                  reference_scan_bases_per_assembly=1000)
    (cfg / 'genome.json').write_text(json.dumps(config))
    cc = json.loads((cfg / 'context.json').read_text()); cc['random_iterations'] = 2
    (cfg / 'context.json').write_text(json.dumps(cc))
    data = root / 'data'; data.mkdir(); cache = root / 'cache'
    primers = [{'primer': 'P1', 'sequence': 'CAGGCCCTTC'},
               {'primer': 'P2', 'sequence': 'AGTCAGCCAC'},
               {'primer': 'P3', 'sequence': 'AATCGGGCTG'}]
    random_set = {'set_id': 'random_test', 'status': 'accepted',
                  'P1': 'CCAGCTCGTC', 'P2': 'ATCGAGCCAC', 'P3': 'ATCGACGGTG'}
    write_tsv(data / 'historical_primers.tsv', primers)
    write_tsv(data / 'random_primer_sequences.tsv', [random_set])
    seqs = [p['sequence'] for p in primers] + [random_set[p['primer']] for p in primers]
    scanner = SiteScanner(seqs)
    target_sets = ['historical_guthrie'] * 3 + ['random_test'] * 3
    arrays, inv, locations = [], [], []
    rng = np.random.default_rng(5874)
    for i in range(10):
        acc = f'GCA_{900000000+i}.1'
        sequence = rng.choice(list('ACGT'), 25000)
        for pi, primer in enumerate(seqs):
            start = pi * 4000 + 80
            size = 180 + (i * 101 + pi * 203) % 1300
            sequence[start:start+10] = list(primer)
            sequence[start+size-10:start+size] = list(revcomp(primer))
            if (i + pi) % 3:
                start += 1900; size = 150 + (i * 73 + pi * 111) % 900
                sequence[start:start+10] = list(primer)
                sequence[start+size-10:start+size] = list(revcomp(primer))
        fasta = root / (acc + '.fna'); fasta.write_text('>test_contig\n' + ''.join(sequence) + '\n')
        arr, meta, _ = scan_assembly(fasta, scanner, target_sets, config, cache, lambda *x: None)
        arrays.append(arr)
        inv.append({'assembly': acc, 'expected_total_bases': len(sequence), 'expected_sequence_records': 1,
                    'biosample': f'test_sample_{i}', 'strain': f'test_strain_{i}', 'isolate': '', 'submitter': ''})
        locations.append({'assembly': acc, 'fasta': fasta.name, 'fasta_sha256': sha256(fasta), 'gff': '', 'gff_sha256': ''})
    gff = root / 'test.gff'; gff.write_text('##gff-version 3\n##sequence-region test_contig 1 25000\ntest_contig\ttest\tgene\t100\t800\t.\t+\t.\tID=gene1\n')
    locations[0].update(gff=gff.name, gff_sha256=sha256(gff))
    write_tsv(data / 'inventory.tsv', inv); write_tsv(root / 'manifest.tsv', locations)
    indices = {p['primer']: n for n, p in enumerate(primers)}
    matrices = {}
    for setting in config['candidate_settings']:
        sizes = select_sizes(arrays, indices, 'both', 100, setting['max_amplicon_bp'])
        features, matrix = pooled_matrix(sizes, setting['relative_tolerance'])
        matrices[setting['id']] = matrix
        columns = [f'band_{j}' for j in range(len(features))]
        write_tsv(data / ('bands_' + setting['id'] + '.tsv'),
                  [{'assembly_accession': r['assembly'], **dict(zip(columns, matrix[n].tolist()))} for n,r in enumerate(inv)])
    distances = jaccard(matrices[config['reported_setting']])
    pairs = [{'assembly_a': inv[i]['assembly'], 'assembly_b': inv[j]['assembly'],
              'mash_distance': 0.001 + (abs(i-j)+0.1*(i+j))/100,
              'auxiliary_distance': 0.002 + (abs(i-j)+0.12*(i+j))/100,
              'rapd_distance': distances[i,j]} for i in range(10) for j in range(i+1,10)]
    write_tsv(data / 'distance_pairs.tsv', pairs)
    write_tsv(root / 'mash.tsv', [{k:r[k] for k in ['assembly_a','assembly_b','mash_distance']} for r in pairs])
    write_tsv(data / 'legacy_bridge.tsv', [{'legacy_band': 'test_band', 'legacy_kb': 0.5}])
    write_tsv(data / 'random_primer_previous_summary.tsv',
              [{'set_id':'historical','set_type':'legacy','amplicons_total':1,'all_pairs_spearman_r':0.1},
               {'set_id':'random_test','set_type':'random','amplicons_total':1,'all_pairs_spearman_r':0.1}])
    ext_cache = root / 'external'; ext_cache.mkdir(); source = ext_cache / 'test.zip'; samplezip(source)
    bc = json.loads((cfg / 'benchmark.json').read_text())
    bc['collections'] = [{'id':'test_collection','file':'test.zip','md5':hashlib.md5(source.read_bytes()).hexdigest(),
                          'expected_image_pairs':2,'role':'software_test_only'}]
    (cfg / 'benchmark.json').write_text(json.dumps(bc))
    inputs = root / 'inputs.json'
    inputs.write_text(json.dumps({'genome_manifest':'manifest.tsv','genome_data':'data',
                                 'scan_cache':'cache','mash_table':'mash.tsv','external_cache':'external'}))
    return inputs, cfg


class WorkflowTests(unittest.TestCase):
    def test_relative_input_resolution(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            self.assertEqual(resolve_path('file.dat', p), (p/'file.dat').resolve())

    def test_output_inside_repository_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main(['detection','--output-dir',str(ROOT/'outputs_that_must_not_be_created')])

    def test_existing_output_rejected(self):
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stderr(io.StringIO()):
            (Path(d)/'keep.txt').write_text('keep')
            with self.assertRaises(SystemExit): main(['detection','--output-dir',d])
            self.assertEqual((Path(d)/'keep.txt').read_text(),'keep')

    def test_failure_status_and_diagnostics(self):
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            out=Path(d)/'output'
            rc=main(['genome','--output-dir',str(out)])
            self.assertNotEqual(rc,0)
            self.assertEqual(json.loads((out/'STATUS.json').read_text())['status'],'FAILED')
            self.assertTrue((out/'error.log').is_file())
            self.assertTrue((out/'SHA256SUMS.txt').is_file())

    def test_complete_artificial_workflow_and_readonly_inputs(self):
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
            root=Path(d); inputs,cfg=fixture(root)
            prior={p.relative_to(root).as_posix():sha256(p) for p in root.iterdir() if p.is_file()}
            out=root/'out'
            with patch('legacy_marker_rescue.benchmark.acquire.urlopen', side_effect=AssertionError('No network in fixture')):
                code=main(['all','--inputs',str(inputs),'--config-dir',str(cfg),'--output-dir',str(out)])
            self.assertEqual(code,0,(out/'error.log').read_text() if (out/'error.log').exists() else '')
            status=json.loads((out/'STATUS.json').read_text())
            self.assertEqual(status['status'],'COMPLETE')
            self.assertEqual(set(status['stages']),{'genome','context','detection','external'})
            self.assertEqual(status['stages']['genome']['assemblies'],10)
            self.assertEqual(status['stages']['context']['assemblies_with_gene_annotation'],1)
            self.assertEqual(status['stages']['external']['images'],2)
            self.assertEqual(prior,{p.relative_to(root).as_posix():sha256(p) for p in root.iterdir() if p.is_file()})
            from legacy_marker_rescue.genome.common import read_tsv
            scans=read_tsv(out/'genome/genome_scan_inventory.tsv')
            self.assertTrue(all(x['cache_reused']=='True' for x in scans))
            groups=read_tsv(out/'context/annotation_inventory.tsv')
            self.assertEqual(sum(x['annotation_status']=='MISSING_GFF' for x in groups),9)

if __name__=='__main__': unittest.main()
