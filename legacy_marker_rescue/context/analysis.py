"""Historical-product selection and interval-level context calculations."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import numpy as np
from .io_utils import read_json, write_json, read_tsv, write_tsv, sha256, write_rows_gzip
from .binning import cluster_sizes
from .metrics import metric_names, interval_metrics, iteration_rng, sample_start, summarize_intervals
EXPECTED_COLUMNS = ['target_index','contig_index','start_0based','end_exclusive','size_bp','orientation_0_inward_1_outward','left_mismatches','right_mismatches']
def load_selected(source: Path, config: dict):
    inventory = read_tsv(source / 'assembly_groups.tsv')
    names = [r['assembly'] for r in inventory]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate assembly rows')
    grouped = defaultdict(list)
    for r in inventory:
        grouped[r['analysis_group']].append(r)
    for group, rows in grouped.items():
        if sum((r['representative'] == '1' for r in rows)) != 1:
            raise ValueError(f'Expected one representative per metadata group: {group}')
    primers = read_tsv(source / 'primer_sequences_used.tsv')
    targets = {i: row for i, row in enumerate(primers) if row['set_id'] == config['historical_set_id']}
    if not targets:
        raise ValueError('Historical primer set missing')
    if len({r['primer'] for r in targets.values()}) != len(targets):
        raise ValueError('Duplicate historical primer label')
    parent_config = read_json(source / 'analysis_config_used.json')
    checks = [('minimum_amplicon_bp', 'min_amplicon_bp'), ('maximum_amplicon_bp', 'max_amplicon_bp'), ('maximum_mismatches_per_binding_site', 'max_mismatches_per_site')]
    for now, old in checks:
        if config[now] != parent_config[old]:
            raise ValueError(f'Follow-up setting differs from parent: {now}')
    setting = next((r for r in parent_config['candidate_settings'] if r['id'] == parent_config['reported_setting']))
    if setting['relative_tolerance'] != config['relative_band_tolerance']:
        raise ValueError('Band tolerance differs from parent')
    if config['orientation'] != 'inward':
        raise ValueError('This follow-up requires inward products')
    selected, metadata, raw_count = ({}, {}, 0)
    for inv in inventory:
        acc = inv['assembly']
        meta = read_json(source / 'candidate_products' / (acc + '.json'))
        array_path = source / 'candidate_products' / (acc + '.npz')
        if sha256(array_path) != meta['npz_sha256']:
            raise ValueError('Candidate-product hash mismatch')
        if meta['amplicon_columns'] != EXPECTED_COLUMNS:
            raise ValueError('Unrecognized candidate-product coordinate schema')
        with np.load(array_path, allow_pickle=False) as z:
            arr = z['amplicons'].copy()
        if arr.ndim != 2 or arr.shape[1] != 8 or (not np.issubdtype(arr.dtype, np.integer)):
            raise ValueError('Unexpected candidate array structure')
        if len(meta['contig_names']) != len(meta['contig_lengths']) or len(set(meta['contig_names'])) != len(meta['contig_names']):
            raise ValueError('Invalid contig metadata')
        mask = np.isin(arr[:, 0], list(targets)) & (arr[:, 5] == 0) & (arr[:, 4] >= config['minimum_amplicon_bp']) & (arr[:, 4] <= config['maximum_amplicon_bp'])
        subset = arr[mask]
        if len(subset) != len(np.unique(subset, axis=0)):
            raise ValueError('Duplicate selected candidate rows')
        products = []
        for target, ci, start, end, size, orientation, lm, rm in sorted(subset.tolist()):
            if not 0 <= ci < len(meta['contig_names']):
                raise ValueError('Invalid contig index')
            if not 0 <= start < end <= meta['contig_lengths'][ci] or end - start != size:
                raise ValueError('Candidate interval length or coordinate mismatch')
            if not 0 <= lm <= config['maximum_mismatches_per_binding_site'] or not 0 <= rm <= config['maximum_mismatches_per_binding_site']:
                raise ValueError('Mismatch bound violated')
            primer = targets[target]
            identifier = f"{acc}|{primer['primer']}|{ci}|{start}|{end}"
            products.append({'product_id': identifier, 'assembly': acc, 'primer': primer['primer'], 'primer_sequence': primer['sequence'], 'contig': meta['contig_names'][ci], 'start_0based': start, 'end_exclusive': end, 'size_bp': size, 'left_mismatches': lm, 'right_mismatches': rm, 'orientation': 'inward', 'analysis_group': inv['analysis_group'], 'is_representative': int(inv['representative'])})
        selected[acc] = products
        metadata[acc] = meta
        raw_count += len(products)
    previous = [r for r in read_tsv(source / 'product_analysis' / 'primer_set_metrics.tsv') if r['set_id'] == config['historical_set_id'] and r['mode'] == 'inward_min100']
    if len(previous) != 1:
        raise ValueError('Parent inward model summary missing or duplicated')
    if raw_count != int(previous[0]['products']):
        raise ValueError('Product total does not reproduce parent output')
    feature_rows, blocks = ([], [])
    for primer in sorted((r['primer'] for r in targets.values())):
        values = [p['size_bp'] for rows in selected.values() for p in rows if p['primer'] == primer]
        mapping, centers, members = cluster_sizes(values, config['relative_band_tolerance'])
        block = np.zeros((len(names), len(centers)), dtype=np.int8)
        for i, acc in enumerate(names):
            for product in selected[acc]:
                if product['primer'] == primer:
                    block[i, mapping[product['size_bp']]] = 1
        for j, center in enumerate(centers):
            count = int(block[:, j].sum())
            key = f'{primer}|bin{j + 1:03d}'
            feature_rows.append({'band_id': key, 'primer': primer, 'center_bp': center, 'presence_count_all_assemblies': count, 'assemblies': len(names), 'prevalence_all_assemblies': count / len(names), 'class': 'fixed' if count == len(names) else 'variable', 'member_sizes_bp': ';'.join((str(v) for v in members[j]))})
            for rows in selected.values():
                for product in rows:
                    if product['primer'] == primer and mapping[product['size_bp']] == j:
                        product.update(band_id=key, band_center_bp=center, band_presence_count_all_assemblies=count, band_prevalence_all_assemblies=count / len(names), band_class='fixed' if count == len(names) else 'variable')
        blocks.append(block)
    matrix = np.concatenate(blocks, axis=1)
    reference_path = source / 'product_analysis' / 'band_matrices' / 'inward_min100' / 'historical_guthrie.npz'
    with np.load(reference_path, allow_pickle=False) as z:
        key = 'presence'
        if key not in z.files:
            raise ValueError('Parent band-matrix key is missing')
        if z['ids'].tolist() != names:
            raise ValueError('Parent matrix assembly order differs')
        if not np.array_equal(matrix, z[key]):
            raise ValueError('Rebuilt inward band matrix differs from parent')
    if matrix.shape[1] != int(previous[0]['band_bins']):
        raise ValueError('Band count differs from parent')
    return (inventory, targets, selected, metadata, feature_rows)

def sequence_site_audit(seqs, products):
    complement = str.maketrans('ACGT', 'TGCA')
    for p in products:
        seq = seqs[p['contig']]
        primer = p['primer_sequence']
        rc = primer.translate(complement)[::-1]
        left = seq[p['start_0based']:p['start_0based'] + len(primer)]
        right = seq[p['end_exclusive'] - len(primer):p['end_exclusive']]
        if set(left + right) - set('ACGT'):
            raise ValueError('Ambiguous nucleotide at a cached primer-binding site')
        if sum((a != b for a, b in zip(left, primer))) != p['left_mismatches']:
            raise ValueError('Cached left-site mismatch does not agree with FASTA')
        if sum((a != b for a, b in zip(right, rc))) != p['right_mismatches']:
            raise ValueError('Cached right-site mismatch does not agree with FASTA')

def cache_is_valid(cache, signature):
    try:
        state = read_json(cache / 'CACHE.json')
        if state['signature'] != signature:
            return False
        return all(((cache / rel).is_file() and sha256(cache / rel) == h for rel, h in state['files'].items()))
    except (OSError, ValueError, KeyError):
        return False

def build_assembly(cache, acc, products, seqs, indexes, coverage, feature_rows, annotation, config, signature, audit):
    cache.mkdir(parents=True, exist_ok=True)
    names = metric_names(config)
    observed = []
    for product in products:
        row = {**product, 'random_iteration': 0, 'annotation_coverage': coverage.get(product['contig'], annotation['annotation_status']), 'contig_length_bp': len(seqs[product['contig']])}
        row.update(interval_metrics(seqs[product['contig']], product['start_0based'], product['end_exclusive'], indexes.get(product['contig']), config))
        observed.append(row)
    write_rows_gzip(cache / 'observed_intervals.tsv.gz', observed)
    summaries = summarize_intervals(observed, acc, 0, names)
    randoms = []
    for iteration in range(1, int(config['random_iterations']) + 1):
        rng = iteration_rng(config['random_seed'], acc, iteration)
        batch = []
        for product in products:
            seq = seqs[product['contig']]
            start = sample_start(len(seq), product['size_bp'], rng)
            end = start + product['size_bp']
            row = {'product_id': product['product_id'], 'assembly': acc, 'random_iteration': iteration, 'primer': product['primer'], 'contig': product['contig'], 'start_0based': start, 'end_exclusive': end, 'contig_length_bp': len(seq), 'annotation_coverage': coverage.get(product['contig'], annotation['annotation_status'])}
            row.update(interval_metrics(seq, start, end, indexes.get(product['contig']), config))
            batch.append(row)
        randoms.extend(batch)
        summaries.extend(summarize_intervals(batch, acc, iteration, names))
    write_rows_gzip(cache / 'matched_random_intervals.tsv.gz', randoms)
    write_tsv(cache / 'iteration_metrics.tsv', summaries)
    write_rows_gzip(cache / 'annotation_intervals_used.tsv.gz', feature_rows, ['contig', 'start_0based', 'end_exclusive', 'feature_id'])
    coverage_rows = [{'contig': name, 'length_bp': len(seq), 'annotation_coverage': coverage.get(name, annotation['annotation_status'])} for name, seq in seqs.items()]
    write_tsv(cache / 'contig_coverage.tsv', coverage_rows)
    audit = {**audit, 'annotation': annotation, 'observed_products': len(observed), 'matched_random_intervals': len(randoms), 'primer_site_checks': len(products) * 2, 'gene_context_observed_intervals': sum((r['overlaps_annotated_gene'] is not None for r in observed))}
    write_json(cache / 'assembly_audit.json', audit)
    files = {p.name: sha256(p) for p in cache.iterdir() if p.is_file() and p.name != 'CACHE.json'}
    write_json(cache / 'CACHE.json', {'signature': signature, 'files': files})
