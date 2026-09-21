"""Genome simulation, distance-matrix inference, and primer-set comparisons."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from .common import read_tsv, write_tsv, write_json, metadata_groups, sha256
from .sequence import SiteScanner, scan_assembly
from .validation import cached_analysis, analyze_products, split_sample_analysis
from .mash_refresh import refresh


def run(manifest: Path, data: Path, config: dict, out: Path, cache: Path,
        mash_reference: Path | None = None) -> dict:
    """Run the study analysis from FASTA files with optional verified scan caching.

    An explicitly supplied Mash table replaces executable invocation, never silently.
    The inventory fixes row order and biological grouping before distances are tested.
    """
    out.mkdir(parents=True, exist_ok=False)
    inventory = read_tsv(data / 'inventory.tsv')
    locations = read_tsv(manifest)
    paths = {}
    location_by_id = {}
    for row in locations:
        acc = row['assembly']
        if acc in paths: raise ValueError('Duplicate assembly in FASTA manifest')
        p = Path(row['fasta']).expanduser()
        p = (p if p.is_absolute() else manifest.parent / p).resolve()
        if not p.is_file(): raise FileNotFoundError(f'FASTA not found for {acc}: {p}')
        expected = row.get('fasta_sha256', '').strip()
        if expected and sha256(p) != expected: raise ValueError(f'FASTA hash mismatch: {acc}')
        paths[acc] = p
        location_by_id[acc] = row
    ids = [r['assembly'] for r in inventory]
    if len(ids) != len(set(ids)) or set(ids) != set(paths):
        raise ValueError('FASTA manifest and analysis inventory must have the same unique assemblies')
    groups, reps = metadata_groups(inventory)
    write_tsv(out / 'assembly_groups.tsv', groups)
    write_json(out / 'analysis_config_used.json', config)
    print('GENOME: distance matrices and assembly-label permutations', flush=True)
    cached = cached_analysis(data, out / 'baseline_statistics', inventory, reps, config, print)
    historical = read_tsv(data / 'historical_primers.tsv')
    randoms = read_tsv(data / 'random_primer_sequences.tsv')
    targets = [{'set_id': 'historical_guthrie', 'primer': r['primer'], 'sequence': r['sequence']} for r in historical]
    for row in randoms:
        if row.get('status') != 'accepted': raise ValueError('Primer set is not an accepted fixed input')
        for primer in historical:
            name, reference = primer['primer'], primer['sequence']
            seq = row[name]
            if sum(b in 'GC' for b in seq) != sum(b in 'GC' for b in reference):
                raise ValueError('Random primer does not match the specified GC count')
            targets.append({'set_id': row['set_id'], 'primer': name, 'sequence': seq})
    write_tsv(out / 'primer_sequences_used.tsv', targets)
    scanner = SiteScanner([r['sequence'] for r in targets], config['max_mismatches_per_site'], config['sequence_scan_chunk_bp'])
    target_sets = [r['set_id'] for r in targets]
    products = out / 'candidate_products'
    products.mkdir()
    arrays, records = [], []
    for i, row in enumerate(inventory, 1):
        acc = row['assembly']
        print(f'GENOME {i}/{len(ids)}: {acc}', flush=True)
        arr, meta, reused = scan_assembly(paths[acc], scanner, target_sets, config, cache, print)
        if int(row['expected_total_bases']) != meta['total_bases'] or int(row['expected_sequence_records']) != meta['records']:
            raise ValueError(f'FASTA inventory differs for {acc}')
        arrays.append(arr)
        path = products / (acc + '.npz')
        np.savez_compressed(path, amplicons=arr)
        meta = {**meta, 'npz_sha256': sha256(path)}
        write_json(products / (acc + '.json'), meta)
        records.append({'assembly': acc, 'fasta_sha256': meta['fasta_sha256'],
                        'sequence_records': meta['records'], 'total_bases': meta['total_bases'],
                        'candidate_intervals': len(arr), 'cache_reused': reused,
                        'scalar_site_check_bases': meta['reference_checked_bases']})
        write_tsv(out / 'genome_scan_inventory.tsv', records)
    if mash_reference:
        rows = read_tsv(mash_reference)
        positions = {a: i for i, a in enumerate(ids)}
        distance = np.full((len(ids), len(ids)), np.nan)
        np.fill_diagonal(distance, 0.)
        seen = set()
        for row in rows:
            i, j = positions[row['assembly_a']], positions[row['assembly_b']]
            value = float(row['mash_distance'])
            key = tuple(sorted((i, j)))
            if i == j or key in seen or not np.isfinite(value) or value < 0:
                raise ValueError('Invalid or duplicate supplied Mash pair')
            seen.add(key)
            distance[i, j] = distance[j, i] = value
        if not np.isfinite(distance).all(): raise ValueError('Incomplete supplied Mash matrix')
        mash_status = {'status': 'SUPPLIED_DISTANCE_TABLE', 'sha256': sha256(mash_reference)}
    else:
        distance, mash_status = refresh(paths, ids, out / 'mash', print)
        if distance is None:
            raise RuntimeError('Mash was not available. Install Mash natively or in default WSL, or explicitly supply a distance table.')
    cached['layers']['mash_distance'] = distance
    write_tsv(out / 'mash_pairs.tsv', [{'assembly_a': ids[i], 'assembly_b': ids[j], 'mash_distance': distance[i, j]}
                                     for i in range(len(ids)) for j in range(i + 1, len(ids))])
    write_json(out / 'mash_status.json', mash_status)
    print('GENOME: primer sets, product geometry, and size-compatibility comparison', flush=True)
    sizes, product_summary = analyze_products(arrays, targets, inventory, reps, data, cached, config, out / 'product_analysis', print)
    representative_sizes = {k: [v[i] for i in reps] for k, v in sizes.items()}
    heldout = split_sample_analysis(representative_sizes, distance[np.ix_(reps, reps)],
                                   [ids[i] for i in reps], config, out / 'internal_holdout', print)
    status = {'status': 'COMPLETE', 'assemblies': len(ids), 'metadata_representatives': len(reps),
              'product_analysis': product_summary, 'holdout': heldout, 'mash': mash_status}
    write_json(out / 'STATUS.json', status)
    return status
