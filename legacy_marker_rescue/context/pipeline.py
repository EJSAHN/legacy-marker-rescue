"""Sequence and annotation context for inward-facing historical products."""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
from .analysis import load_selected, sequence_site_audit, build_assembly
from .io_utils import read_tsv, write_tsv, read_json, write_json, sha256, fasta_records, read_rows_gzip
from .annotation import parse_gff
from .metrics import aggregate_comparisons, metric_names


def run(genome_result: Path, manifest: Path, config: dict, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=False)
    inventory, targets, selected, metadata, bands = load_selected(genome_result, config)
    locations = {r['assembly']: r for r in read_tsv(manifest)}
    if len(locations) != len(inventory): raise ValueError('Context input inventory differs from genome results')
    write_json(out / 'context_config_used.json', config)
    write_tsv(out / 'historical_band_definitions.tsv', bands)
    write_tsv(out / 'historical_inward_products.tsv', [p for r in inventory for p in selected[r['assembly']]])
    summaries, observations, annotations, issues = [], [], [], []
    for index, inv in enumerate(inventory, 1):
        acc = inv['assembly']; row = locations[acc]
        fasta = Path(row['fasta']).expanduser()
        fasta = (fasta if fasta.is_absolute() else manifest.parent / fasta).resolve()
        if sha256(fasta) != metadata[acc]['fasta_sha256']: raise ValueError(f'FASTA changed: {acc}')
        seqs = dict(fasta_records(fasta))
        if list(seqs) != metadata[acc]['contig_names'] or [len(s) for s in seqs.values()] != metadata[acc]['contig_lengths']:
            raise ValueError(f'FASTA contig order or length differs: {acc}')
        sequence_site_audit(seqs, selected[acc])
        gff_text = row.get('gff', '').strip()
        if gff_text and gff_text != 'NA':
            gff = Path(gff_text).expanduser()
            gff = (gff if gff.is_absolute() else manifest.parent / gff).resolve()
            gff_sha = sha256(gff)
            if row.get('gff_sha256') and row['gff_sha256'] != gff_sha:
                raise ValueError(f'GFF changed: {acc}')
            indexes, coverage, features, annotation = parse_gff(gff, {k: len(v) for k, v in seqs.items()}, acc, tuple(config['feature_types']))
            resolution = 'EXPLICIT_MANIFEST'
        else:
            gff = None; gff_sha = None; indexes = {}; features = []
            coverage = {k: 'missing_gff' for k in seqs}
            annotation = {'annotation_status': 'MISSING_GFF', 'errors': [], 'warnings': [], 'gene_covered_contigs': 0}
            resolution = 'MISSING_GFF'
        if annotation['annotation_status'] != 'USABLE':
            issues.append({'assembly': acc, 'issue': annotation['annotation_status'],
                           'consequence': 'Gene metrics remain missing; sequence metrics are retained'})
        audit = {'assembly': acc, 'analysis_group': inv['analysis_group'], 'is_representative': int(inv['representative']),
                 'fasta_sha256': metadata[acc]['fasta_sha256'], 'gff_sha256': gff_sha, 'gff_resolution': resolution}
        signature = hashlib.sha256(json.dumps({'assembly': acc, 'fasta_sha256': audit['fasta_sha256'],
                       'gff_sha256': gff_sha, 'config': config}, sort_keys=True).encode()).hexdigest()
        dest = out / 'context_by_assembly' / acc
        build_assembly(dest, acc, selected[acc], seqs, indexes, coverage, features, annotation, config, signature, audit)
        record = read_json(dest / 'assembly_audit.json')
        annotations.append({'assembly': acc, 'analysis_group': inv['analysis_group'], 'is_representative': inv['representative'],
                            'fasta_sha256': audit['fasta_sha256'], 'gff_sha256': gff_sha,
                            'annotation_status': annotation['annotation_status'],
                            'gene_covered_contigs': annotation['gene_covered_contigs'],
                            'observed_products': record['observed_products'],
                            'gene_context_observed_intervals': record['gene_context_observed_intervals'],
                            'matched_random_intervals': record['matched_random_intervals']})
        summaries.extend(read_tsv(dest / 'iteration_metrics.tsv'))
        observations.extend(read_rows_gzip(dest / 'observed_intervals.tsv.gz'))
        print(f'CONTEXT {index}/{len(inventory)}: {acc}; {len(selected[acc])} products; annotation={annotation["annotation_status"]}', flush=True)
        del seqs
    representatives = {r['assembly'] for r in inventory if r['representative'] == '1'}
    comparison = aggregate_comparisons(summaries, representatives, int(config['random_iterations']), metric_names(config))
    write_tsv(out / 'assembly_iteration_metrics.tsv', summaries)
    write_tsv(out / 'context_comparison.tsv', comparison)
    write_tsv(out / 'historical_observed_context.tsv', observations)
    write_tsv(out / 'annotation_inventory.tsv', annotations)
    write_tsv(out / 'annotation_issues.tsv', issues, ['assembly', 'issue', 'consequence'])
    status = {'status': 'COMPLETE', 'assemblies': len(inventory), 'metadata_representatives': len(representatives),
              'observed_products': sum(len(v) for v in selected.values()),
              'matched_random_intervals': sum(r['matched_random_intervals'] for r in annotations),
              'assemblies_with_gene_annotation': sum(r['annotation_status'] == 'USABLE' for r in annotations)}
    if status['matched_random_intervals'] != status['observed_products'] * int(config['random_iterations']):
        raise ValueError('Paired interval totals differ')
    write_json(out / 'STATUS.json', status)
    return status
