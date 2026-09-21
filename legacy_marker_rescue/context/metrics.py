"""Paired interval metrics and summaries; no independent-pair significance tests."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import math
import random
import statistics
import numpy as np
BASE_METRICS = ['size_bp', 'distance_to_contig_end_bp', 'gc_fraction_all_bases', 'gc_fraction_acgt', 'n_fraction', 'ambiguous_fraction', 'shannon_entropy_acgt', 'longest_acgt_run_bp', 'overlaps_annotated_gene', 'nearest_gene_distance_bp', 'gene_count_within_flank']

def metric_names(config):
    return BASE_METRICS + ['within_' + str(w) + 'bp_of_contig_end' for w in config['contig_end_windows_bp']]

def sequence_metrics(sequence: str):
    seq = sequence.upper()
    n = len(seq)
    if not n:
        raise ValueError('Empty interval sequence')
    counts = [seq.count(c) for c in 'ACGT']
    canonical = sum(counts)
    ent = -sum((c / canonical * math.log2(c / canonical) for c in counts if c)) if canonical else None
    longest, current, last = (0, 0, None)
    for c in seq:
        if c in 'ACGT':
            current = current + 1 if c == last else 1
            longest = max(longest, current)
            last = c
        else:
            current, last = (0, None)
    return {'gc_fraction_all_bases': (counts[1] + counts[2]) / n, 'gc_fraction_acgt': (counts[1] + counts[2]) / canonical if canonical else None, 'n_fraction': seq.count('N') / n, 'ambiguous_fraction': (n - canonical) / n, 'shannon_entropy_acgt': ent, 'longest_acgt_run_bp': longest}

def interval_metrics(sequence: str, start: int, end: int, index, config):
    if not 0 <= start < end <= len(sequence):
        raise ValueError('Interval outside FASTA record')
    distance = min(start, len(sequence) - end)
    row = {'size_bp': end - start, 'distance_to_contig_end_bp': distance}
    row.update(sequence_metrics(sequence[start:end]))
    for window in config['contig_end_windows_bp']:
        row['within_' + str(window) + 'bp_of_contig_end'] = int(distance <= window)
    if index is None:
        row.update(overlaps_annotated_gene=None, nearest_gene_distance_bp=None, gene_count_within_flank=None)
    else:
        row.update(index.query(start, end, int(config['gene_flank_bp'])))
    return row

def iteration_rng(seed: int, assembly: str, iteration: int):
    h = hashlib.sha256(f'{int(seed)}|{assembly}|{int(iteration)}'.encode('ascii')).digest()
    return random.Random(int.from_bytes(h, 'big'))

def sample_start(length: int, size: int, rng) -> int:
    if size <= 0 or length < size:
        raise ValueError('Invalid matched interval length')
    return rng.randrange(length - size + 1)

def as_number(value):
    if value in (None, '', 'NA'):
        return None
    try:
        out = float(value)
    except (ValueError, TypeError):
        raise ValueError(f'Non-numeric metric: {value!r}')
    if not math.isfinite(out):
        raise ValueError('Nonfinite number must be encoded as NA')
    return out

def summarize_intervals(rows: list[dict], assembly: str, iteration: int, names):
    out = []
    for name in names:
        values = [as_number(r.get(name)) for r in rows]
        values = [v for v in values if v is not None]
        total = math.fsum(values) if values else None
        out.append({'assembly': assembly, 'random_iteration': iteration, 'metric': name, 'intervals': len(rows), 'nonmissing_intervals': len(values), 'sum': total, 'mean': total / len(values) if values else None})
    return out

def aggregate_comparisons(rows: list[dict], representatives: set[str], n_iterations: int, names):
    by_key = {}
    assemblies = set()
    for r in rows:
        acc, it, metric = (r['assembly'], int(r['random_iteration']), r['metric'])
        key = (acc, it, metric)
        if key in by_key:
            raise ValueError('Duplicate assembly/iteration/metric summary')
        by_key[key] = r
        assemblies.add(acc)
    output = []
    scopes = [('all_assemblies_interval_weighted', assemblies, False), ('representatives_interval_weighted', representatives, False), ('representatives_equal_assembly_weight', representatives, True)]
    for scope, selected, equal in scopes:
        for metric in names:
            valid = []
            for acc in sorted(selected):
                obs = by_key[acc, 0, metric]
                obs_count = int(obs['nonmissing_intervals'])
                for it in range(1, n_iterations + 1):
                    rand = by_key[acc, it, metric]
                    if int(rand['intervals']) != int(obs['intervals']):
                        raise ValueError('Observed/random interval count mismatch')
                    if metric not in ('gc_fraction_acgt', 'shannon_entropy_acgt') and int(rand['nonmissing_intervals']) != obs_count:
                        raise ValueError('Paired metric denominator differs across background draws')
                if obs_count and all((int(by_key[acc, it, metric]['nonmissing_intervals']) > 0 for it in range(1, n_iterations + 1))):
                    valid.append(acc)

            def combined(it):
                if not valid:
                    return None
                if equal:
                    return statistics.mean((as_number(by_key[a, it, metric]['mean']) for a in valid))
                numerator = math.fsum((as_number(by_key[a, it, metric]['sum']) for a in valid))
                denominator = sum((int(by_key[a, it, metric]['nonmissing_intervals']) for a in valid))
                return numerator / denominator
            observed = combined(0)
            draws = [combined(it) for it in range(1, n_iterations + 1)] if valid else []
            mean = statistics.mean(draws) if draws else None
            output.append({'scope': scope, 'metric': metric, 'assemblies_selected': len(selected), 'assemblies_with_metric': len(valid), 'observed_intervals_selected': sum((int(by_key[a, 0, metric]['intervals']) for a in selected)), 'nonmissing_observed_intervals': sum((int(by_key[a, 0, metric]['nonmissing_intervals']) for a in valid)), 'random_nonmissing_interval_min': min((sum((int(by_key[a, it, metric]['nonmissing_intervals']) for a in valid)) for it in range(1, n_iterations + 1))), 'random_nonmissing_interval_max': max((sum((int(by_key[a, it, metric]['nonmissing_intervals']) for a in valid)) for it in range(1, n_iterations + 1))), 'random_iterations': n_iterations, 'observed': observed, 'random_mean': mean, 'observed_minus_random_mean': observed - mean if mean is not None else None, 'observed_divided_by_random_mean': observed / mean if mean not in (None, 0) else None, 'random_iteration_sd': statistics.stdev(draws) if len(draws) > 1 else None, 'random_iteration_q025': float(np.quantile(draws, 0.025)) if draws else None, 'random_iteration_q975': float(np.quantile(draws, 0.975)) if draws else None, 'range_interpretation': 'background-draw variation, not a confidence interval', 'p_value': None, 'inference': 'descriptive matched-background comparison; no interval-level significance test'})
    return output
