"""Matrix-label permutation and descriptive correlations; no independent-edge tests."""
from __future__ import annotations
import math
import numpy as np
from scipy.stats import rankdata
from .numerics import sum_products, row_products, euclidean_norm
from .common import accession, read_tsv, write_tsv, write_json

def z_vector(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all() or len(x) < 3:
        return None
    y = x - x.mean()
    norm = euclidean_norm(y)
    return y / norm if norm > np.finfo(float).eps else None

def correlation(x, y, method='spearman'):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or not np.isfinite(x).all() or (not np.isfinite(y).all()):
        return float('nan')
    if method == 'spearman':
        x = rankdata(x, method='average')
        y = rankdata(y, method='average')
    elif method != 'pearson':
        raise ValueError('Unknown correlation method')
    x = z_vector(x)
    y = z_vector(y)
    return float(np.clip(sum_products(x, y), -1, 1)) if x is not None and y is not None else float('nan')

def jaccard(binary):
    a = np.asarray(binary)
    if a.ndim != 2 or not np.isin(a, [0, 1]).all():
        raise ValueError('Expected a binary matrix')
    a = a.astype(np.int64)
    intersection = np.einsum('ik,jk->ij', a, a, optimize=False)
    counts = a.sum(axis=1)
    union = counts[:, None] + counts[None, :] - intersection
    out = np.full(union.shape, np.nan, dtype=float)
    np.divide(intersection, union, out=out, where=union > 0)
    out = 1 - out
    np.fill_diagonal(out, 0.0)
    return out

def load_band_matrix(path, assemblies):
    rows = read_tsv(path)
    cols = [k for k in rows[0] if k != 'assembly_accession']
    rr = {accession(r['assembly_accession']): r for r in rows}
    if len(rr) != len(rows) or set(rr) != set(assemblies):
        raise ValueError(f'Assembly mismatch in {path.name}')
    x = np.array([[int(rr[a][c]) for c in cols] for a in assemblies], dtype=np.int8)
    if not np.isin(x, [0, 1]).all():
        raise ValueError('Nonbinary band matrix')
    return (cols, x)

def load_distance_layers(path, assemblies):
    rows = read_tsv(path)
    n = len(assemblies)
    pos = {a: i for i, a in enumerate(assemblies)}
    seen = set()
    matrices = {name: np.zeros((n, n)) for name in ('mash_distance', 'auxiliary_distance', 'rapd_distance')}
    for row in rows:
        a = accession(row['assembly_a'])
        b = accession(row['assembly_b'])
        if a not in pos or b not in pos or a == b:
            raise ValueError('Unexpected assembly or self-pair')
        i, j = (pos[a], pos[b])
        key = tuple(sorted((i, j)))
        if key in seen:
            raise ValueError('Duplicated pair record')
        seen.add(key)
        for name, d in matrices.items():
            v = float(row[name])
            if not np.isfinite(v) or v < 0:
                raise ValueError('Invalid distance')
            d[i, j] = d[j, i] = v
    if len(seen) != n * (n - 1) // 2:
        raise ValueError('Incomplete pair matrix')
    return matrices

def make_permutations(n, b, seed, blocks=None):
    rng = np.random.default_rng(seed)
    out = np.tile(np.arange(n), (b, 1))
    groups = [np.arange(n)] if blocks is None else [np.flatnonzero(np.asarray(blocks) == g) for g in sorted(set(blocks))]
    for i in range(b):
        for ids in groups:
            out[i, ids] = rng.permutation(ids)
    return out

def row_z(a):
    a = np.asarray(a, dtype=float)
    center = a - a.mean(axis=1, keepdims=True)
    norm = euclidean_norm(center, axis=1, keepdims=True)
    return np.divide(center, norm, out=np.full_like(center, np.nan), where=norm > np.finfo(float).eps)

def permutation_family(distances, target, permutations):
    """Same entity-label permutation applied jointly to rows/columns of target for all settings."""
    distances = np.asarray(distances, dtype=float)
    n = target.shape[0]
    i, j = np.triu_indices(n, 1)
    if distances.ndim != 3 or distances.shape[1:] != target.shape:
        raise ValueError('Matrix dimension mismatch')
    if not np.isfinite(distances).all() or not np.isfinite(target).all():
        raise ValueError('Incomplete matrix in permutation test')
    if not np.allclose(target, target.T) or not np.allclose(distances, distances.transpose(0, 2, 1)):
        raise ValueError('Distance matrices must be symmetric')
    x = distances[:, i, j]
    y = target[i, j]
    yperm = target[permutations[:, i], permutations[:, j]]
    results = {}
    for method in ('spearman', 'pearson'):
        xx = rankdata(x, axis=1, method='average') if method == 'spearman' else x
        yy = rankdata(y, method='average') if method == 'spearman' else y
        yp = rankdata(yperm, axis=1, method='average') if method == 'spearman' else yperm
        xz = row_z(xx)
        yz = z_vector(yy)
        if yz is None:
            raise ValueError('Genome distance has no variation')
        obs = row_products(xz, yz[None, :])[:, 0]
        perm = row_products(row_z(yp), xz)
        valid = np.isfinite(obs)
        maxnull = np.nanmax(np.abs(perm[:, valid]), axis=1) if valid.any() else np.full(len(perm), np.nan)
        p = []
        padj = []
        for k, r in enumerate(obs):
            if not valid[k]:
                p.append(np.nan)
                padj.append(np.nan)
                continue
            p.append((1 + np.count_nonzero(np.abs(perm[:, k]) >= abs(r) - 1e-12)) / (len(perm) + 1))
            padj.append((1 + np.count_nonzero(maxnull >= abs(r) - 1e-12)) / (len(perm) + 1))
        results[method] = {'observed': obs, 'permuted': perm, 'p_two_sided': np.array(p), 'p_search_adjusted': np.array(padj)}
    return results

def write_family(out, settings, results, scope, units, blocks_desc='unrestricted entity-label exchangeability assumed'):
    rows = []
    for method, res in results.items():
        for k, setting in enumerate(settings):
            rows.append({'scope': scope, 'setting': setting, 'method': method, 'assemblies_or_representatives': units, 'pairs': units * (units - 1) // 2, 'correlation': res['observed'][k], 'permutations': len(res['permuted']), 'p_two_sided': res['p_two_sided'][k], 'p_max_absolute_over_settings': res['p_search_adjusted'][k], 'assumption': blocks_desc, 'interpretation': 'exploratory matrix association; not external validation'})
    write_tsv(out / 'matrix_permutation_tests.tsv', rows)
    np.savez_compressed(out / 'permutation_statistics.npz', **{m: r['permuted'] for m, r in results.items()})
    return rows

def describe_grid(distances, layers, ids, setting_ids):
    n = len(ids)
    ix = np.triu_indices(n, 1)
    summary = []
    loo = []
    for name, d in zip(setting_ids, distances):
        for layer in ('mash_distance', 'auxiliary_distance'):
            for method in ('pearson', 'spearman'):
                r = correlation(d[ix], layers[layer][ix], method)
                summary.append({'setting': name, 'reference': layer, 'method': method, 'assemblies': n, 'pairs': len(ix[0]), 'correlation': r})
        for k, a in enumerate(ids):
            mask = (ix[0] != k) & (ix[1] != k)
            loo.append({'setting': name, 'omitted_assembly': a, 'pairs': int(mask.sum()), 'pearson': correlation(d[ix][mask], layers['mash_distance'][ix][mask], 'pearson'), 'spearman': correlation(d[ix][mask], layers['mash_distance'][ix][mask], 'spearman')})
    return (summary, loo)
