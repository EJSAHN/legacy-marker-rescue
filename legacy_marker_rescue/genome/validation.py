"""Validation analyses on cached and regenerated RAPD products."""
from __future__ import annotations
import json, math, re
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from .numerics import row_products
from .common import read_tsv, write_tsv, write_json
from .stats import correlation, jaccard, load_band_matrix, load_distance_layers, make_permutations, permutation_family, write_family, describe_grid, row_z
from .sequence import pooled_matrix, training_frozen_matrix, fixed_log_matrix, select_sizes

def cached_analysis(baseline, out, inventory, reps, config, log):
    out.mkdir(parents=True, exist_ok=True)
    ids = [r['assembly'] for r in inventory]
    settings = [c['id'] for c in config['candidate_settings']]
    layers = load_distance_layers(baseline / 'distance_pairs.tsv', ids)
    mats = []
    bits = []
    features = []
    for setting in settings:
        ff, bb = load_band_matrix(baseline / f'bands_{setting}.tsv', ids)
        features.append(ff)
        bits.append(bb)
        mats.append(jaccard(bb))
    mats = np.array(mats)
    ix = np.triu_indices(len(ids), 1)
    selected = settings.index(config['reported_setting'])
    error = float(np.max(np.abs(mats[selected][ix] - layers['rapd_distance'][ix])))
    if error > 5e-08:
        raise ValueError(f'Cached binary profiles disagree with submitted pairwise RAPD distances ({error})')
    summaries, loo = describe_grid(mats, layers, ids, settings)
    write_tsv(out / 'descriptive_grid.tsv', summaries)
    write_tsv(out / 'leave_one_assembly_out.tsv', loo)
    status = {'assemblies': len(ids), 'pairs': len(ix[0]), 'biological_label_groups': len(reps), 'cached_matrix_pairwise_max_abs_difference': error, 'reported_setting_pearson': correlation(mats[selected][ix], layers['mash_distance'][ix], 'pearson'), 'reported_setting_spearman': correlation(mats[selected][ix], layers['mash_distance'][ix], 'spearman'), 'complete_historical_search_known': False, 'search_family_scope': 'five retained permissive settings only; undocumented past exploration is not adjusted', 'historical_parameters_already_inspected': True}
    submats = mats[:, reps][:, :, reps]
    target = layers['mash_distance'][np.ix_(reps, reps)]
    p = make_permutations(len(reps), config['permutations'], config['permutation_seed'])
    log(f'  matrix-label permutations: {len(p):,}; candidate settings: {len(settings)}; sample representatives: {len(reps)}')
    result = permutation_family(submats, target, p)
    tests = write_family(out / 'representative_label_permutations', settings, result, 'metadata-defined sample representatives', len(reps))
    np.savez_compressed(out / 'representative_label_permutations' / 'label_permutations.npz', permutation=p, representatives=reps)
    blocks = [inventory[i].get('submitter') or f'UNKNOWN_SINGLETON_{i}' for i in reps]
    if any((blocks.count(b) > 1 for b in blocks)):
        bp = make_permutations(len(reps), config['permutations'], config['permutation_seed'] + 1, blocks)
        restricted = permutation_family(submats, target, bp)
        write_family(out / 'within_submitter_permutations', settings, restricted, 'within documented submitter', len(reps), 'exchangeable within submitter blocks; unknown provenance remains a singleton')
        write_tsv(out / 'submitter_blocks.tsv', [{'assembly': ids[i], 'block': b} for i, b in zip(reps, blocks)])
    aux = layers['auxiliary_distance'][np.ix_(reps, reps)]
    auxresult = permutation_family(submats, aux, p)['spearman']
    mashresult = result['spearman']
    chosen = int(np.nanargmax(auxresult['observed']))
    permchosen = np.nanargmax(auxresult['permuted'], axis=1)
    permstat = mashresult['permuted'][np.arange(len(p)), permchosen]
    observed = mashresult['observed'][chosen]
    write_json(out / 'auxiliary_selection_diagnostic.json', {'selection_rule': 'maximize training-free full-sample auxiliary Spearman among the five recorded settings', 'scope': 'retrospective rule sensitivity, not reconstruction of the entire original selection process', 'chosen_setting': settings[chosen], 'mash_correlation_of_chosen_setting': observed, 'p_two_sided_under_joint_target_label_permutation': (1 + int(np.count_nonzero(np.abs(permstat) >= abs(observed) - 1e-12))) / (len(p) + 1)})
    write_json(out / 'cached_summary.json', status)
    return {'ids': ids, 'layers': layers, 'matrices': mats, 'bits': bits, 'features': features, 'settings': settings, 'summary': status}

def feature_summary(features, mat):
    n = len(mat)
    counts = mat.sum(axis=0)
    p = counts / n if n else np.array([])
    total = mat.shape[1]
    categories = {'present_in_all': counts == n, 'variable_across_assemblies': (counts > 0) & (counts < n), 'prevalence_at_least_95pct': p >= 0.95, 'rare_private_at_most_5pct': (p > 0) & (p <= 0.05), 'intermediate_5_to_80pct': (p > 0.05) & (p < 0.8), 'near_core_80_to_95pct': (p >= 0.8) & (p < 0.95)}
    out = {'band_bins': total, 'mean_bands_per_assembly': float(mat.sum(axis=1).mean())}
    for name, mask in categories.items():
        out[name + '_count'] = int(mask.sum())
        out[name + '_fraction'] = float(mask.mean()) if total else np.nan
    return out

def write_binary(path, features, mat, ids):
    names = [f'{p}|bin_{i + 1:03d}_{c:g}bp' for i, (p, c) in enumerate(features)]
    write_tsv(path, [{'assembly': a, **{name: int(v) for name, v in zip(names, mat[i])}} for i, a in enumerate(ids)], ['assembly'] + names)

def bridge_matches(features, mat, legacy_rows, tolerance):
    sizes = np.array([float(c) / 1000 for p, c in features])
    counts = mat.sum(axis=0)
    n = mat.shape[0]
    result = []
    for old in legacy_rows:
        size = float(old['legacy_kb'])
        matches = np.abs(sizes - size) <= tolerance + 1e-12
        variable = matches & (counts > 0) & (counts < n)
        fixed = matches & (counts == n)
        result.append({'legacy_band': old['legacy_band'], 'legacy_size_kb': size, 'tolerance_kb': tolerance, 'modern_matches': int(matches.sum()), 'variable_matches': int(variable.sum()), 'fixed_matches': int(fixed.sum()), 'class': 'variable_match' if variable.any() else 'fixed_only' if matches.any() else 'no_match'})
    return result

def bridge_summary(matches):
    return {'legacy_bins': len(matches), 'any_match': sum((r['modern_matches'] > 0 for r in matches)), 'variable_match': sum((r['variable_matches'] > 0 for r in matches)), 'fixed_only': sum((r['class'] == 'fixed_only' for r in matches)), 'no_match': sum((r['class'] == 'no_match' for r in matches))}

def analyze_products(amps, targets, inventory, reps, baseline, cached, config, out, log):
    out.mkdir(parents=True, exist_ok=True)
    ids = cached['ids']
    n = len(ids)
    ix = np.triu_indices(n, 1)
    layers = cached['layers']
    legacy_rows = read_tsv(baseline / 'legacy_bridge.tsv')
    setids = list(dict.fromkeys((t['set_id'] for t in targets)))
    maps = {s: {t['primer']: i for i, t in enumerate(targets) if t['set_id'] == s} for s in setids}
    historical = setids[0]
    oldsummary = {r['set_id']: r for r in read_tsv(baseline / 'random_primer_previous_summary.tsv')}
    oldhistorical = next((r for r in oldsummary.values() if r['set_type'] == 'legacy'))
    metrics = []
    bridge_rows = []
    convention = []
    new_historical_by_setting = []
    historical_sizes = {}
    for si, sid in enumerate(setids):
        log(f'  band matrices and Mash comparison [{si + 1}/{len(setids)}]: {sid}')
        modes = [('legacy_both_min100', 'both', 100), ('legacy_both_min11', 'both', 11), ('inward_min100', 'inward', config['min_amplicon_bp'])]
        for mode, geometry, minbp in modes:
            sr = select_sizes(amps, maps[sid], geometry, minbp, config['max_amplicon_bp'])
            count = sum((len(ss) for rec in sr for ss in rec.values()))
            if count > config['max_amplicons_per_set']:
                raise RuntimeError(f'Fixed set {sid} exceeds resource cap; not replaced')
            feats, mat = pooled_matrix(sr, 0.05)
            d = jaccard(mat)
            desc = feature_summary(feats, mat)
            row = {'set_id': sid, 'set_type': 'historical' if sid == historical else 'GC_matched_random', 'mode': mode, 'assemblies': n, 'products': count, **desc}
            for method in ('pearson', 'spearman'):
                row['mash_' + method] = correlation(d[ix], layers['mash_distance'][ix], method)
                row['auxiliary_' + method] = correlation(d[ix], layers['auxiliary_distance'][ix], method)
                ii = np.triu_indices(len(reps), 1)
                row['representative_mash_' + method] = correlation(d[np.ix_(reps, reps)][ii], layers['mash_distance'][np.ix_(reps, reps)][ii], method)
            metrics.append(row)
            write_binary(out / 'band_matrices' / mode / (sid + '.tsv'), feats, mat, ids)
            np.savez_compressed(out / 'band_matrices' / mode / (sid + '.npz'), presence=mat, distances=d, ids=np.array(ids), size_bp=np.array([x[1] for x in feats]), primer=np.array([x[0] for x in feats]))
            if mode == 'legacy_both_min100' and sid == historical:
                delta = float(np.nanmax(np.abs(d - cached['matrices'][cached['settings'].index(config['reported_setting'])])))
                old = oldhistorical
                convention.append({'set_id': sid, 'comparison': 'regenerated main convention vs archived main', 'previous_products': int(float(old['amplicons_total'])), 'regenerated_products': count, 'product_difference': count - int(float(old['amplicons_total'])), 'max_abs_pair_distance_difference': delta, 'note': 'full cohort bins are re-fitted; archived 22-row profiles may retain pre-filter bin boundaries'})
            if mode == 'legacy_both_min11' and sid != historical and (sid in oldsummary):
                old = oldsummary[sid]
                convention.append({'set_id': sid, 'comparison': 'old random convention vs archived random', 'previous_products': int(float(old['amplicons_total'])), 'regenerated_products': count, 'product_difference': count - int(float(old['amplicons_total'])), 'auxiliary_spearman_difference': row['auxiliary_spearman'] - float(old['all_pairs_spearman_r'])})
            if mode == 'inward_min100':
                for tol in config['bridge_tolerances_kb']:
                    matches = bridge_matches(feats, mat, legacy_rows, tol)
                    bridge_rows.append({'set_id': sid, 'set_type': row['set_type'], 'mode': mode, 'tolerance_kb': tol, **bridge_summary(matches)})
                    if sid == historical:
                        write_tsv(out / 'bridge' / f'historical_matches_tol_{tol:g}.tsv', matches)
                if sid == historical:
                    for setting in config['candidate_settings']:
                        filtered = [{p: [x for x in ss if x <= setting['max_amplicon_bp']] for p, ss in rec.items()} for rec in sr]
                        ff, bb = pooled_matrix(filtered, setting['relative_tolerance'])
                        historical_sizes[setting['id']] = filtered
                        new_historical_by_setting.append(jaccard(bb))
                        write_binary(out / 'inward_sensitivity_matrices' / (setting['id'] + '.tsv'), ff, bb, ids)
    write_tsv(out / 'primer_set_metrics.tsv', metrics)
    write_tsv(out / 'legacy_convention_reproduction.tsv', convention)
    write_tsv(out / 'bridge' / 'all_primer_sets.tsv', bridge_rows)
    tests = []
    for mode in sorted({r['mode'] for r in metrics}):
        hist = next((r for r in metrics if r['set_id'] == historical and r['mode'] == mode))
        rr = [r for r in metrics if r['set_type'] == 'GC_matched_random' and r['mode'] == mode]
        for metric in ('mash_spearman', 'mash_pearson', 'representative_mash_spearman', 'products', 'band_bins', 'variable_across_assemblies_fraction', 'prevalence_at_least_95pct_fraction'):
            observed = float(hist[metric])
            vals = np.array([r[metric] for r in rr], dtype=float)
            valid = np.isfinite(vals)
            if valid.sum() != len(vals) or not np.isfinite(observed):
                pupper = plower = np.nan
            else:
                pupper = (1 + np.count_nonzero(vals >= observed - 1e-12)) / (len(vals) + 1)
                plower = (1 + np.count_nonzero(vals <= observed + 1e-12)) / (len(vals) + 1)
            tests.append({'mode': mode, 'metric': metric, 'historical_value': observed, 'random_mean': np.nanmean(vals), 'random_sd': np.nanstd(vals, ddof=1), 'random_sets': len(vals), 'finite_random_sets': int(valid.sum()), 'empirical_upper': pupper, 'empirical_lower': plower, 'interpretation': 'diagnostic rank in fixed GC-matched set collection; not equivalence proof'})
    write_tsv(out / 'random_primer_mash_comparison.tsv', tests)
    btests = []
    for tol in config['bridge_tolerances_kb']:
        hist = next((r for r in bridge_rows if r['set_id'] == historical and r['tolerance_kb'] == tol))
        random = [r for r in bridge_rows if r['set_id'] != historical and r['tolerance_kb'] == tol]
        for metric in ('any_match', 'variable_match'):
            values = np.array([r[metric] for r in random])
            btests.append({'tolerance_kb': tol, 'metric': metric, 'historical': hist[metric], 'random_mean': float(values.mean()), 'random_sets': len(values), 'empirical_upper': (1 + int((values >= hist[metric]).sum())) / (len(values) + 1), 'interpretation': 'size compatibility null; tolerance search is descriptive, not independent hypothesis tests'})
    write_tsv(out / 'bridge' / 'random_primer_size_compatibility.tsv', btests)
    mats = np.array(new_historical_by_setting)
    rep_mats = mats[:, reps][:, :, reps]
    perms = make_permutations(len(reps), config['permutations'], config['permutation_seed'])
    result = permutation_family(rep_mats, layers['mash_distance'][np.ix_(reps, reps)], perms)
    write_family(out / 'inward_matrix_permutations', cached['settings'], result, 'inward-only products on sample representatives', len(reps))
    summary = {'historical_set_id': historical, 'primer_sets': len(setids), 'recorded_random_sets': len(setids) - 1, 'legacy_product_count_disagreements': sum((r.get('product_difference', 0) != 0 for r in convention)), 'primary_product_model': 'inward-only; min 100 bp; max 3000 bp; up to one mismatch per site', 'old_main_model': 'both F_to_RC and RC_to_F; min 100 bp', 'old_random_model': 'both F_to_RC and RC_to_F; min 11 bp implied by site order', 'interpretation': 'Inward-facing products are primary; other geometries are diagnostic comparisons.'}
    write_json(out / 'product_analysis_summary.json', summary)
    return (historical_sizes, summary)

def split_sample_analysis(size_by_setting, target, ids, config, out, log):
    """Retrospective internal holdout with train-only binning and joint-label null.

All splits are locked before looking at correlations. Never split pairwise edges.
Statistics across splits are descriptive; a single joint permutation tests the
median holdout statistic, recomputing the setting choice within every split.
"""
    out.mkdir(parents=True, exist_ok=True)
    n = len(ids)
    settings = config['candidate_settings']
    testn = max(4, int(round(n * config['holdout_fraction'])))
    if n - testn < 5:
        raise ValueError('Too few independent representatives for split evaluation')
    rng = np.random.default_rng(config['split_seed'])
    splits = []
    seen = set()
    while len(splits) < config['holdout_splits']:
        test = tuple(sorted((int(x) for x in rng.choice(n, testn, replace=False))))
        if test in seen:
            continue
        seen.add(test)
        train = tuple((i for i in range(n) if i not in test))
        splits.append((train, test))
        if len(seen) >= math.comb(n, testn):
            break
    permutations = np.vstack([np.arange(n), make_permutations(n, config['holdout_permutations'], config['permutation_seed'] + 10)])
    write_tsv(out / 'split_membership.tsv', [{'split': s + 1, 'assembly': ids[i], 'partition': 'train' if i in tr else 'test'} for s, (tr, te) in enumerate(splits) for i in range(n)])
    np.savez_compressed(out / 'joint_entity_permutations.npz', permutations=permutations)
    results = []
    coverage = []
    selection_rows = []
    aggregate = {}
    for representation in ('training_frozen_greedy', 'data_independent_log_grid'):
        score_matrix = np.full((len(splits), len(permutations)), np.nan)
        for s, (tr, te) in enumerate(splits):
            if s % 10 == 0:
                log(f'  split-sample {representation}: {s + 1}/{len(splits)}')
            tr = np.array(tr)
            te = np.array(te)
            it, jt = np.triu_indices(len(tr), 1)
            iv, jv = np.triu_indices(len(te), 1)
            train_vectors = []
            test_vectors = []
            test_unassigned = []
            for setting in settings:
                sr = size_by_setting[setting['id']]
                if representation == 'training_frozen_greedy':
                    x, audit = training_frozen_matrix(sr, tr, setting['relative_tolerance'])
                    denom = sum((r['total_products'] for r in audit if r['row'] in te))
                    unused = sum((r['unassigned_products'] for r in audit if r['row'] in te))
                    coverage.append({'split': s + 1, 'setting': setting['id'], 'test_products': denom, 'unassigned_test_products': unused, 'unassigned_fraction': unused / denom if denom else np.nan})
                    test_unassigned.append(unused / denom if denom else np.nan)
                else:
                    x = fixed_log_matrix(sr, setting['relative_tolerance'], config['min_amplicon_bp'])
                    test_unassigned.append(0.0)
                d = jaccard(x)
                train_vectors.append(d[tr[it], tr[jt]])
                test_vectors.append(d[te[iv], te[jv]])
            xtrain = np.array(train_vectors)
            xtest = np.array(test_vectors)
            ztrain = row_z(rankdata(xtrain, axis=1, method='average'))
            ztest = row_z(rankdata(xtest, axis=1, method='average'))
            ztest_pearson = row_z(xtest)
            ytrain = target[permutations[:, tr[it]], permutations[:, tr[jt]]]
            ytest = target[permutations[:, te[iv]], permutations[:, te[jv]]]
            ts = row_products(row_z(rankdata(ytrain, axis=1, method='average')), ztrain)
            vs = row_products(row_z(rankdata(ytest, axis=1, method='average')), ztest)
            vp = row_products(row_z(ytest), ztest_pearson)
            ts[~np.isfinite(ts)] = -np.inf
            chosen = np.argmax(ts, axis=1)
            valid = np.max(ts, axis=1) > -np.inf
            score_matrix[s, valid] = vs[np.arange(len(permutations))[valid], chosen[valid]]
            k = int(chosen[0])
            obs_test = vs[0, k] if valid[0] else np.nan
            results.append({'representation': representation, 'split': s + 1, 'train_assemblies': len(tr), 'test_assemblies': len(te), 'train_pairs': len(it), 'test_pairs': len(iv), 'selected_setting': settings[k]['id'] if valid[0] else 'UNDEFINED', 'training_spearman': ts[0, k] if valid[0] else np.nan, 'test_spearman': obs_test, 'test_pearson': vp[0, k] if valid[0] else np.nan, 'unassigned_test_product_fraction': test_unassigned[k] if valid[0] else np.nan})
            for k, setting in enumerate(settings):
                selection_rows.append({'representation': representation, 'split': s + 1, 'setting': setting['id'], 'training_spearman': ts[0, k] if np.isfinite(ts[0, k]) else np.nan, 'test_spearman': vs[0, k], 'selected': int(k == chosen[0] and valid[0])})
        valid_splits = np.isfinite(score_matrix).all(axis=1)
        if valid_splits.all():
            medians = np.median(score_matrix, axis=0)
            p = (1 + int(np.count_nonzero(np.abs(medians[1:]) >= abs(medians[0]) - 1e-12))) / len(permutations)
        else:
            medians = np.nanmedian(score_matrix, axis=0)
            p = np.nan
        observed = score_matrix[:, 0]
        finite = observed[np.isfinite(observed)]
        aggregate[representation] = {'splits': len(splits), 'finite_observed_splits': len(finite), 'median_test_spearman': float(np.median(finite)) if len(finite) else None, 'descriptive_5th_95th_percentiles': np.quantile(finite, [0.05, 0.95]).tolist() if len(finite) else [], 'fraction_positive_descriptive': float((finite > 0).mean()) if len(finite) else None, 'joint_label_permutation_p_two_sided': p, 'undefined_comparison_policy': 'p-value withheld if any split or permutation is undefined', 'scope': 'retrospective internal resampling; dependent splits; not external cohort validation'}
        np.savez_compressed(out / (representation + '_all_null_statistics.npz'), split_test_correlations=score_matrix, median_statistics=medians)
    write_tsv(out / 'heldout_results.tsv', results)
    write_tsv(out / 'all_candidate_scores.tsv', selection_rows)
    write_tsv(out / 'heldout_binning_coverage.tsv', coverage)
    write_json(out / 'heldout_summary.json', aggregate)
    return aggregate
