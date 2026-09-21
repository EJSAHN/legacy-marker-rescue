"""Deterministic single-primer site scanning and explicit strand geometry.

Coordinates refer to the stored FASTA strand. An F sequence upstream of its
reverse complement produces inward-facing extension. RC upstream of F is
retained ONLY as an archived-model audit, not as an inward-facing PCR product.
"""
from __future__ import annotations
import bisect, hashlib, json, time
from collections import defaultdict
from pathlib import Path
import numpy as np
from .common import read_fasta, sha256, write_json
BASES = 'ACGT'
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')

def revcomp(seq):
    return seq.upper().translate(COMPLEMENT)[::-1]

def variants(primer, max_mismatches=1):
    if max_mismatches not in (0, 1):
        raise ValueError('Only exact or one-mismatch scanning is supported')
    primer = primer.upper()
    if len(primer) != 10 or set(primer) - set(BASES):
        raise ValueError('Each RAPD primer must be an unambiguous 10-mer')
    result = {primer: 0}
    if max_mismatches:
        for i, ch in enumerate(primer):
            for b in BASES:
                if b != ch:
                    result[primer[:i] + b + primer[i + 1:]] = 1
    return result

def encode(seq):
    out = 0
    for ch in seq:
        out = out << 2 | BASES.index(ch)
    return out

class SiteScanner:
    """One sequence pass for all fixed historical and random primers."""

    def __init__(self, primers, max_mismatches=1, chunk_bp=300000):
        self.primers = list(primers)
        self.chunk_bp = chunk_bp
        self.k = 10
        self.variants = []
        self.lookup = defaultdict(list)
        for pi, primer in enumerate(primers):
            f = variants(primer, max_mismatches)
            r = {revcomp(v): mm for v, mm in f.items()}
            self.variants.append((f, r))
            for side, vs in enumerate((f, r)):
                for v, mm in vs.items():
                    self.lookup[encode(v)].append((pi, side, mm))
        self.present = np.zeros(1 << 2 * self.k, dtype=bool)
        self.present[list(self.lookup)] = True
        self.trans = np.full(256, 4, dtype=np.uint8)
        for i, ch in enumerate(BASES):
            self.trans[ord(ch)] = i
            self.trans[ord(ch.lower())] = i

    def scan(self, seq):
        hits = [[[], []] for _ in self.primers]
        n = len(seq) - self.k + 1
        if n <= 0:
            return hits
        raw = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        for offset in range(0, n, self.chunk_bp):
            length = min(self.chunk_bp, n - offset)
            bases = self.trans[raw[offset:offset + length + self.k - 1]]
            codes = np.zeros(length, dtype=np.uint32)
            valid = np.ones(length, dtype=bool)
            for j in range(self.k):
                x = bases[j:j + length]
                codes = codes << 2 | x & 3
                valid &= x < 4
            for pos in np.flatnonzero(valid & self.present[codes]):
                for pi, side, mm in self.lookup[int(codes[pos])]:
                    hits[pi][side].append((offset + int(pos), mm))
        return hits

    def reference_scan(self, seq):
        """Independent scalar implementation used to check all target/orientation hits."""
        hits = [[[], []] for _ in self.primers]
        for pi, (f, r) in enumerate(self.variants):
            for side, vs in enumerate((f, r)):
                for v, mm in vs.items():
                    pos = seq.find(v)
                    while pos >= 0:
                        hits[pi][side].append((pos, mm))
                        pos = seq.find(v, pos + 1)
                hits[pi][side].sort()
        return hits

def pair_positions(left, right, k=10, max_bp=3000, min_bp=11):
    """Products from left and right sites; caller explicitly chooses orientation."""
    rpos = [p for p, m in right]
    for start, lmm in left:
        lo = bisect.bisect_right(rpos, start)
        lo = max(lo, bisect.bisect_left(rpos, start + min_bp - k))
        hi = bisect.bisect_right(rpos, start + max_bp - k)
        for rp, rmm in right[lo:hi]:
            end = rp + k
            yield (start, end, end - start, lmm, rmm)

def scan_assembly(path: Path, scanner: SiteScanner, target_set, config, cache_dir: Path, log):
    fingerprint = sha256(path)
    signature = hashlib.sha256(json.dumps({'engine': 'onepass_10mer_v1', 'primers': scanner.primers, 'max_bp': config['max_amplicon_bp'], 'max_mismatch': config['max_mismatches_per_site'], 'max_sites': config['max_sites_per_primer_contig'], 'max_amps': config['max_amplicons_per_set_assembly'], 'reference_bp': config['reference_scan_bases_per_assembly']}, sort_keys=True).encode()).hexdigest()
    key = hashlib.sha256((fingerprint + signature).encode()).hexdigest()
    npz = cache_dir / (key + '.npz')
    meta = cache_dir / (key + '.json')
    if npz.exists() and meta.exists():
        info = json.loads(meta.read_text())
        if info.get('npz_sha256') == sha256(npz):
            with np.load(npz, allow_pickle=False) as z:
                amps = z['amplicons'].copy()
            return (amps, info, True)
    records = []
    names = []
    lengths = []
    set_counts = defaultdict(int)
    site_totals = np.zeros((len(scanner.primers), 2), dtype=np.int64)
    checked = 0
    last = time.monotonic()
    seen = 0
    for ci, (name, seq) in enumerate(read_fasta(path)):
        names.append(name)
        lengths.append(len(seq))
        seen += len(seq)
        if len(seq) < 20:
            continue
        hits = scanner.scan(seq)
        remaining = config['reference_scan_bases_per_assembly'] - checked
        if remaining >= 20:
            part = seq[:min(len(seq), remaining)]
            if scanner.scan(part) != scanner.reference_scan(part):
                raise RuntimeError('Vectorized site scanner differs from scalar reference')
            checked += len(part)
        for pi, sides in enumerate(hits):
            for side, h in enumerate(sides):
                if len(h) > config['max_sites_per_primer_contig']:
                    raise RuntimeError(f'Binding-site resource cap exceeded; no candidate replaced ({name}, target {pi})')
                site_totals[pi, side] += len(h)
            for orientation, (left, right) in enumerate((sides, sides[::-1])):
                for start, end, size, lmm, rmm in pair_positions(left, right, max_bp=config['max_amplicon_bp']):
                    sid = target_set[pi]
                    set_counts[sid] += 1
                    if set_counts[sid] > config['max_amplicons_per_set_assembly']:
                        raise RuntimeError(f'Product resource cap exceeded for {sid}; no candidate replaced')
                    records.append((pi, ci, start, end, size, orientation, lmm, rmm))
        if time.monotonic() - last > 15:
            log(f'  scanned {ci + 1} records, {seen:,} bases, {len(records):,} candidate intervals')
            last = time.monotonic()
    if not names:
        raise ValueError('FASTA contains no records')
    amps = np.array(records, dtype=np.int64).reshape((-1, 8))
    cache_dir.mkdir(parents=True, exist_ok=True)
    tmp = npz.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, amplicons=amps)
    tmp.replace(npz)
    info = {'fasta_sha256': fingerprint, 'engine_signature': signature, 'npz_sha256': sha256(npz), 'records': len(names), 'total_bases': sum(lengths), 'contig_names': names, 'contig_lengths': lengths, 'site_totals': site_totals.tolist(), 'reference_checked_bases': checked, 'amplicon_columns': ['target_index', 'contig_index', 'start_0based', 'end_exclusive', 'size_bp', 'orientation_0_inward_1_outward', 'left_mismatches', 'right_mismatches']}
    write_json(meta, info)
    return (amps, info, False)

def select_sizes(amplicons_by_assembly, target_indices, mode, min_bp, max_bp):
    output = []
    for arr in amplicons_by_assembly:
        per = {}
        for p, pi in target_indices.items():
            mask = (arr[:, 0] == pi) & (arr[:, 4] >= min_bp) & (arr[:, 4] <= max_bp)
            if mode == 'inward':
                mask &= arr[:, 5] == 0
            elif mode != 'both':
                raise ValueError('Unknown orientation policy')
            per[p] = arr[mask, 4].astype(int).tolist()
        output.append(per)
    return output

def cluster_sizes(sizes, tolerance):
    """Archived greedy bin rule based on unique integer sizes (not counts)."""
    bins = []
    total = 0
    for size in sorted(set((int(s) for s in sizes))):
        center = total / len(bins[-1]) if bins else 0
        if bins and abs(size - center) / max(center, 1) <= tolerance:
            bins[-1].append(size)
            total += size
        else:
            bins.append([size])
            total = size
    mapping = {}
    centers = []
    for i, members in enumerate(bins):
        centers.append(float(round(sum(members) / len(members))))
        for s in members:
            mapping[s] = i
    return (mapping, centers, bins)

def pooled_matrix(size_records, tolerance):
    primers = sorted(set((p for rec in size_records for p in rec)))
    features = []
    blocks = []
    for primer in primers:
        mapping, centers, _ = cluster_sizes([s for rec in size_records for s in rec.get(primer, [])], tolerance)
        block = np.zeros((len(size_records), len(centers)), dtype=np.int8)
        for i, rec in enumerate(size_records):
            for s in rec.get(primer, []):
                block[i, mapping[s]] = 1
        blocks.append(block)
        features.extend(((primer, c) for c in centers))
    return (features, np.concatenate(blocks, axis=1) if blocks else np.zeros((len(size_records), 0), dtype=np.int8))

def training_frozen_matrix(size_records, train_indices, tolerance):
    """Fit greedy bins on training sizes only; test-only unassigned products remain explicit."""
    n = len(size_records)
    primers = sorted(set((p for i in train_indices for p in size_records[i])))
    blocks = []
    audit = []
    for primer in primers:
        mapping, centers, bins = cluster_sizes([s for i in train_indices for s in size_records[i].get(primer, [])], tolerance)
        mu = np.array([sum(b) / len(b) for b in bins])
        radii = np.array([max(max((abs(s - c) for s in b)), tolerance * c) for b, c in zip(bins, mu)])
        block = np.zeros((n, len(bins)), dtype=np.int8)
        for i, rec in enumerate(size_records):
            values = rec.get(primer, [])
            unassigned = 0
            for s in values:
                if s in mapping:
                    idx = mapping[s]
                elif len(mu):
                    valid = np.flatnonzero(np.abs(mu - s) <= radii + 1e-12)
                    idx = int(valid[np.argmin(np.abs(mu[valid] - s))]) if len(valid) else None
                else:
                    idx = None
                if idx is None:
                    unassigned += 1
                else:
                    block[i, idx] = 1
            audit.append({'row': i, 'primer': primer, 'total_products': len(values), 'unassigned_products': unassigned})
        blocks.append(block)
    return (np.concatenate(blocks, axis=1) if blocks else np.zeros((n, 0), dtype=np.int8), audit)

def fixed_log_matrix(size_records, tolerance, anchor=100):
    """Data-independent log-size grid; a separate sensitivity representation, not archived bins."""
    step = (1 + tolerance) / (1 - tolerance)
    keys = []
    for rec in size_records:
        keys.append({(p, int(np.floor(np.log(s / anchor) / np.log(step) + 1e-12))) for p, ss in rec.items() for s in ss})
    universe = sorted(set().union(*keys))
    pos = {k: j for j, k in enumerate(universe)}
    mat = np.zeros((len(keys), len(universe)), dtype=np.int8)
    for i, ks in enumerate(keys):
        for k in ks:
            mat[i, pos[k]] = 1
    return mat
