"""Recalculate Mash distances when an existing native or WSL installation is available."""
from __future__ import annotations
import os, shutil, subprocess
from pathlib import Path
import numpy as np
from .common import accession, write_tsv, write_json, sha256

def write_mash_file_list(listing: Path, file_paths):
    """Write one unquoted path per LF-terminated line on every host OS.

    The list is consumed by Mash, including Linux Mash launched through WSL.
    Binary UTF-8 output deliberately avoids Windows text-mode CRLF translation.
    Spaces and non-ASCII characters in paths are preserved.
    """
    lines = [str(value) for value in file_paths]
    if not lines:
        raise ValueError('Mash input file list is empty')
    for value in lines:
        if not value or any((char in value for char in ('\r', '\n', '\x00'))):
            raise ValueError('Mash input path is empty or contains a line/control character')
    payload = ('\n'.join(lines) + '\n').encode('utf-8')
    if payload.startswith(b'\xef\xbb\xbf'):
        raise ValueError('Mash input file list must not start with a UTF-8 BOM')
    listing.write_bytes(payload)
    stored = listing.read_bytes()
    if stored != payload or b'\r' in stored or stored.count(b'\n') != len(lines):
        raise RuntimeError('Mash input file-list byte verification failed')
    return {'entries': len(lines), 'line_ending': 'LF', 'encoding': 'UTF-8 without BOM', 'carriage_return_bytes': 0, 'file_list_sha256': sha256(listing)}

def refresh(paths, ids, out, log):
    out.mkdir(parents=True, exist_ok=True)
    executable = shutil.which('mash')
    prefix = []
    mode = 'native'
    if executable:
        prefix = [executable]
    elif os.name == 'nt' and shutil.which('wsl.exe'):
        try:
            probe = subprocess.run(['wsl.exe', '--exec', 'mash', '--version'], capture_output=True, text=True, timeout=20)
            if probe.returncode == 0:
                prefix = ['wsl.exe', '--exec', 'mash']
                mode = 'WSL'
        except (OSError, subprocess.TimeoutExpired):
            pass
    if not prefix:
        return (None, {'status': 'NOT_RUN_NO_EXISTING_MASH', 'reason': 'No native or default-WSL Mash found. Cached Mash distances were retained explicitly; no installation attempted.'})

    def native_path(p):
        p = p.resolve()
        if mode != 'WSL':
            return str(p)
        proc = subprocess.run(['wsl.exe', '--exec', 'wslpath', '-a', str(p)], capture_output=True, text=True, timeout=20)
        if proc.returncode:
            raise RuntimeError('WSL could not translate a project path')
        return proc.stdout.strip()
    version = subprocess.run(prefix + ['--version'], capture_output=True, text=True, timeout=20)
    if version.returncode:
        raise RuntimeError('Existing Mash executable failed its version check')
    listing = out / 'genome_list.txt'
    list_audit = write_mash_file_list(listing, [native_path(paths[a]) for a in ids])
    list_audit['execution_mode'] = mode
    write_json(out / 'file_list_validation.json', list_audit)
    log(f'  Mash input list verified: {len(ids)} paths; LF only; no carriage-return bytes')
    sketchprefix = out / 'genomes'
    sketch = str(sketchprefix) + '.msh'
    cmd = prefix + ['sketch', '-k', '21', '-s', '10000', '-p', '1', '-l', native_path(listing), '-o', native_path(sketchprefix)]
    log('  recalculating Mash sketches from the exact local FASTA files (k=21, sketch=10000)')
    with (out / 'mash_console.txt').open('w', encoding='utf-8') as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if p.returncode:
            raise RuntimeError(f'Mash sketch failed (exit {p.returncode}); see mash_console.txt')
    distfile = out / 'all_pairs_raw.tsv'
    with distfile.open('w', encoding='utf-8') as f, (out / 'mash_console.txt').open('a', encoding='utf-8') as e:
        p = subprocess.run(prefix + ['dist', native_path(Path(sketch)), native_path(Path(sketch))], stdout=f, stderr=e)
        if p.returncode:
            raise RuntimeError(f'Mash dist failed (exit {p.returncode})')
    pos = {a: i for i, a in enumerate(ids)}
    dist = np.full((len(ids), len(ids)), np.nan)
    for line in distfile.read_text().splitlines():
        r = line.split('\t')
        if len(r) < 3:
            raise ValueError('Incomplete Mash output line')
        a = accession(r[0])
        b = accession(r[1])
        v = float(r[2])
        i = pos[a]
        j = pos[b]
        if np.isfinite(dist[i, j]) and abs(dist[i, j] - v) > 1e-12:
            raise ValueError('Conflicting Mash duplicate row')
        dist[i, j] = v
    if not np.isfinite(dist).all() or not np.allclose(dist, dist.T) or (not np.allclose(np.diag(dist), 0)):
        raise ValueError('Refreshed Mash output is incomplete or asymmetric')
    write_tsv(out / 'refreshed_mash_pairs.tsv', [{'assembly_a': ids[i], 'assembly_b': ids[j], 'mash_distance': dist[i, j]} for i in range(len(ids)) for j in range(i + 1, len(ids))])
    distfile.unlink()
    listing.unlink()
    Path(sketch).unlink(missing_ok=True)
    return (dist, {'status': 'COMPLETE', 'mash_version': version.stdout.strip() or version.stderr.strip(), 'execution_mode': mode, 'k': 21, 'sketch_size': 10000, 'threads': 1, 'reference_for_regenerated_analysis': 'freshly recalculated Mash'})
