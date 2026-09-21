"""Tabular input, sequence records, and assembly identity groups."""
from __future__ import annotations
import csv, gzip, hashlib, json, math, re
from pathlib import Path
import numpy as np
ACCESSION = re.compile('GC[AF]_\\d+\\.\\d+')

def accession(value: str) -> str:
    found = ACCESSION.findall(str(value))
    if len(set(found)) != 1:
        raise ValueError(f'Expected one assembly accession in identifier: {value!r}')
    return found[0]

def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'Empty or duplicate column headings in {path.name}')
        return list(reader)

def scalar(value):
    if value is None:
        return ''
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return format(float(value), '.12g') if math.isfinite(float(value)) else 'NA'
    return value

def write_tsv(path: Path, rows: list[dict], columns=None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or list(dict.fromkeys((k for row in rows for k in row)))
    if not columns:
        columns = ['status']
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: scalar(row.get(k)) for k in columns})

def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value

def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_fasta(path: Path):
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii', errors='strict') as f:
        name = None
        chunks = []
        seen = set()
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield (name, ''.join(chunks).upper())
                name = line[1:].split()[0]
                chunks = []
                if name in seen:
                    raise ValueError(f'Duplicate FASTA record ID: {name}')
                seen.add(name)
            else:
                if name is None:
                    raise ValueError('Sequence occurs before first FASTA header')
                chunks.append(line)
        if name is not None:
            yield (name, ''.join(chunks).upper())

def metadata_groups(inventory: list[dict]):
    """Unite exact BioSample IDs and normalized explicit strain/isolate identifiers."""
    parent = list(range(len(inventory)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        parent[find(b)] = find(a)
    known = {}
    for i, row in enumerate(inventory):
        tokens = []
        if row.get('biosample'):
            tokens.append(('biosample', row['biosample'].strip().upper()))
        for col in ('strain', 'isolate'):
            v = re.sub('[^A-Z0-9]', '', row.get(col, '').upper())
            if v and v not in {'UNKNOWN', 'NA', 'NONE', 'NOTAPPLICABLE', 'NOTPROVIDED'}:
                tokens.append(('biological_label', v))
        for token in tokens:
            if token in known:
                union(i, known[token])
            else:
                known[token] = i
    groups = {}
    for i in range(len(inventory)):
        groups.setdefault(find(i), []).append(i)
    group_id = {}
    reps = []
    for ids in sorted(groups.values(), key=lambda ids: min((inventory[i]['assembly'] for i in ids))):
        representative = min(ids, key=lambda i: inventory[i]['assembly'])
        reps.append(representative)
        for i in ids:
            group_id[i] = inventory[representative]['assembly']
    rows = []
    for i, row in enumerate(inventory):
        rows.append({**row, 'analysis_group': group_id[i], 'representative': int(i in reps), 'grouping_basis': 'explicit BioSample or strain/isolate; no identity inferred from band similarity'})
    return (rows, np.array(reps, dtype=int))
