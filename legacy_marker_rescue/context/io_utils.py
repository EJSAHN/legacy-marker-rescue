"""Strict, deterministic I/O for a non-destructive follow-up stage."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import zipfile

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_json(path: Path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def write_json(path: Path, value) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + '\n', encoding='utf-8')

def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f'Missing or duplicate headers: {path.name}')
        rows = list(reader)
        if any((None in r or any((v is None for v in r.values())) for r in rows)):
            raise ValueError(f'Ragged TSV: {path.name}')
        return rows

def scalar(v):
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return format(v, '.17g') if math.isfinite(v) else 'NA'
    if isinstance(v, bool):
        return int(v)
    return v

def write_tsv(path: Path, rows: list[dict], columns=None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or list(dict.fromkeys((k for r in rows for k in r))) or ['status']
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='raise', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: scalar(r.get(k)) for k in columns})

def write_rows_gzip(path: Path, rows: list[dict], columns=None) -> None:
    path = Path(path)
    columns = columns or list(dict.fromkeys((k for r in rows for k in r))) or ['status']
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='raise', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: scalar(r.get(k)) for k in columns})

def read_rows_gzip(path: Path) -> list[dict]:
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def fasta_records(path: Path):
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii', errors='strict') as f:
        name, chunks, seen = (None, [], set())
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield (name, ''.join(chunks).upper())
                name = line[1:].split()[0]
                if not name or name in seen:
                    raise ValueError('Empty or duplicated FASTA identifier')
                seen.add(name)
                chunks = []
            else:
                if name is None:
                    raise ValueError('FASTA sequence before the first header')
                if re.search('[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]', line):
                    raise ValueError(f'Unexpected FASTA character in record {name}')
                chunks.append(line)
        if name is not None:
            yield (name, ''.join(chunks).upper())

def relative_to_root(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())

def checked_relative(root: Path, rel: str) -> Path:
    path = (root / rel).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Path escapes the specified root')
    return path
