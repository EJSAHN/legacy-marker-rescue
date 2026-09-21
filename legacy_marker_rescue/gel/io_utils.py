"""Strict, portable file IO for detector inputs and derived tables."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def write_json(path: Path, data: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, ensure_ascii=True, allow_nan=False)
        handle.write('\n')

def write_table(path: Path, rows: list[dict], columns: list[str] | None=None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(dict.fromkeys((k for row in rows for k in row)))
    if not columns:
        columns = ['status']
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        for row in rows:
            clean = {}
            for key, val in row.items():
                if val is None:
                    val = 'NA'
                elif isinstance(val, (dict, list, tuple)):
                    val = json.dumps(val, ensure_ascii=True, allow_nan=False, separators=(',', ':'))
                elif isinstance(val, str):
                    val = val.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')
                clean[key] = val
            writer.writerow(clean)

def read_table(path: Path) -> list[dict]:
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if any((None in r or any((v is None for v in r.values())) for r in rows)):
        raise ValueError(f'Ragged table: {Path(path).name}')
    return rows
