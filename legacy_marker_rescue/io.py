"""Validated file access and content-addressed output records."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> Any:
    with Path(path).open(encoding='utf-8-sig') as handle:
        return json.load(handle)


def write_json(path: Path, value: Any) -> None:
    import numpy as np
    def convert(x):
        if isinstance(x, dict): return {str(k): convert(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [convert(v) for v in x]
        if isinstance(x, Path): return str(x)
        if isinstance(x, np.ndarray): return convert(x.tolist())
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, (float, np.floating)):
            return float(x) if math.isfinite(float(x)) else None
        return x
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(convert(value), ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def read_table(path: Path) -> list[dict]:
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields) or any(not f for f in fields):
            raise ValueError(f'Invalid table header: {path.name}')
        rows = list(reader)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError(f'Ragged table: {path.name}')
        return rows


def write_manifest(root: Path) -> dict[str, str]:
    root = Path(root)
    entries = {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob('*'))
               if p.is_file() and p.name != 'SHA256SUMS.txt' and '__pycache__' not in p.parts}
    (root / 'SHA256SUMS.txt').write_text(''.join(f'{h}  {rel}\n' for rel, h in entries.items()), encoding='utf-8')
    return entries


def verify_manifest(root: Path) -> dict[str, str]:
    root = Path(root).resolve()
    manifest = root / 'SHA256SUMS.txt'
    entries = {}
    for line in manifest.read_text(encoding='utf-8-sig').splitlines():
        if not line.strip(): continue
        checksum, name = line.split('  ', 1)
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or '\\' in name or ':' in name:
            raise ValueError(f'Unsafe manifest entry: {name}')
        if name in entries: raise ValueError(f'Duplicate manifest entry: {name}')
        path = (root / relative).resolve()
        if root not in path.parents or not path.is_file() or digest(path) != checksum:
            raise ValueError(f'Content mismatch: {name}')
        entries[name] = checksum
    if not entries: raise ValueError('Empty content manifest')
    return entries


def resolve_path(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()
