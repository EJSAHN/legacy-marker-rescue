from __future__ import annotations
from pathlib import Path
import csv


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def iter_fasta_records(path: str | Path):
    name = None
    seq = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()
