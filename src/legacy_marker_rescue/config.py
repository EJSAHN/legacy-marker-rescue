from __future__ import annotations
from pathlib import Path
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError(f"Configuration file did not parse as a mapping: {path}")
    return cfg


def resolve_path(path: str | Path, root: str | Path = ".") -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return Path(root) / path
