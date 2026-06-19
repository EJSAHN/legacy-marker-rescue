from __future__ import annotations
import random
from .sequence import gc_fraction

DNA = "ACGT"


def random_10mer(gc_target: float, tolerance: float = 0.10, rng: random.Random | None = None) -> str:
    rng = rng or random.Random()
    for _ in range(10000):
        seq = "".join(rng.choice(DNA) for _ in range(10))
        if abs(gc_fraction(seq) - gc_target) <= tolerance:
            return seq
    raise RuntimeError("Could not generate GC-matched 10-mer")


def random_primer_set(reference_primers: dict[str, str], rng: random.Random | None = None) -> dict[str, str]:
    rng = rng or random.Random()
    return {name: random_10mer(gc_fraction(seq), rng=rng) for name, seq in reference_primers.items()}
