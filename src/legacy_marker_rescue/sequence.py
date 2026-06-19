from __future__ import annotations

DNA_ALPHABET = "ACGT"
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1].upper()


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Sequences must have equal length")
    return sum(x != y for x, y in zip(a.upper(), b.upper()))


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    bases = sum(1 for x in seq if x in "ACGT")
    return 0.0 if bases == 0 else sum(1 for x in seq if x in "GC") / bases
