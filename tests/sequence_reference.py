"""Archived implementations used only for regression tests, not as the revised engine."""
import bisect,itertools,re
from typing import Dict,List,Tuple
BASES="ACGT"

def generate_variants(primer: str, max_mismatch: int) -> Dict[str, int]:
    """Return sequence -> minimum mismatch count for all variants up to max_mismatch."""
    primer = primer.upper()
    variants: Dict[str, int] = {primer: 0}
    if max_mismatch <= 0:
        return variants
    k = len(primer)
    for m in range(1, max_mismatch + 1):
        for positions in itertools.combinations(range(k), m):
            choices = []
            for pos in positions:
                choices.append([b for b in BASES if b != primer[pos]])
            for replacements in itertools.product(*choices):
                seq = list(primer)
                for pos, base in zip(positions, replacements):
                    seq[pos] = base
                s = "".join(seq)
                variants[s] = min(variants.get(s, m), m)
    return variants

def compile_variant_regex(variants: Dict[str, int]) -> re.Pattern:
    # All alternatives have equal length, but sorting makes exact primer preferred in group output.
    alts = sorted(variants, key=lambda x: (variants[x], x))
    return re.compile("(?=(" + "|".join(re.escape(a) for a in alts) + "))")

def find_sites(seq: str, pattern: re.Pattern, variants: Dict[str, int]) -> List[Tuple[int, int, str]]:
    """Return list of (position, mismatch_count, matched_sequence)."""
    out = []
    for m in pattern.finditer(seq):
        hit = m.group(1)
        if hit in variants:
            out.append((m.start(), variants[hit], hit))
    return out

def pair_sites(
    fwd_sites: List[Tuple[int, int, str]],
    rev_sites: List[Tuple[int, int, str]],
    k: int,
    min_amplicon: int,
    max_amplicon: int,
) -> List[Tuple[int, int, int, str, int, int, str, str]]:
    """Pair forward-primer and reverse-complement-primer sites that face inward."""
    amps: List[Tuple[int, int, int, str, int, int, str, str]] = []
    fwd_sites = sorted(fwd_sites, key=lambda x: x[0])
    rev_sites = sorted(rev_sites, key=lambda x: x[0])
    fpos = [x[0] for x in fwd_sites]
    rpos = [x[0] for x in rev_sites]

    # F site upstream, RC site downstream.
    for i, fm, fseq in fwd_sites:
        j = bisect.bisect_right(rpos, i)
        while j < len(rev_sites):
            r, rm, rseq = rev_sites[j]
            size = r + k - i
            if size > max_amplicon:
                break
            if size >= min_amplicon:
                amps.append((i, r + k, size, "F_to_RC", fm, rm, fseq, rseq))
            j += 1

    # RC site upstream, F site downstream.
    for i, rm, rseq in rev_sites:
        j = bisect.bisect_right(fpos, i)
        while j < len(fwd_sites):
            r, fm, fseq = fwd_sites[j]
            size = r + k - i
            if size > max_amplicon:
                break
            if size >= min_amplicon:
                amps.append((i, r + k, size, "RC_to_F", rm, fm, rseq, fseq))
            j += 1

    return sorted(set(amps), key=lambda x: (x[0], x[1], x[3], x[4], x[5]))

def cluster_band_sizes(sizes: List[int], rel_tol: float) -> Dict[int, str]:
    sizes_sorted = sorted(set(int(x) for x in sizes))
    bins: List[List[int]] = []
    for s in sizes_sorted:
        if not bins:
            bins.append([s])
            continue
        center = sum(bins[-1]) / len(bins[-1])
        if abs(s - center) / max(center, 1.0) <= rel_tol:
            bins[-1].append(s)
        else:
            bins.append([s])
    mapping: Dict[int, str] = {}
    for idx, b in enumerate(bins, start=1):
        center = int(round(sum(b) / len(b)))
        label = f"B{idx:03d}_{center}bp"
        for s in b:
            mapping[s] = label
    return mapping
