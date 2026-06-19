from __future__ import annotations
from pathlib import Path
from .io import iter_fasta_records
from .sequence import reverse_complement, hamming_distance


def find_binding_sites(sequence: str, primer: str, max_mismatch: int) -> list[int]:
    primer = primer.upper()
    k = len(primer)
    seq = sequence.upper()
    sites = []
    for i in range(0, len(seq) - k + 1):
        window = seq[i:i+k]
        if "N" in window:
            continue
        if hamming_distance(window, primer) <= max_mismatch:
            sites.append(i)
    return sites


def simulate_rapd_for_fasta(fasta_path: str | Path, primers: dict[str, str], max_mismatch: int, max_amplicon_bp: int) -> list[dict]:
    rows = []
    assembly = Path(fasta_path).name
    for contig, seq in iter_fasta_records(fasta_path):
        for primer_name, primer in primers.items():
            fwd = find_binding_sites(seq, primer, max_mismatch)
            rev = find_binding_sites(seq, reverse_complement(primer), max_mismatch)
            for left in fwd:
                for right in rev:
                    if right <= left:
                        continue
                    size = right + len(primer) - left
                    if 0 < size <= max_amplicon_bp:
                        rows.append({
                            "assembly": assembly,
                            "contig": contig,
                            "primer": primer_name,
                            "start": left + 1,
                            "end": right + len(primer),
                            "amplicon_size_bp": size,
                        })
    return rows


def band_bin(size_bp: int, relative_tolerance: float) -> int:
    if size_bp <= 0:
        return 0
    # Log-scale relative binning; neighboring bins are approximately relative_tolerance apart.
    import math
    base = math.log1p(relative_tolerance)
    return int(round(math.log(size_bp) / base))


def build_band_matrix(amplicons: list[dict], relative_tolerance: float) -> tuple[list[dict], list[dict]]:
    assemblies = sorted({r["assembly"] for r in amplicons})
    bins = {}
    for r in amplicons:
        b = band_bin(int(r["amplicon_size_bp"]), relative_tolerance)
        key = f"{r['primer']}|B{b}"
        bins.setdefault(key, []).append(int(r["amplicon_size_bp"]))
    band_names = sorted(bins)
    matrix_rows = []
    for asm in assemblies:
        present = {f"{r['primer']}|B{band_bin(int(r['amplicon_size_bp']), relative_tolerance)}" for r in amplicons if r["assembly"] == asm}
        row = {"assembly": asm}
        for band in band_names:
            row[band] = 1 if band in present else 0
        matrix_rows.append(row)
    summary_rows = []
    for band in band_names:
        sizes = bins[band]
        summary_rows.append({
            "band": band,
            "mean_size_bp": round(sum(sizes) / len(sizes), 2),
            "n_amplicons": len(sizes),
        })
    return matrix_rows, summary_rows
