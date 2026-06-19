from __future__ import annotations
import math


def jaccard_distance(a: list[int], b: list[int]) -> float:
    intersection = sum(1 for x, y in zip(a, b) if x == 1 and y == 1)
    union = sum(1 for x, y in zip(a, b) if x == 1 or y == 1)
    return 0.0 if union == 0 else 1.0 - intersection / union


def pairwise_band_distances(matrix_rows: list[dict]) -> list[dict]:
    bands = [k for k in matrix_rows[0].keys() if k != "assembly"] if matrix_rows else []
    rows = []
    for i in range(len(matrix_rows)):
        for j in range(i + 1, len(matrix_rows)):
            a = [int(matrix_rows[i][b]) for b in bands]
            b = [int(matrix_rows[j][b]) for b in bands]
            rows.append({
                "assembly_a": matrix_rows[i]["assembly"],
                "assembly_b": matrix_rows[j]["assembly"],
                "rapd_distance": jaccard_distance(a, b),
            })
    return rows


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n == 0:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a-mx)*(b-my) for a,b in zip(x,y))
    denx = math.sqrt(sum((a-mx)**2 for a in x))
    deny = math.sqrt(sum((b-my)**2 for b in y))
    return float("nan") if denx == 0 or deny == 0 else num / (denx * deny)


def rank(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j+1]] == values[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rank(x), rank(y))
