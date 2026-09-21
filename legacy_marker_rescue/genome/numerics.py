"""Small-array reductions without dispatch to a BLAS matrix-product backend.

The calculations are ordinary float64 sums of products. This changes neither
correlation definitions nor permutation order, candidate sets, or decision rules.
"""
from __future__ import annotations
import numpy as np

def sum_products(left, right) -> float:
    """Inner product of equal-length real vectors, using NumPy reductions."""
    a = np.asarray(left, dtype=np.float64)
    b = np.asarray(right, dtype=np.float64)
    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError('Expected equally sized one-dimensional vectors')
    return float(np.sum(a * b, dtype=np.float64))

def row_products(left, right) -> np.ndarray:
    """All pairwise row inner products, without calling dot or matmul.

Only one left-shaped temporary is created at a time. This avoids both a
BLAS dependency for these small products and a large three-dimensional array.
"""
    a = np.asarray(left, dtype=np.float64)
    b = np.asarray(right, dtype=np.float64)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError('Row-product matrices must have the same column count')
    result = np.empty((a.shape[0], b.shape[0]), dtype=np.float64)
    for column in range(b.shape[0]):
        result[:, column] = np.sum(a * b[column], axis=1, dtype=np.float64)
    return result

def euclidean_norm(values, axis=None, keepdims=False):
    """Euclidean norm by elementwise squares and reduction."""
    a = np.asarray(values, dtype=np.float64)
    return np.sqrt(np.sum(a * a, axis=axis, keepdims=keepdims, dtype=np.float64))
